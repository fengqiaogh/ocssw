/*******************************************************************************
 *
 * NAME: CldMask.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute cloud masks for a given DDProcess object class.
 *
 *  Created on: April 3, 2021
 *      Author: Sam Anderson, DT
 *
 *  Modified:
 *
 *******************************************************************************/

#include <math.h>
#include <DDAlgorithm.h>
#include <CldMask.h>

using namespace std;

/**************************************************************************
 * NAME: CldMaskLand()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

CldMaskLand::CldMaskLand() {
}

/**************************************************************************
 * NAME: CldMaskLand()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

CldMaskLand::CldMaskLand(DDAlgorithm* process) {
	p_ = process;
    mask_cirrus_ = 0;
    snowmask_ = 0;
}

/**************************************************************************
 * NAME: ~CldMaskLand()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

CldMaskLand::~CldMaskLand() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * cloud mask operations.
 *
 *************************************************************************/

int CldMaskLand::initialize() {
    return DT_SUCCESS;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Virtual function computes cloud mask over land using spatial
 *      variability of 0.47 (>0.01) and 1.38 um (>0.007) reflectance as well
 *      as absolute  value of 1.38 um > 0.1
 *
 * INPUT PARAMETERS:
 *
 *    ISWATH     Number of pixels at 1 km resolution along scan
 *    ILINE      Number of pixels at 1 km resolution against scan
 *    Refl_3     Reflectance at 0.47 um
 *    Refl_26    Reflectance at 1.38 um
 *
 * OUTPUT PARAMETERS:
 *
 *    CLDMSK_1KM      Cloud mask at 1 km resolution
 *    CldMsk_500      Cloud mask at 500m resolution
 *    CldMsk_250      Cloud mask at 250m resolution
 *
 *************************************************************************/

int CldMaskLand::compute( unsigned char& mask )
{
    int status = DT_SUCCESS;

    mask = 1;

    if (p_->rfl_[(size_t)rhot_band::W410] < filltest) {
        mask = ubytefill;
        return status;
    }
    float m03_avg     = 0;
    float m03_stddev  = 0;
    float m09_avg     = 0;
    float m09_stddev  = 0;

// Gather 3x3 grid
    float rfl[NUM_RFL_BANDS][3][3];
    memset(rfl,0,sizeof(rfl));
    for (size_t iw=0; iw<NUM_RFL_BANDS; iw++) {
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
                rfl[iw][il][ip] = p_->rfla_[iw][il][ip];
            }
        }
    }

// M03 W488
    size_t cnt = 0;
    for (size_t il=0; il<=2; il++) {
        for (size_t ip=0; ip<=2; ip++) {
            if (rfl[(size_t)rhot_band::W490][il][ip] > filltest) {
                m03_avg += rfl[(size_t)rhot_band::W490][il][ip];
                cnt++;
            }
        }
    }
    if (cnt >= 2) {
        m03_avg /= cnt;
        cnt = 0;
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
                if (rfl[(size_t)rhot_band::W490][il][ip] > filltest) {
                    m03_stddev += pow((rfl[(size_t)rhot_band::W490][il][ip]-m03_avg),2);
                    cnt++;
                }
            }
        }
        m03_stddev = sqrt(m03_stddev/(cnt-1));
    } else {
        m03_stddev = floatfill;
    }
// M09 W1380
    cnt = 0;
    for (size_t il=0; il<=2; il++) {
        for (size_t ip=0; ip<=2; ip++) {
            if (rfl[(size_t)rhot_band::W1380][il][ip] > filltest) {
                m09_avg += rfl[(size_t)rhot_band::W1380][il][ip];
                cnt++;
            }
        }
    }
    if (cnt >= 2) {
        m09_avg /= cnt;
        cnt = 0;
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
                if (rfl[(size_t)rhot_band::W1380][il][ip] > filltest) {
                    m09_stddev += pow((rfl[(size_t)rhot_band::W1380][il][ip]-m09_avg),2);
                    cnt++;
                }
            }
        }
        m09_stddev = sqrt(m09_stddev/(cnt-1));
    } else {
        m09_stddev = floatfill;
    }

    if ((m03_stddev > filltest && (m03_stddev > THRHLD470_STD)) ||
            (m09_stddev > filltest && (m09_stddev > THRHLD1380_STD))) {
             mask = 0;
    } else {
        mask = 1;
    }
//   -- perform the check on the 1.38 um band (M09)
    if (rfl[(size_t)rhot_band::W1380][1][1] > THRHLD1380) {
        mask = 0;
    }
//   -- perform check on 488 nm band (M03)
    if (rfl[(size_t)rhot_band::W490][1][1] >  THRHLD470) {
        mask = 0;
    }

    compute_snowmask();

    return status;
}



/**************************************************************************
 * NAME: CldMaskOcean()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

CldMaskOcean::CldMaskOcean() {
}

/**************************************************************************
 * NAME: CldMaskOcean()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

CldMaskOcean::CldMaskOcean(DDAlgorithm* proc) {
    p_ = proc;
}

/**************************************************************************
 * NAME: ~CldMaskOcean()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

CldMaskOcean::~CldMaskOcean() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * cloud mask operations.
 *
 *************************************************************************/

int CldMaskOcean::initialize() {
    return DT_SUCCESS;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Compute cloud mask. This is a port of the Deep Blue cloud
 * mask with modifications for compatibility with Dark Target reflectance
 * units and inversion of the mask bit (here 0 = cloud, 1 = no cloud)
 *
 *************************************************************************/

int CldMaskOcean::compute( unsigned char& mask )
{
    int status = DT_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W410] < filltest) {
        mask = ubytefill;
        return status;
    }
    float lat = p_->lat_;
    float solz = p_->solz_;
    float m01_avg     = 0;
    float m01_stddev  = 0;
    float m08_avg     = 0;
    float m08_stddev  = 0;
    mask = 1;

// Convert to units of normalize radiance
    float rfl[NUM_RFL_BANDS][3][3];
    memset(rfl,0,sizeof(rfl));
    double cossza = cos(solz*DEG2RAD);
    for (size_t iw=0; iw<NUM_RFL_BANDS; iw++) {
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
            	if (p_->rfla_[iw][il][ip] < filltest) {
            		rfl[iw][il][ip] = floatfill;
            	} else {
            		rfl[iw][il][ip] = p_->rfla_[iw][il][ip]*cossza/M_PI;
            	}
            }
        }
    }

// M01 W412
    size_t cnt = 0;
    for (size_t il=0; il<=2; il++) {
        for (size_t ip=0; ip<=2; ip++) {
            if (rfl[(size_t)rhot_band::W410][il][ip] > filltest) {
                m01_avg += rfl[(size_t)rhot_band::W410][il][ip];
                cnt++;
            }
        }
    }
    if (cnt >= 2) {
        m01_avg /= cnt;
        cnt = 0;
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
                if (rfl[(size_t)rhot_band::W410][il][ip] > filltest) {
                    m01_stddev += pow((rfl[(size_t)rhot_band::W410][il][ip]-m01_avg),2);
                    cnt++;
                }
            }
        }
        m01_stddev = sqrt(m01_stddev/(cnt-1));
    } else {
        m01_stddev = floatfill;
    }
// M08 W1240
    cnt = 0;
    for (size_t il=0; il<=2; il++) {
        for (size_t ip=0; ip<=2; ip++) {
            if (rfl[(size_t)rhot_band::W1240][il][ip] > filltest) {
                m08_avg += rfl[(size_t)rhot_band::W1240][il][ip];
                cnt++;
            }
        }
    }
    if (cnt >= 2) {
        m08_avg /= cnt;
        cnt = 0;
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
                if (rfl[(size_t)rhot_band::W1240][il][ip] > filltest) {
                    m08_stddev += pow((rfl[(size_t)rhot_band::W1240][il][ip]-m08_avg),2);
                    cnt++;
                }
            }
        }
        m08_stddev = sqrt(m08_stddev/(cnt-1));
    } else {
        m08_stddev = floatfill;
    }
    float cosza = cos(solz*DEG2RAD);
    if (lat > 65.0) {
        if ((m01_stddev > filltest &&
            (m01_stddev > M01_STDV_THOLD*cosza)) ||
            (m08_stddev > filltest &&
            (m08_stddev > M08_HILAT_STDV_THOLD*cosza))) {
             mask = 0;
        } else {
            mask = 1;
        }
    } else {
        if ((m01_stddev > filltest &&
             m01_stddev > M01_STDV_THOLD*cosza) ||
            (m08_stddev > filltest &&
             m08_stddev > M08_STDV_THOLD*cosza)) {
             mask = 0;
        } else {
            mask = 1;
        }
    }
//   -- perform the check on the 1.38 um band (M09)
    if (rfl[(size_t)rhot_band::W1380][1][1] > M09_THOLD*cosza) {
        mask = 0;
    }
//   -- perform check on 488 nm band (M03)
    if (rfl[(size_t)rhot_band::W490][1][1] >  M03_THOLD*cosza) {
        mask = 0;
    }

    return status;
}


