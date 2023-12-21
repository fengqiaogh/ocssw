/*******************************************************************************
 *
 * NAME: DtMask.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute cloud masks for a given DDProcess object class.
 *
 *  Created on: November 3, 2016
 *      Author: Sam Anderson, DT
 *
 *  Modified:
 *
 *******************************************************************************/

#include "darktarget/DtMask.h"

#include <math.h>

#include <DDProcess.h>
#include "darktarget/DtAlgorithm.h"

using namespace std;

/**************************************************************************
 * NAME: DtCloudMaskLand()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtCloudMaskLand::DtCloudMaskLand() {
}

/**************************************************************************
 * NAME: DtCloudMaskLand()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtCloudMaskLand::DtCloudMaskLand(DtAlgorithm* process) {
	p_ = process;
    mask_cirrus_ = 0;
    snowmask_ = 0;
}

/**************************************************************************
 * NAME: ~DtCloudMaskLand()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DtCloudMaskLand::~DtCloudMaskLand() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * cloud mask operations.
 *
 *************************************************************************/

int DtCloudMaskLand::initialize() {
    return DTDB_SUCCESS;
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

int DtCloudMaskLand::compute( unsigned char& mask )
{
    int status = DTDB_SUCCESS;

    mask = 1;

    if (p_->rfl_[(size_t)rhot_band::W410] < 0) {
        mask = DFILL_UBYTE;
        return status;
    }
    float m03_avg     = 0;
    float m03_stddev  = 0;
    float m09_avg     = 0;
    float m09_stddev  = 0;

// Gather 3x3 grid
    float rfl[NTWL][3][3];
    memset(rfl,0,sizeof(rfl));
    for (size_t iw=0; iw<NTWL; iw++) {
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
            if (rfl[(size_t)rhot_band::W490][il][ip] > 0) {
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
                if (rfl[(size_t)rhot_band::W490][il][ip] > 0) {
                    m03_stddev += pow((rfl[(size_t)rhot_band::W490][il][ip]-m03_avg),2);
                    cnt++;
                }
            }
        }
        m03_stddev = sqrt(m03_stddev/(cnt-1));
    } else {
        m03_stddev = DFILL_FLOAT;
    }
// M09 W1380
    cnt = 0;
    for (size_t il=0; il<=2; il++) {
        for (size_t ip=0; ip<=2; ip++) {
            if (rfl[(size_t)rhot_band::W1380][il][ip] > 0) {
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
                if (rfl[(size_t)rhot_band::W1380][il][ip] > 0) {
                    m09_stddev += pow((rfl[(size_t)rhot_band::W1380][il][ip]-m09_avg),2);
                    cnt++;
                }
            }
        }
        m09_stddev = sqrt(m09_stddev/(cnt-1));
    } else {
        m09_stddev = DFILL_FLOAT;
    }

    if ((m03_stddev > 0 && (m03_stddev > THRHLD470_STD)) ||
            (m09_stddev > 0 && (m09_stddev > THRHLD1380_STD))) {
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
 * NAME: compute_snowmask()
 *
 * DESCRIPTION: Compute transmission corrections and snow mask
 *
 *************************************************************************/

int DtCloudMaskLand::compute_snowmask()
{
    int status = DTDB_SUCCESS;

// For Snow masking
   snowmask_ = 1;
/*
// Compute Temperature(convert from radiance to temperature 11.um channel)
// Derive constants
    float c1 = 2.0*PlancksConstant*(SpeedOfLight*SpeedOfLight);
    float c2 = (PlancksConstant*SpeedOfLight)/BoltzConstant;
//  convert wavelength to meters
    float w_meter = (1.0e-6*WAV2);
//  Compute Rong_rond ratio
    float ratio = 0;
    float r865 = p_->rfl_[(size_t)rhot_band::W865];
    float r1240 = p_->rfl_[(size_t)rhot_band::W1240];
    float r11000 = p_->rfl_[(size_t)rhot_band::W11000];

    if((r865 > 0) && (r1240 > 0)) {
        ratio = (r865 - r1240) / (r865 + r1240);
    }
    else {
        ratio = 0.0;
    }
    float temp = 285;
    if (p_->instrument_ != sensor::POCI) {
        temp = c2 / (w_meter*log(c1/(1.0e+6 * r11000 * pow(w_meter,5)+1.0)));
    }
//  Set snow mask as snowy pixel based on Ratio and temp.
    if ((ratio > 0.01) && (temp < 285)) {
        snowmask_ = 0;
    }
*/
    return status;
}

/**************************************************************************
 * NAME: DtCloudMaskOcean()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtCloudMaskOcean::DtCloudMaskOcean() {
}

/**************************************************************************
 * NAME: DtCloudMaskOcean()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtCloudMaskOcean::DtCloudMaskOcean(DtAlgorithm* proc) {
    p_ = proc;
}

/**************************************************************************
 * NAME: ~DtCloudMaskOcean()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DtCloudMaskOcean::~DtCloudMaskOcean() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * cloud mask operations.
 *
 *************************************************************************/

int DtCloudMaskOcean::initialize() {
    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Compute cloud mask. This is a port of the Deep Blue cloud
 * mask with modifications for compatibility with Dark Target reflectance
 * units and inversion of the mask bit (here 0 = cloud, 1 = no cloud)
 *
 *************************************************************************/

int DtCloudMaskOcean::compute( unsigned char& mask )
{
    int status = DTDB_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W410] < 0) {
        mask = DFILL_UBYTE;
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
    float rfl[NTWL][3][3];
    memset(rfl,0,sizeof(rfl));
    double cossza = cos(solz*DEGtoRAD);
    for (size_t iw=0; iw<NTWL; iw++) {
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
            	if (p_->rfla_[iw][il][ip] < 0) {
            		rfl[iw][il][ip] = DFILL_FLOAT;
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
            if (rfl[(size_t)rhot_band::W410][il][ip] > 0) {
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
                if (rfl[(size_t)rhot_band::W410][il][ip] > 0) {
                    m01_stddev += pow((rfl[(size_t)rhot_band::W410][il][ip]-m01_avg),2);
                    cnt++;
                }
            }
        }
        m01_stddev = sqrt(m01_stddev/(cnt-1));
    } else {
        m01_stddev = DFILL_FLOAT;
    }
// M08 W1240
    cnt = 0;
    for (size_t il=0; il<=2; il++) {
        for (size_t ip=0; ip<=2; ip++) {
            if (rfl[(size_t)rhot_band::W1240][il][ip] > 0) {
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
                if (rfl[(size_t)rhot_band::W1240][il][ip] > 0) {
                    m08_stddev += pow((rfl[(size_t)rhot_band::W1240][il][ip]-m08_avg),2);
                    cnt++;
                }
            }
        }
        m08_stddev = sqrt(m08_stddev/(cnt-1));
    } else {
        m08_stddev = DFILL_FLOAT;
    }
    float cosza = cos(solz*DEGtoRAD);
    if (lat > 65.0) {
        if ((m01_stddev > 0 &&
            (m01_stddev > M01_STDV_THOLD*cosza)) ||
            (m08_stddev > 0 &&
            (m08_stddev > M08_HILAT_STDV_THOLD*cosza))) {
             mask = 0;
        } else {
            mask = 1;
        }
    } else {
        if ((m01_stddev > 0 &&
             m01_stddev > M01_STDV_THOLD*cosza) ||
            (m08_stddev > 0 &&
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

/**************************************************************************
 * NAME: DtSedimentMask()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtSedimentMask::DtSedimentMask() {
}

/**************************************************************************
 * NAME: DtSedimentMask()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtSedimentMask::DtSedimentMask(DtAlgorithm* proc) {
    p_ = proc;
}

/**************************************************************************
 * NAME: ~DtSedimentMask()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DtSedimentMask::~DtSedimentMask() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * cloud mask operations.
 *
 *************************************************************************/

int DtSedimentMask::initialize() {
    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Make sediment array
 *
 *************************************************************************/

int DtSedimentMask::compute(short& mask) {
    int status = DTDB_SUCCESS;

    float x[NWAV];
    float y[NWAV];
    float sig[NWAV] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    bool bUseWav[NWAV];
    bUseWav[D488] = (p_->rfl_[(size_t)rhot_band::W490] > 0) ? true : false;
    bUseWav[D550] = false;
    bUseWav[D670] = false;
    bUseWav[D865] = false;
//  bUseWav[D1240] = (p_->rfl_[W124] > 0) ? true : false;
    bUseWav[D1240] = false;  // conform to bug in fortran code
    bUseWav[D1610] = (p_->rfl_[(size_t)rhot_band::W1610] > 0) ? true : false;
    bUseWav[D2250] = (p_->rfl_[(size_t)rhot_band::W2250] > 0) ? true : false;
    float logRefl[NWAV];
    logRefl[D488] = log(p_->rfl_[(size_t)rhot_band::W490]);
    logRefl[D550] = 0;
    logRefl[D670] = 0;
    logRefl[D865] = 0;
    logRefl[D1240] = log(p_->rfl_[(size_t)rhot_band::W1240]);
    logRefl[D1610] = log(p_->rfl_[(size_t)rhot_band::W1610]);
    logRefl[D2250] = log(p_->rfl_[(size_t)rhot_band::W2250]);

    DtAlgOcean* pdto = static_cast<DtAlgOcean*> (p_);
    int numWaves = 0;
    for (int iWav = 0; iWav < NWAV; iWav++) {
        if (bUseWav[iWav]) {
            x[numWaves] = log(pdto->lut_.WAVE[iWav]);
            y[numWaves] = logRefl[iWav];
            numWaves++;
        }
    }
    float acoef = 0.0;
    float bcoef = 0.0;
    float del_wav55 = 0.0;
    if (numWaves > 2) {
        pdto->fit_line(x, y, sig, numWaves, acoef, bcoef);
        del_wav55 = p_->rfl_[(size_t)rhot_band::W550]
                - exp(acoef + bcoef * log(pdto->lut_.WAVE[D550]));
    } else
        del_wav55 = 0.0;

    mask = 1;
    float r488 = p_->rfl_[(size_t)rhot_band::W490];
    float r2250 = p_->rfl_[(size_t)rhot_band::W2250];
//  mask1 true it is sediment+smoke+dust set mask to sediments
    if ((r2250 < 0.10) && (del_wav55 >= 0.015)) {
        mask = 0;
//  mask2 true it is  smoke+dust and no sediments set previous sediment mask
//  to no sediment because it smoke or dust.
        if ((r488 >= 0.25) && (del_wav55 >= 0.015)) {
            mask = 1;
        }
    }

    return status;
}



