/*******************************************************************************
 *
 * NAME: DbMask.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute cloud masks for a given granule of data.
 *
 *  Created on: July 15, 2018
 *      Author: Sam Anderson, DB
 *
 *  Modified:
 *
 *******************************************************************************/

#include "deepblue/DbMask.h"

#include <math.h>

#include <DDProcess.h>
#include "deepblue/DbAlgorithm.h"

using namespace std;

/**************************************************************************
 * NAME: DbCloudMaskLand()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbCloudMaskLand::DbCloudMaskLand() {
    p_ = 0;
}

DbCloudMaskLand::DbCloudMaskLand(DbAlgLand* proc)
{
    p_ = proc;
}

DbCloudMaskLand::~DbCloudMaskLand() {
}

/**************************************************************************
 * NAME: DbCloudMaskLand::compute_1()
 *
 * DESCRIPTION: Initial set of cloud mask filters applied to modis units
 *
 *************************************************************************/

int DbCloudMaskLand::compute_1( const size_t iy, const size_t ix, short& mask,
        short& snow1, short& snow2 )
{
    int status = DTDB_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W410] < 0) {
        mask = DFILL_SHORT;
        snow1 = DFILL_SHORT;
        snow2 = DFILL_SHORT;
        return status;
    }
    mask = 0;
    snow1 = 0;
    snow2 = 0;

    float lat = p_->lat_;
    float lon = p_->lon_;

    // Convert raw reflectances to modis units
    float rfl[NTWL];
    if (p_->rfl_[(size_t)rhot_band::W410] > 0) {
        for (size_t ib=0; ib<NTWL; ib++) {
            rfl[ib] = p_->rfl_[ib];
        }
    } else {
        for (size_t ib=0; ib<NTWL; ib++) {
            rfl[ib] = DFILL_FLOAT;
        }
    }

    int ilat = (int)((lat + 90.)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = (int)((lon + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    int gzflg = p_->gz_lut_->GEOZONE_FLAG[ilat][ilon];
    // Define pixel range
    size_t ix1 = (ix==0) ? 0 : ix-1;
    size_t ix2 = min(ix+1,p_->pixels_-1);
    size_t iy1 = (iy==0) ? 0 : iy-1;
    size_t iy2 = min(iy+1,p_->lines_-1);
    // Compute Rfl minmax ratios
    float rflmin = 999.0;
    float rflmax = -999.0;
    for (size_t cx=ix1; cx<=ix2; cx++) {
		for (size_t cy=iy1; cy<=iy2; cy++) {
			 float rf = p_->rfla_[(size_t)rhot_band::W410][cy-iy1][cx-ix1];
             if (rf > 0) {
                 rflmin = (rf<rflmin) ? rf : rflmin;
                 rflmax = (rf>rflmax) ? rf : rflmax;
             }
        }
    }
    if (rflmax<1.E-8 || rflmin<1.E-8 ) {
        mask = DFILL_SHORT;
        snow1 = DFILL_SHORT;
        snow2 = DFILL_SHORT;
        return DTDB_SUCCESS;
    }
    float minmax_rfl = rflmax/rflmin + DbAlgLand::delta;
    // Begin cloud tests
    if (minmax_rfl > 1.2) {
        mask = 1;
    }
    if (p_->sr670_ >= 0.08 && minmax_rfl > 1.15) {
        mask = 1;
    }
// Over bright surfaces, apply RR1.38/0.66 and BTD11-12 threshold !JH
// test addition - 6/7/2011
     if (p_->sr670_ >= 0.08 ) {
         if (rfl[(size_t)rhot_band::W1380]  > 0.018 && p_->pwv_ >= 0.9) {
             mask = 1;  // Thin Cirrus Screening
         }
         if (gzflg == 24) {// Taklimakan
             if (rfl[(size_t)rhot_band::W1380] > 0.04 && p_->pwv_ < 0.9 &&
                     rfl[(size_t)rhot_band::W670] < 0.45) {
                 mask = 1;
             }
         } else if (rfl[(size_t)rhot_band::W1380] > 0.018 && p_->pwv_ < 0.9 &&
                     rfl[(size_t)rhot_band::W670] < 0.55) {
                 mask = 1;
         }
     }
// Over dark surfaces
     if (p_->sr670_ > 0.0 && p_->sr670_ < 0.08 ) {
         float dd = rfl[(size_t)rhot_band::W670]/rfl[(size_t)rhot_band::W410];
         if (p_->rfl_[(size_t)rhot_band::W2250] > 0.36 && dd < 0.95) {
             mask = 1;
         }
         if (rfl[(size_t)rhot_band::W1380] > 0.018 &&
             p_->pwv_ > 0.4) {
             mask = 1;
         }
     }
// High AMF
     if (lat > 60.0 && p_->amf_ > 7.0) {
         if (rfl[(size_t)rhot_band::W1380] > 0.018 &&
             p_->pwv_ > 0.9) {
             mask = 1;
         }
         if (rfl[(size_t)rhot_band::W490] > 0.4) {
             mask = 1;
         }
     }
     if ((gzflg == 5 || gzflg == 1 || gzflg == 26 || gzflg == 27) ||
     (lon < 15.0 &&  ((lon >-20.0 && lon < 55.0) && gzflg == -999))) {
         if (p_->sr670_ > 0.0 && p_->sr670_ < 0.08 ) {
             if (rfl[(size_t)rhot_band::W490] > 0.15) {
                 if (rfl[(size_t)rhot_band::W2250] > 0.16) {
                     mask = 1;
                 }
             }
         }
     }

// -- try to detect snow cover and skip pixel.
// -- source: http://modis-snow-ice.gsfc.nasa.gov/?c=atbdt=atbd
     float ndsi = DFILL_FLOAT;
     if (rfl[(size_t)rhot_band::W550] > 0 && rfl[(size_t)rhot_band::W2250] > 0) {
         ndsi = (rfl[(size_t)rhot_band::W550] - rfl[(size_t)rhot_band::W1610]) /
                 (rfl[(size_t)rhot_band::W550] + rfl[(size_t)rhot_band::W1610]);
     } else {
         mask = 1;
     }
     if (ndsi > 0.35 && rfl[(size_t)rhot_band::W865] > 0.11 &&
         rfl[(size_t)rhot_band::W550] > 0.1) {
         snow1 = 1;
     }
     if (ndsi > -0.2 && rfl[(size_t)rhot_band::W865] > 0.11) {
         snow2 = 1;
     }
     if (snow1 == 1) {
         mask = 1;
     }

     return status;
}

/**************************************************************************
 * NAME: DbCloudMaskLand::compute_2()
 *
 * DESCRIPTION: Initial set of cloud mask filters applied to modis units
 *
 *************************************************************************/


int DbCloudMaskLand::compute_2( const size_t iy, const size_t ix,
                                short& mask, const short snow2 )
{
    int status = DTDB_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W410] < 0) {
        mask = DFILL_SHORT;
        return status;
    }
    float lat = p_->lat_;
    float lon = p_->lon_;
    float solz = p_->solz_;
    int ilat = (int)((lat + 90.)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = (int)((lon + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    int gzflg = p_->gz_lut_->GEOZONE_FLAG[ilat][ilon];
    int sfcstd = p_->sp_lut_->SURFACE_ELEVATION[ilat][ilon];
    int month = p_->month_;

    if (p_->ler412_ > 50.0) {
        mask = 1;
    }
    if (p_->sr670_ < 0.1) {
        if (p_->rfl_[(size_t)rhot_band::W2250] > 0.05 &&
            p_->ler412_ > 12.0) {
            mask = 1;
        }
    }
    if (p_->sr670_ < 0.08) {
        if (p_->ler412_ > 40.0) {
                mask = 1;
        }
        if (p_->ler412_ > 20.0 && solz > 72.0) {
            mask = 1;
        }
    }
// Define pixel range
    size_t ix1 = (ix==0) ? 0 : ix-1;
    size_t ix2 = min(ix+1,p_->pixels_-1);
    size_t iy1 = (iy==0) ? 0 : iy-1;
    size_t iy2 = min(iy+1,p_->lines_-1);
// Compute Ler minmax ratios
    float lermin = 999.0;
    float lermax = -999.0;
    for (size_t cx=ix1; cx<=ix2; cx++) {
        for (size_t cy=iy1; cy<=iy2; cy++) {
//            float lr = p_->ler412_a_[cy][cx];
//            if (lr > 0) {
//                lermin = (lr<lermin) ? lr : lermin;
//                lermax = (lr>lermax) ? lr : lermax;
//            }
            lermin = 1;
            lermax = 1;
        }
    }

    if (lermax<1.E-8 || lermin<1.E-8 ) {
        mask = DFILL_FLOAT;
        return DTDB_SUCCESS;
    }
    float minmax_ler = lermax/lermin + DbAlgLand::delta;
//-- N. America
    if (snow2 == 1 && gzflg == 13 &&
        p_->ndvi_ < 0.35 &&
        sfcstd <= 50 && minmax_ler > 1.20) {
        mask = 1;
    }
//-- higher ndvi and minmax thresholds for spring at high latitudes
    if (lat > 45 && month >= 2 && month <= 6 &&
        snow2 == 1 && gzflg == 13 && p_->ndvi_ < 0.45 &&
        sfcstd <= 50 && minmax_ler > 1.30) {
        mask = 1;
    }
//-- Rest of the world
    if (snow2 == 1 && gzflg != 13 && p_->ndvi_ < 0.35 &&
        sfcstd <= 50 && minmax_ler > 1.40) {
        mask = 1;
    }
// --higher ndvi threshold and lower minmax for spring at high latitudes
    if (lat > 50 && month >= 2 && month <= 6 &&  snow2 == 1 &&
        gzflg != 13 && p_->ndvi_ < 0.45 &&
        sfcstd <= 50 && minmax_ler > 1.30) {
        mask = 1;
    }
//-- Mountainous regions
    if (snow2 == 1 && sfcstd > 50 && minmax_ler > 1.80) {
        mask = 1;
    }
// -- detect anomalous D* values in large dust storm over Taklimakan
// -- region and resetthe cloud mask to 0.
    if (gzflg == 24) {
        if (p_->ler412_ > 0 &&
            p_->ler488_ > 0 &&
            p_->ler670_ > 0) {
            if (p_->ler488_/p_->ler670_ < 0.65 &&
                p_->ler412_ > 18.0) {
                mask = 0;
            }
        }
    }

    return status;
}

/**************************************************************************
 * NAME: DbSmokeMask()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbSmokeMask::DbSmokeMask() {
    p_ = 0;
}

DbSmokeMask::DbSmokeMask(DbAlgLand* proc) {
    p_ = proc;
}

DbSmokeMask::~DbSmokeMask() {
}

/**************************************************************************
 * NAME: DbSmokeMask::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/

int DbSmokeMask::compute( const size_t iy, const size_t ix, short& mask  )
{
    int status = DTDB_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W1380] < 0) {
        mask = DFILL_SHORT;
        return status;
    }
    mask = 0;

    size_t ix1 = (ix==0) ? 0 : ix-1;
    size_t ix2 = min(ix+1,p_->pixels_-1);
    size_t iy1 = (iy==0) ? 0 : iy-1;
    size_t iy2 = min(iy+1,p_->lines_-1);
    int cnt2250_1 = 0;
    int cnt2250_2 = 0;

    float tler[3][3][3];
    float tsr[3][3][3];
    float td[3][3][3];
    memset(&tler[0][0][0],0,sizeof(tler));
    memset(&tsr[0][0][0],0,sizeof(tsr));
    memset(&td[0][0][0],0,sizeof(td));
    for (size_t cy=iy1; cy<=iy2; cy++) {
        for (size_t cx=ix1; cx<=ix2; cx++) {
//            tler[0][cy-iy1][cx-ix1] = p_->ler412_a_[cy][cx];
//            tler[1][cy-iy1][cx-ix1] = p_->ler488_a_[cy][cx];
//            tler[2][cy-iy1][cx-ix1] = p_->ler670_a_[cy][cx];
//            tsr[0][cy-iy1][cx-ix1] = p_->sr412_a_[cy][cx];
//            tsr[1][cy-iy1][cx-ix1] = p_->sr488_a_[cy][cx];
//            tsr[2][cy-iy1][cx-ix1] = p_->sr670_a_[cy][cx];
//            td[0][cy-iy1][cx-ix1] = p_->rfla_[(size_t)rhot_band::W1380][cy][cx];
//            td[1][cy-iy1][cx-ix1] = p_->rfla_[(size_t)rhot_band::W2250][cy][cx];
//            td[2][cy-iy1][cx-ix1] = p_->rfla_[(size_t)rhot_band::W11000][cy][cx];

            if (p_->rfla_[(size_t)rhot_band::W2250][cy-iy1][cx-ix1] < 0.035) cnt2250_1++;
            if (p_->rfla_[(size_t)rhot_band::W2250][cy-iy1][cx-ix1] < 0.065) cnt2250_2++;
        }
    }
    float rat488670 = p_->ler488_ / p_->ler670_;
    float rat412488 = p_->ler412_ / p_->ler488_;
    if (p_->ler412_ > 12.0) {
      if ((p_->sr670_ < 0.08 && cnt2250_1 == 9) ||
          (p_->sr670_ >= 0.08 && cnt2250_2 == 9)) {
        if (p_->rfl_[(size_t)rhot_band::W2250] < 0.025) {
          mask = 1;
        } else {
          if (p_->ler412_ > 20.0 &&
              p_->rfl_[(size_t)rhot_band::W2250] < 0.04) {
              mask = 1;
          } else {
            if (p_->ler488_ > 0.0 && p_->ler670_ > 0.0) {
              rat488670 = p_->ler488_ / p_->ler670_;
              if (rat488670 > 0.88) {
                  mask = 1;
              } else {
                if (rat488670 > 0.7 && rat488670 < 0.88 &&
                    p_->rfl_[(size_t)rhot_band::W2250] < 0.04 &&
                    p_->rfl_[(size_t)rhot_band::W1380] > 0.0015) {
                    mask = 1;
                }
              }
            }
          }
        }
      }
    }
    if (p_->sr670_ >= 0.12 || p_->ler670_ > 50.0) {
        mask = 0;
    }
    if (rat412488/rat488670 > 0.94) {
        mask = 0;
    }

    return status;
}

/**************************************************************************
 * NAME: DbPyrocbMask()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbPyrocbMask::DbPyrocbMask() {
    p_ = 0;
}

DbPyrocbMask::DbPyrocbMask(DbAlgLand* proc) {
    p_ = proc;
}

DbPyrocbMask::~DbPyrocbMask() {
}

/**************************************************************************
 * NAME: DbPyrocbMask::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/

int DbPyrocbMask::compute( const size_t iy, const size_t ix, short& mask )
{
    int status = DTDB_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W1380] < 0) {
        mask = DFILL_SHORT;
        return status;
    }
    mask= 0;
    float rat412488 = p_->ler412_ / p_->ler488_;
    if (p_->rfl_[(size_t)rhot_band::W1380] > 0.06  &&
        p_->rfl_[(size_t)rhot_band::W2250] > 0.2   &&
        p_->ler670_   > 24.0  &&
        rat412488 < 0.7) {
        mask = 1;
    }

    return status;
}

/**************************************************************************
 * NAME: DbHighAltSmokeMask()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbHighAltSmokeMask::DbHighAltSmokeMask() {
    p_ = 0;
}

DbHighAltSmokeMask::DbHighAltSmokeMask(DbAlgLand* proc) {
    p_ = proc;
}

DbHighAltSmokeMask::~DbHighAltSmokeMask() {
}

/**************************************************************************
 * NAME: DbHighAltSmokeMask::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/

int DbHighAltSmokeMask::compute( const size_t iy, const size_t ix, short& mask )
{
    int status = DTDB_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W1380] < 0) {
        mask = DFILL_SHORT;
        return status;
    }
    mask= 0;
    float lat = p_->lat_;
    float lon = p_->lon_;

    int ilat = (int)((lat + 90.)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = (int)((lon + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    int gzflg = p_->gz_lut_->GEOZONE_FLAG[ilat][ilon];

    size_t ix1 = (ix==0) ? 0 : ix-1;
    size_t ix2 = min(ix+1,p_->pixels_-1);
    size_t iy1 = (iy==0) ? 0 : iy-1;
    size_t iy2 = min(iy+1,p_->lines_-1);
    int cntW2250 = 0;
    float tler[3][3][3];
    float tsr[3][3][3];
    float td[3][3][3];
    memset(&tler[0][0][0],0,sizeof(tler));
    memset(&tsr[0][0][0],0,sizeof(tsr));
    memset(&td[0][0][0],0,sizeof(td));
    for (size_t cy=iy1; cy<=iy2; cy++) {
        for (size_t cx=ix1; cx<=ix2; cx++) {
//            tler[0][cy-iy1][cx-ix1] = p_->ler412_a_[cy][cx];
//            tler[1][cy-iy1][cx-ix1] = p_->ler488_a_[cy][cx];
//            tler[2][cy-iy1][cx-ix1] = p_->ler670_a_[cy][cx];
//            tsr[0][cy-iy1][cx-ix1] = p_->sr412_a_[cy][cx];
//            tsr[1][cy-iy1][cx-ix1] = p_->sr488_a_[cy][cx];
//            tsr[2][cy-iy1][cx-ix1] = p_->sr670_a_[cy][cx];
//            td[0][cy-iy1][cx-ix1] = p_->rfl_[(size_t)rhot_band::W1380][cy][cx];
//            td[1][cy-iy1][cx-ix1] = p_->rfl_[(size_t)rhot_band::W2250][cy][cx];
//            td[2][cy-iy1][cx-ix1] = p_->rfl_[(size_t)rhot_band::W11000][cy][cx];
            if (p_->rfla_[(size_t)rhot_band::W2250][cy-iy1][cx-ix1] < 0.25) cntW2250++;
        }
    }
    float rat488670 = p_->ler488_ / p_->ler670_;
    float rat412488 = p_->ler412_ / p_->ler488_;

    if (p_->ler412_ > 12.0) {
        if ((p_->sr670_ < 0.08 && cntW2250 == 9) ||
            (p_->sr670_ >= 0.08 && cntW2250 == 9)) {
            if (p_->ler488_ > 18.0 &&
                (rat488670 > 0.7 && rat488670 < 0.88) &&
                p_->rfl_[(size_t)rhot_band::W2250] < 0.22 &&
                p_->rfl_[(size_t)rhot_band::W1380] > 0.0032) {
                mask = 1;
            }
            if (p_->ler488_ > 18.0 && (rat488670 > 0.7 && rat488670 < 0.88) &&
                    p_->rfl_[(size_t)rhot_band::W2250] < 0.22 && rat412488 < 0.8) {
                mask = 1;
            }
            if (p_->ler488_ > 18.0 && p_->rfl_[(size_t)rhot_band::W1380] > 0.0015 &&
            		rat412488 < 0.8) {
                mask = 1;
            }
            if (gzflg == 31) {
                if (p_->ler488_ > 18.0 &&
                    (rat488670 > 0.9 && rat488670 < 1.06) &&
                    (rat412488 > 0.7 && rat412488 < 0.88) &&
                    p_->rfl_[(size_t)rhot_band::W2250] < 0.22 &&
                    p_->rfl_[(size_t)rhot_band::W1380] > 0.0032) {
                    mask = 1;
                }
            }
        }
    }
    if (p_->sr670_ >= 0.12 || p_->rfl_[(size_t)rhot_band::W2250] > 0.06 ||
        p_->ler670_ > 50.0 || p_->rfl_[(size_t)rhot_band::W1380] < 0.0005 ||
        p_->rfl_[(size_t)rhot_band::W1380] > 0.013) {
        mask = 0;
    }
    if (rat412488/rat488670 > 0.94) {
        mask = 0;
    }
    return status;
}

/**************************************************************************
 * NAME: DbCloudMaskOcean()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbCloudMaskOcean::DbCloudMaskOcean() {
}

DbCloudMaskOcean::DbCloudMaskOcean(DbAlgOcean* proc) {
    p_ = proc;
}

DbCloudMaskOcean::~DbCloudMaskOcean() {
}

/**************************************************************************
 * NAME: DbCloudMaskOcean::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/


int DbCloudMaskOcean::compute( short& mask )
{
    int status = DTDB_SUCCESS;

    if (p_->rfl_[(size_t)rhot_band::W410] < 0) {
        mask = DFILL_SHORT;
        return status;
    }
    float lat = p_->lat_;
    float solz = p_->solz_;
    float m01_avg     = 0;
    float m01_stddev  = 0;
    float m08_avg     = 0;
    float m08_stddev  = 0;
    mask = 0;

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
             mask = 1;
        } else {
            mask = 0;
        }
    } else {
        if ((m01_stddev > 0 &&
             m01_stddev > M01_STDV_THOLD*cosza) ||
            (m08_stddev > 0 &&
             m08_stddev > M08_STDV_THOLD*cosza)) {
             mask = 1;
        } else {
            mask = 0;
        }
    }
//   -- perform the check on the 1.38 um band (M09)
    if (rfl[(size_t)rhot_band::W1380][1][1] > M09_THOLD*cosza) {
        mask = 1;
    }
//   -- perform check on 488 nm band (M03)
    if (rfl[(size_t)rhot_band::W490][1][1] >  M03_THOLD*cosza) {
        mask = 1;
    }

    return status;
}

