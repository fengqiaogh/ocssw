/*******************************************************************************
 *
 * NAME: DbMaskOCI.cpp
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

#include "deepblue/DbMaskOCI.h"

#include <math.h>

#include <Granule.h>
#include "deepblue/DbProcess.h"

using namespace std;

/**************************************************************************
 * NAME: DbCloudMaskLandOCI()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbCloudMaskLandOCI::DbCloudMaskLandOCI() {
    p_ = 0;
}

DbCloudMaskLandOCI::DbCloudMaskLandOCI(Granule* granule, DbProcessLand* proc) :
             Mask(granule)
{
    p_ = proc;
}

DbCloudMaskLandOCI::~DbCloudMaskLandOCI() {
}

/**************************************************************************
 * NAME: DbCloudMaskLandOCI::compute_1()
 *
 * DESCRIPTION: Initial set of cloud mask filters applied to modis units
 *
 *************************************************************************/

int DbCloudMaskLandOCI::compute_1( const int iy, const int ix, short& mask,
        short& snow1, short& snow2 )
{
    int status = DT_SUCCESS;

    if (p_->rfl_[W412] < filltest) {
        mask = shortfill;
        snow1 = shortfill;
        snow2 = shortfill;
        return status;
     }
     mask = 0;
     snow1 = 0;
     snow2 = 0;

     float lat = g_->in_->latitude;
     float lon = g_->in_->longitude;

 // Convert raw reflectances to modis units
     float to_modis  = PI;
     float rfl[NUM_RFL_BANDS];
     if (p_->rfl_[W412] > filltest) {
         for (int ib=0; ib<NUM_RFL_BANDS; ib++) {
             rfl[ib] = p_->rfl_[ib]*to_modis;
         }
     } else {
         for (int ib=0; ib<NUM_RFL_BANDS; ib++) {
             rfl[ib] = floatfill;
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
     int ix1 = max(ix-1,0);
     int ix2 = min(ix+1,g_->pixels_-1);
     int iy1 = max(iy-1,0);
     int iy2 = min(iy+1,g_->lines_-1);
// Compute Rfl minmax ratios
     float rflmin = 999.0;
     float rflmax = -999.0;
     for (int cx=ix1; cx<=ix2; cx++) {
         for (int cy=iy1; cy<=iy2; cy++) {
             float rf = p_->rfl_[W412][cy][cx];
             if (rf > filltest) {
                 rflmin = (rf<rflmin) ? rf : rflmin;
                 rflmax = (rf>rflmax) ? rf : rflmax;
             }
         }
     }
     if (rflmax<1.E-8 || rflmin<1.E-8 ) {
         mask = shortfill;
         snow1 = shortfill;
         snow2 = shortfill;
         return DT_SUCCESS;
     }
     float minmax_rfl = rflmax/rflmin + DbProcessLand::delta;
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
//         if ((p_->dstar_ < 1.12 && solz < 68.0) ||
//             (p_->dstar_ < 1.2 && solz  >  68.0)) {
             if (rfl[W1380]  > 0.018 &&
                p_->pwv_ >= 0.9) {
                mask = 1;  // Thin Cirrus Screening
             }
             if (gzflg == 24) {// Taklimakan
                 if (rfl[W1380] > 0.04 &&
                     p_->pwv_ < 0.9 &&
//                     p_->btd11_ > -2.5 &&
                     rfl[W670] < 0.45) {
                     mask = 1;
                 }
//                 if (p_->rfl_[W11000] < 265.0) {
//                     mask = 1;   // Cloud Screening
//                 }
             } else {   // Everywhere else
                 if (rfl[W1380] > 0.018 &&
                     p_->pwv_ < 0.9 &&
//                     p_->btd11_ > -1.0 &&
                     rfl[W670] < 0.55) {
                     mask = 1;
                 }
//                 if (p_->rfl_[W11000] < 270.0) {
//                     mask = 1;   // Cloud Screening
//                 }
             }
//             if (p_->sr670_ >= 0.16) {
//                 if (p_->rfl_[W11000] >= 270.0 &&
//                     p_->rfl_[W11000] < 281.0 &&
//                     p_->btd11_ > -0.5) {
//                     mask = 1;   // Cloud Edge Contamination Removal
//                 }
//             }
//         }
     }
// Over dark surfaces
     if (p_->sr670_ > 0.0 && p_->sr670_ < 0.08 ) {
//         if ((p_->dstar_ < 1.12 && solz < 68.0) ||
//            (p_->dstar_ < 1.2 && solz > 68.0)) {
             float dd = rfl[W670]/rfl[W412];
             if (p_->rfl_[W2250] > 0.36 && dd < 0.95) {
                 mask = 1;
             }
             if (rfl[W1380] > 0.018 &&
                 p_->pwv_ > 0.4) {
                 mask = 1;
             }
//             if (p_->rfl_[W11000] < 270.0) {
//                 mask = 1;
//             }
//         }
     }
// High AMF
     if (lat > 60.0 && p_->amf_ > 7.0) {
         if (rfl[W1380] > 0.018 &&
             p_->pwv_ > 0.9) {
             mask = 1;
         }
         if (rfl[W488] > 0.4) {
             mask = 1;
         }
     }
// -- over *really* bright surfaces, except Arabian Peninsula (6 <= gzflg <= 11)
//     if (gzflg < 6 || gzflg > 11) {
//         if (p_->sr670_ > 0.25) {
//             if (p_->dstar_ < 0.85 && p_->pwv_ < 2.3) {
//               mask = 1;
//             }
//         }
//     }
     if ((gzflg == 5 || gzflg == 1 || gzflg == 26 || gzflg == 27) ||
     (lon < 15.0 &&  ((lon >-20.0 && lon < 55.0) && gzflg == -999))) {
//         if (p_->sr670_ < 0.08 && p_->dstar_ < 1.12) {
         if (p_->sr670_ > 0.0 && p_->sr670_ < 0.08) {
             if (rfl[W488] > 0.15) {
                 if (rfl[W2250] > 0.16) {
                     mask = 1;
                 }
             }
         }
//         if (p_->sr670_ > 0.08 && p_->dstar_ < 1.12) {
//         if (p_->sr670_ > 0.08) {
//             if (rfl[W488] > 0.25) {
//                 if (p_->btd11_ > 3.0) {
//                     mask = 1;
//                 }
//             }
//         }
     }

// -- only apply these cloud tests over region 5 (N. Africa) and southern Africa.
// -- developed using MYD021KM.A2011171.A12[25|30]*.hdf as test case.
// -- over *really* bright surfaces, except Arabian Peninsula (6 <= gzflg <= 11)
//     if (gzflg < 6 || gzflg > 11) {
//         if (p_->sr670_ > 0.25) {
//             if (p_->dstar_ < 0.85 && p_->pwv_ < 2.3) {
//             if (p_->pwv_ < 2.3) {
//                 mask = 1;
//             }
//         }
//     }
//     if ((gzflg == 5  || gzflg ==  1 || gzflg == 26 || gzflg == 27) ||
//          (lat < 15.0 && (lon >-20.0 && lon < 55.0) && gzflg == shortfill)) {
//
//         if (p_->sr670_ < 0.08 && p_->dstar_ < 1.12) {
//         if (p_->sr670_ < 0.08) {
//             if (rfl[W488] > 0.15) {
//                 if (rfl[W2250] > 0.16) {
//                     mask = 1;
//                 }
//             }
//         }
//         if (p_->sr670_ > 0.08 && p_->dstar_ < 1.12) {
//         if (p_->sr670_ > 0.08) {
//             if (rfl[W488] > 0.25) {
//                 if (p_->btd11_ > 3.0) {
//                     mask = 1;
//                 }
//             }
//         }
//     }
// -- try to detect snow cover and skip pixel.
// -- source: http://modis-snow-ice.gsfc.nasa.gov/?c=atbdt=atbd
     float ndsi = floatfill;
     if (rfl[W550] > filltest && rfl[W2250] > filltest) {
         ndsi = (rfl[W550] - rfl[W1610]) /
                 (rfl[W550] + rfl[W1610]);
     } else {
         mask = 1;
     }
     if (ndsi > 0.35 && rfl[W865] > 0.11 && rfl[W550] > 0.1 ){
//         && p_->rfl_[W11000] < 286.0) {
         snow1 = 1;
     }
//     if (ndsi > -0.2 && rfl[W865] > 0.11 &&
//         p_->rfl_[W11000] < 286.0) {
     if (ndsi > -0.2 && rfl[W865] > 0.11) {
         snow2 = 1;
     }
//     if ((p_->dstar_ >= 1.06 && solz < 68.0) &&
//         (gzflg < 1 || (gzflg > 5 && gzflg != 26 && gzflg != 27
//          && gzflg != 24))) {
//         mask = 0;
//     }
     if (snow1 == 1) {
         mask = 1;
     }

     return status;
}

/**************************************************************************
 * NAME: DbCloudMaskLandOCI::compute_2()
 *
 * DESCRIPTION: Initial set of cloud mask filters applied to modis units
 *
 *************************************************************************/

int DbCloudMaskLandOCI::compute_2( const int iy, const int ix,
                                short& mask, const short snow2 )
{
    int status = DT_SUCCESS;

    if (p_->rfl_[W412] < filltest) {
        mask = shortfill;
        return status;
    }
    float lat = g_->in_->latitude;
    float lon = g_->in_->longitude;
    float solz = g_->in_->solar_zenith;
    int ilat = (int)((lat + 90.)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = (int)((lon + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    int gzflg = p_->gz_lut_->GEOZONE_FLAG[ilat][ilon];
    int sfcstd = p_->sp_lut_->SURFACE_ELEVATION[ilat][ilon];
    int month = g_->in_->start_month;

    if (p_->ler412_ > 50.0) {
        mask = 1;
    }
    if (p_->sr670_ < 0.1) {
        if (p_->rfl_[W2250] > 0.05 &&
            p_->ler412_ > 12.0) {
            mask = 1;
        }
    }
    if (p_->sr670_ < 0.08) {
//        if (p_->dstar_ < 1.12) {
//            if (p_->rfl_[W11000] > 270 &&
//                p_->rfl_[W11000] < 274 &&
            if (p_->ler412_ > 40.0) {
                mask = 1;
            }
            if (solz > 72.0) {
//                if (p_->rfl_[W11000] > 270 &&
//                    p_->rfl_[W11000] < 275 &&
                if (p_->ler412_ > 20.0) {
                    mask = 1;
                }
            }
//        }
    }
// Define pixel range
    int iy1 = max(iy-1,0);
    int iy2 = min(iy+1,g_->lines_-1);
    int ix1 = max(ix-1,0);
    int ix2 = min(ix+1,g_->pixels_-1);
// Compute Ler minmax ratios
    float lermin = 999.0;
    float lermax = -999.0;
    for (int cx=ix1; cx<=ix2; cx++) {
        for (int cy=iy1; cy<=iy2; cy++) {
//            float lr = p_->ler412_a_[cy][cx];
//            if (lr > filltest) {
//                lermin = (lr<lermin) ? lr : lermin;
//                lermax = (lr>lermax) ? lr : lermax;
//            }
            lermin = 1;
            lermax = 1;
        }
    }

    if (lermax<1.E-8 || lermin<1.E-8 ) {
        mask = floatfill;
        return DT_SUCCESS;
    }
    float minmax_ler = lermax/lermin + DbProcessLand::delta;
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
        if (p_->ler412_ > filltest &&
            p_->ler488_ > filltest &&
            p_->ler670_ > filltest) {
//            if (p_->rfl_[W12000] < 280.0 &&
            if (p_->ler488_/p_->ler670_ < 0.65 &&
                p_->ler412_ > 18.0) {
                mask = 0;
            }
        }
    }

    return status;
}

/**************************************************************************
 * NAME: DbSmokeMaskOCI()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbSmokeMaskOCI::DbSmokeMaskOCI() {
    p_ = 0;
}

DbSmokeMaskOCI::DbSmokeMaskOCI(Granule* granule, DbProcessLand* proc) :
        Mask(granule) {
    p_ = proc;
}

DbSmokeMaskOCI::~DbSmokeMaskOCI() {
}

/**************************************************************************
 * NAME: DbSmokeMaskOCI::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/

int DbSmokeMaskOCI::compute( const int iy, const int ix, short& mask  )
{
    int status = DT_SUCCESS;

    if (p_->rfl_[W1380] < filltest) {
        mask = shortfill;
        return status;
    }
    mask = 0;

    int ix1 = max(ix-1,0);
    int ix2 = min(ix+1,g_->pixels_-1);
    int iy1 = max(iy-1,0);
    int iy2 = min(iy+1,g_->lines_-1);
//    int cntbt11 = 0;
    int cnt2250_1 = 0;
    int cnt2250_2 = 0;

    float tler[3][3][3];
    float tsr[3][3][3];
    float td[3][3][3];
    memset(&tler[0][0][0],0,sizeof(tler));
    memset(&tsr[0][0][0],0,sizeof(tsr));
    memset(&td[0][0][0],0,sizeof(td));
    for (int cy=iy1; cy<=iy2; cy++) {
        for (int cx=ix1; cx<=ix2; cx++) {
//            tler[0][cy-iy1][cx-ix1] = p_->ler412_a_[cy][cx];
//            tler[1][cy-iy1][cx-ix1] = p_->ler488_a_[cy][cx];
//            tler[2][cy-iy1][cx-ix1] = p_->ler670_a_[cy][cx];
//            tsr[0][cy-iy1][cx-ix1] = p_->sr412_a_[cy][cx];
//            tsr[1][cy-iy1][cx-ix1] = p_->sr488_a_[cy][cx];
//            tsr[2][cy-iy1][cx-ix1] = p_->sr670_a_[cy][cx];
//            td[0][cy-iy1][cx-ix1] = p_->rfl_[W1380][cy][cx];
//            td[1][cy-iy1][cx-ix1] = p_->rfl_[W2250][cy][cx];
//            td[2][cy-iy1][cx-ix1] = p_->rfl_[W11000][cy][cx];

//            if (p_->rfl_[W11000][cy][cx] > 286.0) cntbt11++;
            if (p_->rfl_[W2250][cy][cx] < 0.035) cnt2250_1++;
            if (p_->rfl_[W2250][cy][cx] < 0.065) cnt2250_2++;
        }
    }
    float rat488670 = p_->ler488_ / p_->ler670_;
    float rat412488 = p_->ler412_ / p_->ler488_;
//    if (p_->ler412_ > 12.0 && cntbt11 == 9) {
    if (p_->ler412_ > 12.0) {
      if ((p_->sr670_ < 0.08 && cnt2250_1 == 9) ||
          (p_->sr670_ >= 0.08 && cnt2250_2 == 9)) {
        if (p_->rfl_[W2250] < 0.025) {
          mask = 1;
        } else {
          if (p_->ler412_ > 20.0 &&
              p_->rfl_[W2250] < 0.04) {
              mask = 1;
          } else {
            if (p_->ler488_ > 0.0 && p_->ler670_ > 0.0) {
              rat488670 = p_->ler488_ / p_->ler670_;
              if (rat488670 > 0.88) {
                  mask = 1;
              } else {
                if (rat488670 > 0.7 && rat488670 < 0.88 &&
                    p_->rfl_[W2250] < 0.04 &&
                    p_->rfl_[W1380] > 0.0015) {
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
 * NAME: DbPyrocbMaskOCI()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbPyrocbMaskOCI::DbPyrocbMaskOCI() {
    p_ = 0;
}

DbPyrocbMaskOCI::DbPyrocbMaskOCI(Granule* granule, DbProcessLand* proc) :
        Mask(granule) {
    p_ = proc;
}

DbPyrocbMaskOCI::~DbPyrocbMaskOCI() {
}

/**************************************************************************
 * NAME: DbPyrocbMaskOCI::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/

int DbPyrocbMaskOCI::compute( const int iy, const int ix, short& mask )
{
    int status = DT_SUCCESS;

    if (p_->rfl_[W1380] < filltest) {
        mask = shortfill;
        return status;
    }
    mask= 0;
    float rat412488 = p_->ler412_ / p_->ler488_;
    if (p_->rfl_[W1380] > 0.06  &&
        p_->rfl_[W2250] > 0.2   &&
//        p_->rfl_[W11000] < 268.0 &&
        p_->ler670_   > 24.0  &&
        rat412488 < 0.7) {
        mask = 1;
    }

    return status;
}

/**************************************************************************
 * NAME: DbHighAltSmokeMaskOCI()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbHighAltSmokeMaskOCI::DbHighAltSmokeMaskOCI() {
    p_ = 0;
}

DbHighAltSmokeMaskOCI::DbHighAltSmokeMaskOCI(Granule* granule, DbProcessLand* proc):
        Mask(granule) {
    p_ = proc;
}

DbHighAltSmokeMaskOCI::~DbHighAltSmokeMaskOCI() {
}

/**************************************************************************
 * NAME: DbHighAltSmokeMaskOCI::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/

int DbHighAltSmokeMaskOCI::compute( const int iy, const int ix, short& mask )
{
    int status = DT_SUCCESS;

    if (p_->rfl_[W1380] < filltest) {
        mask = shortfill;
        return status;
    }
    mask= 0;
    float lat = g_->in_->latitude;
    float lon = g_->in_->longitude;

    int ilat = (int)((lat + 90.)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = (int)((lon + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    int gzflg = p_->gz_lut_->GEOZONE_FLAG[ilat][ilon];

    int ix1 = max(ix-1,0);
    int ix2 = min(ix+1,g_->pixels_-1);
    int iy1 = max(iy-1,0);
    int iy2 = min(iy+1,g_->lines_-1);
    int cntW2250 = 0;
//    int cntbt11 = 0;
    float tler[3][3][3];
    float tsr[3][3][3];
    float td[3][3][3];
    memset(&tler[0][0][0],0,sizeof(tler));
    memset(&tsr[0][0][0],0,sizeof(tsr));
    memset(&td[0][0][0],0,sizeof(td));
    for (int cy=iy1; cy<=iy2; cy++) {
        for (int cx=ix1; cx<=ix2; cx++) {
//            tler[0][cy-iy1][cx-ix1] = p_->ler412_a_[cy][cx];
//            tler[1][cy-iy1][cx-ix1] = p_->ler488_a_[cy][cx];
//            tler[2][cy-iy1][cx-ix1] = p_->ler670_a_[cy][cx];
//            tsr[0][cy-iy1][cx-ix1] = p_->sr412_a_[cy][cx];
//            tsr[1][cy-iy1][cx-ix1] = p_->sr488_a_[cy][cx];
//            tsr[2][cy-iy1][cx-ix1] = p_->sr670_a_[cy][cx];
//            td[0][cy-iy1][cx-ix1] = p_->rfl_[W1380][cy][cx];
//            td[1][cy-iy1][cx-ix1] = p_->rfl_[W2250][cy][cx];
//            td[2][cy-iy1][cx-ix1] = p_->rfl_[W11000][cy][cx];
//            if (p_->rfl_[W11000][cy][cx] > 280.0) cntbt11++;
            if (p_->rfl_[W2250][cy][cx] < 0.25) cntW2250++;
        }
    }
    float rat488670 = p_->ler488_ / p_->ler670_;
    float rat412488 = p_->ler412_ / p_->ler488_;

//    if (p_->ler412_ > 12.0 && cntbt11 == 9) {
    if (p_->ler412_ > 12.0) {
        if ((p_->sr670_ < 0.08 && cntW2250 == 9) ||
            (p_->sr670_ >= 0.08 && cntW2250 == 9)) {
            if (p_->ler488_ > 18.0 &&
                (rat488670 > 0.7 && rat488670 < 0.88) &&
                p_->rfl_[W2250] < 0.22 &&
                p_->rfl_[W1380] > 0.0032) {
                mask = 1;
            }
            if (p_->ler488_ > 18.0 &&
                (rat488670 > 0.7 && rat488670 < 0.88) &&
                p_->rfl_[W2250] < 0.22 && rat412488 < 0.8) {
                mask = 1;
            }
            if (p_->ler488_ > 18.0 &&
                p_->rfl_[W1380] > 0.0015 && rat412488 < 0.8) {
//                p_->rfl_[W11000] > 268.0 && rat412488 < 0.8) {
                mask = 1;
            }
            if (gzflg == 31) {
                if (p_->ler488_ > 18.0 &&
                    (rat488670 > 0.9 && rat488670 < 1.06) &&
                    (rat412488 > 0.7 && rat412488 < 0.88) &&
                    p_->rfl_[W2250] < 0.22 &&
                    p_->rfl_[W1380] > 0.0032) {
                    mask = 1;
                }
            }
        }
    }
    if (p_->sr670_ >= 0.12 || p_->rfl_[W2250] > 0.06 ||
        p_->ler670_ > 50.0 || p_->rfl_[W1380] < 0.0005 ||
        p_->rfl_[W1380] > 0.013) {
        mask = 0;
    }
    if (rat412488/rat488670 > 0.94) {
        mask = 0;
    }
    return status;
}

/**************************************************************************
 * NAME: DbCloudMaskOceanOCI()
 *
 * DESCRIPTION: Class Constructor / Destructor
 *
 *************************************************************************/

DbCloudMaskOceanOCI::DbCloudMaskOceanOCI() {
    p_ = 0;
}

DbCloudMaskOceanOCI::DbCloudMaskOceanOCI(Granule* granule, DbProcessOcean* proc) :
        Mask(granule) {
    p_ = proc;
}

DbCloudMaskOceanOCI::~DbCloudMaskOceanOCI() {
}

/**************************************************************************
 * NAME: DbCloudMaskOceanOCI::compute()
 *
 * DESCRIPTION:
 *
 *************************************************************************/

int DbCloudMaskOceanOCI::compute( const int iy, const int ix, short& mask )
{
    int status = DT_SUCCESS;

    if (p_->rfl_[W412] < filltest) {
        mask = shortfill;
        return DT_SUCCESS;
    }
    float lat = g_->in_->latitude;
    float solz = g_->in_->solar_zenith;
    float m01_avg     = 0;
    float m01_stddev  = 0;
    float m08_avg     = 0;
    float m08_stddev  = 0;
    mask = 0;

    int ipmin = max(ix-1,0);
    int ipmax = min(ix+1, g_->pixels_-1);
    int ilmin = max(iy-1,0);
    int ilmax = min(iy+1, g_->lines_-1);
// M01 W412
    int cnt = 0;
    for (int il=ilmin; il<=ilmax; il++) {
        for (int ip=ipmin; ip<=ipmax; ip++) {
            if (p_->rfl_[W412][il][ip] > filltest &&
                    g_->in_->land_water[il][ip] == 0) {
                m01_avg += p_->rfl_[W412][il][ip];
                cnt++;
            }
        }
    }
    if (cnt >= 2) {
        m01_avg /= cnt;
        cnt = 0;
        for (int il=ilmin; il<=ilmax; il++) {
            for (int ip=ipmin; ip<=ipmax; ip++) {
                if (p_->rfl_[W412][il][ip] > filltest &&
                        g_->in_->land_water[il][ip] == 0) {
                    m01_stddev += pow((p_->rfl_[W412][il][ip]-m01_avg),2);
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
    for (int il=ilmin; il<=ilmax; il++) {
        for (int ip=ipmin; ip<=ipmax; ip++) {
            if (p_->rfl_[W1240][il][ip] > filltest &&
                    g_->in_->land_water[il][ip] == 0) {
                m08_avg += p_->rfl_[W1240][il][ip];
                cnt++;
            }
        }
    }
    if (cnt >= 2) {
        m08_avg /= cnt;
        cnt = 0;
        for (int il=ilmin; il<=ilmax; il++) {
            for (int ip=ipmin; ip<=ipmax; ip++) {
                if (p_->rfl_[W1240][il][ip] > filltest &&
                        g_->in_->land_water[il][ip] == 0) {
                    m08_stddev += pow((p_->rfl_[W1240][il][ip]-m08_avg),2);
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
             mask = 1;
        }
    } else {
        if ((m01_stddev > filltest &&
             m01_stddev > M01_STDV_THOLD*cosza) ||
            (m08_stddev > filltest &&
             m08_stddev > M08_STDV_THOLD*cosza)) {
             mask = 1;
        }
    }
//   -- perform the check on the 1.38 um band (M09)
    if (p_->rfl_[W1380] > M09_THOLD*cosza) {
        mask = 1;
    }
//   -- perform check on 488 nm band (M03)
    if (p_->rfl_[W488] >  M03_THOLD*cosza) {
        mask = 1;
    }

    return status;
}



