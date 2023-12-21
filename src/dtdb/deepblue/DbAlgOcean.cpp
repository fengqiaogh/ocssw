/*******************************************************************************
 *
 * NAME: DbAlgOcean.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute processings for a given DDProcess object class.
 *
 *  Created on: October 13, 2018
 *      Author: Sam Anderson, DT
 *
 *  Modified:
 *
 *******************************************************************************/

#include <math.h>
#include <fstream>
#include <math.h>
#include <levmar.h>

#include <DDProcess.h>
#include <DDOptions.h>
#include "deepblue/DbLutNetcdf.h"
#include "deepblue/DbMask.h"
#include "deepblue/DbAlgorithm.h"

#define real __minpack_real__

using namespace std;

/**
*  Function to be minimized in cminpack lmdif1
*/
void fcn(const int* m, const int* n, const double* x,
            double* fvec, double* fdif, int* iflag);

/**
*  Function to be minimized in levmar
*/
void fcn_lm(double *fvec, double *x, const int m, const int n, void *adata);

float rfl[NOWL];
float plut[NOWL][NFMF3][NAOT1];
int num_fmf;
int num_aot;


/**************************************************************************
 * NAME: DbAlgOcean()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DbAlgOcean::DbAlgOcean()
{
}

/**************************************************************************
 * NAME: ~DbAlgOcean()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DbAlgOcean::~DbAlgOcean()
{
    delete bath_lut_;
    delete chl_lut_;
    delete oaLut_[DUST];
    delete oaLut_[FINE];
    delete oaLut_[MARI];
    delete oaLut_[MIX];
    delete cm_;
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * processing operations.
 *
 *************************************************************************/

int DbAlgOcean::initialize( map<string, ddata*> imap )
{
	int status = DTDB_SUCCESS;

    status = DbAlgorithm::initialize( imap );
    if (status != DTDB_SUCCESS) {
        std::cerr << "DbAlgLand:: Base class initialization failure" << std::endl;
        return status;
    }
    status = initialize_LUT_data();
    if (status != DTDB_SUCCESS) {
        std::cerr << "DbAlgOcean:: LUT initialization failure" << std::endl;
        return status;
    }
	cm_ = new DbCloudMaskOcean(this);

	return status;
}

/**************************************************************************
 * NAME: initialize_LUT_data()
 *
 * DESCRIPTION: Read all calibration LUTs from binary files and create the
 * radiance and brightness temperature output LUTs.
 *
 *************************************************************************/

int DbAlgOcean::initialize_LUT_data()
{
	int status = DTDB_SUCCESS;

// Load LUTs
    DbLutNetcdf* lutgen = new DbLutNetcdf();
    bath_lut_ = new dbBathymetryLUT;
    memset(bath_lut_, 0, sizeof(dbBathymetryLUT));
    status = lutgen->read_bathymetry_lut(bath_lut_);
    chl_lut_ = new dbChlLUT();
    memset(chl_lut_, 0, sizeof(dbChlLUT));
    status = lutgen->read_chl_lut(chl_lut_);
    oaLut_[DUST] = new dbOceanAerosolLUMA;
    status = lutgen->read_ocean_aero_lut(oaLut_[DUST], LUT_OCEAN_AEROSOL_DUST);
    oaLut_[FINE] = new dbOceanAerosolLUMA;
    status = lutgen->read_ocean_aero_lut(oaLut_[FINE], LUT_OCEAN_AEROSOL_FINE);
    oaLut_[MARI] = new dbOceanAerosolLUMA;
    status = lutgen->read_ocean_aero_lut(oaLut_[MARI], LUT_OCEAN_AEROSOL_MARI);
    oaLut_[MIX] = new dbOceanAerosolLUMA;
    status = lutgen->read_ocean_aero_lut(oaLut_[MIX], LUT_OCEAN_AEROSOL_MIX);
    delete lutgen;

	return status;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Virtual function computes processing.
 *
 *************************************************************************/

map<string, ddata*> DbAlgOcean::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap)
{
    int status = get_inputs(start, count, imap);
	size_t iy = start[0];
	size_t ix = start[1];

// -- Skip pixels with invalid angles
    if (((solz_ > threshsolz_) || (solz_ < 0.0)) && bmasksolz_) {
    	l2_flags_ |= (unsigned int) flags::HISOLZEN;
    	return set_fills();
    }
    if (((senz_ < 0.0) || (senz_ > threshsenz_)) && bmasksenz_) {
    	l2_flags_ |= (unsigned int) flags::HISATZEN;
    	return set_fills();
    }
    if ((raa_ < 0.0) || (raa_ > 180.0) || (ws_ < 0.0) || (ws_ > 99.0)) {
    	return set_fills();
    }
// -- Skip pixels with invalid reflectances
    if ((rfl_[(size_t)rhot_band::W490] < 0) || (rfl_[(size_t)rhot_band::W550] < 0) || (rfl_[(size_t)rhot_band::W670] < 0) ||
        (rfl_[(size_t)rhot_band::W865] < 0)) {
//        std::cerr << "DbAlgOcean:: Invalid reflectances at "<< iy << ":" << ix << std::endl;
    	l2_flags_ |= (int) flags::BOWTIEDEL;
    	return set_fills();
    }
	// -- skip bad pixels (cloud, deletion pixels, etc)
	mask_cm_ = cloud_mask_;
	if (mask_cm_ == DFILL_UBYTE) {
		cm_->compute( mask_cm_ );
		cloud_mask_ = mask_cm_;
	}
	if (bmaskcloud_ && cloud_mask_) {
		l2_flags_ |= (int) flags::CLDICE;
		return set_fills();
	}
// -- skip glint-affected regions.
	compute_scatter_angle(scatter_angle_);
	status = compute_glint_refl(glint_angle_);
	if (status != 0) {
		std::cerr << "DbAlgOcean:: Failed to calculate glint angle" << std::endl;
		return set_fills();
	}
// Old threshold was 0.003. This recovers area and doesn't introduce problems.
	if (bmaskglint_) {
//		if (glint > 0.005) {
//        std::cerr << "DbAlgOcean:: Glint exceeds threshold" << std::endl;
		if (glint_refl_ > threshglint_) {
			l2_flags_ |= (int) flags::HIGLINT;
			return set_fills();
		}
	}

// convert from VIIRS to IOF units
    double cossza = cos(solz_*DEGtoRAD);
    for ( int i=0; i<NTWL; i++) {
    	if (rfl_[i] > DFILL_TEST) {
    		rfl_[i] *= (cossza/M_PI);
    	}
        for (size_t il=0; il<=2; il++) {
            for (size_t ip=0; ip<=2; ip++) {
            	if (rfla_[i][il][ip] > DFILL_TEST) {
            		rfla_[i][il][ip] *= (cossza/M_PI);
            	}
            }
        }
    }

    set_fill_out();
    OMODEL_ENUM min_model = NMODEL;

    if (bgascorrect_) {
		compute_gas_correction();
		for (size_t ib=0; ib<NTWL; ib++) {
	    	if (rfl_[ib] > DFILL_TEST) {
	    		rfl_[ib] *= gasc_[ib];
	    	}
		}
    }

    int ilat = min((int)((lat_ + 90.)*10.0 + 0.5),1800-1);
    int ilon = min((int)((lon_ + 180.0)*10.0 + 0.5),3600-1);
    chl_  = chl_lut_->log_chl[ilat][ilon][month_-1];
    int xind = min((int)((lon_ + 180.0)*60.0 + 0.5), 21600);
    int yind = min((int)((lat_ + 90.0)*60.0 + 0.5), 10800);
    bathy_ = bath_lut_->z[yind][xind];

    for (int i=0; i<NDBMDL+1; i++) {
        memset(&oOut_[i], 0, sizeof(osOut));
    }

// -- run ocean retrieval for this pixel. Check for status < 0 as 0 is a successful
// -- retrieval while 1 means the pixel was too turbid -- not an error.
    oOut_[DUST].model_flag = DUST;
    status = run_inversion(iy, ix, oaLut_[DUST], &oOut_[DUST]);
    if (status < 0) {
       std::cerr << "DbAlgOcean:: Coarse aerosol inversion failed" << std::endl;
    }
    oOut_[FINE].model_flag = FINE;
    status = run_inversion(iy, ix, oaLut_[FINE], &oOut_[FINE]);
    if (status < 0) {
       std::cerr << "DbAlgOcean:: Fine aerosol inversion failed" << std::endl;
    }
    oOut_[MARI].model_flag = MARI;
    status = run_inversion(iy, ix, oaLut_[MARI], &oOut_[MARI]);
    if (status < 0) {
       std::cerr << "DbAlgOcean:: Maritime aerosol inversion failed" << std::endl;
    }
    oOut_[MIX].model_flag = MIX;
    status = run_inversion(iy, ix, oaLut_[MIX], &oOut_[MIX]);
    if (status < 0) {
       std::cerr << "DbAlgOcean:: Mixed aerosol inversion failed" << std::endl;
    }

    if (oOut_[DUST].ss > 0) {
        if ((oOut_[DUST].ss < oOut_[FINE].ss || oOut_[FINE].ss < 0) &&
            (oOut_[DUST].ss < oOut_[MARI].ss || oOut_[MARI].ss < 0) &&
            (oOut_[DUST].ss < oOut_[MIX].ss || oOut_[MIX].ss < 0)) {
            min_model = DUST;
        }
    }
    if (oOut_[FINE].ss > 0) {
        if ((oOut_[FINE].ss < oOut_[DUST].ss || oOut_[DUST].ss < 0) &&
            (oOut_[FINE].ss < oOut_[MARI].ss || oOut_[MARI].ss < 0) &&
            (oOut_[FINE].ss < oOut_[MIX].ss || oOut_[MIX].ss < 0)) {
            min_model = FINE;
        }
    }
    if (oOut_[MARI].ss > 0) {
        if ((oOut_[MARI].ss < oOut_[DUST].ss || oOut_[DUST].ss < 0) &&
            (oOut_[MARI].ss < oOut_[FINE].ss || oOut_[FINE].ss < 0) &&
            (oOut_[MARI].ss < oOut_[MIX].ss || oOut_[MIX].ss < 0)) {
            min_model = MARI;
        }
    }
    if (oOut_[MIX].ss > 0) {
        if ((oOut_[MIX].ss < oOut_[DUST].ss || oOut_[DUST].ss < 0) &&
            (oOut_[MIX].ss < oOut_[FINE].ss || oOut_[FINE].ss < 0) &&
            (oOut_[MIX].ss < oOut_[MARI].ss || oOut_[MARI].ss < 0)) {
            min_model = MIX;
        }
    }

    if (min_model < NMODEL) {
        oOut_[BEST].aot550 = oOut_[min_model].aot550;
        oOut_[BEST].fmf = oOut_[min_model].fmf;
        oOut_[BEST].ae = oOut_[min_model].ae;
        oOut_[BEST].ss = oOut_[min_model].ss;
        oOut_[BEST].chl = oOut_[min_model].chl;
        oOut_[BEST].alg_flag = oOut_[min_model].alg_flag;
        oOut_[BEST].model_flag = oOut_[min_model].model_flag;
        for (int i=0; i<NOWL; i++) {
            oOut_[BEST].aot[i] = oOut_[min_model].aot[i];
        }
    }

    // write to generic outputs

    qual_flag_ = (oOut_[BEST].aot550 < DFILL_TEST) ? 0 : 3;
	if (oOut_[BEST].model_flag == 1 ) {
		aerosol_type_ = 0;
    }else if (oOut_[BEST].model_flag == 2 ) {
    	aerosol_type_ = 7;
    }else if (oOut_[BEST].model_flag == 3 ) {
    	aerosol_type_ = 6;
    }else if (oOut_[BEST].model_flag == 4 ) {
    	aerosol_type_ = 5;
    } else {
    	aerosol_type_ = 8;
    }
    error_flag_ = DFILL_SHORT;

    scatter_ang_ = scatter_angle_;
    glint_ang_ = glint_angle_;
    sse_ = oOut_[BEST].ss;
    fmf_ = oOut_[BEST].fmf;
    aot_550_ = oOut_[BEST].aot550;
    ae1_ = oOut_[BEST].ae;
    ae2_ = DFILL_FLOAT;
    ndv_ = DFILL_FLOAT;
    chlor_ = chl_;
    aot_[0] = DFILL_FLOAT;
    for ( int ib=0; ib<NOWL; ib++ ) {
    	aot_[ib+1] = oOut_[BEST].aot[ib];
    }
    for ( int ib=0; ib<NLWL+1; ib++ ) {
        sr_[ib] = DFILL_FLOAT;
        ssa_[ib] = DFILL_FLOAT;
    }
    for (int ib = 0; ib < DB_RFL_BANDS; ib++) {
        rfl_[ib] *= (M_PI/cossza);
    }

    return set_outputs();
}

/**************************************************************************
 * NAME: run_inversion()
 *
 * DESCRIPTION: Run inversion to extract tau estimates.
 *
 **************************************************************************
 * Library Function
 *
 * NAME: lmdif1
 *
 *     the purpose of lmdif1 is to minimize the sum of the squares of
 *     m nonlinear functions in n variables by a modification of the
 *     levenberg-marquardt algorithm. this is done by using the more
 *     general least-squares solver lmdif. the user must provide a
 *     subroutine which calculates the functions. the jacobian is
 *     then calculated by a forward-difference approximation.
 *
 *     the subroutine statement is
 *
 *       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
 *
 *     where
 *
 *       fcn is the name of the user-supplied subroutine which
 *         calculates the functions. fcn must be declared
 *         in an external statement in the user calling
 *         program, and should be written as follows.
 *
 *         subroutine fcn(m,n,x,fvec,iflag)
 *         integer m,n,iflag
 *         double precision x(n),fvec(m)
 *         ----------
 *         calculate the functions at x and
 *         return this vector in fvec.
 *         ----------
 *         return
 *         end
 *
 *         the value of iflag should not be changed by fcn unless
 *         the user wants to terminate execution of lmdif1.
 *         in this case set iflag to a negative integer.
 *
 *       m is a positive integer input variable set to the number
 *         of functions.
 *
 *       n is a positive integer input variable set to the number
 *         of variables. n must not exceed m.
 *
 *       x is an array of length n. on input x must contain
 *         an initial estimate of the solution vector. on output x
 *         contains the final estimate of the solution vector.
 *
 *       fvec is an output array of length m which contains
 *         the functions evaluated at the output x.
 *
 *       tol is a nonnegative input variable. termination occurs
 *         when the algorithm estimates either that the relative
 *         error in the sum of squares is at most tol or that
 *         the relative error between x and the solution is at
 *         most tol.
 *
 *       info is an integer output variable. if the user has
 *         terminated execution, info is set to the (negative)
 *         value of iflag. see description of fcn. otherwise,
 *         info is set as follows.
 *
 *         info = 0  improper input parameters.
 *
 *         info = 1  algorithm estimates that the relative error
 *                   in the sum of squares is at most tol.
 *
 *         info = 2  algorithm estimates that the relative error
 *                   between x and the solution is at most tol.
 *
 *         info = 3  conditions for info = 1 and info = 2 both hold.
 *
 *         info = 4  fvec is orthogonal to the columns of the
 *                   jacobian to machine precision.
 *
 *         info = 5  number of calls to fcn has reached or
 *                   exceeded 200*(n+1).
 *
 *         info = 6  tol is too small. no further reduction in
 *                   the sum of squares is possible.
 *
 *         info = 7  tol is too small. no further improvement in
 *                   the approximate solution x is possible.
 *
 *       iwa is an integer work array of length n.
 *
 *       wa is a work array of length lwa.
 *
 *       lwa is a positive integer input variable not less than
 *         m*n+5*n+m.
 *
 *     subprograms called
 *
 *       user-supplied ...... fcn
 *
 *       minpack-supplied ... lmdif
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 *
 *************************************************************************/

int DbAlgOcean::run_inversion(size_t iy, size_t ix, dbOceanAerosolLUMA* lut,
                                  osOut* iout )
{
    int status = DTDB_SUCCESS;

    // -- find index for relative azimuth angle for linear interpolation
    int rbeg = locate(NVRAA, &lut->raa[0], raa_, status);
    if (status != 0) {
      std::cerr << "DbAlgOcean::run_inversion:: Failed to locate raa index" << std::endl;
      status = -1;
      return status;
    }
    float raa_frac = (raa_ - lut->raa[rbeg]) / (lut->raa[rbeg+1]-lut->raa[rbeg]);

// -- find indices for sza for linear interpolation
    int mbeg = locate(NVSZA, &lut->sza[0], solz_, status);
    if (status != 0) {
      std::cerr << "DbAlgOcean::run_inversion:: Failed to locate sza index" << std::endl;
      return status;
    }
    float sza_frac = (solz_-lut->sza[mbeg]) / (lut->sza[mbeg+1]-lut->sza[mbeg]);

// -- find indices for vza for linear interpolation
    int nbeg = locate(NVVZA, &lut->vza[0], senz_, status);
    if (status != 0) {
      std::cerr << "DbAlgOcean::run_inversion:: Failed to locate vza index" << std::endl;
      return status;
    }
    float vza_frac = (senz_-lut->vza[nbeg]) / (lut->vza[nbeg+1]-lut->vza[nbeg]);

// -- find indices for windspeed
    float ws_frac = 0.0;
    int wbeg = locate(NWS, &lut->wspd[0], ws_, status);
    if (status != 0) {
      switch (status) {
        case -1:
          if (ws_ < lut->wspd[0]) {
             wbeg = 0;
             ws_frac = 0.0;
          } else {
             std::cerr << "DbAlgOcean::run_inversion:: Failed windspeed lookup" << std::endl;
             return status;
          }
          break;
        case 1:
          if (ws_ > lut->wspd[lut->nwspd-1]) {
             wbeg = lut->nwspd - 2;
             ws_frac = 1.0;
          } else {
             std::cerr << "DbAlgOcean::run_inversion:: Failed windspeed lookup" << std::endl;
             return status;
          }
          break;
        case -2:
           std::cerr << "DbAlgOcean::run_inversion:: Windspeeds not monotonic" << std::endl;
           return status;
           break;
        default:
           std::cerr << "DbAlgOcean::run_inversion:: Unexpected status returned from locate" << std::endl;
           return status;
           break;
      }
    } else {
      ws_frac = (ws_ - lut->wspd[wbeg])/(lut->wspd[wbeg+1]-lut->wspd[wbeg]);
    }
// -- find indices for chl
    float chl_frac = 0.0;
    int cbeg = locate(NCHL, &lut->chl[0], chl_, status);
    if (status != 0) {
      switch (status) {
        case (-1):
          if (chl_ < lut->chl[0]) {
            cbeg = 0;
            chl_frac = 0.0;
          } else {
            std::cerr << "DbAlgOcean::run_inversion:: Failed CHL lookup" << std::endl;
            return status;
          }
          break;
        case (1):
          if (chl_ > lut->chl[lut->nchl-1]) {
            cbeg = lut->nchl - 2;
            chl_frac = 1.0;
          } else {
            std::cerr << "DbAlgOcean::run_inversion:: Failed CHL lookup" << std::endl;
            return status;
          }
          break;
        case (-2):
          std::cerr << "DbAlgOcean::run_inversion:: CHL not monotonic" << std::endl;
          return status;
          break;
        default:
          std::cerr << "DbAlgOcean::run_inversion:: Unexpected status returned from locate" << std::endl;
          return status;
          break;
      }
    } else {
      chl_frac = (chl_- lut->chl[cbeg])/(lut->chl[cbeg+1]-lut->chl[cbeg]);
    }

// Note: in here, the yXXX quantities are the interpolated I/F.
// the xXXX quantities are the value of parameter (e.g. angle) we are interpolating too.
// Those are not actually needed, but calculate them so you can print as a sanity check
// that you are interpolating where you think you are.
    float ysza[2];
    float yvza[2];
    float ywspd[2];
    float ychl[2];
    num_fmf = lut->nfmf;
    num_aot = lut->naot;
    for (size_t iwv=0; iwv<lut->nwave; iwv++) {        // wavelength
      for (size_t ifmf=0; ifmf<lut->nfmf; ifmf++) {    // fine mode fraction
        for (size_t iaot=0; iaot<lut->naot; iaot++) {  // AOT
          for (size_t iws=0; iws<=1; iws++) {          // WSPD
            for (size_t isza=0; isza<=1; isza++) {     // SZA
              for (size_t ivza=0; ivza<=1; ivza++) {   // VZA
                switch (iwv) {
// M03
                  case (0):
                    for (size_t ichl=0; ichl<=1; ichl++) {
                       ychl[ichl] = lut->m03[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg][nbeg+ivza][mbeg+isza]*(1.0-raa_frac) +
                                     lut->m03[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg+1][nbeg+ivza][mbeg+isza]*raa_frac;
                    }
                    yvza[ivza] = ychl[0]*(1.0-chl_frac) + ychl[1]*chl_frac;
                    break;
// M04
                  case (1):
                    for (size_t ichl=0; ichl<=1; ichl++) {
                       ychl[ichl] = lut->m04[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg][nbeg+ivza][mbeg+isza]*(1.0-raa_frac) +
                                     lut->m04[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg+1][nbeg+ivza][mbeg+isza]*raa_frac;
                    }
                    yvza[ivza] = ychl[0]*(1.0-chl_frac) + ychl[1]*chl_frac;
                    break;
// M05
                  case (2):
                    for (size_t ichl=0; ichl<=1; ichl++) {
                       ychl[ichl] = lut->m05[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg][nbeg+ivza][mbeg+isza]*(1.0-raa_frac) +
                                     lut->m05[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg+1][nbeg+ivza][mbeg+isza]*raa_frac;
                    }
                    yvza[ivza] = ychl[0]*(1.0-chl_frac) + ychl[1]*chl_frac;
                    break;
// M07
                  case (3):
                    for (size_t ichl=0; ichl<=1; ichl++) {
                       ychl[ichl] = lut->m07[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg][nbeg+ivza][mbeg+isza]*(1.0-raa_frac) +
                                     lut->m07[cbeg+ichl][wbeg+iws][ifmf][iaot][rbeg+1][nbeg+ivza][mbeg+isza]*raa_frac;
                    }
                    yvza[ivza] = ychl[0]*(1.0-chl_frac) + ychl[1]*chl_frac;
                    break;
// M08
                  case (4):
                    yvza[ivza] = lut->m08[wbeg+iws][ifmf][iaot][rbeg][nbeg+ivza][mbeg+isza]*(1.0-raa_frac) +
                                  lut->m08[wbeg+iws][ifmf][iaot][rbeg+1][nbeg+ivza][mbeg+isza]*raa_frac;
                    break;
// M10
                  case (5):
                    yvza[ivza] = lut->m10[wbeg+iws][ifmf][iaot][rbeg][nbeg+ivza][mbeg+isza]*(1.0-raa_frac) +
                                  lut->m10[wbeg+iws][ifmf][iaot][rbeg+1][nbeg+ivza][mbeg+isza]*raa_frac;
                    break;
// M11
                  case (6):
                    yvza[ivza] = lut->m11[wbeg+iws][ifmf][iaot][rbeg][nbeg+ivza][mbeg+isza]*(1.0-raa_frac) +
                                  lut->m11[wbeg+iws][ifmf][iaot][rbeg+1][nbeg+ivza][mbeg+isza]*raa_frac;
                     break;
                  default:
                     std::cerr << "DbAlgOcean::run_inversion:: Unexpected wavelength" << std::endl;
                     return status;
                     break;
                 }
              }
              ysza[isza] = yvza[0]*(1.0-vza_frac) + yvza[1]*vza_frac;
            }
            ywspd[iws] = ysza[0]*(1.0-sza_frac) + ysza[1]*sza_frac;
          }
          plut[iwv][ifmf][iaot] = ywspd[0]*(1.0-ws_frac) + ywspd[1]*ws_frac;
        }
      }
    }
    rfl[O488] = rfl_[(size_t)rhot_band::W490];
    rfl[O550] = rfl_[(size_t)rhot_band::W550];
    rfl[O670] = rfl_[(size_t)rhot_band::W670];
    rfl[O865] = rfl_[(size_t)rhot_band::W865];
    rfl[O1240] = rfl_[(size_t)rhot_band::W1240];
    rfl[O1610] = rfl_[(size_t)rhot_band::W1610];
    rfl[O2250] = rfl_[(size_t)rhot_band::W2250];

/*
//-- calculate turbid water residual. If too high, skip becayse too turbid.
//-- Otherwise, run turbid IR retrieval below.
    float turb = calc_turbid_residual(solz_, rfl_[(size_t)rhot_band::W490], rfl_[(size_t)rhot_band::W1240],
                                 rfl_[(size_t)rhot_band::W1610], rfl_[(size_t)rhot_band::W2250], rfl_[(size_t)rhot_band::W550], status);
    if (status != 0) {
       std::cerr << "DbAlgOcean::run_inversion:: Turbid water residual calculation failure" << std::endl;
       return status;
    }

//-- convert reduced look up table to differences and feed to
//-- Leverberg-Marquardt (lmdif1).
//-- For more info on lmdif1: http://www.netlib.org/minpack/lmdif1.f
//-- If over turbid water, switch to 4 IR band retrieval or skip.
    short nbands = 0;
    float rfactor = cos(solz_*DEGtoRAD)/M_PI;
    if (turb >= 0.1*rfactor) {
//      std::cerr << "DbAlgOcean::run_inversion:: Turbid water exceeds threshold." << std::endl;
      status = 1;
      return status;
    } else if ((turb > 0.011*rfactor && turb < 0.1*rfactor) ||
               (bathy_ > -20) || (pow(10,chl_) >= 3)) {
// mildly turbid water or shallow water (that is not too turbid), try IR retrieval.
// bathy > -x means water with a depth < x m. Note this also catches elevated inland water.
       for (size_t i=0; i<lut->naot; i++) {
          for (size_t j=0; j<lut->nfmf; j++) {
              plut[O488][j][i] = 0.0;
              plut[O550][j][i] = 0.0;
              plut[O670][j][i] = 0.0;
              plut[O865][j][i] = (plut[O865][j][i] - rfl_[(size_t)rhot_band::W865]) /
                      (0.08*(rfl_[(size_t)rhot_band::W865] + 0.00001));
              plut[O1240][j][i] = (plut[O1240][j][i] - rfl_[(size_t)rhot_band::W1240]) /
                      (0.07*(rfl_[(size_t)rhot_band::W1240] + 0.00001));
              plut[O1610][j][i] = (plut[O1610][j][i] - rfl_[(size_t)rhot_band::W1610]) /
                      (0.06*(rfl_[(size_t)rhot_band::W1610] + 0.00001));
              plut[O2250][j][i] = (plut[O2250][j][i] - rfl_[(size_t)rhot_band::W2250]) /
                      (0.10*(rfl_[(size_t)rhot_band::W2250] + 0.00001));
          }
       }
       nbands = 4;
       iout->alg_flag = 1;
    } else {     // clear water, do full retrieval.
        for (size_t i=0; i<lut->naot; i++) {
           for (size_t j=0; j<lut->nfmf; j++) {
              plut[O488][j][i] = (plut[O488][j][i] - rfl_[(size_t)rhot_band::W490]) /
                      (0.06*(rfl_[(size_t)rhot_band::W490] + 0.00001));
              plut[O550][j][i] = (plut[O550][j][i] - rfl_[(size_t)rhot_band::W550]) /
                      (0.06*(rfl_[(size_t)rhot_band::W550] + 0.00001));
              plut[O670][j][i] = (plut[O670][j][i] - rfl_[(size_t)rhot_band::W670]) /
                      (0.04*(rfl_[(size_t)rhot_band::W670] + 0.00001));
              plut[O865][j][i] = (plut[O865][j][i] - rfl_[(size_t)rhot_band::W865]) /
                      (0.04*(rfl_[(size_t)rhot_band::W865] + 0.00001));
              plut[O1240][j][i] = (plut[O1240][j][i] - rfl_[(size_t)rhot_band::W1240]) /
                      (0.07*(rfl_[(size_t)rhot_band::W1240] + 0.00001));
              plut[O1610][j][i] = (plut[O1610][j][i] - rfl_[(size_t)rhot_band::W1610]) /
                      (0.06*(rfl_[(size_t)rhot_band::W1610] + 0.00001));
              plut[O2250][j][i] = (plut[O2250][j][i] - rfl_[(size_t)rhot_band::W2250]) /
                      (0.10*(rfl_[(size_t)rhot_band::W2250] + 0.00001));
           }
        }
        nbands = 7;
        iout->alg_flag = 0;
    }

// -- calculate sum of squares of node points, find min and use as starting point
    double node_ss[NFMF3][NAOT1];
    double fvec[NOWL];
    double mins[2] = {0.0,0.0};
    float min_ss = 999999;
    memset(node_ss, 0, NFMF3*NAOT1*sizeof(double));
    memset(fvec, 0, NOWL*sizeof(double));

    for (size_t i=0; i<lut->naot; i++) {
       for (size_t j=0; j<lut->nfmf; j++) {
           float sum = 0;
           for (size_t k=0; k<lut->nwave; k++) {
               sum += plut[k][j][i]*plut[k][j][i];
           }
           node_ss[j][i] = sum / (nbands - 2);
           if (node_ss[j][i] < min_ss) {
               min_ss = node_ss[j][i];
               mins[0] = j;
               mins[1] = i;
           }
       }
    }

    const double tol = 1e-6;
    int   iwa[2];
    const int  lwa = 32;
    double wa[32];

// -- mins[1] = AOT index, mins[0] = FMAF index
    __minpack_func__(lmdif1)(fcn, &m, &n, &mins[0], &fvec[0], &tol,
                              &status, &iwa[0], &wa[0], &lwa);

    if (status == 0 || status > 4) { // || status > 3)  lmdif1 failed!
        iout->aot550  = DFILL_FLOAT;
        for (int j=0; j<NOWL; j++) {
            iout->aot[j] = DFILL_FLOAT;
        }
        iout->fmf = DFILL_FLOAT;
        iout->ae  = DFILL_FLOAT;
        iout->chl  = DFILL_FLOAT;
        iout->alg_flag  = DFILL_SHORT;
        iout->model_flag  = DFILL_SHORT;
        status  = -1;
        return status;
    }
*/
    double info[LM_INFO_SZ];
    double opts[LM_OPTS_SZ];

    /*Optimisation control parameters*/
    opts[0] = 1E-3; //1E-3 LM_INIT_MU; /* scale factor for initial mu */
    opts[1] = 1E-10; //1E-10 LM_INIT_MU;  convergence in gradient
    opts[2] = 1E-5; //1E-5 relative error desired in the approximate solution
    opts[3] = 1E-7; //1E-7 relative error desired in the sum of squares
    opts[4] = 1E-6; //1E-6  LM_DIFF_DELTA;  /*step used in difference approximation to the Jacobian*/

    double lowerBounds[2] = {0.0, 0.0}; //Lower bounds for two parameters
    double upperBounds[2] = {8.0, 14.0}; //Upper bounds for two parameters
    int maxit = 1000;
    double fvec[NOWL] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double fdif[NOWL] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double mins[2] = {0.0,1.0};
    const int m = 2;
    const int n = NOWL;

    /*Run optimization*/
    /*levmar routine (dlevmar), box contstraints (bc), numerical partials (diff)*/
    //For further details, see source code in: lmbc_core.c //
    status = dlevmar_bc_dif(fcn_lm, &mins[0], &fvec[0], m, n,
            lowerBounds, upperBounds, NULL, maxit, opts, info, NULL,
            NULL, NULL);

    fcn( &m, &n, &mins[0], &fvec[0], &fdif[0], 0 );

    if (status < 0) { // dlevmar_dif failed!
        return status;
    }

// -- adjust indices to maximum table values if lmdif1 was successful.
    if (mins[0] < 0) mins[0] = 0;
    if (mins[0] >= lut->nfmf-1) mins[0] =
            (double)lut->nfmf-1.001;   // needs to be just under
    if (mins[1] < 0) mins[1] = 0;      // max node value to be
    if (mins[1] >= lut->naot-1) mins[1] =
            (double)lut->naot-1.001;   // correctly interpolated.

// -- calculate the sum of the squares of the residuals of all bands
//    normalized by nbands - 2.
    iout->ss = 0.0;
    for (int i=0; i<NOWL; i++) {
        iout->ss += fdif[i]*fdif[i];
    }
    iout->ss /= (NOWL - 2);

 //   -- convert indices to AOT, FMAF, and AE values, linear interpolation used.
    int j1 = floor(mins[0]);
    int j2 = j1 + 1;
    int i1 = floor(mins[1]);
    int i2 = i1 + 1;

    iout->aot550 = lut->aot550[i1] + (mins[1]-i1)*(lut->aot550[i2]-lut->aot550[i1]);

    iout->aot[O488] = (lut->aot[O488][j1][i1]*(i2-mins[1])*(j2-mins[0])) +
                      (lut->aot[O488][j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
                      (lut->aot[O488][j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
                      (lut->aot[O488][j2][i2]*(mins[1]-i1)*(mins[0]-j1));
    iout->aot[O550] = (lut->aot[O550][j1][i1]*(i2-mins[1])*(j2-mins[0])) +
                      (lut->aot[O550][j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
                      (lut->aot[O550][j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
                      (lut->aot[O550][j2][i2]*(mins[1]-i1)*(mins[0]-j1));
    iout->aot[O670] = (lut->aot[O670][j1][i1]*(i2-mins[1])*(j2-mins[0])) +
                      (lut->aot[O670][j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
                      (lut->aot[O670][j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
                      (lut->aot[O670][j2][i2]*(mins[1]-i1)*(mins[0]-j1));
    iout->aot[O865] = (lut->aot[O865][j1][i1]*(i2-mins[1])*(j2-mins[0])) +
                      (lut->aot[O865][j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
                      (lut->aot[O865][j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
                      (lut->aot[O865][j2][i2]*(mins[1]-i1)*(mins[0]-j1));
    iout->aot[O1240] = (lut->aot[O1240][j1][i1]*(i2-mins[1])*(j2-mins[0])) +
                       (lut->aot[O1240][j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
                       (lut->aot[O1240][j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
                       (lut->aot[O1240][j2][i2]*(mins[1]-i1)*(mins[0]-j1));
    iout->aot[O1610] = (lut->aot[O1610][j1][i1]*(i2-mins[1])*(j2-mins[0])) +
                       (lut->aot[O1610][j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
                       (lut->aot[O1610][j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
                       (lut->aot[O1610][j2][i2]*(mins[1]-i1)*(mins[0]-j1));
    iout->aot[O2250] = (lut->aot[O2250][j1][i1]*(i2-mins[1])*(j2-mins[0])) +
                       (lut->aot[O2250][j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
                       (lut->aot[O2250][j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
                       (lut->aot[O2250][j2][i2]*(mins[1]-i1)*(mins[0]-j1));

    iout->ae = (lut->ae[j1][i1]*(i2-mins[1])*(j2-mins[0])) +
               (lut->ae[j1][i2]*(mins[1]-i1)*(j2-mins[0])) +
               (lut->ae[j2][i1]*(i2-mins[1])*(mins[0]-j1)) +
               (lut->ae[j2][i2]*(mins[1]-i1)*(mins[0]-j1));

    iout->fmf = lut->fmf[j1]+(mins[0]-j1)*(lut->fmf[j2]-lut->fmf[j1]);

    iout->chl = chl_;

    return status;
}


 /**************************************************************************
  * NAME: fcn_lm()
  *
  * DESCRIPTION: function for levmar
  *
  *************************************************************************/

 void fcn_lm(double *x, double *fvec, int m, int n, void *adata)
 {
     int flag;
     double fdif[NOWL] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
     string func_name = "fcn_lm";

     fcn( &n, &m, x, fvec, fdif, &flag );
 }


 /**************************************************************************
  * NAME: fcn()
  *
  * DESCRIPTION: function for lmdif1, returns lut values for indices
  * specified in x
  *
  *************************************************************************/

 void fcn(const int* m, const int* n, const double* x,
          double* fvec, double* fdif, int* iflag)
 {
    string func_name = "fcn";

// -- x values may be out of bounds for the lut.
// -- limit to size of plut
    double tx[2];
    tx[0] = x[0];
    tx[1] = x[1];
    int j = floor(x[0]);
    int i = floor(x[1]);
    if (i < 0) {
        i = 0;
//        tx[1] = 0.0;
    }
    if (j < 0) {
        j = 0;
        tx[0] = 0.0;
    }
    if (i > num_aot-2) {
        i = num_aot-2;
        tx[1] = num_aot-1;
    }
    if (j > num_fmf-2) {
        j = num_fmf-2;
        tx[0] = num_fmf-1;
    }

// -- bilinear interpolation of reduced look up table to indices x[0], x[1].
    float f[NOWL];
    f[O488] = (plut[O488][j][i])*(i+1-tx[1])*(j+1-tx[0]) +
              (plut[O488][j][i+1])*(tx[1]-i)*(j+1-tx[0]) +
              (plut[O488][j+1][i])*(i+1-tx[1])*(tx[0]-j) +
              (plut[O488][j+1][i+1])*(tx[1]-i)*(tx[0]-j);
    f[O550] = (plut[O550][j][i])*(i+1-tx[1])*(j+1-tx[0]) +
              (plut[O550][j][i+1])*(tx[1]-i)*(j+1-tx[0]) +
              (plut[O550][j+1][i])*(i+1-tx[1])*(tx[0]-j) +
              (plut[O550][j+1][i+1])*(tx[1]-i)*(tx[0]-j);
    f[O670] = (plut[O670][j][i])*(i+1-tx[1])*(j+1-tx[0]) +
              (plut[O670][j][i+1])*(tx[1]-i)*(j+1-tx[0]) +
              (plut[O670][j+1][i])*(i+1-tx[1])*(tx[0]-j) +
              (plut[O670][j+1][i+1])*(tx[1]-i)*(tx[0]-j);
    f[O865] = (plut[O865][j][i])*(i+1-tx[1])*(j+1-tx[0]) +
              (plut[O865][j][i+1])*(tx[1]-i)*(j+1-tx[0]) +
              (plut[O865][j+1][i])*(i+1-tx[1])*(tx[0]-j) +
              (plut[O865][j+1][i+1])*(tx[1]-i)*(tx[0]-j);
    f[O1240] = (plut[O1240][j][i])*(i+1-tx[1])*(j+1-tx[0]) +
              (plut[O1240][j][i+1])*(tx[1]-i)*(j+1-tx[0]) +
              (plut[O1240][j+1][i])*(i+1-tx[1])*(tx[0]-j) +
              (plut[O1240][j+1][i+1])*(tx[1]-i)*(tx[0]-j);
    f[O1610] = (plut[O1610][j][i])*(i+1-tx[1])*(j+1-tx[0]) +
              (plut[O1610][j][i+1])*(tx[1]-i)*(j+1-tx[0]) +
              (plut[O1610][j+1][i])*(i+1-tx[1])*(tx[0]-j) +
              (plut[O1610][j+1][i+1])*(tx[1]-i)*(tx[0]-j);
    f[O2250] = (plut[O2250][j][i])*(i+1-tx[1])*(j+1-tx[0]) +
              (plut[O2250][j][i+1])*(tx[1]-i)*(j+1-tx[0]) +
              (plut[O2250][j+1][i])*(i+1-tx[1])*(tx[0]-j) +
              (plut[O2250][j+1][i+1])*(tx[1]-i)*(tx[0]-j);

    float fd[NOWL];
    fd[O488] = f[O488] - rfl[O488];
    fd[O550] = f[O550] - rfl[O550];
    fd[O670] = f[O670] - rfl[O670];
    fd[O865] = f[O865] - rfl[O865];
    fd[O1240] = f[O1240] - rfl[O1240];
    fd[O1610] = f[O1610] - rfl[O1610];
    fd[O2250] = f[O2250] - rfl[O2250];

    fvec[O488] = fd[O488] / (0.06*(rfl[O488] + 0.00001));
    fvec[O550] = fd[O550] / (0.06*(rfl[O550] + 0.00001));
    fvec[O670] = fd[O670] / (0.04*(rfl[O670] + 0.00001));
    fvec[O865] = fd[O865] / (0.04*(rfl[O865] + 0.00001));
    fvec[O1240] = fd[O1240] / (0.07*(rfl[O1240] + 0.00001));
    fvec[O1610] = fd[O1610] / (0.06*(rfl[O1610] + 0.00001));
    fvec[O2250] = fd[O2250] / (0.10*(rfl[O2250] + 0.00001));

    fdif[O488] = fd[O488] / rfl[O488];
    fdif[O550] = fd[O550] / rfl[O550];
    fdif[O670] = fd[O670] / rfl[O670];
    fdif[O865] = fd[O865] / rfl[O865];
    fdif[O1240] = fd[O1240] / rfl[O1240];
    fdif[O1610] = fd[O1610] / rfl[O1610];
    fdif[O2250] = fd[O2250] / rfl[O2250];

    return;
}

/**************************************************************************
 * NAME: calc_glint_refl()
 *
 * DESCRIPTION: Calculate reflectance from glint region.
 *
 *************************************************************************/

float DbAlgOcean::calc_glint_refl(size_t iy, size_t ix, int& status)
{
    double cc = DEGtoRAD;
    float nr = 1.341;

    double zx = (sin(senz_*cc)*sin(raa_*cc))/(cos(solz_*cc)+cos(senz_*cc));
    double zy = (sin(solz_*cc)-sin(senz_*cc)*cos(raa_*cc))/(cos(solz_*cc)+cos(senz_*cc));
    double sigx = sqrt(0.003+0.00192*ws_);
    double sigy = sqrt(0.00316*ws_);
    double zeta = zx / sigx;
    double eta  = zy / sigy;
    double p = (1.0/(2.0*M_PI*sigx*sigy))*exp(-0.5*(zeta*zeta + eta*eta));
    double costwoomega = cos(senz_*cc)*cos(solz_*cc) -
                         sin(senz_*cc)*sin(solz_*cc)*cos(raa_*cc);
    double cosbeta = (cos(solz_*cc)+cos(senz_*cc))/(sqrt(2.0 + 2.0*costwoomega));
    double w = 0.5 * acos(costwoomega);
    double wp = asin(sin(w)/nr);
    double a1 = sin(w - wp);
    double b1 = sin(w+wp);
    double c1 = tan(w-wp);
    double d1 = tan(w+wp);
    double R = 0.5*((a1*a1)/(b1*b1)+(c1*c1)/(d1*d1));
//    double rgl = PI * p * R/(4*cos(solz_*cc)*cos(senz_*cc)*pow(cosbeta,4.0));
    float glint_refl = p*R/(4*cos(senz_*cc)*pow(cosbeta,4.0));

    return glint_refl;
}

/**************************************************************************
 * NAME: calc_turbid_residual()
 *
 * DESCRIPTION: Calculate turbid water residual.
 *
 *************************************************************************/

float DbAlgOcean::calc_turbid_residual(float sza, float r488,
        float r1240, float r1600, float r2250, float r550, int& status )
{

    float x[4] = {0.486, 1.238, 1.601, 2.257};
    float y[4] = {r488, r1240, r1600, r2250};
    for (int i=0; i<4; i++) {
        x[i] = log10(x[i]);
        y[i] = log10(y[i]);
    }
    float r[2];
    status = linfit(4, x, y, r);
    if (status != 0) {
        std::cerr << "DbAlgOcean::calc_turbid_residual:: Linfit error" << std::endl;
        return DFILL_FLOAT;
    }
    float est550 = r[0] + log10(0.551)*r[1];
    est550 = pow(10.0, est550);
//   -- calculate residual and convert to standard reflectance units via cos(sza)/pi.
    float res   = r550 - est550;

    return res;
}

/**************************************************************************
 * NAME: linfit()
 *
 * DESCRIPTION: Compute linear fit for data in x and y arrays.
 *
 *************************************************************************/

int DbAlgOcean::linfit(int size, float x[], float y[], float r[])
{
    float sx  = 0.0;
    float sy  = 0.0;
    float sxy = 0.0;
    float sxx = 0.0;
    float syy = 0.0;
    for (int i=0; i<size; i++) {
      sx = sx + x[i];
      sy = sy + y[i];
      sxx = sxx + (x[i] * x[i]);
      sxy = sxy + (x[i] * y[i]);
      syy = syy + (y[i] * y[i]);
    }

    r[1] = ((size*sxy) - (sx*sy))/((size*sxx)-(sx*sx));
    r[0] = (sy/size)-(r[1]*sx/size);

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: set_fill_out()
 *
 * DESCRIPTION: Set output fill values.
 *
 *************************************************************************/

int DbAlgOcean::set_fill_out()
{
    int status = DTDB_SUCCESS;

    for (int i=0; i<NDBMDL+1; i++) {
        oOut_[i].aot550 = DFILL_FLOAT;
        oOut_[i].fmf = DFILL_FLOAT;
        oOut_[i].ae = DFILL_FLOAT;
        oOut_[i].ss = DFILL_FLOAT;
        oOut_[i].chl = DFILL_FLOAT;
        oOut_[i].alg_flag = DFILL_SHORT;
        oOut_[i].model_flag = DFILL_SHORT;
        for (int j=0; j<NOWL; j++) {
            oOut_[i].aot[j] = DFILL_FLOAT;
        }
    }
    mask_cm_ = DFILL_SHORT;

    return status;
}
