/*******************************************************************************
 *
 * NAME: DtAlgOcean.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute processings for a given DDProcess object class.
 *
 *  Created on: November 3, 2016
 *      Author: Sam Anderson, DT
 *
 *  Modified:
 *
 *******************************************************************************/

#include <math.h>
#include <cstring>
#include <fstream>

#include <DDProcess.h>
#include "darktarget/DtLutNetcdf.h"
#include "darktarget/DtMask.h"
#include "darktarget/DtAlgorithm.h"

using namespace std;

/*******************************************************************************
*
*  DESCRIPTION:
*  The MODIS ocean aerosol product consists of aerosol optical thickness
*  and size parameters estimates derived on a 10x10 (1-km) pixel spatial
*  array.  The measured radiances in a wide spectral range (0.47-2.13
*  microns) are inverted assuming a bi-lognormal size distribution.
*  The volume and the mean particle size for each log-normal mode are
*  determined.  When fully developed, the aerosol ocean algorithm with
*  use the seven MODIS bands: 1, 2, 3, 4, 5, 6, and 7.
*
*  INPUT PARAMETERS:   NONE
*
*  OUTPUT PARAMETERS:
*  SDS_ref_STD     Standard deviation of reflectances at 7 bands
*  SDSTAU_best     Optical thickness for best solution
*  SDSTAUS_best    Optical thickness contribution small particles for best solution
*  SDSTAUB_best    Optical thickness contribution large particles for best solution
*  SDSTAU_average  Optical thickness for best solution
*  SDSTAUS_average Optical thickness contribution small particles for best solution
*  SDSTAUB_average Optical thickness contribution large particles for best solution
*  SDS_Least_error         Least square error estimated
*  SDS_small_weighting     Small mode weighting factor
*  SDS_sol_INDX_small      Solution Number index small particles
*  SDS_sol_INDX_large      Solution Number index large particles
*  SDSASSY_best      Asymmetry_Factor at 7 bands for best solution
*  SDSASSY_average   Asymmetry_Factor at 7 bands for average solution
*  SDSBAcK_best      Backscattering ratio at 7 bands of best solution
*  SDSBAcK_average   Backscattering ratio at 7 bands of average solution
*  SDS_effrad         Effective_Radius at 0.55 micron of both solutions
*  SDS_AOT_model Ratio of optical depth of small mode vs effective optical depth at 550
*  SDS_RefF_best     Normalized reflected_flux at 7 bands of best solution
*  SDS_RefF_average  Normalized reflected_flux at 7 bands of average solution
*  SDS_TranF_best    Normalized Transmitted_flux at 7 bands of best solution
*  SDS_TranF_average Normalized Transmitted_flux at 7 bands of average solution
*  SDS_SCAT_ANGLE_OCEAN Scattering angle ocean
*  SDS_QCONTROL         Quality control SDS array
*  SDS_NUMPIXELS        Number of Pixels used for 0.55 micron
*  SDS_ccn              cloud_Fraction in percentage
*  SDS_mass_conc      Mass concentration
*  SDS_angs_coeff1      Angstrom Exponent for 0.550 and 0.865 miron
*  SDS_angs_coeff2      Angstrom exponent for 0.865 and 2.130 micron
*
*  REVISION HISTORY:
*  $Log: Process_ocean_V2.f,v $
*  Revision 2.1  1996/05/28  15:34:05  vlin
*  Updated to write out one scancube at a time
*
*  TEAM-UNIQUE HEADER:
*
*  *REFERENcES AND CREDITS:
*
*  WRITTEN BY:
*  SHANA MATTOO                E-MAIL:mattoo@climate.gsfc.nasa.gov
*  APPLIED RESEARCH CORPORATION          PHONE:  (301) 614-6214
*  NASA/GODDARD SPAcE FLIGHT cENTER      FAX:    (301) 614-6307
*  CODE 913
*  GREENBELT, MD 20771
*
*  DESIGN NOTES:
*
*  Internals Variables:
*
*  Reflectances of Wavelengths used modis,ch1-ch7  identified by REFW*
*
*      REFW670
*      REFW865
*      REFW488
*      REFW550
*      REFW1240
*      REFW1610
*      REFW2250
*      CLDMSK     cloud mask 0=cloudy,1=clear
*      MTHET0     Measured solar zenith angle in deg.
*      MTHET      Measured viewangle from ground in deg.
*      MPHI       Measured azimuth  in deg.
*      NWAV       Number of wavelengths.
*      NAOT       Number of optical thicknesses.
*      NTH0       Number of solar zenith angles.
*      NTHET      Number of view angles.
*      NPHI       Number of azimuth angles.
*      PHc,JPHI   Azimuth angles.
*      THET       View angles.
*      THET0      Solar zenith angles.
*      TVALUE     Transmission factor.
*      AINT       Reflectance from look-up table.
*      TAUA       optical thicknesses.
*      Numdata    Number of input data sets.
*
*************************************************************************/

/**************************************************************************
 * NAME: DtAlgOcean()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtAlgOcean::DtAlgOcean()
{
	cm_ = nullptr;
    memset(refl_rayl_, 0, NWAV*sizeof(float));
    memset(refl_big_, 0, NWAV*NUM_CASES_BIG*NAOT*sizeof(float));
    memset(refl_small_, 0, NWAV*NUM_CASES_SMALL*NAOT*sizeof(float));
    memset(albedo_R_small_tau_, 0, NWAV*NUM_CASES_SMALL*NAOT*sizeof(float));
    memset(albedo_R_big_tau_, 0, NWAV*NUM_CASES_BIG*NAOT*sizeof(float));
    memset(albedo_T_small_tau_, 0, NWAV*NUM_CASES_SMALL*NAOT*sizeof(float));
    memset(albedo_T_big_tau_, 0, NWAV*NUM_CASES_BIG*NAOT*sizeof(float));
    memset(refl_, 0, NWAV*sizeof(float));
    memset(sdev_, 0, NWAV*sizeof(float));
    memset(numData_, 0, NWAV*sizeof(int));
    memset(good_pixels_, 0, NWAV*sizeof(int));
    memset(tau_, 0, NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(tau_small_, 0, NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(tau_big_, 0, NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(backscatter_, 0, NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(assym_, 0, NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(refl_flux_, 0, NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(trans_flux_, 0, NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(angstrom_exp_, 0, NUM_STATS*NUMCASES*NUMCASEB*sizeof(float));
    memset(solution_index_, 0, NUM_SIZES*NUMCASES*NUMCASEB*sizeof(float));
    memset(ccn_, 0, NUMCASES*NUMCASEB*sizeof(float));
    memset(eff_radius_, 0, NUMCASES*NUMCASEB*sizeof(float));
    memset(eff_variance_, 0, NUMCASES*NUMCASEB*sizeof(float));
    memset(mass_con_ocean_, 0, NUMCASES*NUMCASEB*sizeof(float));
    memset(xmin_, 0, 2*NUMCASES*NUMCASEB*sizeof(float));
    memset(funmin_, 0, 2*NUMCASES*NUMCASEB*sizeof(float));
    memset(tau_X55_, 0, 2*NUMCASES*NUMCASEB*sizeof(float));
    memset(error_, 0, 2*NWAV*NUMCASES*NUMCASEB*sizeof(float));
    memset(refl_avg_, 0, NWAV*sizeof(float));
    memset(sdev_avg_, 0, NWAV*sizeof(float));
    memset(tau_avg_, 0, NUM_STATS*NWAV*sizeof(float));
    memset(tau_small_avg_, 0, NUM_STATS*NWAV*sizeof(float));
    memset(tau_big_avg_, 0, NUM_STATS*NWAV*sizeof(float));
    memset(backscatter_avg_, 0, NUM_STATS*NWAV*sizeof(float));
    memset(assym_avg_, 0, NUM_STATS*NWAV*sizeof(float));
    memset(refl_flux_avg_, 0, NUM_STATS*NWAV*sizeof(float));
    memset(trans_flux_avg_, 0, NUM_STATS*NWAV*sizeof(float));
    memset(angstrom_exp_avg_, 0, NUM_STATS*NUM_SIZES*sizeof(float));
    memset(solution_index_avg_, 0, NUM_STATS*NUM_SIZES*sizeof(float));
    memset(ccn_avg_, 0, NUM_STATS*sizeof(float));
    memset(eff_radius_avg_, 0, NUM_STATS*sizeof(float));
    memset(eff_variance_avg_, 0, NUM_STATS*sizeof(float));
    memset(mass_con_ocean_avg_, 0, NUM_STATS*sizeof(float));
    memset(xmin_avg_, 0, NUM_STATS*sizeof(float));
    memset(funmin_avg_, 0, NUM_STATS*sizeof(float));
    memset(tau_X55_avg_, 0, NUM_STATS*sizeof(float));
    qcontrol_= 0;
    qcontrol_exclude_= 0;
    qcontrol_special_= 0;
    qcontrol_cirrus_= 0;
    quality_dust_flag_glint_= 0;
    quality_dust_flag_off_glint_= 0;
    memset(quality_to_pass_, 0, 2*sizeof(short));
    memset(quality_flag_, 0, 12*sizeof(short));
    memset(sds_refl_, 0, NWAV*sizeof(float));
    memset(sds_refl_sdev_, 0, NWAV*sizeof(float));
    memset(sds_numPixels_, 0, NWAV*sizeof(short));
    memset(sds_tau_best_, 0, NWAV*sizeof(float));
    memset(sds_tau_avg_, 0, NWAV*sizeof(float));
    memset(sds_tau_big_best_, 0, NWAV*sizeof(float));
    memset(sds_tau_big_avg_, 0, NWAV*sizeof(float));
    memset(sds_tau_small_best_, 0, NWAV*sizeof(float));
    memset(sds_tau_small_avg_, 0, NWAV*sizeof(float));
    memset(sds_assy_best_, 0, NWAV*sizeof(float));
    memset(sds_assy_avg_, 0, NWAV*sizeof(float));
    memset(sds_back_best_, 0, NWAV*sizeof(float));
    memset(sds_back_avg_, 0, NWAV*sizeof(float));
    memset(sds_reff_best_, 0, NWAV*sizeof(float));
    memset(sds_reff_avg_, 0, NWAV*sizeof(float));
    memset(sds_tranf_best_, 0, NWAV*sizeof(float));
    memset(sds_tranf_avg_, 0, NWAV*sizeof(float));
    sds_Tau_Land_Ocean_= 0;
    sds_Tau_Land_Ocean_img_= 0;
    memset(sds_Small_Weighting_, 0, NUM_SIZES*sizeof(float));
    memset(sds_Least_Error_, 0, NUM_SIZES*sizeof(float));
    memset(sds_tau_X55_, 0, NUM_SIZES*sizeof(float));
    memset(sds_EffRad_, 0, NUM_SIZES*sizeof(float));
    memset(sds_EffVar_, 0, NUM_SIZES*sizeof(float));
    memset(sds_Sol_Index_Small_, 0, NUM_SIZES*sizeof(short));
    memset(sds_Sol_Index_Large_, 0, NUM_SIZES*sizeof(short));
    memset(sds_Angs_Coeff1_, 0, NUM_SIZES*sizeof(float));
    memset(sds_Angs_Coeff2_, 0, NUM_SIZES*sizeof(float));
    memset(sds_Mass_Conc_, 0, NUM_SIZES*sizeof(float));
    memset(sds_CCN_, 0, NUM_SIZES*sizeof(float));
    memset(sds_AOT_model_, 0, (NUMCASES+NUMCASEB)*sizeof(float));
    sds_land_ocean_quality_= 0;
}

/**************************************************************************
 * NAME: ~DtAlgOcean()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DtAlgOcean::~DtAlgOcean()
{
    if (cm_ != 0) {
        delete cm_;
    }
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * processing operations.
 *
 *************************************************************************/

int DtAlgOcean::initialize( map<string, ddata*> imap )
{
	int status = DTDB_SUCCESS;

    qcontrol_special_ = 0;

    DtAlgorithm::initialize( imap );

	DtLutNetcdf* lutgen = new DtLutNetcdf();
	status = lutgen->read_ocean_aerosol_lut( lut_ );
	delete lutgen;

    cm_ = new DtCloudMaskOcean(this);

	return status;
}

/**************************************************************************
 * NAME: process()
 *
 * DESCRIPTION: Dark target ocean algorithm process.
 *
 *************************************************************************/

map<string, ddata*> DtAlgOcean::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap)
{
    get_inputs(start, count, imap);

    set_fill_out();
    if (rfl_[(size_t)rhot_band::W550]<0.0) {
    	l2_flags_ |= (unsigned int) flags::BOWTIEDEL;
    	return set_fills();
    }
    if (((solz_ > threshsolz_) || (solz_ < 0.0)) && bmasksolz_) {
    	l2_flags_ |= (unsigned int) flags::HISOLZEN;
    	return set_fills();
    }
    if (((senz_ < 0.0) || (senz_ > threshsenz_)) && bmasksenz_) {
    	l2_flags_ |= (unsigned int) flags::HISATZEN;
    	return set_fills();
    }
    if ((raa_ < 0.0) || (raa_ > MAXMPHI) || (ws_ < 0.0) || (ws_ > 99.0)) {
    	return set_fills();
    }
	rfld_[D412] = rfl_[(size_t)rhot_band::W410];
	rfld_[D488] = rfl_[(size_t)rhot_band::W490];
	rfld_[D550] = rfl_[(size_t)rhot_band::W550];
	rfld_[D670] = rfl_[(size_t)rhot_band::W670];
	rfld_[D865] = rfl_[(size_t)rhot_band::W865];
	rfld_[D1240] = rfl_[(size_t)rhot_band::W1240];
	rfld_[D1380] = rfl_[(size_t)rhot_band::W1380];
	rfld_[D1610] = rfl_[(size_t)rhot_band::W1610];
	rfld_[D2250] = rfl_[(size_t)rhot_band::W2250];
	if (bgascorrect_) {
		compute_gas_correction();
		rfld_[D412] *= gasc_[D412];
		rfld_[D488] *= gasc_[D488];
		rfld_[D550] *= gasc_[D550];
		rfld_[D670] *= gasc_[D670];
		rfld_[D865] *= gasc_[D865];
		rfld_[D1240] *= gasc_[D1240];
		rfld_[D1610] *= gasc_[D1610];
		rfld_[D2250] *= gasc_[D2250];
	}

	if (cloud_mask_ == DFILL_UBYTE) {
		cm_->compute(cmask_);
		cloud_mask_ = (cmask_+1)%2;
	} else {
		cmask_ = (cloud_mask_+1)%2;
	}
	if (!bmaskcloud_) cmask_ = 1;
	if (bmaskcloud_ && cloud_mask_) {
		l2_flags_ |= (unsigned int) flags::CLDICE;
		return set_fills();
	}

	compute_glint_refl();
	index_geometry( solz_, senz_, raa_ );
	interpolate_rayleigh();
	compute_avg_refl();
	store_reflectance();

	if((!bmaskglint_) ||
	   ((refl_[D550] > 0) && (refl_[D670] > 0) &&
	   (refl_[D865] > 0) && (refl_[D2250] > 0) &&
	   (qcontrol_special_== 0 || qcontrol_special_== 2))) {

		interpolate_angles();

		int iSol = 0;
		for (int iSmall=0; iSmall<NUMCASES; iSmall++) {
			for (int iBig=0; iBig<NUMCASEB; iBig++) {

//                  compute_minimum_baseline(iBig, iSmall, iSol);
				compute_minimum(iBig, iSmall, iSol);
				compute_tau_flux(iBig, iSmall, iSol);
				solution_index_[ISMALL][iSol] = iSmall;
				solution_index_[IBIG][iSol] = iBig;
				iSol++;
			}
		}

		average_output();
		store_output();

		qcontrol_ = (tau_avg_[1][0] > MAXTAU) ? -29 : qcontrol_;
		qcontrol_ = (tau_avg_[1][0] < 0) ? 9 : qcontrol_;
		qcontrol_ = (tau_avg_[1][0] <= MINTAU) ? -28 : qcontrol_;
	} else {
		qcontrol_ = -21;
		return set_fills();
	}

    assign_quality();

    // write to generic outputs

    qual_flag_ = quality_to_pass_[0];
	aerosol_type_ = (NUMCASEB-1)*sds_Sol_Index_Large_[1] +
    		sds_Sol_Index_Small_[1];
    error_flag_ = abs(quality_to_pass_[1]);

    scatter_ang_ = 180.0 - glint_angle_;
    glint_ang_ = glint_angle_;
    sse_ = sds_Least_Error_[1];
    fmf_ = sds_Small_Weighting_[1];
    aot_550_ = sds_tau_X55_[1];
    ae1_ = sds_Angs_Coeff1_[1];
    ae2_ = sds_Angs_Coeff2_[1];
    ndv_ = DFILL_FLOAT;
    chlor_ = DFILL_FLOAT;
    aot_[0] = DFILL_FLOAT;
    for ( int iWav=0; iWav<NOWL; iWav++ ) {
    	aot_[iWav+1] = sds_tau_avg_[iWav];
    }
    for ( int iWav=0; iWav<NLWL+1; iWav++ ) {
        sr_[iWav] = DFILL_FLOAT;
        ssa_[iWav] = DFILL_FLOAT;
    }

    return set_outputs();
}

/**************************************************************************
 * NAME: index_geometry()
 *
 * DESCRIPTION: Reset LUT indices based on measured geometry.
 *
 *************************************************************************/

int DtAlgOcean::index_geometry( float sza, float azim, float phi)
{
	int status = DTDB_SUCCESS;

	float del = (sza <= lut_.THET0[4]) ? 12.0 : 6.0;
	SZA_0_ = (sza <= lut_.THET0[4]) ? (int)(sza/del) : (int)(sza/del)-4;
	SZA_1_ = SZA_0_+1;
	THE_0_ = (int)(azim/6.0);
	THE_1_ = THE_0_+1;
	PHI_0_ = (int)(phi/12.0);
	PHI_1_ = PHI_0_+1;

	return status;
}

/**************************************************************************
 * NAME: interpolate_rayleigh()
 *
 * DESCRIPTION: Subroutine interpolate_rayleigh interpolates the lookup
 *              Rayleigh reflectances to the measured geometry.
 *
 *************************************************************************/

int DtAlgOcean::interpolate_rayleigh()
{
    int status = DTDB_SUCCESS;

    int   lutindex = 0;
    while ((wind_[++lutindex] < ws_) &&
           (lutindex < (WIND_LUT_ENTRIES-1))) {}
    short iws1 = lutindex-1;
    short iws2 = lutindex;
    float ws_lut[WIND_LUT_ENTRIES];
    memset(ws_lut, 0, WIND_LUT_ENTRIES*sizeof(float));
    float phi_out[2] = {0.0,0.0};
    float the_out[2] = {0.0,0.0};
    float sza_out[2] = {0.0,0.0};
    for (int iWav=0; iWav<NWAV; iWav++) {
        int numWS = 0;
        for (int iWS=iws1; iWS <= iws2; iWS++ ) {
            ws_lut[numWS] = wind_[iWS];
            int numSza = 0;
            for (int iSza=SZA_0_; iSza<=SZA_1_; iSza++) {
                int numThe = 0;
                for (int iThe=THE_0_; iThe<=THE_1_; iThe++) {
                    int numPhi = 0;
                    for (int iPhi=PHI_0_; iPhi<=PHI_1_; iPhi++) {
                        numPhi++;
                    }
                    interp_extrap(numPhi, raa_, &lut_.PHC[PHI_0_],
                            &lut_.REF_RAYALL[iWav][iWS][iSza][iThe][PHI_0_],
                            phi_out[numThe]);
                    numThe++;
                }
                interp_extrap(numThe, senz_, &lut_.THET[THE_0_],
                        phi_out, the_out[numSza]);
                numSza++;
            }
            interp_extrap(numSza, solz_, &lut_.THET0[SZA_0_],
                    the_out, sza_out[numWS]);
            numWS++;
        }
        interp_extrap(numWS, ws_, ws_lut, sza_out, refl_rayl_[iWav]);

        refl_rayl_[iWav] *= M_PI/cos(DEGtoRAD*solz_);
    }

    return status;
}

/**************************************************************************
 * NAME: interpolate_angles()
 *
 * DESCRIPTION: Subroutine interpolates the lookup reflectances to the
 *              measured geometry.
 *
 *************************************************************************/

int DtAlgOcean::interpolate_angles()
{
    int status = DTDB_SUCCESS;

    int   lutindex = 0;
    while ((wind_[++lutindex] < ws_) &&
           (lutindex < (WIND_LUT_ENTRIES-1))) {}
    short iws1 = lutindex-1;
    short iws2 = lutindex;
    float ws_lut[WIND_LUT_ENTRIES];
    memset(ws_lut, 0, WIND_LUT_ENTRIES*sizeof(float));

    float phi_out[2][2] = {{0.0,0.0},{0.0,0.0}};
    float the_out[2][2] = {{0.0,0.0},{0.0,0.0}};
    float sza_out[2][2] = {{0.0,0.0},{0.0,0.0}};
    float ws_out = 0.0;

    for (int iWav=0; iWav<NWAV; iWav++) {
        //  INTERPOLATE FOR SMALL CASES
        for (int iCase=0; iCase<NUMCASES; iCase++) {
            for (int iTau=0; iTau<NAOT; iTau++) {
                int numWS = 0;
                for (int iWS=iws1; iWS <= iws2; iWS++) {
                    ws_lut[numWS] = wind_[iWS];
                    int numSza = 0;
                    for (int iSza=SZA_0_; iSza<=SZA_1_; iSza++) {
                        numSza++;
                    }
                    interp_extrap(numSza, solz_, &lut_.THET0[SZA_0_],
                        &lut_.ALBEDO_R_SMALL[iWav][iCase][iTau][iWS][SZA_0_],
                        sza_out[0][numWS]);
                    interp_extrap(numSza, solz_, &lut_.THET0[SZA_0_],
                        &lut_.ALBEDO_T_SMALL[iWav][iCase][iTau][iWS][SZA_0_],
                        sza_out[1][numWS]);
                    numWS++;
                }
                interp_extrap(numWS, ws_, ws_lut,
                        sza_out[0], albedo_R_small_tau_[iWav][iCase][iTau]);
                interp_extrap(numWS, ws_, ws_lut,
                        sza_out[1], albedo_T_small_tau_[iWav][iCase][iTau]);
            }
        }
        //  INTERPOLATE FOR BIG CASES
        for (int iCase=0; iCase<NUMCASEB; iCase++) {
            for (int iTau=0; iTau<NAOT; iTau++) {
                int numWS = 0;
                for (int iWS=iws1; iWS <= iws2; iWS++) {
                    ws_lut[numWS] = wind_[iWS];
                    int numSza = 0;
                    for (int iSza=SZA_0_; iSza<=SZA_1_; iSza++) {
                        numSza++;
                    }
                    interp_extrap(numSza, solz_, &lut_.THET0[SZA_0_],
                        &lut_.ALBEDO_R_BIG[iWav][iCase][iTau][iWS][SZA_0_],
                        sza_out[0][numWS]);
                    interp_extrap(numSza, solz_, &lut_.THET0[SZA_0_],
                        &lut_.ALBEDO_T_BIG[iWav][iCase][iTau][iWS][SZA_0_],
                        sza_out[1][numWS]);
                    numWS++;
                }
                interp_extrap(numWS, ws_, ws_lut,
                        sza_out[0], albedo_R_big_tau_[iWav][iCase][iTau]);
                interp_extrap(numWS, ws_, ws_lut,
                        sza_out[1], albedo_T_big_tau_[iWav][iCase][iTau]);
            }
        }
    }

    for (int iWav = 0; iWav < NWAV; iWav++) {
        for (int iCase = 0; iCase < NUMCASES; iCase++) {
            for (int iTau = 1; iTau < NAOT; iTau++) {
                int numWS = 0;
                for (int iWS = iws1; iWS <= iws2; iWS++) {
                    ws_lut[numWS] = wind_[iWS];
                    int numSza = 0;
                    for (int iSza = SZA_0_; iSza <= SZA_1_; iSza++) {
                        int numThe = 0;
                        for (int iThe = THE_0_; iThe <= THE_1_; iThe++) {
                            int numPhi = 0;
                            for (int iPhi = PHI_0_; iPhi <= PHI_1_; iPhi++) {
                                numPhi++;
                            }
                            interp_extrap(numPhi, raa_, &lut_.PHC[PHI_0_],
                                    &lut_.AINTS[iWav][iCase][iTau][iWS][iSza][iThe][PHI_0_],
                                    phi_out[ISMALL][numThe]);
                            numThe++;
                        }
                        interp_extrap(numThe, senz_, &lut_.THET[THE_0_],
                                phi_out[ISMALL], the_out[ISMALL][numSza]);
                        numSza++;
                    }
                    interp_extrap(numSza, solz_, &lut_.THET0[SZA_0_],
                            the_out[ISMALL], sza_out[ISMALL][numWS]);
                    numWS++;
                }
                interp_extrap(numWS, ws_, ws_lut, sza_out[ISMALL],
                        ws_out);
                refl_small_[iWav][iCase][iTau] = ws_out * M_PI
                        / cos( solz_*DEGtoRAD);
             }
            refl_small_[iWav][iCase][0] = refl_rayl_[iWav];
        }
        for (int iCase = 0; iCase < NUMCASEB; iCase++) {
            for (int iTau = 1; iTau < NAOT; iTau++) {
                int numWS = 0;
                for (int iWS = iws1; iWS <= iws2; iWS++) {
                    ws_lut[numWS] = wind_[iWS];
                    int numSza = 0;
                    for (int iSza = SZA_0_; iSza <= SZA_1_; iSza++) {
                        int numThe = 0;
                        for (int iThe = THE_0_; iThe <= THE_1_; iThe++) {
                            int numPhi = 0;
                            for (int iPhi = PHI_0_; iPhi <= PHI_1_; iPhi++) {
                                numPhi++;
                            }
                            interp_extrap(numPhi, raa_, &lut_.PHC[PHI_0_],
                            &lut_.AINTB[iWav][iCase][iTau][iWS][iSza][iThe][PHI_0_],
                            phi_out[IBIG][numThe]);
                            numThe++;
                        }
                        interp_extrap(numThe, senz_, &lut_.THET[THE_0_],
                                phi_out[IBIG], the_out[IBIG][numSza]);
                        numSza++;
                    }
                    interp_extrap(numSza, solz_, &lut_.THET0[SZA_0_],
                            the_out[IBIG], sza_out[IBIG][numWS]);
                    numWS++;
                }
                interp_extrap(numWS, ws_, ws_lut, sza_out[IBIG],
                        ws_out);
                refl_big_[iWav][iCase][iTau] = ws_out * M_PI
                        / cos( solz_*DEGtoRAD);
            }
            refl_big_[iWav][iCase][0] = refl_rayl_[iWav];
        }
    }

    return status;
}

/**************************************************************************
 * NAME: compute_tau_flux()
 *
 * DESCRIPTION: This subroutine computes CCN, assymetry factor,
 *              backscattering ratio, effective radius, effective
 *              variance, and optical thicknesses for large and small
 *              mode.
 *
 *************************************************************************/

int DtAlgOcean::compute_tau_flux( int iBig, int iSmall, int iSol )
{
	int status = DTDB_SUCCESS;

    int   lutindex = 0;
    while ((wind_[++lutindex] < ws_) &&
           (lutindex < (WIND_LUT_ENTRIES-1))) {}
    short iws1 = lutindex-1;
    short iws2 = lutindex;
    int   numWS = 0;

//  Since we see the glint effect around tau of 0.5  no retrievals will be done
//  less than tau of 0.5 in glint area. It is set to negative value so that later on it is
//  filled with fill values
    if((tau_X55_[0][iSol] <= 0.7) && (quality_dust_flag_glint_ == 1)) {
    	tau_X55_[0][iSol] = -0.02;
    }
//  If tau at 0.55um is greater than maxtau(ie 5.00) taux55 is set to maxtau for
//  glint and off glint heavy dust episodes
    if(((tau_X55_[0][iSol] > MAXTAU) && (quality_dust_flag_glint_== 1)) ||
    		((tau_X55_[0][iSol]  > MAXTAU) &&
    		(quality_dust_flag_off_glint_== 1))) {
        tau_X55_[0][iSol]= MAXTAU;
    }
//  Compute normalized extinction coefficient
    float tau_X55_small = tau_X55_[0][iSol]*xmin_[0][iSol];
    float tau_X55_big = tau_X55_[0][iSol]- tau_X55_small;

    float ws_lut[WIND_LUT_ENTRIES];
    float rsmall[WIND_LUT_ENTRIES];
    float rbig[WIND_LUT_ENTRIES];
    memset(ws_lut, 0, WIND_LUT_ENTRIES*sizeof(float));
    memset(rsmall, 0, WIND_LUT_ENTRIES*sizeof(float));
    memset(rbig, 0, WIND_LUT_ENTRIES*sizeof(float));

    for (int iWav=0; iWav<NWAV; iWav++) {
        numWS = 0;
        for (int iWS=iws1; iWS <= iws2; iWS++) {
            ws_lut[numWS] = wind_[iWS];
            rsmall[numWS] =
              lut_.EXTSMALL[iWav][iSmall][iWS]/lut_.EXTSMALL[D550][iSmall][iWS];
            rbig[numWS] =
              lut_.EXTBIG[iWav][iBig][iWS]/lut_.EXTBIG[D550][iBig][iWS];
            numWS++;
        }
        float extsmall = 0.0;
        interp_extrap(numWS,ws_,ws_lut,rsmall,extsmall);
        float extbig = 0.0;
        interp_extrap(numWS,ws_,ws_lut,rbig,extbig);

        tau_[iWav][iSol] = tau_X55_small*extsmall + tau_X55_big*extbig;
        tau_small_[iWav][iSol] = tau_X55_small*extsmall;
        tau_big_[iWav][iSol] = tau_X55_big*extbig;
    }

    if ((tau_[D550][iSol] > 0) && (tau_[D865][iSol] > 0)) {
        float tlog = log(tau_[D550][iSol]/tau_[D865][iSol]) /
                log(lut_.WAVE[D865]/lut_.WAVE[D550]);
        angstrom_exp_[0][iSol] = tlog;
    }
    if ((tau_[D865][iSol] > 0) && (tau_[D2250][iSol] > 0)) {
        float tlog = log(tau_[D865][iSol]/tau_[D2250][iSol]) /
                log(lut_.WAVE[D2250]/lut_.WAVE[D865]);
        angstrom_exp_[1][iSol] = tlog;
    }

//  Interpolate for large and small cases for total optical depth at 0.55
    float albedo_R_small[NWAV][NUM_CASES_SMALL];
    float albedo_R_big[NWAV][NUM_CASES_BIG];
    float albedo_T_small[NWAV][NUM_CASES_SMALL];
    float albedo_T_big[NWAV][NUM_CASES_BIG];
    memset(albedo_R_small, 0, NWAV*NUM_CASES_SMALL*sizeof(float));
    memset(albedo_R_big, 0, NWAV*NUM_CASES_BIG*sizeof(float));
    memset(albedo_T_small, 0, NWAV*NUM_CASES_SMALL*sizeof(float));
    memset(albedo_T_big, 0, NWAV*NUM_CASES_BIG*sizeof(float));

    for (int iWav=0; iWav<NWAV; iWav++) {
        interp_extrap(NAOT,tau_[D550][iSol],lut_.TAUAS[D550][iSmall],
                albedo_R_small_tau_[iWav][iSmall],albedo_R_small[iWav][iSmall]);
        interp_extrap(NAOT,tau_[D550][iSol],lut_.TAUAS[D550][iSmall],
                albedo_T_small_tau_[iWav][iSmall],albedo_T_small[iWav][iSmall]);
        interp_extrap(NAOT,tau_[D550][iSol],lut_.TAUAB[D550][iBig],
                albedo_R_big_tau_[iWav][iBig],albedo_R_big[iWav][iBig]);
        interp_extrap(NAOT,tau_[D550][iSol],lut_.TAUAB[D550][iBig],
                albedo_T_big_tau_[iWav][iBig],albedo_T_big[iWav][iBig]);
    }

//  Compute CCN, asymmetry factor, backscattering ratio,
//  effective radius, effective variance, and optical thicknesses for
//  large and small mode.
    float wout = 0.0;
    interp_extrap(numWS, ws_, &ws_lut[iws1],
            &lut_.EXTSMALL[D550][iSmall][iws1], wout);
    float nSmall = tau_X55_small/wout;
    interp_extrap(numWS, ws_, &ws_lut[iws1],
            &lut_.EXTBIG[D550][iBig][iws1], wout);
    float nBig = tau_X55_big/wout;
    interp_extrap(numWS, ws_, &ws_lut[iws1],
            &lut_.CCNSMALL[iSmall][iws1], wout);
    ccn_[iSol] = nSmall* wout;

    eff_radius_[iSol] = 0.0;
    eff_variance_[iSol] = 0.0;
    mass_con_ocean_[iSol] = 0.0;

    if (nSmall > 0.0 ||nBig > 0.0 ) {
        float m[NUM_SIZES][NUM_MOMENTS];
        memset(m, 0, NUM_SIZES*NUM_MOMENTS*sizeof(float));
        for ( int iM=1; iM < NUM_MOMENTS; iM++) {
            interp_extrap(numWS, ws_, &ws_lut[iws1],
                    &lut_.MOMENTSSMALL[iSmall][iM][iws1], m[ISMALL][iM]);
            interp_extrap(numWS, ws_, &ws_lut[iws1],
                    &lut_.MOMENTSBIG[iBig][iM][iws1], m[IBIG][iM]);
        }
        eff_radius_[iSol] = (nSmall*m[ISMALL][2]+ nBig*m[IBIG][2]) /
                (nSmall*m[ISMALL][1]+ nBig*m[IBIG][1]);
        mass_con_ocean_[iSol] = nSmall*m[ISMALL][2]+ nBig*m[IBIG][2];
        eff_variance_[iSol] = (nSmall*m[ISMALL][3]+ nBig*m[IBIG][3])/
                (pow(eff_radius_[iSol],2)*(nSmall*m[ISMALL][1]+nBig*4))-1.0;
    }

    for (int iWav=0; iWav<NWAV; iWav++) {
        float al[NUM_SIZES] = {0.0,0.0};
        float bs[NUM_SIZES] = {0.0,0.0};
        float as[NUM_SIZES] = {0.0,0.0};
        interp_extrap(numWS, ws_, &ws_lut[iws1],
                &lut_.ALBEDOSMALL[iWav][iSmall][iws1], al[ISMALL]);
        interp_extrap(numWS, ws_, &ws_lut[iws1],
                &lut_.ALBEDOBIG[iWav][iBig][iws1], al[IBIG]);
        interp_extrap(numWS, ws_, &ws_lut[iws1],
                &lut_.BACKSCTTSMALL[iWav][iSmall][iws1], bs[ISMALL]);
        interp_extrap(numWS, ws_, &ws_lut[iws1],
                &lut_.BACKSCTTBIG[iWav][iBig][iws1], bs[IBIG]);
        interp_extrap(numWS, ws_, &ws_lut[iws1],
                &lut_.ASSYMSMALL[iWav][iSmall][iws1], as[ISMALL]);
        interp_extrap(numWS, ws_, &ws_lut[iws1],
                &lut_.ASSYMBIG[iWav][iBig][iws1], as[IBIG]);
//  compute backscattering ratio
        if ((tau_small_[iWav][iSol] > 0.0) ||
                (tau_big_[iWav][iSol] > 0.0)) {
            backscatter_[iWav][iSol] =
                (al[ISMALL]*tau_small_[iWav][iSol]*bs[ISMALL] +
                 al[IBIG]*tau_big_[iWav][iSol]*bs[IBIG]) /
                (al[ISMALL]*tau_small_[iWav][iSol] +
                 al[IBIG]*tau_big_[iWav][iSol]);
            assym_[iWav][iSol] =
                (al[ISMALL]*tau_small_[iWav][iSol]*as[ISMALL] +
                 al[IBIG]*tau_big_[iWav][iSol]*as[IBIG]) /
                (al[ISMALL]*tau_small_[iWav][iSol] +
                 al[IBIG]*tau_big_[iWav][iSol]);
            refl_flux_[iWav][iSol] =
                (tau_small_[iWav][iSol]* albedo_R_small[iWav][iSmall]
                + tau_big_[iWav][iSol]*albedo_R_big[iWav][iBig])/
                (tau_small_[iWav][iSol] + tau_big_[iWav][iSol]);
// subtract rayleigh and apply angle correction because it is divided in LUT
            refl_flux_[iWav][iSol] =
                (refl_flux_[iWav][iSol] - albedo_R_small_tau_[iWav][0][0]) *
                    cos(solz_*DEGtoRAD);
// compute transmittance
            trans_flux_[iWav][iSol]=
                (tau_small_[iWav][iSol]*albedo_T_small[iWav][iSmall]
                + tau_big_[iWav][iSol]*albedo_T_big[iWav][iBig])/
                (tau_small_[iWav][iSol]+tau_big_[iWav][iSol]);
        }
    }

    return status;
}

/**************************************************************************
 * NAME: compute_avg_refl()
 *
 * DESCRIPTION: This subroutine averages the reflectances for each
 *              10*10 pixel  square and finds standard deviation for
 *              averaged Reflectances
 *
 *************************************************************************/

int DtAlgOcean::compute_avg_refl()
{
	int status = DTDB_SUCCESS;

    float var[NWAV];
    float dvar[NWAV];
    float refl_interm[NWAV][GRIDX*GRIDY];
    float array_interm[NWAV][GRIDX*GRIDY];
    memset(var, 0, NWAV*sizeof(float));
    memset(dvar, 0, NWAV*sizeof(float));
    memset(refl_interm, 0, NWAV*GRIDX*GRIDY*sizeof(float));
    memset(array_interm, 0, NWAV*GRIDX*GRIDY*sizeof(float));

    qcontrol_exclude_ = 0;
    qcontrol_special_ = 0;
    quality_dust_flag_glint_ = 0;
    quality_dust_flag_off_glint_ = 0;
    for (int iWav=0; iWav<NWAV; iWav++) {
    	good_pixels_[iWav] = 0;
    }

    qcontrol_cirrus_ = 0;

    for (int iWav = 0; iWav<NWAV; iWav++) {
        numData_[iWav]= 0;
        for (int ix = 0; ix < GRIDX*GRIDY; ix++) {
            refl_interm[iWav][ix] = DFILL_FLOAT;
        }
    }

    DtSedimentMask* sm = new DtSedimentMask(this);
    short smask = 1;

    for ( size_t iy=0; iy<GRIDY; iy++ ) {
        for ( size_t ix = 0; ix < GRIDX; ix++ ) {

            sm->compute(smask);

            if( rfld_[D1240] > 0) {
                if ((cmask_ == 0) &&
                    (rfld_[D670] > 1.50*refl_rayl_[D670]) &&
                    ((rfld_[D1380] / rfld_[D1240] > 0.10) &&
                    (rfld_[D1380] / rfld_[D1240] < 0.30)) &&
                    ((rfld_[D1380] > 0.01) && (rfld_[D1380] <= 0.03))) {
                    qcontrol_cirrus_ = 1;
                }
            }
            if (cmask_*smask > 0) {
                refl_interm[D488][numData_[D488]] = rfld_[D488];
                numData_[D488] += 1;
                refl_interm[D550][numData_[D550]] = rfld_[D550];
                numData_[D550] += 1;
                refl_interm[D670][numData_[D670]] = rfld_[D670];
                numData_[D670] += 1;
                refl_interm[D865][numData_[D865]] = rfld_[D865];
                numData_[D865] += 1;
                refl_interm[D1240][numData_[D1240]] = rfld_[D1240];
                numData_[D1240] += 1;
                refl_interm[D1610][numData_[D1610]] = rfld_[D1610];
                numData_[D1610] += 1;
                refl_interm[D2250][numData_[D2250]] = rfld_[D2250];
                numData_[D2250] += 1;
            }
        }
    }
    delete sm;

//  Sort in ascending order and get the index for wavelength at 0.865
	int index[GRIDX*GRIDY];

    sort_index(numData_[D865], &refl_interm[D865][0], index);

//  Reject 1/4 of brightest and darkest cloud-free pixels to eliminate
//  the residual shadows and sub_cloud pixel based on wavelength of 865 nm
	int n10 = numData_[D865]/4;
	int n40 = (numData_[D865]-n10);
//    if(n10 == 0)  n10 = 1;
//    if(n40 == 0)  n40 = 1;
    int q1_pixels = (int) LINE*LINE*0.025;
    int q2_pixels = (int) LINE*LINE*0.050;
    int q3_pixels = (int) LINE*LINE*0.075;
//  Throw away 1/4 of brightest & darkest pixels in all wavelengths
    if ((n40-n10) > q1_pixels) {
    	for (int iWav=0; iWav < NWAV; iWav++) {
    		for (int ix=n10; ix < n40; ix++) {
    	        if ((refl_interm[iWav][index[ix]] > 0.0) &&
    	        	(refl_interm[iWav][index[ix]]<=1.3)) {
    	        	array_interm[iWav][good_pixels_[iWav]] =
    	        	        refl_interm[iWav][index[ix]];
    	        	good_pixels_[iWav] += 1;
    	        }
    		}
    	}
		for (int iWav=0; iWav < NWAV; iWav++) {
			if (good_pixels_[iWav] > 0) {
				float refl_val = 0.0;
				float std_val = 0.0;
				mean_std(good_pixels_[iWav], array_interm[iWav],
				         refl_val, std_val);
				refl_[iWav] = refl_val;
				sdev_[iWav] = std_val;
				var[iWav] = std_val/refl_val;
			}
			else {
				refl_[iWav] = DFILL_FLOAT;
				sdev_[iWav] = DFILL_FLOAT;
				var[iWav] = DFILL_FLOAT;
			}
		}
		dvar[D670] = -(1.0/(log(lut_.WAVE[D670]/lut_.WAVE[D865])))*
						log((1.+ var[D670])/(1.+var[D865]));
		dvar[D1240] = -(1./(log(lut_.WAVE[D1240]/lut_.WAVE[D865])))*
						log((1.+ var[D1240])/(1.+var[D865]));
		dvar[D1610] = -(1./(log(lut_.WAVE[D1610]/lut_.WAVE[D865])))*
						log((1.+ var[D1610])/(1.+var[D865]));
		dvar[D2250] = -(1./(log(lut_.WAVE[D2250]/lut_.WAVE[D865])))*
						log((1.+ var[D2250])/(1.+var[D865]));

//  if GLINT_ANGLE <= glint threshold store the
//  reflectances standard deviation and number of pixels only.
		//  Apply glint mask
        if ((glint_angle_ > 0.75*GLINT_ANGLE_THRESHOLD) &&
                (glint_angle_ <= GLINT_ANGLE_THRESHOLD)) {
            qcontrol_ = 10;
        }
		if( glint_refl_< threshglint_  &&
				((refl_[D488] > 0  && refl_[D670]> 0) &&
				(refl_[D488]/refl_[D670] <= 0.75))) {
			 quality_dust_flag_off_glint_ = 1;
		}
		if( glint_refl_>= threshglint_  &&
				((refl_[D488] > 0  && refl_[D670] > 0) &&
				(refl_[D488]/refl_[D670] <= 0.95))) {
			 quality_dust_flag_glint_ = 1;
		}
		if( glint_refl_ < threshglint_  ||
				quality_dust_flag_off_glint_ == 1 ||
				quality_dust_flag_glint_ == 1) {
//   If there is valid data in 865 nm channel and at least one channel has valid data
			int total_good = good_pixels_[D550]+good_pixels_[D670]+
					good_pixels_[D1240]+good_pixels_[D1610]+good_pixels_[D2250];
			if((good_pixels_[D865] > q1_pixels) &&
				  (total_good > q3_pixels)) {
//   If  Reflectance at wave 865 nm is greater than rhomin1
//   the aerosol content is  enough to process
				if (refl_[D865] > (1.10*refl_rayl_[D865])) {
//   Quality control flags are set for the number of cloud free pixels used
//   qcontrol_=1 means that box is less than 10% cloud free.
					numData_[D550] = good_pixels_[D865];
					if ( good_pixels_[D865] < q2_pixels) {
						qcontrol_=1;
					}
					if ( good_pixels_[D865] > q2_pixels) {
						qcontrol_ = 0;
					}
//    If reflectance at wave 865 nm is less than rhomin2, then
//    the aerosol content is not enough to derive size-distribution, but
//    is enough to derive optical thickness ... process
					if ( refl_[D865] < (1.50*refl_rayl_[D865])) {
						qcontrol_= 3;
						qcontrol_special_= 2;
					}
//  Test if channel 1. 24  is good to process
					if (refl_[D1240]*cos(solz_*DEGtoRAD) < 3.600E-04) {
						qcontrol_exclude_= 6;
					}
//  Test if channel 2.13 and 1.64 um are good to process
					if (refl_[D1610]*cos(solz_*DEGtoRAD) < 3.600E-04) {
						qcontrol_exclude_= 3;
					}
					if (refl_[D2250]*cos(solz_*DEGtoRAD) < 3.110E-04) {
						qcontrol_exclude_= 4;
					}
					if (refl_[D2250]*cos(solz_*DEGtoRAD) < 3.110E-04 &&
							refl_[D1610]*cos(solz_*DEGtoRAD)  < 3.600E-04) {
						qcontrol_exclude_= 5;
					}
//  Process only for valid data
					qcontrol_ = qcontrol_exclude_;
					if ( qcontrol_ > 0) {
//  The aerosol type and content are variable
						if (var[D865] > 0.05 &&
						(fabs(dvar[D670]) > 0.15 ||
						 fabs(dvar[D1610]) > 0.15 ||
						 fabs(dvar[D2250]) > 0.15)) {
							qcontrol_= 6;
						}
//  The aerosol content are variable not the aerosol content
						if (var[D865] > 0.05 &&
						(fabs(dvar[D670]) < 0.15 ||
						 fabs(dvar[D1610]) < 0.15 ||
						 fabs(dvar[D2250]) < 0.15)) {
							qcontrol_= 7;
						}
					}
				}
				else {
//   If  Reflectance at wavelength 865 nm is less than rhomin1
//   the aerosol content is  not enough to process
					qcontrol_= 17;
					qcontrol_special_= 1;
				}
			}
			else {
//  If there is no valid data in all channels
				qcontrol_ = -20;
				qcontrol_special_ = -2;
			}
		}
		else {
//  Glint and only Reflectances, sdev and number of pixels will be stored
			qcontrol_special_= 3;
			numData_[D550] = good_pixels_[D865];
			qcontrol_= 11;
	    	l2_flags_ |= (unsigned int) flags::HIGLINT;
		}
    }
	else {
// else for number of cloud-free pixels in box is too cloudy
// If cloud free pixels are less than 10 no processing is performed .
		qcontrol_= -22;
		qcontrol_special_= -2;
	}
    if (qcontrol_ < 0) {
         refl_[D488] = DFILL_FLOAT;
         refl_[D550] = DFILL_FLOAT;
         refl_[D670] = DFILL_FLOAT;
         refl_[D865] = DFILL_FLOAT;
         refl_[D1240] = DFILL_FLOAT;
         refl_[D1610] = DFILL_FLOAT;
         refl_[D2250] = DFILL_FLOAT;
    }

    return status;
}

/**************************************************************************
 * NAME: store_results()
 *
 * This subroutine stores all the output variables to be written to output file.
 *
 ***************************************************************************/

int DtAlgOcean::store_reflectance()
{
	int status = DTDB_SUCCESS;

	float maxrefl = 1.2;
	float minrefl = 0.0;
	for (int iWav=0; iWav < NWAV; iWav++) {
		if (refl_[iWav] >= minrefl) {
		    if (refl_[iWav] <= maxrefl) {
                sds_refl_[iWav] = refl_[iWav];
                sds_refl_sdev_[iWav] = sdev_[iWav];
		    } else {
                sds_refl_[iWav] = maxrefl;
                sds_refl_sdev_[iWav] = sdev_[iWav];
		    }
		} else {
            sds_refl_[iWav] = minrefl;
            sds_refl_sdev_[iWav] = 0.0;
		}
    	if (good_pixels_[iWav] >= 0) {
    		sds_numPixels_[iWav] = good_pixels_[iWav];
    	} else {
    	    sds_numPixels_[iWav] = 0;
    	}
	}

	return status;
}

/****************************************************************************
 * NAME: store_output()
 *
 * This subroutine stores all the output variables to be written to output file.
 *
 ***************************************************************************/

int DtAlgOcean::store_output()
{
	int status = DTDB_SUCCESS;

	for ( int iWav=0; iWav<NWAV; iWav++ ) {
		if ((tau_avg_[IBEST][iWav] > -0.1) /* && (tau_avg_[IBEST][iWav] < 5.0) */) {
			sds_tau_best_[iWav] = tau_avg_[IBEST][iWav];
		}
		if ((tau_avg_[IAVG][iWav] > -0.1) /* && (tau_avg_[IBEST][iWav] < 5.0) */) {
			sds_tau_avg_[iWav] = tau_avg_[IAVG][iWav];
		}
		if (tau_small_avg_[IBEST][iWav] > 0.0) {
			sds_tau_small_best_[iWav] = tau_small_avg_[IBEST][iWav];
		}
		if (tau_small_avg_[IAVG][iWav] > 0.0) {
			sds_tau_small_avg_[iWav] = tau_small_avg_[IAVG][iWav];
		}
		if (tau_big_avg_[IBEST][iWav] > 0.0) {
			sds_tau_big_best_[iWav] = tau_big_avg_[IBEST][iWav];
		}
		if (tau_big_avg_[IAVG][iWav] > 0.0) {
			sds_tau_big_avg_[iWav] = tau_big_avg_[IAVG][iWav];
		}
		if (assym_avg_[IBEST][iWav] > 0.0) {
			sds_assy_best_[iWav] = assym_avg_[IBEST][iWav];
		}
		if (assym_avg_[IAVG][iWav] > 0.0) {
			sds_assy_avg_[iWav] = assym_avg_[IAVG][iWav];
		}
		if (backscatter_avg_[IBEST][iWav] > 0.0) {
			sds_back_best_[iWav] = backscatter_avg_[IBEST][iWav];
		}
		if (backscatter_avg_[IAVG][iWav] > 0.0) {
			sds_back_avg_[iWav] = backscatter_avg_[IAVG][iWav];
		}
		if (refl_flux_avg_[IBEST][iWav] > 0.0) {
			sds_reff_best_[iWav] = refl_flux_avg_[IBEST][iWav];
		}
		if (refl_flux_avg_[IAVG][iWav] > 0.0) {
			sds_reff_avg_[iWav] = refl_flux_avg_[IAVG][iWav];
		}
		if (trans_flux_avg_[IBEST][iWav] > 0.0) {
			sds_tranf_best_[iWav] = trans_flux_avg_[IBEST][iWav];
		}
		if (trans_flux_avg_[IAVG][iWav] > 0.0) {
			sds_tranf_avg_[iWav] = trans_flux_avg_[IAVG][iWav];
		}
	}
    for ( int iS=0; iS<NUM_STATS; iS++ ) {
		if (mass_con_ocean_avg_[iS] > 0.0) {
			sds_Mass_Conc_[iS] = mass_con_ocean_avg_[iS];
		}
		if (eff_radius_avg_[iS] > 0.0) {
			sds_EffRad_[iS] = eff_radius_avg_[iS];
		}
		if (ccn_avg_[iS] > 0.0) {
			sds_CCN_[iS] = ccn_avg_[iS];
		}
		if (angstrom_exp_avg_[iS][0] > -1.0) {
			sds_Angs_Coeff1_[iS] = angstrom_exp_avg_[iS][0];
		}
		if (angstrom_exp_avg_[iS][1] > -1.0) {
			sds_Angs_Coeff2_[iS] = angstrom_exp_avg_[iS][1];
		}
        if (funmin_[iS][0] >= 0.0) {
            sds_Least_Error_[iS] = funmin_avg_[iS];
        }
        if (tau_X55_[iS][0] >= 0.0) {
            sds_tau_X55_[iS] = tau_X55_avg_[iS];
        }
		if (xmin_[iS][0] >= 0.0) {
			sds_Small_Weighting_[iS] = xmin_avg_[iS];
		}
		if (solution_index_avg_[iS][0] >= 0.0) {
			sds_Sol_Index_Small_[iS] = (short) solution_index_avg_[iS][0];
		}
		if (solution_index_avg_[iS][1] >= 0.0) {
			sds_Sol_Index_Large_[iS] = (short) solution_index_avg_[iS][1];
		}
    }
    sds_land_ocean_quality_ = quality_to_pass_[0];

	return status;
}

/****************************************************************************
 * NAME: average_output()
 *
 * Subroutine sorts according to minimum error and averages all variables of output
 *
 ***************************************************************************/

int DtAlgOcean::average_output()
{
	int status = DTDB_SUCCESS;

	float sdev = 0.0;
	int   index[NUMCASES*NUMCASEB];
    memset(index, 0, NUMCASES*NUMCASEB*sizeof(int));

	sort_index(NUMCASES*NUMCASEB, funmin_[0], index);

	for ( int iWav=0; iWav<NWAV; iWav++) {
		tau_avg_[IBEST][iWav] = tau_[iWav][index[0]];
		tau_small_avg_[IBEST][iWav] = tau_small_[iWav][index[0]];
		tau_big_avg_[IBEST][iWav] = tau_big_[iWav][index[0]];
		backscatter_avg_[IBEST][iWav] = backscatter_[iWav][index[0]];
		assym_avg_[IBEST][iWav] = assym_[iWav][index[0]];
		refl_flux_avg_[IBEST][iWav] = refl_flux_[iWav][index[0]];
		trans_flux_avg_[IBEST][iWav] = trans_flux_[iWav][index[0]];
	}
	angstrom_exp_avg_[IBEST][0] = angstrom_exp_[0][index[0]];
	angstrom_exp_avg_[IBEST][1] = angstrom_exp_[1][index[0]];
	solution_index_avg_[IBEST][ISMALL] = solution_index_[ISMALL][index[0]];
	solution_index_avg_[IBEST][IBIG] = solution_index_[IBIG][index[0]];
	ccn_avg_[IBEST] = ccn_[index[0]];
    eff_radius_avg_[IBEST] = eff_radius_[index[0]];
    eff_variance_avg_[IBEST] = eff_variance_[index[0]];
	mass_con_ocean_avg_[IBEST] = mass_con_ocean_[index[0]];
	xmin_avg_[IBEST] = xmin_[0][index[0]];
    funmin_avg_[IBEST] = funmin_[0][index[0]];
    tau_X55_avg_[IBEST] = tau_X55_[0][index[0]];

	float weight[NUMCASES*NUMCASEB];
    memset(weight, 0, NUMCASES*NUMCASEB*sizeof(float));
	int num_good = 0;
	for ( int iSol=0; iSol < NUMCASES*NUMCASEB; iSol++ ) {
		if(100.0*funmin_[0][index[iSol]] <= Threshold_LSQ_Error) {
			weight[index[iSol]] = 1.0;
			num_good++;
		}
		else weight[index[iSol]] = 0.0;
	}
	if(num_good < 2) {
		qcontrol_ = 8;
		for ( int iSol=0; iSol<NUMCASES*NUMCASEB; iSol++ )  {
			if (iSol < 3) {
				if(funmin_[0][index[iSol]] > 0.0) {
					weight[index[iSol]] =
					        1.0/pow(100.0*funmin_[0][index[iSol]],2);
				}
			}
			else weight[index[iSol]] = 0.0;
		}
	}
	else qcontrol_ = 0;

	for ( int iWav=0; iWav<NWAV; iWav++) {
		mean_std_weighted( NUMCASES*NUMCASEB, tau_[iWav],
						tau_avg_[IAVG][iWav], sdev, weight );
		mean_std_weighted( NUMCASES*NUMCASEB, tau_small_[iWav],
						tau_small_avg_[IAVG][iWav], sdev, weight );
		mean_std_weighted( NUMCASES*NUMCASEB, tau_big_[iWav],
						tau_big_avg_[IAVG][iWav], sdev, weight );
		mean_std_weighted( NUMCASES*NUMCASEB, backscatter_[iWav],
						backscatter_avg_[IAVG][iWav], sdev, weight );
		mean_std_weighted( NUMCASES*NUMCASEB, assym_[iWav],
						assym_avg_[IAVG][iWav], sdev, weight );
		mean_std_weighted( NUMCASES*NUMCASEB, refl_flux_[iWav],
						refl_flux_avg_[IAVG][iWav], sdev, weight );
		mean_std_weighted( NUMCASES*NUMCASEB, trans_flux_[iWav],
						trans_flux_avg_[IAVG][iWav], sdev, weight );
	}
	mean_std_weighted( NUMCASES*NUMCASEB, angstrom_exp_[0],
			angstrom_exp_avg_[IAVG][0], sdev, weight );
	mean_std_weighted( NUMCASES*NUMCASEB, angstrom_exp_[1],
			angstrom_exp_avg_[IAVG][1], sdev, weight );
	mean_std_weighted( NUMCASES*NUMCASEB, ccn_,
			ccn_avg_[IAVG], sdev, weight );
    mean_std_weighted( NUMCASES*NUMCASEB, eff_radius_,
            eff_radius_avg_[IAVG], sdev, weight );
    mean_std_weighted( NUMCASES*NUMCASEB, eff_variance_,
            eff_variance_avg_[IAVG], sdev, weight );
	mean_std_weighted( NUMCASES*NUMCASEB, mass_con_ocean_,
			mass_con_ocean_avg_[IAVG], sdev, weight );
	mean_std_weighted( NUMCASES*NUMCASEB, xmin_[0],
			xmin_avg_[IAVG], sdev, weight );
    mean_std_weighted( NUMCASES*NUMCASEB, funmin_[0],
            funmin_avg_[IAVG], sdev, weight );
    mean_std_weighted( NUMCASES*NUMCASEB, tau_X55_[0],
            tau_X55_avg_[IAVG], sdev, weight );
	solution_index_avg_[IAVG][ISMALL] = solution_index_[ISMALL][index[0]];
	solution_index_avg_[IAVG][IBIG] = solution_index_[IBIG][index[0]];

	return status;
}

/***************************************************************************
 * NAME: compute_minimum()
 *
 * This function computes a minimum value for aerosol retrieval algorithm.
 * It finds the minimum using a direct least mean square calculation.
 *
 ***************************************************************************/

int DtAlgOcean::compute_minimum(int iBig, int iSmall, int iSol)
{
    int status = DTDB_SUCCESS;

    bool  bUseWav[NWAV];
    //  For summation do not use 1st wavelength
    bUseWav[D488] = false;
    bUseWav[D550] = true;
    bUseWav[D670] = true;
    bUseWav[D865] = true;
    //  if wavelength 1.24 or 1.64 or 2.13 is small  do not use for inversion
    bUseWav[D1240] = (qcontrol_exclude_== 6) ? false : true;
    bUseWav[D1610] = (qcontrol_exclude_== 3 ||
                     qcontrol_exclude_== 5) ? false : true;
    bUseWav[D2250] = (qcontrol_exclude_== 4 ||
                     qcontrol_exclude_== 5) ? false : true;

    //  Compute for tau using radiance value of 0.86 um and tau values of wav550
    float tau_X55[2] = {0.0,0.0};
    interp_extrap(NAOT,refl_[D865],refl_small_[D865][iSmall],
                    lut_.TAUAS[D550][0],tau_X55[ISMALL]);
    interp_extrap(NAOT,refl_[D865],refl_big_[D865][iBig],
                    lut_.TAUAB[D550][0],tau_X55[IBIG]);

    float rbdif[NWAV];
    float sbdif[NWAV];
    float trefl[NWAV][2];
    float scale[NWAV];
    memset(rbdif, 0, NWAV*sizeof(float));
    memset(sbdif, 0, NWAV*sizeof(float));
    memset(trefl, 0, NWAV*2*sizeof(float));
    memset(scale, 0, NWAV*sizeof(float));
    float rb2 = 0;
    float sb2 = 0;
    float rbsb = 0;
    for ( int iWav=0; iWav<NWAV; iWav++ ) {
        float denom = refl_[iWav]-refl_rayl_[iWav] + 0.01;
        denom = (denom < 0.01) ? 0.01 : denom;
        scale[iWav] = 1.0/denom/denom;
        interp_extrap(NAOT,tau_X55[ISMALL],lut_.TAUAS[D550][0],
                refl_small_[iWav][iSmall],trefl[iWav][ISMALL]);
        interp_extrap(NAOT,tau_X55[IBIG],lut_.TAUAB[D550][0],
                refl_big_[iWav][iBig],trefl[iWav][IBIG]);
        if (bUseWav[iWav]) {
            rbdif[iWav] = refl_[iWav] - trefl[iWav][IBIG];
            sbdif[iWav] = trefl[iWav][ISMALL] - trefl[iWav][IBIG];
            rb2 += rbdif[iWav]*rbdif[iWav]*scale[iWav];
            sb2 += sbdif[iWav]*sbdif[iWav]*scale[iWav];
            rbsb += rbdif[iWav]*sbdif[iWav]*scale[iWav];
        }
    }

    float xm = rbsb / sb2;
    xm = (xm > 1) ? 1.0 : xm;
    xm = (xm < 0) ? 0.0 : xm;
    float asum = 0.0;
    float sum_good_pixels = 0.0;
    for ( int iWav=0; iWav<NWAV; iWav++ ) {
    	float mrfl = (xm*trefl[iWav][ISMALL]+(1.0-xm)*trefl[iWav][IBIG]);
        error_[0][iWav][iSol] = (refl_[iWav] - mrfl) / refl_[iWav];
        if (bUseWav[iWav]) {
            asum += error_[0][iWav][iSol]*error_[0][iWav][iSol];
            sum_good_pixels += good_pixels_[iWav];
        }
    }
    xmin_[0][iSol] = xm;
    funmin_[0][iSol] = asum/(NWAV-2);
    tau_X55_[0][iSol] = xm*tau_X55[ISMALL] + (1-xm)*tau_X55[IBIG];

    return status;
}

/***************************************************************************
 * NAME: compute_minimum_baseline()
 *
 * This function computes a minimum value for aerosol retrieval algorithm.
 * It finds the minimum of a function by using an interval halving search.
 *
 ***************************************************************************/

int DtAlgOcean::compute_minimum_baseline(int iBig, int iSmall, int iSol)
{
    int status = DTDB_SUCCESS;

    float iFlag = 0;
    float a[2][5] = {{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0}};
    float x[5] = {0.0, 0.5, 1.0, 0.0, 0.0};
    float f[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    f[0] = fun_tau(x[0], iBig, iSmall, iSol);
	f[1] = fun_tau(x[1], iBig, iSmall, iSol);
	f[2] = fun_tau(x[2], iBig, iSmall, iSol);

    for ( int i=0; i<10; i++ ) {
        x[3] = (x[0]+x[1])/2.0;
        x[4] = (x[1]+x[2])/2.0;
        f[3] = fun_tau(x[3], iBig, iSmall, iSol);
        f[4] = fun_tau(x[4], iBig, iSmall, iSol);
        std::copy(x, x+5, a[0]);
        std::copy(f, f+5, a[1]);

        sort_inplace(5,a[1],a[0]);

        float rmin = a[0][0];
        if( rmin <= x[3]) {
            iFlag++;
            x[2] = x[1];
            x[1] = x[3];
            f[2] = f[1];
            f[1] = f[3];
        }
        if((fabs(rmin-x[1]) < fabs(x[1])*0.000001) && (rmin > x[3])) {
            iFlag++;
            x[0] = x[3];
            x[2] = x[4];
            f[0] = f[3];
            f[2] = f[4];
        }
        if((rmin >= x[4]) &&
        	(fabs(rmin-x[1]) > fabs(x[1])*0.000001) && (rmin > x[3])) {
            iFlag++;
            x[0] = x[1];
            x[1] = x[4];
            f[0] = f[1];
            f[1] = f[4];
        }
    }
    if( iFlag > 1) {

        sort_inplace(3, f, x);

        xmin_[1][iSol] = x[0];
        funmin_[1][iSol] = f[0];
    }
    else {
        xmin_[1][iSol] = DFILL_FLOAT;
        funmin_[1][iSol] = DFILL_FLOAT;
    }

    return status;
}

/***************************************************************************
 * NAME: fun_tau()
 *
 * DESCFRIPTION: Function to be minimized
 *
 ***************************************************************************/

float DtAlgOcean::fun_tau( float xmin, int iBig, int iSmall, int iSol )
{
    float xs[NWAV][NAOT];
    float ys[NWAV][NAOT];
    memset(xs, 0, NWAV*NAOT*sizeof(float));
    memset(ys, 0, NWAV*NAOT*sizeof(float));
//  Compute function to be minimized.
	for ( int iWav=0; iWav<NWAV; iWav++ ) {
		for ( int iTau=0; iTau<NAOT; iTau++ ) {
			xs[iWav][iTau] = lut_.TAUAS[D550][0][iTau];
			ys[iWav][iTau]= xmin*refl_small_[iWav][iSmall][iTau] +
			                (1.-xmin)*refl_big_[iWav][iBig][iTau];
		}
	}

//	Compute for tau using radiance value of 0.86 um and tau values of wav550
    interp_extrap(NAOT,refl_[D865],ys[D865],xs[D550],tau_X55_[1][iSol]);

//  For TAUX55 compute reflectance for all wavelengths.
    float alxwav[NWAV];
    memset(alxwav, 0, NWAV*sizeof(float));
	for ( int iWav=0; iWav<NWAV; iWav++ ) {
		interp_extrap(NAOT,tau_X55_[1][iSol],xs[D550],ys[iWav],alxwav[iWav]);
		error_[1][iWav][iSol] = refl_[iWav] - alxwav[iWav];
	}

	bool  bUseWav[NWAV];
//  For summation do not use 1st wavelength
	bUseWav[D488] = false;
	bUseWav[D550] = true;
	bUseWav[D670] = true;
	bUseWav[D865] = true;
//  if wavelength 1.24 or 1.64 or 2.13 is small  do not use for inversion
	bUseWav[D1240] = (qcontrol_exclude_== 6) ? false : true;
	bUseWav[D1610] = (qcontrol_exclude_== 3 ||
	        qcontrol_exclude_== 5) ? false : true;
	bUseWav[D2250] = (qcontrol_exclude_== 4 ||
	        qcontrol_exclude_== 5) ? false : true;

    float asum = 0.0;
    int   sum_good_pixels = 0;
    for ( int iWav=0; iWav<NWAV; iWav++ ) {
		if (bUseWav[iWav]) {
			float denom = (refl_[iWav]-refl_rayl_[iWav]) + 0.01;
	        if( denom < 0.01) {
	        	denom = 0.01;
	        }
	        asum += pow(error_[1][iWav][iSol]/denom,2.0)*good_pixels_[iWav];
	        sum_good_pixels += good_pixels_[iWav];
		}
	}
	float result = sqrt(asum/((float)sum_good_pixels));

    return result;
}

/***************************************************************************
 * NAME: assign_quality()
 *
 * DESCFRIPTION: This subroutine assigns quality flags and attributes.
 * Product Quality & Retrieval Processing QA flags over ocean
 * Quality of 3 and 4 (average solution) are repeated tion have
 * as best and average solution have same quality.
 *
 * quality_flag_[1]=0     NO retrieval ( NOT useful )
 * quality_flag_[1]=1     retrieval    ( useful)
 *
 * For all non retrieval boxes quality_flag_ is set to zero to indicate
 * not useful and quality_flag_[4] is set to values of qcontrol_ to
 * indicate the cause of Retrieval. quality_flag_[5] is set to 11 for
 * no retrieval.
 *
 * Estimated quality of aerosol parameters
 *	quality_flag_[1]=3     very Good
 *  quality_flag_[1]=2     Good
 *  quality_flag_[1]=1     Marginal
 *  quality_flag_[1]=0     Bad
 *  quality_flag_[3]=3     very Good
 *  quality_flag_[3]=2     Good
 *  quality_flag_[3]=1     Marginal
 *  quality_flag_[3]=0     Bad
 *
 ***************************************************************************/

int DtAlgOcean::assign_quality()
{
	int status = DTDB_SUCCESS;

    quality_to_pass_[0] = 0;
    quality_to_pass_[1] = 0;
    if ( qcontrol_ < 0) {
        quality_flag_[0] = 0;
        quality_flag_[1] = 0;
        quality_flag_[2] = 0;
        quality_flag_[3] = 0;
        quality_flag_[4] = abs(qcontrol_);
        quality_flag_[5] = 15;
    }
    else {
        quality_flag_[0] = 1;
        quality_flag_[2] = 1;
        if( quality_dust_flag_glint_ == 1) {
        	qcontrol_= 12;
        }
        if( qcontrol_cirrus_== 1) {
        	qcontrol_= 13;
        }
        if( quality_dust_flag_off_glint_ == 1) {
        	qcontrol_= 14;
        }
//  if qcontrol_ is 0 ( see doc.) then quality of retrieval is very good
        if( qcontrol_ == 0 ) {
        	quality_flag_[1] = 3;
        	quality_flag_[3] = 3;
        }
//  if qcontrol_ is 7 or 14 ( see doc.) then quality of retrieval is good
        if( qcontrol_ == 7 || qcontrol_ == 14) {
        	quality_flag_[1] = 2;
        	quality_flag_[3] = 2;
        }
//  if qcontrol_ is 1,3,4,6,8 or 10 then quality of retrieval is Average
        if( qcontrol_== 1  || qcontrol_ == 3 || qcontrol_ == 4
            || qcontrol_== 6 || qcontrol_ == 8 || qcontrol_== 10) {
        	quality_flag_[1] = 1;
        	quality_flag_[3] = 1;
        }
//   if qcontrol_ is 2,5,9,12,13 quality of retrieval  is poor
        if( qcontrol_== 2  || qcontrol_== 5 || qcontrol_ == 9
            || qcontrol_ == 12 || qcontrol_== 13 ) {
            quality_flag_[1] = 0;
            quality_flag_[3] = 0;
        }
        quality_flag_[4] = 0;
        quality_flag_[5] = qcontrol_;
        if(qcontrol_ == 17) {
            quality_flag_[1] = 3;
            quality_flag_[3] = 3;
        }
    }
    quality_to_pass_[0] = quality_flag_[1];
    quality_to_pass_[1] = qcontrol_;

    return status;
}

/***************************************************************************
 * NAME: fill_values()
 *
 * DESCFRIPTION: This subroutine assigns fill values to output data
 * where appropriate.
 *
 ***************************************************************************/

int DtAlgOcean::set_fill_out()
{
	int status = DTDB_SUCCESS;

    for ( int iWav=0; iWav<NWAV; iWav++ ) {
        sds_refl_[iWav] = DFILL_FLOAT;
        sds_refl_sdev_[iWav] = DFILL_FLOAT;
        sds_numPixels_[iWav] = DFILL_SHORT;
    }
    for ( int iWav=0; iWav<NWAV; iWav++ ) {
        sds_tau_best_[iWav]  = DFILL_FLOAT;
        sds_tau_avg_[iWav]  = DFILL_FLOAT;
        sds_tau_small_best_[iWav] = DFILL_FLOAT;
        sds_tau_small_avg_[iWav] = DFILL_FLOAT;
        sds_tau_big_best_[iWav] = DFILL_FLOAT;
        sds_tau_big_avg_[iWav] = DFILL_FLOAT;
        sds_assy_best_[iWav] = DFILL_FLOAT;
        sds_assy_avg_[iWav] = DFILL_FLOAT;
        sds_back_best_[iWav] = DFILL_FLOAT;
        sds_back_avg_[iWav] = DFILL_FLOAT;
        sds_reff_best_[iWav] = DFILL_FLOAT;
        sds_reff_avg_[iWav] = DFILL_FLOAT;
        sds_tranf_best_[iWav] = DFILL_FLOAT;
        sds_tranf_avg_[iWav] = DFILL_FLOAT;
    }

	for ( int iCase=0; iCase<NUM_SIZES; iCase++ ) {
		sds_Mass_Conc_[iCase] = DFILL_FLOAT;
		sds_EffRad_[iCase] = DFILL_FLOAT;
		sds_CCN_[iCase] = DFILL_FLOAT;;
		sds_Angs_Coeff1_[iCase] = DFILL_FLOAT;
		sds_Angs_Coeff2_[iCase] = DFILL_FLOAT;
		sds_Least_Error_[iCase] = DFILL_FLOAT;
		sds_Small_Weighting_[iCase] = DFILL_FLOAT;
		sds_Sol_Index_Small_[iCase] = DFILL_FLOAT;
        sds_Sol_Index_Large_[iCase] = DFILL_FLOAT;
        sds_tau_X55_[iCase] = DFILL_FLOAT;
	}

    return status;
}
