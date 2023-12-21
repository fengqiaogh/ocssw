/*******************************************************************************
 *
 * NAME: DtAlgLand.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute processings for a given DtAlgLand object class.
 *
 *  Created on: November 3, 2016
 *      Author: Sam Anderson, DT
 *
 *  Modified:
 *
 *******************************************************************************/

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <fstream>

#include <DDataset.hpp>
#include <DDProcess.h>
#include "darktarget/DtLutNetcdf.h"
#include "darktarget/DtMask.h"
#include "darktarget/DtAlgorithm.h"

using namespace std;

/*
*   This program derives aerosol optical depth (AOD) at 10km x 10km
*   from cloud and water screened 500 meter and 250 meter MODIS measured
*   radiances over land. AOD are derived at 0.553 micron, by inverting
*   radiance data at 0.466,0.644 and 2.113 micron channels. Dark pixels
*   are identified by the 2.113 mid-shortwave-IR (SWIR) *hannel, and
*   sorted so that the brightest (50%) and darkest (20%) pixels are discarded
*   to reduce cloud and other contaminations. Surface reflectance in the
*   the visible channels (0.466 and 0.644 micron) are estimated from
*   the 2.113 channel, using relationships parameterized from vegetation
*   indices in the mid-IR (MVI) channels of 1.24 and 2.113 microns), and
*   also as a function of the solar/surface/satellite scattering angles.
*   Note, that the 2.113 SWIR is not assumed transparent to aerosol, so that
*   the inversion of the two visible and one SWIR channels results in
*   three parameters: 1) spectral AOD normalized to 0.553 micron
*   2) The spectral fine mode AOD (AODfine = AOD * non-dust weighting)
*   and 3) Surface reflectance at 2.113 micron.
*
*   The aerosol models have been derived from the AERONET data base and
*   are considered "dynamic" models of the optical depth.
*   and are fixed as a function of location and season.
*
*   The program assumes that the input data are in reflectance units
*   of PI*L/(F0*cos(theta0)) for all visible and mid-IR channels.
*   The input data are on MODIS spatial resolutions where the pixel size
*   of the 0.66 micron is double the resolution of the 0.466, 1.24 and 2.11
*   mi*ron *hannels
*
*   Note: the program uses input cloud mask file (cloudy=0,clear=1),
*   which is determined by the thresholds of visible reflectance and
*   the visible reflectance ratios.  It requires that for 0.47 micron
*   channel 75% of the subpixels (250 m resolution) to be cloud free
*   and 100% of the subpixels to be water free. For 0.66 micron channel
*   it requires 100% subpixels to be cloud and water free. The selection
*   of 20th to 50th percentile will eliminate the residual clouds.
*   A valid retrieval needs at least 12 remaining pixels.
*
**INPUT PARAMETERS:
*
*   HANDLE_LUTs          Logical unit numbers to open multi-wavelength lookup tables
*   IMONTH               calendar month
*   ISACN                Scan number
*   IDATA                Index of 10x10 km box
*   NUMSQ                Total number of 10x10 km boxes
*   solz_               Solar zenith angle
*   senz_                viewing angle
*   raa_                Relative azimuth angle
*   LAT_CENTER           center latitude of 10x10 km box
*   LON_CENTER           center longitude of 10x10 km box
*   START_500            Starting index for 500 m resolution array
*   END_500              Ending index for 500 m resolution array
*   START_250            Starting index for 250 m resolution array
*   END_250              Ending index for 250 m resolution array
*   start_line            Starting index for 1 km resolution array
*   END_1KM              Ending index for 1 km resolution array
*   W488_SYN             L1B data at 0.466 micron
*   W2250_SYN             L1B data at 2.113 micron
*   W1240_SYN             L1B data at 1.24 micron
*   W1610_SYN             L1B data at 1.64 micron
*   W553_SYN             L1B data at 0.553 micron
*   W670_SYN             L1B data at 0.644 micron
*   W865_SYN             L1B data at 0.865 micron
*   CLDMSK_250           cloud mask in 250 m resolution
*   SET_COUNTER_LAND     counter of land
*
**OUTPUT PARAMETERS:
*
*   QA_Flag_Land         QA flag araay for aerosol over land
*   Success_Ret_Land     Number of successfully retrieved value
*   Fail_Ret_Land        Number of failed retrieved value
*   SDSLAT               HDF SDS array for latitude
*   SDSLON               HDF SDS array for longitude
*   SDS_MTHET0           HDF SDS array for solar zenith angle
*   SDS_MTHET            HDF SDS array for viewing angle
*   SDS_MPHI             HDF SDS array for relative azimuth
*   SDS_Aerosol_Type     HDF SDS array for aerosol type
*   SDS_SCAT_ANGLE_land  HDF SDS array for scattering angle
*   SDS_mass_conc_lan    HDF SDS array for mass concentration
*   SDS_angs_coeff_land  HDF SDS array for angstrom coefficient
*   SDS_CLDFRC_land      HDF SDS array for number of cloudy pixels
*   SDS_dust_weighting   HDF SDS array for dust weighting factor
*   SDS_est_uncer        HDF SDS array for estimated uncertainties
*   SDS_NUMPIXELS_land   HDF SDS array for number of pixels used in retrieval
*   SDSTAU_corrected     HDF SDS array for corrected aerosol optical depth
*   SDSTAU_small_land    HDF SDS array for fine mode aerosol optical depth
*   SDS_ref_land         HDF SDS array for mean apparent reflectance
*   SDS_ref_STD_land     HDF SDS array for std of apparent reflectance
*   SDS_QCONTROL_land    HDF SDS array for quality control
*
**REVISION HISTORY:
* $Log: Process_land_V6.f,v $
* 02/15/2006 levy
* Initial revision
*
**TEAM-UNIQUE HEADER:
*
* Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
* GSFc, Greenbelt, MD
*
*
**DESIGN NOTES:
*
* The following ERROR cODE indicates no retrieval possibly done:
*   Shanacc
*      1     Solar, satellite zenith angle and relative azimuth angles
*            are out of bound of look-up tables
*      2     Threshold for 2.13 um is not met in the entire grid box
*      3     The thresold for number of cloudfree points is not met.
*            (NUMcLDFREE see the parameter) or there is water.
*      4     The computed reflectance for 0.47 or 0.66 channels is out
*            of bounds from look-up table
*      5     Threshold for 2.13 um is not met in option 2 and 3
*      6     Aerosol type can not be determined
*
* The values of aerosol optical thickness corresponding to error code
* Error code       Value of aerosol optical thickness
*
*      1              -3.0
*      2              -4.0
*      3              -5.0
*      4              -6.0
*      5              -7.0
*      6              -8.0
*
*  Aerosol types (FTABLE)
*      1              continental
*      2              Generic : SSA ~ 0.9
*      3              Smoke   : SSA ~ 0.85
*      4              Urban   : SSA ~ 0.95
*      5              Dust
*
* Procedure used in retrieving aerosol optical thickness:
*
*      0 = no dark targets possibly met by the following thresholds
*      1 = threshold used for wave 2.13 um 0.01 - 0.25
*      2 = threshold used for wave 2.13 um 0.25 - 0.40
*
* INTERNALS:
*
*   Subroutines:
*      RNLOOKUP      -  READS THE LOOK-UP TABLES
*      ERROR         -  SETS THE ERROR cODES
*      OUTPUT        -  WRITES THE OUTPUT TO HDF FILE
*      AEROSOL_MAP   -  Reads map of aerosol model as function of season and place
*      INTANGLE_NL_RAY  -  Interpolates LUT Rayleigh parameters to measured geometry
*      Process_land_Rob -  The main code for doing the radiance to aerosol inversion
*      cOMPUTE_ScATTANGLE_LAND - computes scattering angle
*      Statistics_land  -  Statistics for Path Radiance and critical Reflectance
*      FILLVALUE_LAND - Inserts fill values for parameters if errors
*
*
*   Variables:
*     W488_SYN(2*IGRIDX,2*IGRIDY)    Reflectance for wav=0.47um
*     W550_SYN(2*IGRIDX,2*IGRIDY)    Reflectance for wav=0.55um
*     W2250_SYN(2*IGRIDX,2*IGRIDY)    Reflectance for wav=2.13um
*     W1240_SYN(2*IGRIDX,2*IGRIDY)    Reflectance for wav=1.24um
*     W1610_SYN(2*IGRIDX,2*IGRIDY)    Reflectance for wav=1.64um
*     W670_SYN(4*IGRIDX,4*IGRIDY)    Reflectance for wav=0.66um
*     W865_SYN(4*IGRIDX,4*IGRIDY)    Reflectance for wav=0.86um
*     cLDMSK_250(ISWATH_B,ILINE)   cloud mask of 250 m resolution
*                                    (0=cloudy,1=clear)
*     solz_                         solar zenith angle (degree)
*     senz_                          satellite viewangle angle (degree)
*     raa_                          relative azimuth in angle (degree)
*     MHGHT                          Topographic altitude (km)
*     yint644,yint466,slope466,slope644
*                 characteristics of VIS/IR surface reflectance relationship
*
* LAND SDS_ARRAYS.........
*
*    SDS_QCONTROL_land      Quality control SDS array
*    SDS_Aerosol_Type       Index of Aerosol type
*    SDS_SCAT_ANGLE_land    Scattering Angle
*    SDS_angs_coeff_land    Angstrom exponent for 0.47 and 0.67 micron
*    SDS_CLDFRC_land        cloud fraction (%)
*    SDS_dust_weighting     Dust aerosol weighting factor
*    SDS_est_uncer          Uncertainty of optical thickness at 0.47 and 0.66 micron
*    SDS_NUMPIXELS_land     Number of pixels with desired percentile
*    SDSTAU_corrected       corrected optical thickness at 0.47 0.55 and 0.66 micron
*    SDS_ref_land           Mean reflectance at five bands
*    SDS_ref_STD_land       Standard deviation of reflectance at five bands
*    SDS_mass_conc_land     Mass concentration
*/


/**************************************************************************
 * NAME: DtAlgLand()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtAlgLand::DtAlgLand()
{
	cm_= nullptr;
    scatter_angle_= 0;
    memset(mtable_, 0, NLSIZE*sizeof(int));
    memset(refl_ray_nl_, 0, NLWAV*sizeof(float));
    memset(opth_nl_, 0, NLSIZE*NLWAV*NLTAU*sizeof(float));
    memset(int_nl_, 0, NLSIZE*NLWAV*NLTAU*sizeof(float));
    memset(fd_nl_, 0, NLSIZE*NLWAV*NLTAU*sizeof(float));
    memset(t_nl_, 0, NLSIZE*NLWAV*NLTAU*sizeof(float));
    memset(fdt_nl_, 0, NLSIZE*NLWAV*NLTAU*sizeof(float));
    memset(sbar_nl_, 0, NLSIZE*NLWAV*NLTAU*sizeof(float));
    memset(refl_inter_, 0, NLWAV*GRIDX*GRIDY*sizeof(float));
	memset(refl_, 0, NLWAV*sizeof(float));	// mean reflectance
	memset(sdev_, 0, NLWAV*sizeof(float));	// reflectance standard deviation
	error_= 0;
    memset(rho_star_, 0, NLSIZE*NLWAV*NLTAU*sizeof(float));
    memset(rho_star_tot_, 0, NLWAV*NLTAU*sizeof(float));
    memset(rho_S212_, 0, NLSIZE*NLTAU*sizeof(float));
    memset(rho_S212_tot_, 0, NLTAU*sizeof(float));
    memset(errwave_, 0, NLWAV*sizeof(float));
    rho_S466_= 0;
	yint_466_= 0;
	slope_466_= 0;
    rho_S644_= 0;
	yint_644_= 0;
	slope_644_= 0;
	eta_flag_= 0;
    memset(aot_d_, 0, NLWAV*sizeof(float));  // tau corrected
    memset(aot_f_, 0, NLWAV*sizeof(float));  // tau fine (small)
    memset(aot_c_, 0, NLWAV*sizeof(float));  // tau coarse (big)
    memset(rho_sfc_, 0, NLWAV*sizeof(float));  // surface reflectance
    eta_= 0;  // dust weighting
    err644_= 0;  // fitting error
    masscon_= 0;  // mass concentration
    angstrom_= 0;  // angs coefficient
    ndvi_= 0;  // ndvi
    iaer_= 0;  // aerosol type
    memset(good_pixels_, 0, NLWAV*sizeof(int)); // good pixels
	num_pixels_used_= 0;  // = good_pixels_[W659]
	iproc_= 0;
	ifinish_= 0;
	thresh_min_= 0;
	thresh_max_= 0;
	return_quality_cirrus_= 0;
	success_ret_= 0;
	fail_ret_= 0;
	memset(quality_flag_for_joint_, 0, 2*sizeof(short));
	quality_flag_for_retr_= 0;
	qcontrol_special_= 0;
	memset(sds_qcontrol_, 0, QA_LAND*sizeof(short));
	sds_aerosol_type_= 0;
	sds_scat_angle_= 0;
    sds_fitting_error_= 0;
    sds_mass_conc_= 0;
    sds_cloud_fraction_= 0;
	sds_angs_coeff_= 0;
	sds_dust_weighting_= 0;
	sds_ndvi_= 0;
	memset(sds_numpixels_, 0, NLWAV*sizeof(float));
	memset(sds_tau_corrected_, 0, NLWAV*sizeof(float));
	memset(sds_refl_, 0, NLWAV*sizeof(float));
	memset(sds_refl_std_, 0, NLWAV*sizeof(float));
    memset(sds_tau_small_, 0, NLWAV*sizeof(float));
    memset(sds_tau_big_, 0, NLWAV*sizeof(float));
	memset(sds_surface_reflectance_, 0, NLWAV*sizeof(float));

}

/**************************************************************************
 * NAME: ~DtAlgLand()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DtAlgLand::~DtAlgLand()
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

int DtAlgLand::initialize( map<string, ddata*> imap )
{
	int status = DTDB_SUCCESS;

    DtAlgorithm::initialize( imap );

	DtLutNetcdf* lutgen = new DtLutNetcdf();
	status = lutgen->read_land_aerosol_lut( lut_ );
	delete lutgen;

    cm_ = new DtCloudMaskLand(this);

    return status;
}

/**************************************************************************
 * NAME: process()
 *
 * DESCRIPTION: Dark Target land algorithm process.
 *
 *************************************************************************/

map<string, ddata*> DtAlgLand::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap)
{
    get_inputs(start, count, imap);

    l2_flags_ |= (unsigned int) flags::LAND;

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
    rfld_[D412] = rfl_[(size_t)rhot_band::W410];
    rfld_[D488] = rfl_[(size_t)rhot_band::W490];
    rfld_[D550] = rfl_[(size_t)rhot_band::W550];
    rfld_[D670] = rfl_[(size_t)rhot_band::W670];
    rfld_[D865] = rfl_[(size_t)rhot_band::W865];
    rfld_[D1240] = rfl_[(size_t)rhot_band::W1240];
    rfld_[D1610] = rfl_[(size_t)rhot_band::W1610];
    rfld_[D2250] = rfl_[(size_t)rhot_band::W2250];
    if ((raa_ < 0.0) || (raa_ > 180.0) || rfld_[D488] <= 0.0) {
    	return set_fills();
    }

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
		if (cmask_ != DFILL_UBYTE) {
			cloud_mask_ = (cmask_+1)%2;
		}
	} else {
		cmask_ = (cloud_mask_+1)%2;
	};
	if (!bmaskcloud_) cmask_ = 1;
	if (bmaskcloud_ && cloud_mask_) {
		l2_flags_ |= (unsigned int) flags::CLDICE;
		return set_fills();
	}

	index_geometry(solz_, senz_, raa_);
	compute_glint_refl();
	compute_scatter_angle(scatter_angle_);
	interpolate_rayleigh();

	iproc_ = 0;
	ifinish_ = 0;
	thresh_min_ = THR213MIN_1;
	thresh_max_ = THR213MAX_1;

	if ((cmask_ > 0) &&
		rfld_[D488]>0 && rfld_[D670]>0 && rfld_[D865]>0 && rfld_[D2250]>0) {

		ndvi_ = (rfld_[D865]-rfld_[D670])/(rfld_[D865]+rfld_[D670]);
		int ilon = 180 + round(lon_ + DLON);
		int ilat = 90 - round(lat_ + DLAT);

//  Select LUT corresponding to season and location
//  Fine mode LUT (choice of tables #2, #3 or #4)
//  Continental mode LUT = #1;  Coarse mode LUT = #5

		mtable_[ISMALL] = lut_.AEROSOL_ALL[season_][ilat][ilon] + 1;
		mtable_[IBIG] = DTABLE-1;

		index_geometry( solz_, senz_, raa_);
		interpolate_angle();
		interpolate_elevation();
		simulate_toa();
		eta_flag_ = 0;
		retrieve_first();
		if (aot_d_[DL670] > 0) {
			angstrom_ = -1.0 * log(aot_d_[DL488]/aot_d_[DL670])/ log(0.466/0.644);
		}
		iproc_=1;
	}
	iaer_ = mtable_[ISMALL]+1;

    assign_quality();
    store_output();

// write to generic outputs

    qual_flag_ = quality_flag_for_joint_[0];
	aerosol_type_ = sds_aerosol_type_;
    error_flag_ = quality_flag_for_joint_[1];

    scatter_ang_ = 180.0 - glint_angle_;
    glint_ang_ = glint_angle_;
    sse_ = sds_fitting_error_;
    fmf_ = sds_dust_weighting_;
    aot_550_ = sds_tau_corrected_[DL550];
    ae1_ = sds_angs_coeff_;
    ndv_ = sds_ndvi_;
    ae2_ = DFILL_FLOAT;
    chlor_ = DFILL_FLOAT;
    for ( int iWav=0; iWav<NOWL; iWav++ ) {
    	aot_[iWav] = DFILL_FLOAT;
    }
    aot_[(size_t)aot_band::W490] = sds_tau_corrected_[DL488];
    aot_[(size_t)aot_band::W550] = sds_tau_corrected_[DL550];
    aot_[(size_t)aot_band::W670] = sds_tau_corrected_[DL670];
    aot_[(size_t)aot_band::W2250] = sds_tau_corrected_[DL2250];

    sr_[(size_t)srf_band::W410] = DFILL_FLOAT;
    ssa_[(size_t)srf_band::W410] = DFILL_FLOAT;
    sr_[(size_t)srf_band::W490] = sds_surface_reflectance_[DL488];
    ssa_[(size_t)srf_band::W490] = DFILL_FLOAT;
    sr_[(size_t)srf_band::W670] = sds_surface_reflectance_[DL670];
    ssa_[(size_t)srf_band::W670] = DFILL_FLOAT;
    sr_[(size_t)srf_band::W2250] = sds_surface_reflectance_[DL2250];
    ssa_[(size_t)srf_band::W2250] = DFILL_FLOAT;

	vector<float> tba = {490.0,550.0,670.0,2250.0};
	vector<float> yba = {aot_[(size_t)aot_band::W490],aot_[(size_t)aot_band::W550],aot_[(size_t)aot_band::W670],aot_[(size_t)aot_band::W2250]};
	using boost::math::interpolators::barycentric_rational;
	barycentric_rational<float> interp(move(tba), move(yba));
	aot_[(size_t)aot_band::W865] = interp(865.0);
	aot_[(size_t)aot_band::W1240] = interp(1240.0);
	aot_[(size_t)aot_band::W1610] = interp(1610.0);

    return set_outputs();
}

/**************************************************************************
 * NAME: index_geometry()
 *
 * DESCRIPTION: Reset LUT indices based on measured geometry.
 *
 *************************************************************************/

int DtAlgLand::index_geometry( float sza, float azim, float phi)
{
	int status = DTDB_SUCCESS;

    for (int iTh0=0; iTh0< NLTHET0-1; iTh0++) {
        if ((sza >= lut_.THET0_NL[iTh0]) && (sza <= lut_.THET0_NL[iTh0+1])) {
        	SZA_0_= iTh0;
        	SZA_1_= iTh0+1;
        }
    }
    for (int iTh=0; iTh < NLTHE-1; iTh++) {
        if ((azim >= lut_.THE_NL[iTh]) && (azim <= lut_.THE_NL[iTh+1])) {
        	THE_0_= iTh;
        	THE_1_= iTh+1;
        }
    }
    for (int iPhi=0; iPhi < NLPHI-1; iPhi++) {
        if ((phi >= lut_.PHI_NL[iPhi]) && (phi <= lut_.PHI_NL[iPhi+1])) {
        		PHI_0_=iPhi;
				PHI_1_=iPhi+1;
        }
    }

	return status;
}

/**************************************************************************
 * NAME: average()
 *
 * DESCRIPTION: This subroutine processes 10*10 pixel box for cloud detection and
 *    finds the average reflectance for red and blue channels. Surface
 *    Reflectance from wavelength 2.13 . This surface Reflectance and
 *    average Reflectance for red and blue channel are send to lookup
 *    table and optical thickness is derived.
 *
 *************************************************************************/

int DtAlgLand::compute_average(size_t iy, size_t ix)
{
	int status = DTDB_SUCCESS;

	int numCldRed = 0;
	int numCldBlue = 0;
	ifinish_ = 0;
	error_ = 4;

    compute_cloudmask_ndvi(iy, ix, numCldRed, numCldBlue);

    int aindex[GRIDX*GRIDY];
    memset(aindex, 0, GRIDX*GRIDY*sizeof(int));
    float array_interm[NLWAV][GRIDX*GRIDY];
    memset(array_interm, 0, NLWAV*GRIDX*GRIDY*sizeof(float));

//  I213>0 and Icld>0 indicate that the threshold for 2.13 micron
//  is met and the data set is cloud free as defined.
//  index in ascending order and get the index for reflectances at 0.66 um

	 sort_index(numCldBlue,refl_inter_[DL670],aindex);

	 int N20 = (int) ((float)numCldBlue/5.0);
	 int N50 = (int) ((float)numCldBlue/2.0);
//	 if (N20 == 0) N20=1;
	 if ((N50 == 0) && (numCldBlue > 0)) N50=1;
	 int number_pixels = 1;
//  Throw away dark & bright pixels
//  At least 5% after pixels are thrown out
     if((N50-N20) > number_pixels) {
    	 ifinish_ = 1;
    	 for ( int iWav=0; iWav<NLWAV; iWav++ ) {
			good_pixels_[iWav]=0;
			for ( int ix=N20; ix<N50; ix++ ) {
				float ril = refl_inter_[iWav][aindex[ix]];
				if ((ril > 0.0) && (ril <= 1.0)) {
					array_interm[iWav][good_pixels_[iWav]] = ril;
					good_pixels_[iWav] += 1;
				}
			}
    	 }
//  Report  valid pixels used for Blue channel at 500 meter resolution
         num_pixels_used_ = good_pixels_[DL670];
         for ( int iWav = 0; iWav < NLWAV; iWav++ ) {
             if(good_pixels_[iWav] > 0) {
            	 float refl_val = 0;
            	 float sd_val = 0;

            	 mean_std(good_pixels_[iWav],
            			 array_interm[iWav], refl_val, sd_val);

                 refl_[iWav] = refl_val;
                 sdev_[iWav] = sd_val;
             }
             else {
                 refl_[iWav] = DFILL_FLOAT;
                 sdev_[iWav] = DFILL_FLOAT;
             }
         }
//  Check if  3 wavelengths have valid range
         if( !((refl_[DL488] > 0) &&
               (refl_[DL670] > 0) &&
			   (refl_[DL2250] > 0))) {

             error_ = 2;
        	 for ( int iWav=0; iWav<NLWAV; iWav++ ) {
        		 refl_[iWav] = DFILL_FLOAT;
        		 sdev_[iWav] = DFILL_FLOAT;
        	 }
         }
     }
     else {
         error_ = 3;
    	 for ( int iWav=0; iWav<NLWAV; iWav++ ) {
    		 refl_[iWav] = DFILL_FLOAT;
    		 sdev_[iWav] = DFILL_FLOAT;
    		 good_pixels_[iWav] = DFILL_SHORT;
    	 }
     }

     return status;
}

/**************************************************************************
 * NAME: cloudmask_ndvi()
 *
 * DESCRIPTION: This subroutine finds the dark targets using a
 *               threshold in the 2.13 micron channel.it averages the
 *               cloud free pixels in red channel to the 0.5 km resolution.
 *               If too cloudy or too bright at 2.13 it leaves the value
 *               as zero.
 *
 *************************************************************************/

int DtAlgLand::compute_cloudmask_ndvi(size_t iy, size_t ix,
								  	  	  int& numCldRed, int& numCldBlue)
{
	int status = DTDB_SUCCESS;

    numCldBlue = 0;
    numCldRed = 0;
    for (int iWav=0; iWav < NLWAV; iWav++) {
        for (int j=0; j<GRIDX*GRIDY; j++) {
        	refl_inter_[iWav][j] = 0;
        }
    }
//  If snowy pixels set cloud mask to Zero
    if(cm_->snowmask_ == 0) {
        cmask_ = 0;
    }
    numCldRed = 0;
    int noWater = 0;
//  Check if threshold for 2.13 um is met
//  (missing pixels or noisy detectors are ignored)
    if ((rfld_[D2250] <= thresh_max_) &&
        (rfld_[D2250] > thresh_min_) && (cmask_ == 1)) {
        numCldRed++;
        if ((rfld_[D865] > 0.0) && (rfld_[D670] > 0.0)) {
//  NDVI test according to Eric Vermote
            ndvi_ = (rfld_[D865]-rfld_[D670])/(rfld_[D865]+rfld_[D670]);
            if (ndvi_ > 0.1) {
                noWater++;
            }
        }
    }

	return status;
}

/**************************************************************************
 * NAME: interpolate_rayleigh()
 *
 * DESCRIPTION: Subroutine INTANGLE_NL interpolates the lookup
 *               Rayleigh reflectances to the measured geometry.
 *
 *************************************************************************/

int DtAlgLand::interpolate_rayleigh()
{
	int status = DTDB_SUCCESS;

//  Initialize
    float x0[2] = {0.0,0.0};
    float y0[2] = {0.0,0.0};
    float x1[2] = {0.0,0.0};
    float y1[2] = {0.0,0.0};
    float x2[2] = {0.0,0.0};
    float y2[2] = {0.0,0.0};
//  Interpolate Rayleigh
    int iTable = 0;
    int iTau = 0;
    for ( int iWav=0; iWav<NLWAV; iWav++ ) {
        int numWave=0;
        for ( int iTh0=SZA_0_; iTh0<=SZA_1_; iTh0++ ) {
            int numSza=0;
            for ( int iTh=THE_0_; iTh<=THE_1_; iTh++ ) {
                int numThe=0;
                for ( int iPhi=PHI_0_; iPhi<=PHI_1_; iPhi++ ) {
                    x0[numThe]=lut_.PHI_NL[iPhi];
                    y0[numThe]=lut_.INT_NL0[iTable][iWav][iTau][iTh0][iTh][iPhi];
                    numThe++;
                }
                interp_extrap(numThe,raa_,x0,y0,y1[numSza]);
                x1[numSza]=lut_.THE_NL[iTh];
                numSza++;
            }
            interp_extrap(numSza,senz_,x1,y1,y2[numWave]);
            x2[numWave]=lut_.THET0_NL[iTh0];
            numWave++;
        }
        interp_extrap(numWave,solz_,x2,y2,refl_ray_nl_[iWav]);
    }

	return status;
}

/**************************************************************************
 * NAME: interpolate_angle()
 *
 * DESCRIPTION: Subroutine interpolates the lookup reflectances to the
 * 				measured geometry.
 *
 *************************************************************************/

int DtAlgLand::interpolate_angle()
{
	int status = DTDB_SUCCESS;

    float x0[2] = {0.0,0.0};
    float y0[2] = {0.0,0.0};
    float x1[2] = {0.0,0.0};
    float y1[2] = {0.0,0.0};
    float w1[2] = {0.0,0.0};
    float x2[2] = {0.0,0.0};
    float y2[2] = {0.0,0.0};
    float w2[2] = {0.0,0.0};
    float v2[2] = {0.0,0.0};
    float u2[2] = {0.0,0.0};
    for ( int  iSize=0; iSize<NLSIZE; iSize++ ) {
    	for ( int iWav=0; iWav<NLWAV; iWav++ ) {
            for ( int iTau=0; iTau<NLTAU; iTau++ ) {
            	int nTau = 0;
                for ( int iSza=SZA_0_; iSza<=SZA_1_; iSza++ ) {
                	int nTh = 0;
                    for ( int iTh=THE_0_; iTh<=THE_1_; iTh++ ) {
                    	int nPhi = 0;
                        for ( int iPhi=PHI_0_; iPhi<=PHI_1_; iPhi++ ) {
                            x0[nPhi]=lut_.PHI_NL[iPhi];
                            y0[nPhi]=lut_.INT_NL0[mtable_[iSize]][iWav][iTau][iSza][iTh][iPhi];
                            nPhi++;
                        }
                        interp_extrap(nPhi,raa_,x0,y0,y1[nTh]);
                        x1[nTh]=lut_.THE_NL[iTh];
                        w1[nTh]=lut_.T_NL0[mtable_[iSize]][iWav][iTau][iSza][iTh];
                        nTh++;
                    }
					interp_extrap(nTh,senz_,x1,y1,y2[nTau]);
					interp_extrap(nTh,senz_,x1,w1,w2[nTau]);
					x2[nTau]=lut_.THET0_NL[iSza];
					v2[nTau]=lut_.Fd_NL0[mtable_[iSize]][iWav][iTau][iSza];
					u2[nTau]=lut_.SBAR_NL0[mtable_[iSize]][iWav][iTau][iSza];
					nTau++;
                }
                interp_extrap(nTau,solz_,x2,y2,int_nl_[iSize][iWav][iTau]);
                interp_extrap(nTau,solz_,x2,w2,t_nl_[iSize][iWav][iTau]);
                interp_extrap(nTau,solz_,x2,v2,fd_nl_[iSize][iWav][iTau]);
                interp_extrap(nTau,solz_,x2,u2,sbar_nl_[iSize][iWav][iTau]);
                fdt_nl_[iSize][iWav][iTau] =
                		t_nl_[iSize][iWav][iTau]*fd_nl_[iSize][iWav][iTau];
			}
		}
	}

	return status;
}

/**************************************************************************
 * NAME: interpolate_elevation()
 *
 * DESCRIPTION: The subroutine interpolates the lookup
 * reflectances to the target elevation.
 * It interpolates between wavelengths to simulate
 * elevation by a longer wavelength
 *
 *************************************************************************/

int DtAlgLand::interpolate_elevation()
{
	int status = DTDB_SUCCESS;

	float rod_pres[NLWAV];
    float eqwav_nl[NLWAV];
    memset(rod_pres, 0, NLWAV*sizeof(float));
    memset(eqwav_nl, 0, NLWAV*sizeof(float));

//  Estimate surface pressure (hyposometric EQ)
    float pressure = PRESSURE_P0 * exp(-(height_/7.5));
//  Estimate ROD at nominal wavelengths at p0 and at pres
    for ( int iWav=0; iWav<NLWAV; iWav++ ) {
        float expfactor = -4.15 + (0.2 * lut_.WAV_NL[iWav]);
        rod_pres[iWav] = 0.0088*(pressure/PRESSURE_P0) *
        				 pow(lut_.WAV_NL[iWav],expfactor);
//  Estimate equivalent wavelengths for ROD at pressure
        float lambda0 = 0.0;
        float lambda1 = 0.1;
        float lambda2 = 4.0;
        float diff0 = 99.0;
        while (diff0 > 0.00001) {
            lambda0 = (lambda1 + lambda2) / 2.0;
            float fexp0 = -4.15 + 0.2*lambda0;
            float fexp1 = -4.15 + 0.2*lambda1;
			float fexp2 = -4.15 + 0.2*lambda2;
			float ftau0 = 0.0088*pow(lambda0,fexp0);
			float ftau1 = 0.0088*pow(lambda1,fexp1);
			float ftau2 = 0.0088*pow(lambda2,fexp2);
            if ((ftau1 > rod_pres[iWav]) && (ftau2 < rod_pres[iWav])) {
                if (ftau0 > rod_pres[iWav]) {
                    lambda1 = (lambda1 + lambda2)/2.0;
                }
                else {
                    lambda2 = (lambda1 + lambda2)/2.0;
                }
            }
            diff0 = fabs(ftau0 - rod_pres[iWav]);
        }
        eqwav_nl[iWav] = log(lambda0);
    }
// Apply slight modification to optical depth lut -SSA
    for (int iTab=0; iTab<DTABLE; iTab++) {
        for (int iWav=0; iWav<NLUTWAV; iWav++) {
            lut_.OPTH_NL0[iTab][iWav][0] =
                    0.75*lut_.OPTH_NL0[iTab][iWav][1];
        }
    }
//  Interpolate lookup tables to equiv Waves (let's start only with
//  1st two wavelengths until we derive 0.86 table)
    for ( int iSize=0; iSize<NLSIZE; iSize++ ) {
        for ( int iWav=0; iWav<NLWAV; iWav++ ) {

            for ( int iTau=0; iTau<NLTAU; iTau++ ) {
            	float yout = 0.0;
            	float x[2] = {0.0,0.0};
            	float y[2] = {0.0,0.0};
            	float w[2] = {0.0,0.0};
            	float v[2] = {0.0,0.0};
            	float u[2] = {0.0,0.0};
            	float t[2] = {0.0,0.0};
            	float z[2] = {0.0,0.0};
                int iWav1 = (iWav == 3) ? iWav-1 : iWav;
                int iWav2 = (iWav == 3) ? iWav : iWav+1;
                int numWav=0;
                for ( int jWav=iWav1; jWav<=iWav2; jWav++ ) {
 //     Interpolate on log log scale
                    x[numWav]=log(lut_.WAV_NL[jWav]);
                    y[numWav]=log(int_nl_[iSize][jWav][iTau]);
                    w[numWav]=log(fdt_nl_[iSize][jWav][iTau]);
                    u[numWav]=log(fd_nl_[iSize][jWav][iTau]);
                    t[numWav]=log(t_nl_[iSize][jWav][iTau]);
                    z[numWav]=log(sbar_nl_[iSize][jWav][iTau]);
                    v[numWav]=lut_.OPTH_NL0[mtable_[iSize]][jWav][iTau];
                    if (lut_.OPTH_NL0[mtable_[iSize]][jWav][iTau] > 0.) {
                        v[numWav]=log(lut_.OPTH_NL0[mtable_[iSize]][jWav][iTau]);
                    }
                    numWav++;
                }
                interp_extrap(numWav,eqwav_nl[iWav],x,y,yout);
                int_nl_[iSize][iWav][iTau] = exp(yout);
                interp_extrap(numWav,eqwav_nl[iWav],x,w,yout);
                fdt_nl_[iSize][iWav][iTau] = exp(yout);
                interp_extrap(numWav,eqwav_nl[iWav],x,u,yout);
                fd_nl_[iSize][iWav][iTau] = exp(yout);
                interp_extrap(numWav,eqwav_nl[iWav],x,t,yout);
                t_nl_[iSize][iWav][iTau] = exp(yout);
                interp_extrap(numWav,eqwav_nl[iWav],x,z,yout);
                sbar_nl_[iSize][iWav][iTau] = exp(yout);
                interp_extrap(numWav,eqwav_nl[iWav],x,v,yout);
                if (yout == 0.0) {
                    opth_nl_[iSize][iWav][iTau] = yout;
                } else {
                    opth_nl_[iSize][iWav][iTau] = exp(yout);
                }
			}
            refl_ray_nl_[iWav] = int_nl_[iSize][iWav][0];
		}
	}

	return status;
}

/**************************************************************************
 * NAME: retrieve_first()
 *
 * DESCRIPTION: This subroutine retrieves optical thickness, surface
 * reflectance and error parameters by comparing observations and LUT data
 *
 * INPUTS
 *       REFW466L,REFW644L,REFW212L   reflectance
 *       OPTH_NL      LUT based optical thickness
 *       MASSCOEF_NL  LUT based mass concentration coeficients
 *       EXTNORM_NL   LUT based normalized extinction coeficients
 *       RHO_S212     simulated surface reflectance
 *       RHOSTAR      simulated TOA reflectance
 *       yint644,yint466,slope466,slope644
 *          characteristics of VIS/IR surface reflectance relationship
 * OUTPUTS
 *       ERR644	    retrieval fitting errors
 *       AOD        retrieved Optical depths
 *       AODF/C     Fine and coarse mode optical depth
 *       RHOSFC     retrieved surface reflectance
 *       ETA        retrieved fine mode weighting
 *       MASSCON    retrieved mass concentration coefficient
 *       ETA_FLAG   0 if eta within 0.0 - 1.0
 *
 *************************************************************************/

int DtAlgLand::retrieve_first()
{
	int status = DTDB_SUCCESS;

//  Initialize
    float yout=0;
    for ( int iWav=0; iWav<NLWAV; iWav++ ) {
        rho_sfc_[iWav] = DFILL_FLOAT;
        aot_d_[iWav] = DFILL_FLOAT;
    }
//  Now loop through ETA values (kind of like ocean)
    float aot553_temp[NLETA];
    float rhosfc212_temp[NLETA];
    float err644_temp[NLETA];
    float rho_star_temp[NLETA][NLWAV];
    memset(aot553_temp, 0, NLETA*sizeof(float));
    memset(rhosfc212_temp, 0, NLETA*sizeof(float));
    memset(err644_temp, 0, NLETA*sizeof(float));
    memset(rho_star_temp, 0, NLETA*NLWAV*sizeof(float));

    for ( int iEta=0; iEta<NLETA; iEta++ ) {
//  Average, best, QCONTROL, FILLVALUE, OUTPUT, etc
        float eta_temp = -0.2 + (iEta+1) * 0.1;
//  Compute total apparent reflectance
        for ( int iTau=0; iTau<NLTAU; iTau++ ) {
            for ( int iWav=0; iWav<NLWAV; iWav++ ) {
                rho_star_tot_[iWav][iTau] = eta_temp *
                rho_star_[ISMALL][iWav][iTau] +
                ((1-eta_temp) * rho_star_[IBIG][iWav][iTau]);
            }
// 	Compute total surface reflectance
            rho_S212_tot_[iTau] = eta_temp * rho_S212_[ISMALL][iTau] +
                ((1-eta_temp) * rho_S212_[IBIG][iTau]);
        }
//  Interpolate everything based on measured reflectance at 466.

        interp_extrap(NLTAU,rfld_[D488],rho_star_tot_[DL488],
        			opth_nl_[ISMALL][DL550],aot553_temp[iEta]);
        aot553_temp[iEta] = (aot553_temp[iEta]<0) ? 0.0 : aot553_temp[iEta];
        //  Compute apparent reflectance for all wavelengths
        for ( int iWav=0; iWav<NLWAV; iWav++ ) {

            interp_extrap(NLTAU,aot553_temp[iEta],
            		opth_nl_[ISMALL][DL550],
            		rho_star_tot_[iWav],rho_star_temp[iEta][iWav]);
        }
//   COMPUTE Surface Reflectance at 2.1

        interp_extrap(NLTAU,aot553_temp[iEta],
        		opth_nl_[ISMALL][DL550],
        		rho_S212_tot_,rhosfc212_temp[iEta]);

//  Errors
        errwave_[DL488] =
        		fabs(rfld_[D488]-rho_star_temp[iEta][DL488]) / rfld_[D488];
        errwave_[DL670] =
        		fabs(rfld_[D670]-rho_star_temp[iEta][DL670]) / rfld_[D670];
        errwave_[DL2250] =
        		fabs(rfld_[D2250]-rho_star_temp[iEta][DL2250]) / rfld_[D2250];
//  Determine error at 644nm
        err644_temp[iEta] = errwave_[DL670];
    }
//  Sort results in ascending error order
    int nleta_index[NLETA];

	sort_index(NLETA, err644_temp, nleta_index);

//  Best solution
	eta_ = -0.2 + (nleta_index[0]+1) * 0.1;
	err644_ = err644_temp[nleta_index[0]];
	rho_sfc_[DL2250] = rhosfc212_temp[nleta_index[0]];
	rho_sfc_[DL670] = yint_644_ + slope_644_* rho_sfc_[DL2250];
	rho_sfc_[DL488] = yint_466_ + slope_466_*rho_sfc_[DL670];
//  set  out of bounds Eta to 0 or 1
	if( eta_ < 0) {
		eta_= 0.0;
		eta_flag_ = -1;
	}
	if( eta_ > 1) {
		eta_= 1.0;
		eta_flag_ = -1;
	}
//  Compute fine and coarse AOD for mass concentration
	float aot_F553 = aot553_temp[nleta_index[0]] * eta_;
	float aot_C553 = aot553_temp[nleta_index[0]] * (1.-eta_);
//  COMPUTE Mass Concentration
//  Fine

	interp_extrap(NLTAU,aot_F553,opth_nl_[ISMALL][DL550],
			lut_.MASSCOEF_NL0[mtable_[ISMALL]][DL550],yout);

	float massconf = yout*aot_F553;
//  Coarse

	interp_extrap(NLTAU,aot_C553,opth_nl_[IBIG][DL550],
			lut_.MASSCOEF_NL0[mtable_[IBIG]][DL550],yout);

	float massconc = yout*aot_C553;
	masscon_ = massconf + massconc;
//  Compute AOD for fine mode at all wavelengths
	for ( int iWav=0; iWav<NLWAV; iWav++ ) {
		interp_extrap(NLTAU,aot_F553,opth_nl_[ISMALL][DL550],
				lut_.EXTNORM_NL0[mtable_[ISMALL]][iWav],yout);

		aot_f_[iWav] = aot_F553 * yout;
	}
//  Compute AOD for coarse mode at all wavelengths
//  Compute total AOD at all wavelengths
	for ( int iWav=0; iWav<NLWAV; iWav++ ) {
		interp_extrap(NLTAU,aot_C553,opth_nl_[IBIG][DL550],
				lut_.EXTNORM_NL0[mtable_[IBIG]][iWav],yout);

		aot_c_[iWav] = aot_C553 * yout;
		aot_d_[iWav] = aot_f_[iWav] + aot_c_[iWav];
	}

	return status;
}

/**************************************************************************
 * NAME: retrieve_second()
 *
 * DESCRIPTION: This subroutine retrieves optical thickness, surface
 * reflectance and error parameters by comparing
 * observations and LUT data
 *
 * INPUTS: 	REFW466L, REFW644L, REFW212L   reflectance
 *      	OPTH       array of LUT optical thickness
 *      	RHO_S212   simulated surface reflectance
 *      	RHOSTAR    simulated TOA reflectance
 *      	yint644, yint466, slope466, slope644
 *         	characteristics of VIS/IR surface reflectance relationship
 * OUTPUTS:
 *    		ERR644	   retrieval fitting errors
 *    		AOD        retrieved optical depths
 *    		RHOSFC     retrieved surface reflectance
 *    		ETA        retrieved fine mode weighting
 *
 *************************************************************************/

int DtAlgLand::retrieve_second()
{
	int status = DTDB_SUCCESS;

//      Initialize
    float rho_star_temp[NLWAV];
    memset(rho_star_temp, 0, NLWAV*sizeof(float));
    float yout=0;

    for ( int iWav=0; iWav<NLWAV; iWav++ ) {
        rho_sfc_[iWav] = DFILL_FLOAT;
        aot_d_[iWav] = DFILL_FLOAT;
    }
//  Now loop through ETA values (kind of like ocean)
//  Average, best, QCONTROL, FILLVALUE, OUTPUT, etc
//  Compute total apparent reflectance
    for ( int iTau=0; iTau<NLTAU; iTau++ ) {
        for ( int iWav=0; iWav<NLWAV; iWav++ ) {
            rho_star_tot_[iWav][iTau] = rho_star_[ISMALL][iWav][iTau];
        }
//  Compute total surface reflectance
        rho_S212_tot_[iTau] = rho_S212_[ISMALL][iTau];
    }
//  Interpolate everything based on measured reflectance at 466.
//  Compute tau using radiance of 0.455 um and tau wav553
    for ( int iWav=0; iWav<NLWAV; iWav++ ) {

        interp_extrap(NLTAU,rfld_[D488],
        			  rho_star_tot_[DL488],
					  opth_nl_[ISMALL][iWav], aot_d_[iWav]);
    }
//  Compute apparent reflectance for all wavelengths
    for ( int iWav=0; iWav<NLWAV; iWav++ ) {

        interp_extrap(NLTAU,aot_d_[DL550],
        			  opth_nl_[ISMALL][DL550],
					  rho_star_tot_[iWav], rho_star_temp[iWav]);
    }
//  Compute surface reflectance at 2.1

    interp_extrap(NLTAU,aot_d_[DL550],
    			  opth_nl_[ISMALL][DL550],
				  rho_S212_tot_, rho_sfc_[DL2250]);

//  Errors
    errwave_[DL488] =
    		fabs(rfld_[D488] - rho_star_temp[DL488]) / rfld_[D488];
    errwave_[DL670] =
    		fabs(rfld_[D670] - rho_star_temp[DL670]) / rfld_[D670];
    errwave_[DL2250] =
    		fabs(rfld_[D2250] - rho_star_temp[DL2250]) / rfld_[D2250];
	err644_ = errwave_[DL670];

    rho_sfc_[DL550] = DFILL_FLOAT;
    rho_sfc_[DL670] = yint_644_ + slope_644_*rho_sfc_[DL2250];
    rho_sfc_[DL488] = yint_466_ + slope_466_*rho_sfc_[DL670];
//  Compute Mass Concentration
//  Continental only
    interp_extrap(NLTAU,aot_d_[DL550],
    			  opth_nl_[ISMALL][DL550],
    			  lut_.MASSCOEF_NL0[mtable_[ISMALL]][DL550],yout);
    masscon_ = yout*aot_d_[DL550];

	return status;
}

/**************************************************************************
 * NAME: simulate_toa()
 *
 * DESCRIPTION: This subroutine simulates TOA reflectance
 *
 *************************************************************************/

int DtAlgLand::simulate_toa()
{
	int status = DTDB_SUCCESS;

// 	At each tau index, calculate surface reflectance at 2.1 and
//	apparent reflectance at all wavelengths
    for ( int iSize=0; iSize<NLSIZE; iSize++ ) {
        for ( int iTau=0; iTau<NLTAU; iTau++ ) {
        	rho_S212_[iSize][iTau] =
                (int_nl_[iSize][DL2250][iTau]-rfld_[D2250]) /
                (sbar_nl_[iSize][DL2250][iTau] *
                (int_nl_[iSize][DL2250][iTau]-rfld_[D2250]) -
                fdt_nl_[iSize][DL2250][iTau]);
            if (rho_S212_[iSize][iTau] > 1.0) {
            	rho_S212_[iSize][iTau] = DFILL_FLOAT;
            }
//  Estimate surface reflectance
//  Use Red/IR and Blue/Red surface correlations
//  For simple ratios, overwrite (or comment out all lines above)
//  Current VIIRS slope_466 = ratio of M3(.48) / M5(0.67)
//  slope_644 =ratio of M5(0.67) / M11 (2.25)
            const float slpc = 1.0;
            slope_644_ = 0.559*slpc;
            yint_644_ = 0.0;
            slope_466_ = 0.645/slpc;
            yint_466_ = 0.0;
            rho_S644_ = slope_644_*rho_S212_[iSize][iTau] + yint_644_;
            rho_S466_ = slope_466_*rho_S644_ + yint_466_;
//  Compute model differentiated apparent reflectance
            rho_star_[iSize][DL488][iTau] =
                int_nl_[iSize][DL488][iTau] +
                fdt_nl_[iSize][DL488][iTau] * rho_S466_ /
                (1 - sbar_nl_[iSize][DL488][iTau] * rho_S466_);

            rho_star_[iSize][DL550][iTau] = 0.0;
            rho_star_[iSize][DL670][iTau] =
                int_nl_[iSize][DL670][iTau] +
                fdt_nl_[iSize][DL670][iTau] * rho_S644_ /
                (1 - sbar_nl_[iSize][DL670][iTau] * rho_S644_);

            rho_star_[iSize][DL2250][iTau] =
                int_nl_[iSize][DL2250][iTau] +
                fdt_nl_[iSize][DL2250][iTau] *
                rho_S212_[iSize][iTau] /
                (1 - sbar_nl_[iSize][DL2250][iTau] * rho_S212_[iSize][iTau]);
        }
    }

	return status;
}

/**************************************************************************
 * NAME: assign_quality()
 *
 * DESCRIPTION: This subroutine assigns quality variable to be
 * written to output file.
 *
 *************************************************************************/

int DtAlgLand::assign_quality()
{
	int status = DTDB_SUCCESS;

	short quality_land = 0;
    short qa_flag[19];
    memset( qa_flag, 0, sizeof(qa_flag));

	quality_flag_for_joint_[0] = 0;
	quality_flag_for_joint_[1] = 0;
    if (aot_d_[DL550] < -0.10) {
    	error_= 5;
    }
    if (aot_d_[DL550] > 5.00) {
    	error_= 6;
    }
    if (iproc_== 0) {
    	error_= 4;
    }
    if (error_ > 0) {
        quality_flag_for_retr_ = 11;
        fail_ret_+= 1;
        qa_flag[6] = 0;
        qa_flag[7] = 0;
        qa_flag[8] = 0;
        qa_flag[9] = 0;
// set it to 12 to indicate No retrievals
        qa_flag[10] = quality_flag_for_retr_;
        qa_flag[11] = error_;
        qa_flag[14] = 0;
        qcontrol_special_ = 0;
        quality_flag_for_joint_[0] = qa_flag[7];
// adding 10 to set Non ret. flag so that can tell them apart from Ret. Flag
        quality_flag_for_joint_[1] = error_+20;
// fill values
        fill_values();
    }
    else if (error_== 0) {
//   Intilize the quality_flag_for_retr =0  to indicate good quality
        quality_flag_for_retr_ = 0;
// Store SDS arrays for population to HDF file
        success_ret_ += 1;
        qa_flag[6] = 1;
        qa_flag[7] = 3;
// quality_land ==0 is if there is one single pixel of water
        if (quality_land == 0){
            qa_flag[7] = 0;
            quality_flag_for_retr_ = 2;
        }
// If Fitting error is greater than = 0.25 quality is bad
        if (err644_> 0.25) {
            qa_flag[7]  = 0;
            quality_flag_for_retr_ = 4;
        }
// Set Mass concentration & Fine optical depth to zero if tau is -ve.
        if (aot_d_[DL550] < 0.0) {
            qcontrol_special_ = 2;
            quality_flag_for_retr_ = 5;
        }
// If optical depth negative then set angs_coeff_land to fill value
// do not set the flag to zero
		if (!((angstrom_> -1.00) && (angstrom_<= 5.0))) {
			qa_flag[7] = 0;
			quality_flag_for_retr_ = 9;
		}

        float number_pixels = 1;
        int num_Q_pixel1 = (int)(number_pixels *.03);
		int num_Q_pixel2 = (int)(number_pixels *.05);
		int num_Q_pixel3 = (int)(number_pixels *.08);
		int num_Q_pixel4 = (int)(number_pixels *.12);
// set quality for the number of pixels
// 10%
// if( num_pixels_used_ >= 12 && num_pixels_used_ <=20){
        if((num_pixels_used_ > num_Q_pixel1) && (num_pixels_used_ <= num_Q_pixel2)) {
            qa_flag[7]=0;
            quality_flag_for_retr_=6;
        }
// 10-15%
// if( num_pixels_used_ >= 21 && num_pixels_used_ <=30) {
        if((num_pixels_used_ > num_Q_pixel2) && (num_pixels_used_ <= num_Q_pixel3)) {
            qa_flag[7]=1;
            quality_flag_for_retr_=7;
        }
// 15-25%
// if( num_pixels_used_ >=31 && num_pixels_used_ <=50){
        if((num_pixels_used_ > num_Q_pixel3) &&
           (num_pixels_used_ <= num_Q_pixel4)) {
            qa_flag[7]=2;
            quality_flag_for_retr_=8;
        }
// Above 25%
// if( num_pixels_used_ >=51)qa_flag[8]=3
        if( num_pixels_used_ > num_Q_pixel4) {
        	qa_flag[7]=3;
        }
        if (return_quality_cirrus_ == 0){
            quality_flag_for_retr_=3;
            qa_flag[7]=0;
        }

//If Procedure is 2 set quality to report optical depths at 0.47 & 55 only
        if (iproc_>1) {
            qcontrol_special_ = 1;
            quality_flag_for_retr_ = 1;
// If Procedure is 2  Quality is zero
            qa_flag[7] = 0;
        }
// set quality flag 9 & 10 so that Level 3 uses the qulity flag to average
        qa_flag[8] = qa_flag[6];
        qa_flag[9] = qa_flag[7];
// report eta only when optical depth is < 0.2
        if (aot_d_[DL550] < 0.2) {
            quality_flag_for_retr_ = 10;
//            sds_dust_weighting_ = DFILL_FLOAT;
        }
        qa_flag[10] = quality_flag_for_retr_;
        qa_flag[9] = error_;
        quality_flag_for_joint_[0] = qa_flag[7];
        quality_flag_for_joint_[1] = quality_flag_for_retr_;

//  Store QA flags into Quality_Assurance_Land array according to the order
//  of bits in MODIS atmosphere QA plan
		short qa_temp = 0;
		set_byte(qa_flag[6], 0, qa_temp);
		set_byte(qa_flag[7], 1, qa_temp);
		set_byte(qa_flag[8], 4, qa_temp);
		set_byte(qa_flag[9], 5, qa_temp);
		sds_qcontrol_[0] = qa_temp;
		qa_temp=0;
		qa_flag[12] = 0;
		qa_flag[13] = 0;
		set_byte(qa_flag[10], 0, qa_temp);
		set_byte(qa_flag[11], 4, qa_temp);
		sds_qcontrol_[1] = qa_temp;
		qa_temp = 0;
		qa_flag[16] = 3;
		qa_flag[17] = 3;
		set_byte(qa_flag[14], 0, qa_temp);
		set_byte(qa_flag[15], 2, qa_temp);
		set_byte(qa_flag[16], 4, qa_temp);
		set_byte(qa_flag[17], 6, qa_temp);
		sds_qcontrol_[2] = qa_temp;
		qa_temp=0;
		qa_flag[18] = 1;
		set_byte(qa_flag[18], 0, qa_temp);
		sds_qcontrol_[3] = qa_temp;
		sds_qcontrol_[4] = 0;
	//  Re-initialized working variables
		for (int i=0; i<19; i++) {
			qa_flag[i]=0;
		}
    }

    return status;
}


/**************************************************************************
 * NAME: store_output()
 *
 * DESCRIPTION: This subroutine stores all the output variables to be
 * written to output file.
 *
 *************************************************************************/

int DtAlgLand::store_output()
{
	int status = DTDB_SUCCESS;

	// Set Mass concentration & Fine optical depth to zero if tau is -ve.
	if (aot_d_[DL550] < 0.0) {
		masscon_= 0.0;
		for (int iWav=0; iWav<NLWAV; iWav++) {
            aot_f_[iWav] = 0.00;
            aot_c_[iWav] = 0.00;
		}
	}
// If -0.10 <= aot <= -0.05 then set aot_d_[all waves] = -0.05 and Quality = zero
	if ((aot_d_[DL550] >= -0.10) && (aot_d_[DL550] <= -0.05)) {
		if (iproc_ > 1) {
			aot_d_[DL550] = -0.05;
			aot_d_[DL488] = -0.05;
		}
		else {
// Set aot to -0.05, aotFine to 0.00
			for (int iWav=0; iWav<NLWAV; iWav++) {
				aot_d_[iWav] = -0.05;
			}
		}
	}
// If optical depth negative then set angs_coeff_land to fill value
// do not set the flag to zero
	if(qcontrol_special_ == 2){
		sds_angs_coeff_ = DFILL_FLOAT;
	}
	else {
		if ((angstrom_> -1.00) && (angstrom_<= 5.0)) {
			sds_angs_coeff_ = angstrom_;
		}
		else {
			sds_angs_coeff_ = DFILL_FLOAT;
		}
	}
	if (iaer_ >= 0) {
		sds_aerosol_type_ = iaer_;
	}
	if ((scatter_angle_ >= -180) && (scatter_angle_ <= 180)) {
		sds_scat_angle_ = scatter_angle_;
	}
	float aotmax = 5.0;
	if ((aot_d_[DL488] >= -0.1) && (aot_d_[DL488] <= aotmax)) {
		sds_tau_corrected_[DL488] = aot_d_[DL488];
	}
	if ((aot_d_[DL550] >= -0.1) && (aot_d_[DL550] <= aotmax)) {
		sds_tau_corrected_[DL550] = aot_d_[DL550];
	}
	if ((aot_d_[DL670] >= -0.1) && (aot_d_[DL670] <= aotmax)) {
		sds_tau_corrected_[DL670] = aot_d_[DL670];
	}
	if ((aot_d_[DL2250] >= -0.1) && (aot_d_[DL2250] <= aotmax)) {
		sds_tau_corrected_[DL2250] = aot_d_[DL2250];
	}
    sds_tau_small_[DL488] = aot_f_[DL488];
    sds_tau_small_[DL550] = aot_f_[DL550];
    sds_tau_small_[DL670] = aot_f_[DL670];
    sds_tau_small_[DL2250] = aot_f_[DL2250];
    sds_tau_big_[DL488] = aot_c_[DL488];
    sds_tau_big_[DL550] = aot_c_[DL550];
    sds_tau_big_[DL670] = aot_c_[DL670];
    sds_tau_big_[DL2250] = aot_c_[DL2250];
	sds_fitting_error_ = err644_;
	if (cloud_fraction_ >= 0) {
		sds_cloud_fraction_ = cloud_fraction_;
	}
	if (eta_ >= 0) {
		sds_dust_weighting_ = eta_;
	}
	for (int iWav=0; iWav<NLWAV; iWav++) {
		if (good_pixels_[iWav] >= 0) {
			sds_numpixels_[iWav] = good_pixels_[iWav];
		}
	}
	if (masscon_ >= 0) {
		sds_mass_conc_ = masscon_;
	}
	if (ndvi_ >= -1.0 && ndvi_ <= 1.0) {
		sds_ndvi_ = ndvi_;
	}
	if ((rfld_[D488] >= 0.0) && (rfld_[D488] <= 1.2)) {
		sds_refl_[DL488] = rfld_[D488];
		sds_refl_std_[DL488] = sdev_[DL488];
        sds_surface_reflectance_[DL488] = rho_sfc_[DL488];
	}
	if ((rfld_[D550] >= 0.0) && (rfld_[D550] <= 1.2)) {
		sds_refl_[DL550] = rfld_[D550];
		sds_refl_std_[DL550] = sdev_[DL550];
	}
	if ((rfld_[D670] >= 0.0) && (rfld_[D670] <= 1.2)) {
		sds_refl_[DL670] = rfld_[D670];
		sds_refl_std_[DL670] = sdev_[DL670];
        sds_surface_reflectance_[DL670] = rho_sfc_[DL670];
	}
	if ((rfld_[D865] >= 0.0) && (rfld_[D865] <= 1.2)) {
		sds_refl_[D865] = rfld_[D865];
		sds_refl_std_[D865] = sdev_[D865];
	}
	if ((rfld_[D1240] >= 0.0) && (rfld_[D1240] <= 1.2)) {
		sds_refl_[D1240] = rfld_[D1240];
		sds_refl_std_[D1240] = sdev_[D1240];
	}
	if ((rfld_[D1610] >= 0.0) && (rfld_[D1610] <= 1.2)) {
		sds_refl_[D1610] = rfld_[D1610];
		sds_refl_std_[D1610] = sdev_[D1610];
	}
	if ((rfld_[D2250] >= 0.0) && (rfld_[D2250] <= 1.2)) {
		sds_refl_[DL2250] = rfld_[D2250];
		sds_refl_std_[DL2250] = sdev_[DL2250];
        sds_surface_reflectance_[DL2250] = rho_sfc_[DL2250];
	}

	iaer_=DFILL_SHORT;
    iproc_=DFILL_SHORT;
    error_=DFILL_SHORT;
    eta_ = DFILL_FLOAT;
    err644_ = DFILL_FLOAT;
    masscon_ = DFILL_FLOAT;
    angstrom_ = DFILL_FLOAT;
    ndvi_ = DFILL_FLOAT;
    for (int iWav=0; iWav<NLWAV; iWav++) {
        aot_d_[iWav] = DFILL_FLOAT;
        aot_f_[iWav] = DFILL_FLOAT;
        aot_c_[iWav] = DFILL_FLOAT;
        rho_sfc_[iWav] = DFILL_FLOAT;
    }
    for (int iWav=0; iWav<NLWAV; iWav++) {
        rfld_[iWav] = DFILL_FLOAT;
        sdev_[iWav] = DFILL_FLOAT;
    }

    return status;
}

/**
 * NAME: fill_values()
 *
 * This subroutine assigns fill values to output data where appropriate.
 */

int DtAlgLand::fill_values()
{
	int status = DTDB_SUCCESS;

    sds_aerosol_type_=8;
    sds_scat_angle_= DFILL_FLOAT;
    sds_cloud_fraction_= cloud_fraction_;
	sds_dust_weighting_ = DFILL_FLOAT;
    sds_fitting_error_= DFILL_FLOAT;
    sds_mass_conc_ = DFILL_FLOAT;
	sds_angs_coeff_ = DFILL_FLOAT;
	for (int iWav=0; iWav<NLWAV; iWav++) {
        sds_tau_corrected_[iWav]= DFILL_FLOAT;
        sds_tau_small_[iWav]= DFILL_FLOAT;
        sds_tau_big_[iWav]= DFILL_FLOAT;
        sds_surface_reflectance_[iWav] = DFILL_FLOAT;
		sds_numpixels_[iWav]=DFILL_SHORT;
		sds_refl_[iWav] = DFILL_FLOAT;
		sds_refl_std_[iWav] = DFILL_FLOAT;
    }

    return status;
}





