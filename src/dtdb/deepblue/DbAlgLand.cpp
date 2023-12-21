/*******************************************************************************
 *
 * NAME: DbAlgLand.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute processings for a given DbAlgLand object class.
 *
 *  Created on: October 13, 2018
 *      Author: Sam Anderson, DB
 *
 *  Modified:
 *
 *******************************************************************************/

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <fstream>

#include <DDAlgorithm.h>
#include <DDProcess.h>
#include <DDOptions.h>
#include "deepblue/DbLutNetcdf.h"
#include "deepblue/DbMask.h"
#include "deepblue/DbAlgorithm.h"

using namespace std;

extern"C" {
void deepblue_initialize_( char* config_file, char* nc4_file,
		int* lines, int* pixels, float* lat, float* lon,
		int* year, int* month, int* day);
void deepblue_cleanup_();
void find_v_viirs_( float* realbuf, float* tmpvg, float* outbuf, int* flags,
        float* elev, int* ls_flag, float* windsp, float* wv);
};

/**************************************************************************
 * NAME: DbAlgLand()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DbAlgLand::DbAlgLand()
{
}

/**************************************************************************
 * NAME: ~DbAlgLand()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DbAlgLand::~DbAlgLand()
{
    deepblue_cleanup_();

    delete msr_lut_;
    delete vsr_lut_;
    delete scl_lut_;
    delete sp_lut_;
    delete mt_lut_;
    delete gz_lut_;

    delete cm_;
    delete smoke_;
    delete ha_smoke_;
    delete pyrocb_;
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * processing operations.
 *
 *************************************************************************/

int DbAlgLand::initialize( map<string, ddata*> imap )
{
	int status = DTDB_SUCCESS;

    status = DbAlgorithm::initialize( imap );
    if (status != DTDB_SUCCESS) {
        std::cerr << "DbAlgLand:: Base class initialization failure" << std::endl;
        return status;
    }
    status = initialize_LUT_data( imap );
// Compute various arrays, including dstar, ndvi, sca, etc.

// calculate values needed in table interpolation
    for (int j = 0; j < 7; j++) {
        for (int k = j; k <= j + 3; k++) {
            float xdenom = 1.0;
            for (int i = j; i <= j + 3; i++) {
                if (i != k) xdenom = xdenom*(xzlog[k] - xzlog[i]);
            }
            densol_[j][k-j] = xdenom;
        }
    }
    for (int j = 0; j < 5; j++) {
        for (int k = j; k <= j + 3; k++) {
            float xdenom = 1.0;
            for (int i = j; i <= j + 3; i++) {
                if (i != k) xdenom = xdenom*(xlog[k] - xlog[i]);
            }
            denscn_[j][k-j] = xdenom;
        }
    }

    cm_ = new DbCloudMaskLand(this);
    smoke_ = new DbSmokeMask(this);
    pyrocb_ = new DbPyrocbMask(this);
    ha_smoke_ = new DbHighAltSmokeMask(this);

    return status;
}

int DbAlgLand::initialize_LUT_data( map<string, ddata*> imap )
{
	int status = DTDB_SUCCESS;

	DbAlgorithm::initialize_LUT_data( imap);
    string config_file = static_cast<ddstr*>(imap["config_file"])->str;
    int num_lines = static_cast<ddval<int>*>(imap["num_lines"])->val;
    int num_pixels = static_cast<ddval<int>*>(imap["num_pixels"])->val;
    int start_year = static_cast<ddval<int>*>(imap["start_year"])->val;
    int start_month = static_cast<ddval<int>*>(imap["start_month"])->val;
    int start_day = static_cast<ddval<int>*>(imap["start_day"])->val;
    int season = static_cast<ddval<int>*>(imap["season"])->val;
    ddma<float,2>* plat = static_cast<ddma<float,2>*>(imap["latitude"]);
    ddma<float,2>* plon = static_cast<ddma<float,2>*>(imap["longitude"]);

    char nc4_file[255] = "";
    string nc4_str = get_option(INPUT_NC4_LUT);
	if (nc4_str.empty()) {
		nc4_str = get_option(INPUT_DB_NC4_LUT);
	}
    nc4_str.copy(nc4_file, nc4_str.length());
    char cfile[255] = "";
    config_file.copy(cfile, config_file.length());

// Load LUTs
    deepblue_initialize_( cfile, nc4_file, &num_lines,
            &num_pixels, &plat->pts[0][0], &plon->pts[0][0],
            &start_year, &start_month, &start_day );

    DbLutNetcdf* lutgen = new DbLutNetcdf();
    msr_lut_ = new dbModisSurfReflLimited;
    status = lutgen->read_modis_surf_refl_lut(msr_lut_, &ler_start_[0],
            &ler_edge_[0], season, dateline_);
    vsr_lut_ = new dbViirsSurfReflLimited;
    status = lutgen->read_viirs_surf_refl_lut(vsr_lut_, &ler_start_[0],
            &ler_edge_[0], season, dateline_);
    scl_lut_ = new dbSurfCoeffLimited();
    status = lutgen->read_surf_coeff_lut(scl_lut_, &ler_start_[0],
            &ler_edge_[0], season, dateline_);
    sp_lut_ = new dbSurfacePressureLUT();
    status = lutgen->read_surface_pressure_lut(sp_lut_);
    mt_lut_ = new dbTablesLUT();
    status = lutgen->read_tables_lut(mt_lut_);
    gz_lut_ = new dbGeozoneLUT();
    status = lutgen->read_geozone_lut(gz_lut_);
    delete lutgen;

	return status;
}

/**************************************************************************
 * NAME: compute()
 *
 *************************************************************************/

map<string, ddata*> DbAlgLand::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap)
{
    get_inputs(start, count, imap);

    l2_flags_ |= (unsigned int) flags::LAND;

	size_t iy = start[0];
	size_t ix = start[1];

    float realbuf[26];
    float ob[21];
    float tmpvg[7];
    int   flags[4];
    int   ls_flag;
    float cl_flag = DFILL_FLOAT;

    set_fill_out();
    for (size_t i=0; i<21; i++) {
        ob[i] = DFILL_FLOAT;
    }
    mask_cm_ = DFILL_SHORT;
    short snow_1 = DFILL_SHORT;
    short snow_2 = DFILL_SHORT;

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
    if ((rfl_[(int)rhot_band::W490] < 0) || (rfl_[(int)rhot_band::W550] < 0) ||
    		(rfl_[(int)rhot_band::W670] < 0) || (rfl_[(int)rhot_band::W865] < 0)) {
//        std::cerr << "DbAlgOcean:: Invalid reflectances at "<< iy << ":" << ix << std::endl;
    	l2_flags_ |= (unsigned int) flags::BOWTIEDEL;
    	return set_fills();
    }

    compute_glint_refl(glint_angle_);
    compute_scatter_angle(scatter_angle_);

    compute_dstar( iy, ix );
    if (sr670_<=DFILL_TEST) {
    	return set_fills();
    }

    cm_->compute_1( iy, ix, mask_cm_, snow_1, snow_2 );

// at this point convert from VIIRS to IOF units
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
    NC_[NC412] = rfl_[(int)rhot_band::W410];
    NC_[NC488] = rfl_[(int)rhot_band::W490];
    NC_[NC670] = rfl_[(int)rhot_band::W670];

    if (bgascorrect_) {
		compute_gas_correction();
		for (int ib = 0; ib < NTWL; ib++) {
			rfl_[ib] *= gasc_[ib];
		}
    }

    compute_sr( iy, ix );
    if (sr412_<=DFILL_TEST || sr488_<=DFILL_TEST || sr670_<=DFILL_TEST) {
    	return set_fills();
    }
    compute_ler( iy, ix );
    if (ler412_<=DFILL_TEST || ler488_<=DFILL_TEST || ler670_<=DFILL_TEST) {
    	return set_fills();
    }

    cm_->compute_2( iy, ix, mask_cm_, snow_2 );

    smoke_->compute( iy, ix, mask_smoke_ );
    pyrocb_->compute( iy, ix, mask_pyrocb_ );
    ha_smoke_->compute( iy, ix, mask_ha_smoke_ );

    if (mask_smoke_ == 1 || mask_pyrocb_ == 1 || mask_ha_smoke_ == 1) {
        mask_cm_ = 0;
    }

    if (cloud_mask_ == DFILL_UBYTE) {
    	cloud_mask_ = mask_cm_;
    } else {
    	mask_cm_ = cloud_mask_;
    }
    if (!bmaskcloud_) {
        mask_cm_ = 0;
    }
	if (bmaskcloud_ && cloud_mask_) {
		l2_flags_ |= (unsigned int) flags::CLDICE;
		return set_fills();
	}

    if (mask_cm_ == 0 && rfl_[(int)rhot_band::W410] > 0) {

        tmpvg[0] = rfl_[(int)rhot_band::W670];
        tmpvg[1] = rfl_[(int)rhot_band::W490];
        tmpvg[2] = rfl_[(int)rhot_band::W410];
        tmpvg[3] = rfl_[(int)rhot_band::W2250];
        tmpvg[4] = rfl_[(int)rhot_band::W865];
        tmpvg[5] = rfl_[(int)rhot_band::W1240];
        tmpvg[6] = rfl_[(int)rhot_band::W2250];

        realbuf[0]  = lat_;
        realbuf[1]  = lon_;
        realbuf[2]  = solz_;
        realbuf[3]  = senz_;
        realbuf[4]  = raa_;
        realbuf[5]  = ler670_;
        realbuf[6]  = NC_[NC670];   //no
        realbuf[7]  = rfl_[(int)rhot_band::W490];   //new
        realbuf[8]  = NC_[NC488];   //no
        realbuf[9]  = rfl_[(int)rhot_band::W410];   //new
        realbuf[10] = NC_[NC412];   //no
        realbuf[11] = rfl_[(int)rhot_band::W2250];  //old
        realbuf[12] = DFILL_FLOAT;
        realbuf[13] = ndvi_;
        realbuf[14] = rfl_[(int)rhot_band::W12000];
        realbuf[15] = dstar_;
        realbuf[16] = cl_flag;
        realbuf[17] = qdf412_;
        realbuf[18] = qdf488_;
        realbuf[19] = qdf670_;
        realbuf[20] = ps_/1013.25;
        realbuf[21] = ler412_;
        realbuf[22] = ler488_;
        realbuf[23] = DFILL_FLOAT;
        realbuf[24] = DFILL_FLOAT;
        realbuf[25] = DFILL_FLOAT;

        find_v_viirs_( &realbuf[0], &tmpvg[0], &ob[0], &flags[0],
                &height_, &ls_flag, &ws_, &pwv_);
    }

    mask_cm_ = (mask_cm_ >= 0) ? mask_cm_ : DFILL_SHORT;

    float aotmax = (mask_cm_==0) ? 30.0 : 5.0;
    lOut_.aot[OL412] = (ob[0]  > 0.0 && ob[0] <= aotmax) ? ob[0] : DFILL_FLOAT;
    lOut_.aot[OL488] = (ob[1]  > 0.0 && ob[1] <= aotmax) ? ob[1] : DFILL_FLOAT;
    lOut_.aot[OL670] = (ob[2]  > 0.0 && ob[2] <= aotmax) ? ob[2] : DFILL_FLOAT;
    lOut_.ssa[OL412] = (ob[3]  > 0.0 && ob[3] <= aotmax) ? ob[3] : DFILL_FLOAT;
    lOut_.ssa[OL488] = (ob[4]  > 0.0 && ob[4] <= aotmax) ? ob[4] : DFILL_FLOAT;
    lOut_.ssa[OL670] = (ob[5]  > 0.0 && ob[5] <= aotmax) ? ob[5] : DFILL_FLOAT;
    lOut_.aot550    = (ob[6]  > DFILL_TEST) ? ob[6] : DFILL_FLOAT;
    lOut_.ae        = (ob[7]  > DFILL_TEST) ? ob[7] : DFILL_FLOAT;
    lOut_.sr[OL412]  = (ob[8]  > 0.0) ? ob[8]/100.0 : DFILL_FLOAT;
    lOut_.sr[OL488]  = (ob[10] > 0.0) ? ob[10]/100.0 : DFILL_FLOAT;
    lOut_.sr[OL670]  = (ob[11] > 0.0) ? ob[11]/100.0 : DFILL_FLOAT;
    lOut_.sfc_type  = (ob[14] > DFILL_TEST) ? ob[14] : DFILL_SHORT;
    lOut_.alg_flag  = (ob[17] > DFILL_TEST) ? ob[17] : DFILL_SHORT;
    lOut_.ndvi  = (ndvi_ > DFILL_TEST) ? ndvi_ : DFILL_FLOAT;
    lOut_.aerosol_type = 8;

    if (lOut_.aot550 > DFILL_TEST) {
        lOut_.aerosol_type = 5; // mixed, default
    }
    if (lOut_.ae > 1.2 && lOut_.aot550 > 0.4) {
        lOut_.aerosol_type = 4; // non-smoke fine mode
    }
    if (mask_smoke_ == 1  && lOut_.aot550 > DFILL_TEST) {
        lOut_.aerosol_type = 1; // smoke
    }
    if (mask_ha_smoke_ == 1  && lOut_.aot550 > DFILL_TEST) {
        lOut_.aerosol_type = 2; // high altitude smoke
    }
    if (mask_pyrocb_ == 1  && lOut_.aot550 > DFILL_TEST) {
        lOut_.aerosol_type = 3; // pyrocumulonimbus clouds
    }
    if (lOut_.aot550 > DFILL_TEST && lOut_.aot550 < 0.2) {
        lOut_.aerosol_type = 6; // background
    }
    if (dstar_ > 1.1 && lOut_.aot550 > DFILL_TEST) {
        lOut_.aerosol_type = 0; // dust
    }
    if (lOut_.ae < 0.1 && ob[20] < 0.78 && lOut_.aot550 > 0.4) {
        lOut_.aerosol_type = 0; // dust
    }
/*
    if (lOut_.ae < 0.5 && gzone > 0 && gzone < 12 && lOut_.aot550 > 0.2) {
        lOut_.aerosol_type = 0; // dust
    }
*/

// write to generic outputs

    qual_flag_ = (lOut_.aot550 < DFILL_TEST) ? 0 : 3;
    aerosol_type_ = lOut_.aerosol_type;
    error_flag_ = lOut_.alg_flag;

    scatter_ang_ = scatter_angle_;
    glint_ang_ = glint_angle_;
    sse_ = DFILL_FLOAT;
    fmf_ = DFILL_FLOAT;
    aot_550_ = lOut_.aot550;
    ae1_ = lOut_.ae;
    ae2_ = DFILL_FLOAT;
    ndv_ = lOut_.ndvi;
    chlor_ = DFILL_FLOAT;
    for ( int ib=0; ib<NOWL+1; ib++ ) {
    	aot_[ib] = DFILL_FLOAT;
    }
    aot_[(size_t)aot_band::W410] = lOut_.aot[(size_t)srf_band::W410];
    aot_[(size_t)aot_band::W490] = lOut_.aot[(size_t)srf_band::W490];
    aot_[(size_t)aot_band::W550] = lOut_.aot550;
    aot_[(size_t)aot_band::W670] = lOut_.aot[(size_t)srf_band::W670];
    for ( int ib=0; ib<NLWL; ib++ ) {
        sr_[ib] = lOut_.sr[ib];
        ssa_[ib] = lOut_.ssa[ib];;
    }
    sr_[(size_t)srf_band::W2250] = DFILL_FLOAT;
    ssa_[(size_t)srf_band::W2250] = DFILL_FLOAT;
    for (int ib = 0; ib < DB_RFL_BANDS; ib++) {
        rfl_[ib] *= (M_PI/cossza);
    }

    size_t cgood = 0;
 	vector<float> tba, yba;
	if (aot_[(size_t)aot_band::W410]>0.0) {
		tba.push_back(410.0);
		yba.push_back(aot_[(size_t)aot_band::W410]);
		cgood++;
	}
	if (aot_[(size_t)aot_band::W490]>0.0) {
		tba.push_back(490.0);
		yba.push_back(aot_[(size_t)aot_band::W490]);
		cgood++;
	}
	if (aot_[(size_t)aot_band::W550]>0.0) {
		tba.push_back(550.0);
		yba.push_back(aot_[(size_t)aot_band::W550]);
		cgood++;
	}
	if (aot_[(size_t)aot_band::W670]>0.0) {
		tba.push_back(670.0);
		yba.push_back(aot_[(size_t)aot_band::W670]);
		cgood++;
	}
	if (cgood >= 3) {
		using boost::math::interpolators::barycentric_rational;
		barycentric_rational<float> interp(move(tba), move(yba), 2);
		if (aot_[(size_t)aot_band::W670] < 0) {
			aot_[(size_t)aot_band::W670] = interp(670.0);
		}
		aot_[(size_t)aot_band::W865] = interp(865.0);
		aot_[(size_t)aot_band::W1240] = interp(1240.0);
		aot_[(size_t)aot_band::W1610] = interp(1610.0);
		aot_[(size_t)aot_band::W2250] = interp(2250.0);
	}
    return set_outputs();
}

/**************************************************************************
* NAME: compute_dstar()
*
* DESCRIPTION: Compute D* and LER670 required for initial cloud filter
*
*************************************************************************/

int DbAlgLand::compute_dstar( const size_t iy, const size_t ix )
{
    int status = DTDB_SUCCESS;

    if (lat_ < DFILL_TEST || lon_ < DFILL_TEST) {
        return status;
    }
    float psi = acos( cos(solz_*DEGtoRAD)*cos(senz_*DEGtoRAD) -
                      sin(solz_*DEGtoRAD)*sin(senz_*DEGtoRAD) *
                      cos(raa_*DEGtoRAD) );
    sca_  = 180.0 - psi / DEGtoRAD;
    gla_  = psi / DEGtoRAD;
    amf_  = 1.0/cos(solz_*DEGtoRAD) + 1.0/cos(senz_*DEGtoRAD);
    if (rfl_[(int)rhot_band::W865] > DFILL_TEST
             && rfl_[(int)rhot_band::W670] > DFILL_TEST) {
         ndvi_  = (rfl_[(int)rhot_band::W865] - rfl_[(int)rhot_band::W670])
                 / (rfl_[(int)rhot_band::W865] + rfl_[(int)rhot_band::W670]);
    } else {
         ndvi_  = DFILL_FLOAT;
    }
    btd8_  = rfl_[(int)rhot_band::W8550] - rfl_[(int)rhot_band::W11000];
    btd11_  = rfl_[(int)rhot_band::W11000] -rfl_[(int)rhot_band::W12000];
//calculate D* parameter. Skip bad/cloudy pixels
    float A = -0.05;
    float B = 10.0;
    float ratio = (btd11_  - A) / (btd8_  - B);
    dstar_  = (ratio >= 5.0) ? DFILL_FLOAT : exp(ratio);

    if (rfl_[(int)rhot_band::W670] <= 0.0 ||
        rfl_[(int)rhot_band::W490] <= 0.0 ||
        rfl_[(int)rhot_band::W410] <= 0.0) {
        sr670_  = DFILL_FLOAT;
        return status;
    }
    int ilat = floor((lat_ + 90.0)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = floor((lon_ + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    int sx;
    if (dateline_ == 0 || ilon > ler_start_[0]) {
        sx = ilon - ler_start_[0];
    } else {
        sx = ilon + dateline_;
    }
    int sy = ilat - ler_start_[1];
    int nidx;
    if (ndvi_  < NDVI1_CUTOFF) {
        nidx = 0;
    } else if (ndvi_  < NDVI2_CUTOFF) {
        nidx = 1;
    } else {
        nidx = 2;
    }
    float c670[4], sr670;
    if (raa_ < 90.0) {
        for (int i = 0; i < 4; i++) {
            c670[i] = scl_lut_->SC650_FWD_L[nidx][i][sy][sx];
        }
        sr670 = vsr_lut_->SR670_ALL_L[sy][sx]/100.0;
    } else {
        for (int i = 0; i < 4; i++) {
            c670[i] = scl_lut_->SC650_ALL_L[nidx][i][sy][sx];
        }
        sr670 = vsr_lut_->SR670_ALL_L[sy][sx]/100.0;
    }
    float sfref670 = DFILL_FLOAT;
    if (c670[0]>0 || c670[1]>0 || c670[2]>0 || c670[3]>0) {
        sfref670 =
        (c670[0] + sca_*(c670[1] + sca_*(c670[2] + sca_*c670[3])))/100.0;
    }
    sfref670 = (sfref670 < 0.0) ? sr670 : sfref670;
    sr670_  = sfref670;

    return status;
}

/**************************************************************************
* NAME: compute_sr()
*
* DESCRIPTION: Retrieve and compute Lambertian Equivalent Radiance
*
*************************************************************************/

int DbAlgLand::compute_sr( const size_t iy, const size_t ix ) {

    int status = DTDB_SUCCESS;

    if (lat_ <= DFILL_TEST || lon_ <= DFILL_TEST) {
        return status;
    }

    if (rfl_[(int)rhot_band::W670] <= 0.0 ||
        rfl_[(int)rhot_band::W490] <= 0.0 ||
        rfl_[(int)rhot_band::W410] <= 0.0) {
        sr412_ = DFILL_FLOAT;
        sr488_ = DFILL_FLOAT;
        sr670_ = DFILL_FLOAT;
        return status;
    }
    if (rfl_[(int)rhot_band::W865] > DFILL_TEST
             && rfl_[(int)rhot_band::W670] > DFILL_TEST) {
         ndvi_ = (rfl_[(int)rhot_band::W865] - rfl_[(int)rhot_band::W670])
                 / (rfl_[(int)rhot_band::W865] + rfl_[(int)rhot_band::W670]);
    } else {
         ndvi_ = DFILL_FLOAT;
    }
    int ilat = floor((lat_ + 90.0)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = floor((lon_ + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    int sx;
    if (dateline_ == 0 || ilon > ler_start_[0]) {
        sx = ilon - ler_start_[0];
    } else {
        sx = ilon + dateline_;
    }
    int sy = ilat - ler_start_[1];
    int nidx;
    if (ndvi_ < NDVI1_CUTOFF) {
        nidx = 0;
    } else if (ndvi_ < NDVI2_CUTOFF) {
        nidx = 1;
    } else {
        nidx = 2;
    }
    float c412[4], sr412;
    float c488[4], sr488;
    float c670[4], sr670;
    sr412 = vsr_lut_->SR412_ALL_L[sy][sx];
    sr488 = vsr_lut_->SR488_ALL_L[sy][sx];
    sr670 = vsr_lut_->SR670_ALL_L[sy][sx];
    if (sr412<=DFILL_TEST || sr488<=DFILL_TEST || sr670<=DFILL_TEST) {
        sr412_ = DFILL_FLOAT;
        sr488_ = DFILL_FLOAT;
        sr670_ = DFILL_FLOAT;
        return status;
    } else {
        sr412 /= 100.0;
        sr488 /= 100.0;
        sr670 /= 100.0;
    }
    if (raa_ < 90.0) {
        for (int i = 0; i < 4; i++) {
            c412[i] = scl_lut_->SC412_FWD_L[nidx][i][sy][sx];
            c488[i] = scl_lut_->SC470_FWD_L[nidx][i][sy][sx];
            c670[i] = scl_lut_->SC650_FWD_L[nidx][i][sy][sx];
        }
    } else {
        for (int i = 0; i < 4; i++) {
            c412[i] = scl_lut_->SC412_ALL_L[nidx][i][sy][sx];
            c488[i] = scl_lut_->SC470_ALL_L[nidx][i][sy][sx];
            c670[i] = scl_lut_->SC650_ALL_L[nidx][i][sy][sx];
        }
    }
    float sfref412 = DFILL_FLOAT;
    float sfref488 = DFILL_FLOAT;
    float sfref670 = DFILL_FLOAT;
    if (c412[0]>0 || c412[1]>0 || c412[2]>0 || c412[3]>0) {
        sfref412 =
        (c412[0] + sca_*(c412[1] + sca_*(c412[2] + sca_*c412[3])))/100.0;
    }
    if (c488[0]>0 || c488[1]>0 || c488[2]>0 || c488[3]>0) {
        sfref488 =
        (c488[0] + sca_*(c488[1] + sca_*(c488[2] + sca_*c488[3])))/100.0;
    }
    if (c670[0]>0 || c670[1]>0 || c670[2]>0 || c670[3]>0) {
        sfref670 =
        (c670[0] + sca_*(c670[1] + sca_*(c670[2] + sca_*c670[3])))/100.0;
    }
    sfref412 = (sfref412 < 0.0) ? sr412 : sfref412;
    sfref488 = (sfref488 < 0.0) ? sr488 : sfref488;
    sfref670 = (sfref670 < 0.0) ? sr670 : sfref670;
    sr412_ = sfref412;
    sr488_ = sfref488;
    sr670_ = sfref670;

    return status;
}

/**************************************************************************
*NAME: compute_ler()
 *
*DESCRIPTION: Retrieve and compute Lambertian Equivalent Radiance
 *
 *************************************************************************/

int DbAlgLand::compute_ler( const size_t iy, const size_t ix )
{
    int status = DTDB_SUCCESS;

    if (rfl_[(int)rhot_band::W670] <= 0.0 ||
        rfl_[(int)rhot_band::W490] <= 0.0 ||
        rfl_[(int)rhot_band::W410] <= 0.0) {
        ler412_ = DFILL_FLOAT;
        ler488_ = DFILL_FLOAT;
        ler670_ = DFILL_FLOAT;
        qdf412_ = DFILL_FLOAT;
        qdf488_ = DFILL_FLOAT;
        qdf670_ = DFILL_FLOAT;
        return status;
    }

    int ilat = floor((lat_ + 90.0)*10.0);
    if (ilat >= 1800) ilat = 1800-1;
    if (ilat < 0) ilat = 0;
    int ilon = floor((lon_ + 180.0)*10.0);
    if (ilon >= 3600) ilon = 3600-1;
    if (ilon < 0) ilon = 0;
    ilat = floor((90.0 - lat_)*12.0);
    if (ilat >= 2160) ilat = 2160-1;
    if (ilat < 0) ilat = 0;
    ilon = floor((lon_ + 180.0)*12.0);
    if (ilon >= 4320) ilon = 4320-1;
    if (ilon < 0) ilon = 0;

    float tmp0, tmp1;
    compute_pressure(height_, tmp0, pteran_, tmp1);
    double convrt = 0.005729577951;
    ilat = 1;
    if (abs(lat_) > 15.)
        ilat = 2;
    if (abs(lat_) > 60.)
        ilat = 3;
//  unpack thir cloud top pressure, terrain height, surface
//  category, percent cloudiness and snow thickness
//  compute phase factor terms,lagrange coeffs. and index offset
//  for solar and satellite zenith angle interpolations
    float xtemp1 = solz_ / (convrt*10000.);
    float xtemp2 = sin(xtemp1);
    xtemp1 = cos(xtemp1);
    float xzlog1 = 0.0;
    if (xtemp1 > 0.)
        xzlog1 = log(1.0/xtemp1);
//  raa variables
    float ztemp1 = raa_/(convrt*10000.0);
    cphi_ = cos(ztemp1);
    cphir_ = cos(ztemp1*2.0);
//  theta variables
    float ytemp1 = senz_/(convrt*10000.0);
    float ytemp2 = sin(ytemp1);
    ytemp1 = cos(ytemp1);
    float xlog1 = 0.0;
    if (ytemp1 > 0.)
        xlog1 = log(1.0/ytemp1);
//  phase factors
    phs_ = -0.375*xtemp1*xtemp2*ytemp2;
    phsr_ = 2.0*phs_/(3.0*ytemp1*xtemp1*xtemp1);
//  set up index offset value for theta,theta0
    int indsol = 0;
    for (int i = 0; i < 10; i++) {
        indsol = i;
        if (xzlog[i] >= xzlog1) break;
    }
    indsol = indsol - 2;
//  check for range
    if (indsol < 0) indsol = 0;
    if (indsol >= 7) indsol = 7-1;
//  theta
    int indscn = 0;
    for (int i = 0; i < 8; i++) {
        indscn = i;
        if (xlog[i] >= xlog1) break;
    }
    indscn = indscn - 2;
//  check range
    if (indscn < 0) indscn = 0;
    if (indscn >= 5) indscn = 5-1;
//  compute lagrange coeffs for solar zenith angle
    float indmax = indsol + 3;
    int j = 0;
    float xnom = 1.0;
    float cthet0[4];
    for (int k = indsol; k <= indmax; k++) {
        xnom = 1.0;
        for (int i = indsol; i <= indmax; i++) {
            if (i != k) {
                xnom = (xzlog1 - xzlog[i])*xnom;
            }
        }
        cthet0[j] = xnom / densol_[indsol][j];
        j++;
    }
//  compute lagrange coeffs. for theta interp
    indmax = indscn + 3;
    j = 0;
    float ctheta[4];
    for (int k = indscn; k <= indmax; k++) {
        xnom = 1.0;
        for (int i = indscn; i <= indmax; i++) {
            if (i != k) {
                xnom = (xlog1 - xlog[i])*xnom;
            }
        }
        ctheta[j] = xnom / denscn_[indscn][j];
        j++;
    }
//  compute offset into tables for theta0, theta block
    int iofset = indsol*8 + indscn;
//  store products of coeffs
    int l = 0;
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 4; k++) {
            cofs_[l] = ctheta[k]*cthet0[i];
            l++;
        }
    }
    float sfref412 = sr412_;
    float sfref488 = sr488_;
    float sfref670 = sr670_;
    float clref = 0.80;
    float grref = sfref670;
    float qgc670 = 0.0;
    float qcc670 = 0.0;
    float qdif670 = 0.0;
    float qgc670x = 0.0;
    float qcc670x = 0.0;
    float qgc412 = 0.0;
    float qcc412 = 0.0;
    float qdif412 = 0.0;
    float qgc412x = 0.0;
    float qcc412x = 0.0;
    float qgc488 = 0.0;
    float qcc488 = 0.0;
    float qdif488 = 0.0;
    float qgc488x = 0.0;
    float qcc488x = 0.0;
//    int ip670_1[2] = { 0, 80 };
//    int ip670_2[2] = { 0, 1 };
    int ipt_1[2] = { 0, 10400 };
    int ipt_2[2] = { 0, 130 };
    int lamtb1[5] = { 0, 2080, 4160, 6240, 8320 };
    int lamtb2[5] = { 0, 26, 52, 78, 104 };
    float ez670[2], t670[2], sb670[2];
    float ez412[2], t412[2], sb412[2];
    float ez488[2], t488[2], sb488[2];
//   set terrain and cloud weighting fractions
    float pwtlo = (pteran_ - .4) / .6;
    float pwthi = (pcloud_ - .4) / .6;
//  convert 670, 412, and 488 nm channel n values to albedos
    float alb670 = NC_[NC670];
    float alb412 = NC_[NC412];
    float alb488 = NC_[NC488];
    float qsh670 = DFILL_FLOAT;
    float qsh412 = DFILL_FLOAT;
    float qsh488 = DFILL_FLOAT;
    float qsl670 = DFILL_FLOAT;
    float qsl412 = DFILL_FLOAT;
    float qsl488 = DFILL_FLOAT;
//  calculate Rayleigh correction
    for (int ip = 0; ip < 2; ip++) {
//  compute table indices for 670 channel
        int i670_1 = ipt_1[ip] + lamtb1[3] + iofset;
        int i670_2 = ipt_2[ip] + lamtb2[3];
//  compute table indices for 412 channel
        int i412_1 = ipt_1[ip] + lamtb1[0] + iofset;
        int i412_2 = ipt_2[ip] + lamtb2[0];
//  compute table indices for 488 channel
        int i488_1 = ipt_1[ip] + lamtb1[2] + iofset;
        int i488_2 = ipt_2[ip] + lamtb2[2];
//  perform interpolations
//  for 670
        interx(i670_1, i670_2, rhot_band::W670, ez670[ip], t670[ip], sb670[ip]);
//  for 412
        interx(i412_1, i412_2, rhot_band::W410, ez412[ip], t412[ip], sb412[ip]);
//  for 488
        interx(i488_1, i488_2, rhot_band::W490, ez488[ip], t488[ip], sb488[ip]);
//  determine calculated q values
//  670 --
        qgc670 = ez670[ip] + grref*t670[ip] / (1 - grref*sb670[ip]);
        qcc670 = ez670[ip] + clref*t670[ip] / (1 - clref*sb670[ip]);
//  412 --
        qgc412 = ez412[ip] + sfref412*t412[ip] / (1 - sfref412*sb412[ip]);
        qcc412 = ez412[ip] + clref*t412[ip] / (1 - clref*sb412[ip]);
//  488 --
        qgc488 = ez488[ip] + sfref488*t488[ip] / (1 - sfref488*sb488[ip]);
        qcc488 = ez488[ip] + clref*t488[ip] / (1 - clref*sb488[ip]);
        if (ip == 0) {
// for ip = 0 (p = 1 atm) save q values
// and denominator used in ref calculation
            qsl670 = qgc670;
            qsh670 = qcc670;
            qsl412 = qgc412;
            qsh412 = qcc412;
            qsl488 = qgc488;
            qsh488 = qcc488;
        } else {
// for ip = 2 (p = 0.4 atm) interpolate
// 0.4 and 1 atm q values to obtain final q values
            qgc670 = pwtlo*qsl670 + qgc670*(1. - pwtlo);
            qcc670 = pwthi*qsh670 + qcc670*(1. - pwthi);
            qdif670 = qsl670 - qgc670;
// 412 --
            qgc412 = pwtlo*qsl412 + qgc412*(1. - pwtlo);
            qcc412 = pwthi*qsh412 + qcc412*(1. - pwthi);
            qdif412 = qsl412 - qgc412;
// 488 --
            qgc488 = pwtlo*qsl488 + qgc488*(1. - pwtlo);
            qcc488 = pwthi*qsh488 + qcc488*(1. - pwthi);
            qdif488 = qsl488 - qgc488;
        }
    }
//  Calculate reflectivity
    float rh670 = DFILL_FLOAT;
    float rh412 = DFILL_FLOAT;
    float rh488 = DFILL_FLOAT;
    float rl670 = DFILL_FLOAT;
    float rl412 = DFILL_FLOAT;
    float rl488 = DFILL_FLOAT;
    float cl670 = DFILL_FLOAT;
    float cl412 = DFILL_FLOAT;
    float cl488 = DFILL_FLOAT;
    float r670_10 = 0.0;
    float r412_10 = 0.0;
    float r488_10 = 0.0;
    float r670_04 = 0.0;
    float r412_04 = 0.0;
    float r488_04 = 0.0;
    float d670 = 0.0;
    float d412 = 0.0;
    float d488 = 0.0;
    float grr = 0.02;

    for (int ip = 0; ip < 2; ip++) {
// compute table indices for 670 channel
        int i670_1 = ipt_1[ip] + lamtb1[3]  + iofset;
        int i670_2 = ipt_2[ip] + lamtb2[3] ;
// compute table indices for 412 channel
        int i412_1 = ipt_1[ip] + lamtb1[0] + iofset;
        int i412_2 = ipt_2[ip] + lamtb2[0];
// compute table indices for 488 channel
        int i488_1 = ipt_1[ip] + lamtb1[2] + iofset;
        int i488_2 = ipt_2[ip] + lamtb2[2];
// perform interpolations
        interx(i670_1, i670_2, rhot_band::W670, ez670[ip], t670[ip], sb670[ip]);
        interx(i412_1, i412_2, rhot_band::W410, ez412[ip], t412[ip], sb412[ip]);
        interx(i488_1, i488_2, rhot_band::W490, ez488[ip], t488[ip], sb488[ip]);
// determine calculated q values
        qgc670x = ez670[ip] + grr*t670[ip] / (1 - grr*sb670[ip]);
        qcc670x = ez670[ip] + clref*t670[ip] / (1 - clref*sb670[ip]);
        qgc412x = ez412[ip] + grr*t412[ip] / (1 - grr*sb412[ip]);
        qcc412x = ez412[ip] + clref*t412[ip] / (1 - clref*sb412[ip]);
        qgc488x = ez488[ip] + grr*t488[ip] / (1 - grr*sb488[ip]);
        qcc488x = ez488[ip] + clref*t488[ip] / (1 - clref*sb488[ip]);

        if (ip == 0) {
//  for ip = 1 (p = 1 atm) save q values
//  and denominator used in ref calculation
            qsl670 = qgc670x;
            qsh670 = qcc670x;
            qsl412 = qgc412x;
            qsh412 = qcc412x;
            qsl488 = qgc488x;
            qsh488 = qcc488x;
            d670 = alb670 - ez670[ip];
            d412 = alb412 - ez412[ip];
            d488 = alb488 - ez488[ip];
// calculate reflectivity using version 6 method
            if (d670 != 0.) {
                r670_10 = 1. / (t670[ip] / d670 + sb670[ip]);
                r412_10 = 1. / (t412[ip] / d412 + sb412[ip]);
                r488_10 = 1. / (t488[ip] / d488 + sb488[ip]);
            } else {
                r670_10 = 0.0;
                r412_10 = 0.0;
                r488_10 = 0.0;
            }
        } else {
//   for ip = 2 (p = 0.4 atm) interpolate
//   0.4 and 1 atm q values to obtain
//   final q values
            float pwtlo_base = 1.0;
            qgc670x = pwtlo_base*qsl670
                    + qgc670x*(1. - pwtlo_base);
            qcc670x = pwthi*qsh670 + qcc670x*(1. - pwthi);
//   412 --
            qgc412x = pwtlo_base*qsl412
                    + qgc412x*(1. - pwtlo_base);
            qcc412x = pwthi*qsh412 + qcc412x*(1. - pwthi);
//   488 --
            qgc488x = pwtlo_base*qsl488
                    + qgc488x*(1. - pwtlo_base);
            qcc488x = pwthi*qsh488 + qcc488x*(1. - pwthi);
//  determine denominator used in Version 6
//  ref calculation and calculate Version 6
//  reflectivity
            d670 = alb670 - ez670[ip];
            d412 = alb412 - ez412[ip];
            d488 = alb488 - ez488[ip];
            if (d670 != 0.) {
                r670_04 = 1. / (t670[ip] / d670 + sb670[ip]);
                rl670 = pwtlo_base*r670_10 +
                        (1. - pwtlo_base)*r670_04;
                rh670 = pwthi*r670_10 + (1. - pwthi)*r670_04;
//  412 --
                r412_04 = 1. / (t412[ip] / d412 + sb412[ip]);
                rl412 = pwtlo_base*r412_10 +
                        (1. - pwtlo_base)*r412_04;
                rh412 = pwthi*r412_10 + (1. - pwthi)*r412_04;
//  488 --
                r488_04 = 1. / (t488[ip] / d488 + sb488[ip]);
                rl488 = pwtlo_base*r488_10 +
                        (1. - pwtlo_base)*r488_04;
                rh488 = pwthi*r488_10 + (1. - pwthi)*r488_04;
            } else {
                rl670 = r670_10;
                rh670 = r670_10;
//  412 --
                rl412 = r412_10;
                rh412 = r412_10;
//  488 --
                rl488 = r488_10;
                rh488 = r488_10;
            }
        }
    }
//  Calculate cloud fraction
    isnow_ = 0;
    int partial = 1;
    if (qcc670x - qgc670x != 0.0) {
        cl670 = (alb670 - qgc670x) / (qcc670x - qgc670x);
        cl412 = (alb412 - qgc412x) / (qcc412x - qgc412x);
        cl488 = (alb488 - qgc488x) / (qcc488x - qgc488x);
    } else {
        cl670 = 0.0;
    }
    if (partial == 0) {
        cl670 = -1.0;
    } else {
//    if snow/ice probability > 50% recalculate cl670, grref
//    and, if necessary, clref.
        if (isnow_ >= 5) {
            algflg_ = 10;
            if (cl670 > 0.0 && cl670 < 1.0) {
                cl670 = cl670 / 2.0;
                alb670 = (alb670 - cl670*qcc670x) / (1. - cl670);
                for (int ip = 0; ip < 2; ip++) {
                    d670 = alb670 - ez670[ip];
                    if (ip == 0) {
                        if (d670 != 0.0) {
                            r670_10 = 1. / (t670[ip] / d670 + sb670[ip]);
                        } else {
                            r670_10 = 0.0;
                        }
                    } else {
                        if (d670 != 0.) {
                            r670_04 = 1.0 / (t670[ip] / d670 + sb670[ip]);
                            grref = pwtlo*r670_10 +
                                    (1. - pwtlo)*r670_04;
                        } else {
                            grref = r670_10;
                        }
                    }
                }
            } else if (cl670 >= 1.0) {
                cl670 = 0.5;
                grref = r670_10*pwtlo + r670_04*(1. - pwtlo);
                clref = r670_10*pwthi + r670_04*(1. - pwthi);
            } else {
                cl670 = 0.0;
                grref = r670_10*pwtlo + r670_04*(1. - pwtlo);
            }
        } else {
            algflg_ = 0;
        }
    }
//  calculate reflectivity
//  use partial cloud algorithm if 0 < cl670 < 1
//  otherwise use Version 6 calculation
    float ler412 = 0;
    float ler488 = 0;
    float ler670 = 0;
    if (cl670 < 0.0) {
        ler412 = rl412;
        ler488 = rl488;
        ler670 = rl670;
    } else if (cl670 > 1.0) {
        ler412 = rh412;
        ler488 = rh488;
        ler670 = rh670;
    } else {
        ler412 = grr + cl412*(clref - grr);
        ler488 = grr + cl488*(clref - grr);
        ler670 = grr + cl670*(clref - grr);
    }

    ler412_ = ler412*100;
    ler488_ = ler488*100;
    ler670_ = ler670*100;
    qdf412_ = qdif412;
    qdf488_ = qdif488;
    qdf670_ = qdif670;

    return status;
}

/**************************************************************************
 * NAME: compute_stdv()
 *
 * DESCRIPTION: Compute spatial standard deviation.
 *
 *************************************************************************/

int DbAlgLand::compute_stdv() {

    int status = DTDB_SUCCESS;
    return status;
}

/**************************************************************************
* NAME: intero()
*
* DESCRIPTION: Calculates intensity terms used in computation of
*   reflectivity for non 670 nm channels.
*   Performs table interpolation for theta and theta0.
*   Uses precomputed lagrangian coefs.
*   Calculates izero (intensity term arising from atmospheric
*   scattering; also calculates transmittance for reflected
*   radiation).
 *
 *************************************************************************/
int DbAlgLand::interx(int i1, int i2, rhot_band w, float& ozero, float& otr,
        float& osb) {

    int status = DTDB_SUCCESS;

    float inot = 0.0;
    float zone = 0.0;
    float ztwo = 0.0;
    float tran = 0.0;
    int m = 0;
    for (int i = 0; i < 4; i++) {
        int l = i1 + 8*i;
        for (int k = 0; k < 4; k++) {
                inot = inot + mt_lut_->LOGI0[l]*cofs_[m];
                zone = zone + mt_lut_->Z1I0[l]*cofs_[m];
                ztwo = ztwo + mt_lut_->Z2I0[l]*cofs_[m];
                tran = tran + mt_lut_->TI0[l]*cofs_[m];
            m++;
            l++;
        }
    }
//  convert table values into i1, i2, t, and i0
    inot = pow(10.0,inot);
    float ione = zone*phs_*inot;
    float itwo = ztwo*phsr_*phs_*inot;
    tran = tran*inot;
//  compute ezero, t, and sbar
    ozero = inot + ione*cphi_ + itwo*cphir_;
    otr = tran;
        osb = mt_lut_->SB[i2];

    return status;
}

/**************************************************************************
 * NAME: int DbAlgLand::set_fill_out()
 *
 * DESCRIPTION: Set output fill values.
 *
 *************************************************************************/

int DbAlgLand::set_fill_out()
{
    int status = DTDB_SUCCESS;

    for (int iw=0; iw<NLWL; iw++) {
        lOut_.aot[iw] = DFILL_FLOAT;
        lOut_.ssa[iw] = DFILL_FLOAT;
        lOut_.sr[iw] = DFILL_FLOAT;
    }
    lOut_.aot550    = DFILL_FLOAT;
    lOut_.ae        = DFILL_FLOAT;
    lOut_.sfc_type  = DFILL_SHORT;
    lOut_.alg_flag  = DFILL_SHORT;
    lOut_.ndvi  = DFILL_FLOAT;
    lOut_.aerosol_type = 8;
    mask_cm_ = DFILL_SHORT;

    return status;
}


