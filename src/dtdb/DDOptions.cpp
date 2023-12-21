/**************************************************************************
*
* NAME: DDOptions.cpp
*
* DESCRIPTION: Sets up the options defined in the program parameters file
*              and the command line.
*
*  Created on: Aug 25, 2020
*      Author: Sam Anderson, DT
*
**************************************************************************/

#include <cassert>
#include <iostream>
#include <ostream>
#include <cstring>
#include <vector>
#include <libgen.h>

#include <DDOptions.h>

using namespace std;

//---------------------------------------------------------------------------
// Static constants
//---------------------------------------------------------------------------

const string INPUT_PAR               	= "par";
const string DTDB_XML               	= "dtdb_xml";
const string PRODUCT_XML               	= "prod_xml";

const string INPUT_L1B					= "ifile";
const string INPUT_GEO                 	= "geofile";
const string INPUT_MET1               	= "met1";
const string INPUT_MET2               	= "met2";
const string INPUT_MET                	= "pixel_anc_file";
const string INPUT_CLDMASK              = "cloud_mask_file";
const string OUTPUT_NC4					= "ofile";

const string INPUT_SPIX                = "spixl";
const string INPUT_EPIX                = "epixl";
const string INPUT_DPIX                = "dpixl";
const string INPUT_SLINE               = "sline";
const string INPUT_ELINE               = "eline";
const string INPUT_DLINE               = "dline";

const string BOOL_MASK_GLINT	  	   = "maskglint";
const string BOOL_MASK_CLOUD		   = "maskcloud";
const string BOOL_MASK_SOLZ		       = "masksunzen";
const string BOOL_MASK_SENZ			   = "masksatzen";
const string BOOL_GAS_CORRECTION       = "gas_correction";

const string INPUT_THRESH_SOLZ		   = "sunzen";
const string INPUT_THRESH_SENZ         = "satzen";
const string INPUT_THRESH_GLINT		   = "glint_thresh";

const string INPUT_LT_NOISE_SCALE      = "lt_noise_scale";
const string INPUT_LINES_PER_RW        = "lines_per_rw";

const string INPUT_IFFSVM				= "ifile_svm";
const string INPUT_LANDMASK            	= "landmask";
const string INPUT_GAS_CORRECTION		= "lut_gascorr";
const string INPUT_SAT_ID				= "lut_satid";
const string INPUT_GRIB                	= "met";
const string INPUT_OZONE               	= "ozone";
//Dark Target
const string INPUT_OCEAN_BIG1 			= "lut_ocean_big1";
const string INPUT_OCEAN_BIG2 			= "lut_ocean_big2";
const string INPUT_OCEAN_BIG3 			= "lut_ocean_big3";
const string INPUT_OCEAN_BIG4 			= "lut_ocean_big4";
const string INPUT_OCEAN_SMALL1 		= "lut_ocean_small1";
const string INPUT_OCEAN_SMALL2 		= "lut_ocean_small2";
const string INPUT_OCEAN_SMALL3 		= "lut_ocean_small3";
const string INPUT_OCEAN_SMALL4 		= "lut_ocean_small4";
const string INPUT_LAND_W0466 			= "lut_land_w0466";
const string INPUT_LAND_W0554 			= "lut_land_w0554";
const string INPUT_LAND_W0645 			= "lut_land_w0645";
const string INPUT_LAND_W2113 			= "lut_land_w2113";
const string INPUT_LAND_MAP				= "lut_land_aerosol_map";
// Water vapor
const string INPUT_TRANSM_H2O_1        = "lut_transm_h2o_1";
const string INPUT_TRANSM_H2O_2        = "lut_transm_h2o_2";
const string INPUT_TRANSM_H2O_3        = "lut_transm_h2o_3";
const string INPUT_TRANSM_H2O_4        = "lut_transm_h2o_4";
const string INPUT_TRANSM_H2O_5        = "lut_transm_h2o_5";
const string INPUT_TRANSM_H2O_6        = "lut_transm_h2o_6";
const string INPUT_WEIGHT_TABLE        = "lut_weight_table";
const string INPUT_REFL_CH2            = "lut_refl_ch2";
const string INPUT_RATIO_CH19_TO_CH2   = "lut_ratio_ch19_to_ch2";
// Deep Blue
//const string INPUT_MODIS_TABLE         = "lut_modis_table";
const string INPUT_GLOBAL_IGBP         = "lut_global_igbp";
const string INPUT_NVALX_412           = "lut_nvalx_412";
const string INPUT_NVALX_470           = "lut_nvalx_470";
const string INPUT_NVALX_650           = "lut_nvalx_650";
const string INPUT_NVALX21             = "lut_nvalx21";
const string INPUT_DBDT_REGIONS        = "lut_dbdt_regions";
const string INPUT_MODIS_SWIR          = "lut_modis_swir";
const string INPUT_MODIS_SURF_REFL     = "lut_modis_surfdb";
const string INPUT_MODIS_XCAL_412      = "lut_modis_xcal_412";
const string INPUT_MODIS_XCAL_470      = "lut_modis_xcal_470";
const string INPUT_GAIN_412            = "lut_gain_412";
const string INPUT_GAIN_470            = "lut_gain_470";
const string INPUT_ATMOSPHERE          = "lut_atmosphere";

const string INPUT_AERO_LAND_FINE      = "lut_aero_land_fine";
const string INPUT_AERO_LAND_DUST      = "lut_aero_land_dust";
const string INPUT_AERO_OCEAN_DUST     = "lut_aero_ocean_dust";
const string INPUT_AERO_OCEAN_FINE     = "lut_aero_ocean_fine";
const string INPUT_AERO_OCEAN_MARI     = "lut_aero_ocean_mari";
const string INPUT_AERO_OCEAN_MIX      = "lut_aero_ocean_mix";
const string INPUT_BATHYMETRY          = "lut_bathymetry";
const string INPUT_CHL                 = "lut_chl";
const string INPUT_LER_TABLE           = "lut_ler_table";
const string INPUT_LANDCOVER           = "lut_landcover";
const string INPUT_SURFACE_PRESSURE    = "lut_surface_pressure";
const string INPUT_GEOZONE             = "lut_geozone";
const string INPUT_SEASONAL_DESERTS    = "lut_seasonal_deserts";
const string INPUT_BRDF                = "lut_brdf";
const string INPUT_VIIRS_SURF_REFL     = "lut_viirs_surfdb";
const string INPUT_SURF_COEFF          = "lut_surfcoeff";
const string INPUT_SWIR                = "lut_swir";
const string INPUT_VEG_LANDCOVER       = "lut_veg_landcover";
const string INPUT_VEG_21SFC           = "lut_veg_21sfc";
const string INPUT_RAYL_412            = "lut_rayl_412";
const string INPUT_RAYL_470            = "lut_rayl_488";
const string INPUT_RAYL_650            = "lut_rayl_670";
const string INPUT_VIIRS_XCAL_412      = "lut_viirs_xcal_412";
const string INPUT_VIIRS_XCAL_488      = "lut_viirs_xcal_488";
const string INPUT_VIIRS_XCAL_670      = "lut_viirs_xcal_670";

const string INPUT_NC4_LUT	           = "lut";
const string INPUT_DT_NC4_LUT          = "lut_nc4_dt";
const string INPUT_DB_NC4_LUT          = "lut_nc4_db";
const string INPUT_SENSOR_INFO     	   = "sensor_info";
const string LEAP_SEC_PATH 	       	   = "leapsec_file";

const string LUT_GRIB 		           = "GRIB";
const string LUT_GAS_CORRECTION 	   = "GAS_CORRECTION";
const string LUT_LAND_AEROSOL 	       = "LAND_AEROSOL";
const string LUT_OCEAN_AEROSOL         = "OCEAN_AEROSOL";
const string LUT_WATER_VAPOR           = "WATER_VAPOR";
const string LUT_SURFACE_PRESSURE      = "SURFACE_PRESSURE";
const string LUT_MODIS_SWIR            = "MODIS_SWIR_VIS";
const string LUT_MODIS_SURFACE_REFL    = "MODIS_SURFACE_REFLECTANCE";
const string LUT_MODIS_CORRECTIONS     = "MODIS_CORRECTIONS";
const string LUT_NVALX                 = "NVALX";
const string LUT_RAYLEIGH              = "RAYLEGH";

const string LUT_NVALX21               = "VEG_21SFC";
const string LUT_LER_TABLES            = "LER_TABLES";
const string LUT_SWIR                  = "SWIR_VS_VISIBLE";
const string LUT_VIIRS_SURFACE_REFL    = "VIIRS_SURFACE_REFLECTANCE";
const string LUT_SURFACE_COEFF         = "SURFACE_COEFFICIENTS";
const string LUT_LANDCOVER             = "LANDCOVER";
const string LUT_GEOZONE               = "GEOZONE";
const string LUT_OCEAN_AEROSOL_FINE    = "OCEAN_AEROSOL_FINE";
const string LUT_OCEAN_AEROSOL_DUST    = "OCEAN_AEROSOL_DUST";
const string LUT_OCEAN_AEROSOL_MARI    = "OCEAN_AEROSOL_MARITIME";
const string LUT_OCEAN_AEROSOL_MIX     = "OCEAN_AEROSOL_MIXED";
const string LUT_LAND_AEROSOL_FINE     = "LAND_AEROSOL_FINE";
const string LUT_LAND_AEROSOL_DUST     = "LAND_AEROSOL_DUST";
const string LUT_BATHYMETRY            = "BATHYMETRY";
const string LUT_CHL                   = "LOG_CHL";

const string BOOL_SENSOR               = "sensor_band_parameters";
const string BOOL_SCANS                = "scan_line_attributes";
const string BOOL_NAVIGATION           = "navigation_data";
const string BOOL_PROCESS              = "processing_control";
const string BOOL_WAVELENGTHS          = "wavelengths";
const string BOOL_DATA                 = "geophysical_data";
const string BOOL_FLAGS                = "flag_percentages";
const string BOOL_ANCILLARY            = "ancillary";
const string BOOL_GEOLOCATION          = "geolocation";
const string BOOL_OBSERVATIONS         = "observations";
const string BOOL_STATISTICS		   = "statistics";
const string BOOL_ADD_LT_NOISE         = "add_lt_noise";
const string BOOL_SHORT_FORMAT		   = "short_format";

const string BOOL_QUALITY      			= "quality";
const string BOOL_AEROSOL_TYPE 			= "aerosol_type";
const string BOOL_SCATTANG 				= "scattang";
const string BOOL_FMF_550 				= "fmf_550";
const string BOOL_ANGSTROM 				= "angstrom";
const string BOOL_AOT_380 				= "aot_380";
const string BOOL_AOT_410 				= "aot_410";
const string BOOL_AOT_490 				= "aot_490";
const string BOOL_AOT_550 				= "aot_550";
const string BOOL_AOT_670 				= "aot_670";
const string BOOL_AOT_865 				= "aot_865";
const string BOOL_AOT_1240 				= "aot_1240";
const string BOOL_AOT_1610 				= "aot_1610";
const string BOOL_AOT_2250 				= "aot_2250";
const string BOOL_SURF_410 				= "surface_410";
const string BOOL_SURF_490 				= "surface_490";
const string BOOL_SURF_670 				= "surface_670";
const string BOOL_SURF_2250 			= "surface_2250";

const string ALGORITHM                  = "alg";
const string NETCDF_LUT_PATH 			= "lut_nc4_path";

// global variable to hold the command line parameters
static clo_optionList_t* global_optionList = NULL;

void set_optionList(clo_optionList_t* list) {
    global_optionList = list;
}

clo_optionList_t* get_optionList() {
    if(global_optionList == NULL) {
        cerr << "get_optionList: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    return global_optionList;
}

string  get_option(const string& name) {
    if(global_optionList == NULL) {
        cerr << "get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    clo_option_t* option = clo_findOption(global_optionList, name.c_str());
    if(option == NULL || !clo_isOptionSet(option))
        return "";
    string result(clo_getOptionString(option));
    return result;
}

int get_option_int(const string& name) {
    if(global_optionList == NULL) {
        cerr << "get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    return clo_getInt(global_optionList, name.c_str());
}

float get_option_float(const string& name) {
    if(global_optionList == NULL) {
        cerr << "get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    return clo_getFloat(global_optionList, name.c_str());
}

int get_bool(const string& name) {
    if(global_optionList == NULL) {
        cerr << "get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    return clo_getBool(global_optionList, name.c_str());
}

float*  get_option_floats(const string& name, int *count) {
    if(global_optionList == NULL) {
        cerr << "get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    clo_option_t* option = clo_findOption(global_optionList, name.c_str());
    if(option == NULL || !clo_isOptionSet(option))
        return nullptr;
    float* result(clo_getOptionFloats(option, count));
    return result;
}

void add_options(clo_optionList_t* list) {
    if(global_optionList == NULL)
        global_optionList = list;

    const string INPUT_SPIX                = "spixl";
    const string INPUT_EPIX                = "epixl";
    const string INPUT_DPIX                = "dpixl";
    const string INPUT_SLINE               = "sline";
    const string INPUT_ELINE               = "eline";
    const string INPUT_DLINE               = "dline";
    const string INPUT_LT_NOISE_SCALE      = "lt_noise_scale";

    clo_addOption(list, INPUT_SPIX.c_str(), CLO_TYPE_INT, "1", "start pixel number");
    clo_addOption(list, INPUT_EPIX.c_str(), CLO_TYPE_INT, "-1", "end pixel number (-1=the last pixel)");
    clo_addOption(list, INPUT_DPIX.c_str(), CLO_TYPE_INT, "1", "pixel sub-sampling interval");
    clo_addOption(list, INPUT_SLINE.c_str(), CLO_TYPE_INT, "1", "start line number");
    clo_addOption(list, INPUT_ELINE.c_str(), CLO_TYPE_INT, "-1", "end line number (-1=the last line)");
    clo_addOption(list, INPUT_DLINE.c_str(), CLO_TYPE_INT, "1", "line sub-sampling interval");
    clo_addOption(list, INPUT_LINES_PER_RW.c_str(), CLO_TYPE_INT, "1", "lines per read/process/write");
    clo_addOption(list, INPUT_LT_NOISE_SCALE.c_str(), CLO_TYPE_FLOAT, NULL, "noise scaling factors");
    clo_addOption(list, INPUT_THRESH_SOLZ.c_str(), CLO_TYPE_FLOAT, "80", "maximum solar zenith angle");
    clo_addOption(list, INPUT_THRESH_SENZ.c_str(), CLO_TYPE_FLOAT, "90", "maximum sensor zenith angle");
    clo_addOption(list, INPUT_THRESH_GLINT.c_str(), CLO_TYPE_FLOAT, "0.005", "glint threshold");

    clo_addOption(list, DTDB_XML.c_str(), CLO_TYPE_IFILE, NULL, "dtdb xml file name");
    clo_addOption(list, PRODUCT_XML.c_str(), CLO_TYPE_IFILE, NULL, "product xml file name");
    clo_addOption(list, INPUT_L1B.c_str(), CLO_TYPE_IFILE, NULL, "input L1B file name");
    clo_addOption(list, INPUT_GEO.c_str(), CLO_TYPE_IFILE, NULL, "input geolocation file name");
    clo_addOption(list, INPUT_MET1.c_str(), CLO_TYPE_IFILE, NULL, "input ancillary LUT file before");
    clo_addOption(list, INPUT_MET2.c_str(), CLO_TYPE_IFILE, NULL, "input ancillary LUT file after");
    clo_addOption(list, INPUT_MET.c_str(), CLO_TYPE_IFILE, NULL, "input ancillary LUT file");
    clo_addOption(list, INPUT_GRIB.c_str(), CLO_TYPE_IFILE, NULL, "input meteorological LUT file");
    clo_addOption(list, INPUT_OZONE.c_str(), CLO_TYPE_IFILE, NULL, "input ozone LUT file");
    clo_addOption(list, INPUT_ATMOSPHERE.c_str(), CLO_TYPE_IFILE, NULL, "input atmosphere file name");
    clo_addOption(list, INPUT_CLDMASK.c_str(), CLO_TYPE_IFILE, NULL, "input cloud mask file name");
    clo_addOption(list, INPUT_LANDMASK.c_str(), CLO_TYPE_IFILE, NULL, "input land mask LUT file");

    clo_addOption(list, INPUT_GAS_CORRECTION.c_str(), CLO_TYPE_IFILE, NULL, "input gas correction LUT file");
    clo_addOption(list, INPUT_OCEAN_BIG1.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_OCEAN_BIG2.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_OCEAN_BIG3.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_OCEAN_BIG4.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_OCEAN_SMALL1.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_OCEAN_SMALL2.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_OCEAN_SMALL3.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_OCEAN_SMALL4.c_str(), CLO_TYPE_IFILE, NULL, "input ocean aerosol LUT file");
    clo_addOption(list, INPUT_LAND_W0466.c_str(), CLO_TYPE_IFILE, NULL, "input land aerosol LUT file");
    clo_addOption(list, INPUT_LAND_W0554.c_str(), CLO_TYPE_IFILE, NULL, "input land aerosol LUT file");
    clo_addOption(list, INPUT_LAND_W0645.c_str(), CLO_TYPE_IFILE, NULL, "input land aerosol LUT file");
    clo_addOption(list, INPUT_LAND_W2113.c_str(), CLO_TYPE_IFILE, NULL, "input land aerosol LUT file");
    clo_addOption(list, INPUT_LAND_MAP.c_str(), CLO_TYPE_IFILE, NULL, "input land aerosol LUT file");

    clo_addOption(list, INPUT_TRANSM_H2O_1.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_TRANSM_H2O_2.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_TRANSM_H2O_3.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_TRANSM_H2O_4.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_TRANSM_H2O_5.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_TRANSM_H2O_6.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_WEIGHT_TABLE.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_REFL_CH2.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_RATIO_CH19_TO_CH2.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");

    clo_addOption(list, INPUT_SURFACE_PRESSURE.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_NVALX_412.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_NVALX_470.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_NVALX_650.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_SEASONAL_DESERTS.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_DBDT_REGIONS.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_MODIS_SURF_REFL.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_GLOBAL_IGBP.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_NVALX21.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_MODIS_SWIR.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_MODIS_XCAL_412.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_MODIS_XCAL_470.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_GAIN_412.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_GAIN_470.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_RAYL_412.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_RAYL_470.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_RAYL_650.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");

    clo_addOption(list, INPUT_AERO_LAND_FINE.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_AERO_LAND_DUST.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_AERO_OCEAN_DUST.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_AERO_OCEAN_FINE.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_AERO_OCEAN_MARI.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_AERO_OCEAN_MIX.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_BATHYMETRY.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_CHL.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_LER_TABLE.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_SWIR.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_LANDCOVER.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_GEOZONE.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_VIIRS_SURF_REFL.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_SURF_COEFF.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_BRDF.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_VEG_LANDCOVER.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_VEG_21SFC.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_VIIRS_XCAL_412.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_VIIRS_XCAL_488.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");
    clo_addOption(list, INPUT_VIIRS_XCAL_670.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue LUT file");

    clo_addOption(list, INPUT_NC4_LUT.c_str(), CLO_TYPE_IFILE, NULL, "input DTDB NetCDF4 LUT file");
    clo_addOption(list, INPUT_DT_NC4_LUT.c_str(), CLO_TYPE_IFILE, NULL, "input Darktarget NetCDF4 LUT file");
    clo_addOption(list, INPUT_DB_NC4_LUT.c_str(), CLO_TYPE_IFILE, NULL, "input DeepBlue NetCDF4 LUT file");
    clo_addOption(list, INPUT_SENSOR_INFO.c_str(), CLO_TYPE_IFILE, NULL, "OCI sensor info, F0");

    clo_addOption(list, OUTPUT_NC4.c_str(), CLO_TYPE_OFILE, NULL, "output dark target filename");
    clo_addOption(list, NETCDF_LUT_PATH.c_str(), CLO_TYPE_OFILE, NULL, "path to NetCDF4 LUT output directory");
    clo_addOption(list, ALGORITHM.c_str(), CLO_TYPE_STRING, NULL, "Algorithm identifier");

    clo_addOption(list, BOOL_SENSOR.c_str(), CLO_TYPE_BOOL, "True",
            "Report input sensor data");
    clo_addOption(list, BOOL_SCANS.c_str(), CLO_TYPE_BOOL, "True",
            "Report input scan line data");
    clo_addOption(list, BOOL_NAVIGATION.c_str(), CLO_TYPE_BOOL, "True",
            "Report input navigation data");
    clo_addOption(list, BOOL_PROCESS.c_str(), CLO_TYPE_BOOL, "True",
            "Report process data");
    clo_addOption(list, BOOL_WAVELENGTHS.c_str(), CLO_TYPE_BOOL, "True",
            "Report all wavelength sets");
    clo_addOption(list, BOOL_DATA.c_str(), CLO_TYPE_BOOL, "True",
            "Report all geophysical product data");
    clo_addOption(list, BOOL_FLAGS.c_str(), CLO_TYPE_BOOL, "True",
            "Report flag percentages");
    clo_addOption(list, BOOL_ANCILLARY.c_str(), CLO_TYPE_BOOL, "False",
            "Report input ancillary data");
    clo_addOption(list, BOOL_GEOLOCATION.c_str(), CLO_TYPE_BOOL, "False",
            "Report input geolocation data");
    clo_addOption(list, BOOL_OBSERVATIONS.c_str(), CLO_TYPE_BOOL, "False",
            "Report input TOA and surface reflectance data");
    clo_addOption(list, BOOL_STATISTICS.c_str(), CLO_TYPE_BOOL, "False",
            "Report statistical representations");
    clo_addOption(list, BOOL_ADD_LT_NOISE.c_str(), CLO_TYPE_BOOL, "False",
            "Add noise to observed inputs");
    clo_addOption(list, BOOL_GAS_CORRECTION.c_str(), CLO_TYPE_BOOL, "True",
            "Apply gas corrections");
    clo_addOption(list, BOOL_MASK_GLINT.c_str(), CLO_TYPE_BOOL, "True",
            "Apply glint mask");
    clo_addOption(list, BOOL_MASK_CLOUD.c_str(), CLO_TYPE_BOOL, "True",
            "Apply cloud mask");
    clo_addOption(list, BOOL_MASK_SOLZ.c_str(), CLO_TYPE_BOOL, "True",
            "Apply solar zenith angle mask");
    clo_addOption(list, BOOL_MASK_SENZ.c_str(), CLO_TYPE_BOOL, "True",
            "Apply sensor zenith angle mask");
    clo_addOption(list, BOOL_SHORT_FORMAT.c_str(), CLO_TYPE_BOOL, "True",
            "Convert floats to 16 bit short integer");

//output product options
    clo_addOption(list, BOOL_QUALITY.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AEROSOL_TYPE.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_SCATTANG.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_FMF_550.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_ANGSTROM.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_380.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_410.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_490.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_550.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_670.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_865.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_1240.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_1610.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_AOT_2250.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_SURF_410.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_SURF_490.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_SURF_670.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");
    clo_addOption(list, BOOL_SURF_2250.c_str(), CLO_TYPE_BOOL, "False",
            "Product selection option");

}

void copy_options() {
//    clo_option_t* option;
//    clo_optionList_t* list = get_optionList();
}

string get_source() {
    vector<string> sourcesList;

    sourcesList.push_back("met1");
    sourcesList.push_back("met2");
    sourcesList.push_back("lut_gascorr");
    sourcesList.push_back("lut_satid");
    sourcesList.push_back("lut_ocean_big1");
    sourcesList.push_back("lut_ocean_big2");
    sourcesList.push_back("lut_ocean_big3");
    sourcesList.push_back("lut_ocean_big4");
    sourcesList.push_back("lut_ocean_small1");
    sourcesList.push_back("lut_ocean_small2");
    sourcesList.push_back("lut_ocean_small3");
    sourcesList.push_back("lut_ocean_small4");
    sourcesList.push_back("lut_land_w0466");
    sourcesList.push_back("lut_land_w0554");
    sourcesList.push_back("lut_land_w0645");
    sourcesList.push_back("lut_land_w2113");
    sourcesList.push_back("lut_land_aerosol_map");
    sourcesList.push_back("lut_nc4");
    
    string source = basename((char*)get_option(INPUT_IFFSVM).c_str());
    
    for(vector<string>::iterator it = sourcesList.begin(); it < sourcesList.end(); it++) {
        string str = get_option(*it);
        if(!str.empty()) {
            source.append(",");
            source.append(basename((char*)str.c_str()));
        }
    }
    return source;
}

string get_history(int argc, char* argv[]) {
    string history = basename(argv[0]);
    for (int i=1; i<argc; i++) {
        history.append(" ");
        history.append(basename(argv[i]));
    }
    return history;
}

// should the par file processing descend into other par files
static int enableFileDescending = 1;

/** CLO callback function for the "par" option.  Loads the parameter file
 into the list that is stored in option->cb_data */
void par_option_cb(struct clo_option_t *option) {
    if (enableFileDescending)
        clo_readFile((clo_optionList_t*) option->cb_data, option->valStr);
}

/** add all of the accepted command line options to list */
void init_options(clo_optionList_t* list, const char* softwareVersion) {
    char tmpStr[2048];
    clo_option_t* option;

    clo_setVersion2("dtdb", softwareVersion);

    sprintf(tmpStr, "Usage: dtdb argument-list\n\n");

    strcat(tmpStr, "  This program generates a Deep Blue or Dark Target product.\n\n");

    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the command line, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with command line over-riding.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    // add the parameters common to all VIIRS programs
    add_options(list);

    option = clo_addOption(list, "verbose", CLO_TYPE_BOOL, "false", "turn on verbose output");
    clo_addOptionAlias(option, "v");

}


/*
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - load the main program defaults file
    - load the command line (including specified par files)
    - re-load the command line disabling file descending so command
       line arguments will over ride

 */
void read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];
    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        exit(EXIT_FAILURE);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);

    enableFileDescending = 1;

    // load program defaults
    sprintf(tmpStr, "%s/common/dtdb_defaults.par", dataRoot);
    clo_readFile(list, tmpStr);
    clo_setString(list, "par",  tmpStr, "");

    // read all arguments including descending par files
    clo_readArgs(list, argc, argv);

    // enable the dump option
    clo_setEnableDumpOptions(1);

    // if a file was descended, make sure the command args over-rides
    enableFileDescending = 0;
    clo_readArgs(list, argc, argv);
    enableFileDescending = 1;

    copy_options();
}
