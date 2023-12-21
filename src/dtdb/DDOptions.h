
/**************************************************************************
 *
 * NAME: DDOptions.h
 *
 * DESCRIPTION: Header file for declaration of string names and
 * command line options
 *
 *  Created on: Aug 25, 2020
 *      Author: Sam Anderson
 *
 **************************************************************************/

#ifndef _Options_h_
#define _Options_h_

#include <string>
#include <clo.h>

#define VIIRSSCANS  203
#define VDETECTORS  16
#define VLINES      VIIRSSCANS*VDETECTORS
#define VELEMS      3200
#define VNUMSCAN_B  VLINES/VGRIDY
#define VSWATH_B    VELEMS
#define VNUMCELLS_B VELEMS/VGRIDX

#define MODISSCANS  203
#define MDETECTORS  10
#define MLINES      MODISSCANS*MDETECTORS
#define MELEMS      1354
#define MNUMSCAN_B  MLINES/MGRIDY
#define MSWATH_B    MELEMS
#define MNUMCELLS_B MELEMS/MGRIDX

#define PACESCANS   1721
#define PDETECTORS  1
#define PLINES      PACESCANS*PDETECTORS
#define PELEMS      1271
#define PNUMSCAN_B  PLINES/PGRIDY
#define PSWATH_B    PELEMS
#define PNUMCELLS_B PELEMS/PGRIDX
#define NUM_PACE_RSB   60
#define NUM_PACE_SWIR  9

#define NUM_BANDS    	13
#define NUM_DB_BANDS    12
#define DB_RFL_BANDS    9
#define DB_BT_BANDS     3
#define DB_NC_BANDS     3

#define NUM_BT_BANDS    5
#define NUM_EM_BANDS    5
#define NUM_RF_BANDS    11
#define NUM_SW_BANDS    3
#define NUM_DT_QUALS    2
#define NUM_LAND_SOL3   3
#define NUM_LAND_SOL4   4

static constexpr int DD_SUCCESS = 0;
static constexpr int DD_FAIL = 1;

//
// Short name constants
//
extern const std::string INPUT_PAR;
extern const std::string DTDB_XML;
extern const std::string PRODUCT_XML;

extern const std::string INPUT_L1B;
extern const std::string INPUT_GEO;
extern const std::string INPUT_MET1;
extern const std::string INPUT_MET2;
extern const std::string INPUT_MET;
extern const std::string OUTPUT_NC4;

extern const std::string INPUT_SPIX;
extern const std::string INPUT_EPIX;
extern const std::string INPUT_DPIX;
extern const std::string INPUT_SLINE;
extern const std::string INPUT_ELINE;
extern const std::string INPUT_DLINE;

extern const std::string BOOL_MASK_GLINT;
extern const std::string BOOL_MASK_CLOUD;
extern const std::string BOOL_MASK_SOLZ;
extern const std::string BOOL_MASK_SENZ;
extern const std::string BOOL_GAS_CORRECTION;

extern const std::string INPUT_THRESH_SOLZ;
extern const std::string INPUT_THRESH_SENZ;
extern const std::string INPUT_THRESH_GLINT;

extern const std::string INPUT_LT_NOISE_SCALE;
extern const std::string INPUT_LINES_PER_RW;

extern const std::string INPUT_IFFSVM;
extern const std::string INPUT_GRIB;
extern const std::string INPUT_OZONE;
extern const std::string INPUT_LANDMASK;
extern const std::string INPUT_GAS_CORRECTION;
extern const std::string INPUT_SAT_ID;
extern const std::string INPUT_MISSION;

extern const std::string INPUT_OCEAN_BIG1;
extern const std::string INPUT_OCEAN_BIG2;
extern const std::string INPUT_OCEAN_BIG3;
extern const std::string INPUT_OCEAN_BIG4;
extern const std::string INPUT_OCEAN_SMALL1;
extern const std::string INPUT_OCEAN_SMALL2;
extern const std::string INPUT_OCEAN_SMALL3;
extern const std::string INPUT_OCEAN_SMALL4;
extern const std::string INPUT_LAND_W0466;
extern const std::string INPUT_LAND_W0554;
extern const std::string INPUT_LAND_W0645;
extern const std::string INPUT_LAND_W2113;
extern const std::string INPUT_LAND_MAP;

extern const std::string INPUT_TRANSM_H2O_1;
extern const std::string INPUT_TRANSM_H2O_2;
extern const std::string INPUT_TRANSM_H2O_3;
extern const std::string INPUT_TRANSM_H2O_4;
extern const std::string INPUT_TRANSM_H2O_5;
extern const std::string INPUT_TRANSM_H2O_6;
extern const std::string INPUT_WEIGHT_TABLE;
extern const std::string INPUT_SURFACE_PRESSURE;

extern const std::string INPUT_REFL_CH2;
extern const std::string INPUT_RATIO_CH19_TO_CH2;
extern const std::string INPUT_DBDT_REGIONS;
extern const std::string INPUT_GLOBAL_IGBP;
extern const std::string INPUT_NVALX_412;
extern const std::string INPUT_NVALX_470;
extern const std::string INPUT_NVALX_650;
extern const std::string INPUT_NVALX21;
extern const std::string INPUT_RAYL_412;
extern const std::string INPUT_RAYL_470;
extern const std::string INPUT_RAYL_650;
extern const std::string INPUT_MODIS_SWIR;
extern const std::string INPUT_MODIS_SURF_REFL;
extern const std::string INPUT_MODIS_GEOZONE;
extern const std::string INPUT_MODIS_XCAL_412;
extern const std::string INPUT_MODIS_XCAL_470;
extern const std::string INPUT_GAIN_412;
extern const std::string INPUT_GAIN_470;
extern const std::string INPUT_CLDMASK;

extern const std::string INPUT_AERO_LAND_FINE;
extern const std::string INPUT_AERO_LAND_DUST;
extern const std::string INPUT_AERO_OCEAN_DUST;
extern const std::string INPUT_AERO_OCEAN_FINE;
extern const std::string INPUT_AERO_OCEAN_MARI;
extern const std::string INPUT_AERO_OCEAN_MIX;
extern const std::string INPUT_BATHYMETRY;
extern const std::string INPUT_CHL;
extern const std::string INPUT_LER_TABLE;
extern const std::string INPUT_LANDCOVER;
extern const std::string INPUT_GEOZONE;
extern const std::string INPUT_SEASONAL_DESERTS;
extern const std::string INPUT_BRDF;
extern const std::string INPUT_VIIRS_SURF_REFL;
extern const std::string INPUT_SURF_COEFF;
extern const std::string INPUT_SWIR;
extern const std::string INPUT_VEG_LANDCOVER;
extern const std::string INPUT_VEG_21SFC;
extern const std::string INPUT_VIIRS_XCAL_412;
extern const std::string INPUT_VIIRS_XCAL_488;
extern const std::string INPUT_VIIRS_XCAL_670;

extern const std::string INPUT_NC4_LUT;
extern const std::string INPUT_DT_NC4_LUT;
extern const std::string INPUT_DB_NC4_LUT;
extern const std::string INPUT_SENSOR_INFO;

extern const std::string LUT_GRIB ;
extern const std::string LUT_GAS_CORRECTION;
extern const std::string LUT_LAND_AEROSOL;
extern const std::string LUT_WATER_VAPOR;
extern const std::string LUT_LER_TABLES;
extern const std::string LUT_SURFACE_PRESSURE;
extern const std::string LUT_NVALX;
extern const std::string LUT_DESERTS;
extern const std::string LUT_RAYLEIGH;
extern const std::string LUT_LANDCOVER;
extern const std::string LUT_NVALX21;
extern const std::string LUT_OCEAN_AEROSOL_FINE;
extern const std::string LUT_OCEAN_AEROSOL_DUST;
extern const std::string LUT_OCEAN_AEROSOL_MARI;
extern const std::string LUT_OCEAN_AEROSOL_MIX;
extern const std::string LUT_LAND_AEROSOL_FINE;
extern const std::string LUT_LAND_AEROSOL_DUST;
extern const std::string LUT_BATHYMETRY;
extern const std::string LUT_CHL;
extern const std::string LUT_SWIR;
extern const std::string LUT_VIIRS_SURFACE_REFL;
extern const std::string LUT_SURFACE_COEFF;
extern const std::string LUT_GEOZONE;

extern const std::string LUT_OCEAN_AEROSOL;
extern const std::string LUT_MODIS_CORRECTIONS;
extern const std::string LUT_MODIS_SWIR;
extern const std::string LUT_MODIS_SURFACE_REFL;

extern const std::string NETCDF_LUT_PATH;
extern const std::string LEAP_SEC_PATH;

extern const std::string BOOL_SENSOR;
extern const std::string BOOL_SCANS;
extern const std::string BOOL_NAVIGATION;
extern const std::string BOOL_PROCESS;
extern const std::string BOOL_ANCILLARY;
extern const std::string BOOL_GEOLOCATION;
extern const std::string BOOL_OBSERVATIONS;
extern const std::string BOOL_STATISTICS;
extern const std::string BOOL_WAVELENGTHS;
extern const std::string BOOL_DATA;
extern const std::string BOOL_FLAGS;
extern const std::string BOOL_ADD_LT_NOISE;
extern const std::string BOOL_SHORT_FORMAT;

extern const std::string BOOL_L2_FLAGS;
extern const std::string BOOL_QUALITY;
extern const std::string BOOL_AEROSOL_TYPE;
extern const std::string BOOL_SCATTANG;
extern const std::string BOOL_FMF_550;
extern const std::string BOOL_ANGSTROM;
extern const std::string BOOL_AOT_380;
extern const std::string BOOL_AOT_410;
extern const std::string BOOL_AOT_490;
extern const std::string BOOL_AOT_550;
extern const std::string BOOL_AOT_670;
extern const std::string BOOL_AOT_865;
extern const std::string BOOL_AOT_1240;
extern const std::string BOOL_AOT_1610;
extern const std::string BOOL_AOT_2250;
extern const std::string BOOL_SURF_410;
extern const std::string BOOL_SURF_490;
extern const std::string BOOL_SURF_670;
extern const std::string BOOL_SURF_2250;

extern const std::string ALGORITHM;

void set_optionList(clo_optionList_t* list);
clo_optionList_t* get_optionList();
std::string get_option(const std::string& name);
int get_option_int(const std::string& name);
float get_option_float(const std::string& name);
float* get_option_floats(const std::string& name, int *count);
int get_bool(const std::string& name);
std::string get_group(const std::string& group);

void add_options(clo_optionList_t* list);
void copy_options();
std::string get_source();
std::string get_history(int argc, char* argv[]);

void init_options(clo_optionList_t* list, const char* softwareVersion);
void read_options(clo_optionList_t* list, int argc, char* argv[]);

#endif
