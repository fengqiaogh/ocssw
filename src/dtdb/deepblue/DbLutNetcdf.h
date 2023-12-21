/******************************************************************************
 *  NAME: DbLutNetcdf.h
 *
 *  DESCRIPTION: Object class that generates a netCDF4 LUT.
 *
 *  Created on: April 25, 2017
 *      Author: Sam Anderson
 *
 ******************************************************************************/

#ifndef DbLutNetcdf_H_
#define DbLutNetcdf_H_

#include <boost/multi_array.hpp>
#include <netcdf>
#include <DDProcess.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

static const int NDET = 10;          // number of detectors per channel
static const int NAOI = 6;          // number of mirror incidence angles in LUT
static const int NSIDE = 2;          // number of mirror sides
static const int NBANDLUT = 3;       // number of bands in LUT (incl. TDI bands)
// 11 bands = 8,9,10,11,12,13L,13H,14L,14H,15,16
static const int NBANDS = 3; // number of bands in polar. correction from MODIS, band 1,3,8
//   static const int ntimer = 200;     // number of times (~months) for which polarization coeff. were derived(C051), Terra
static const int NTIMER = 138; // number of times for which Terra/MODIS polar. coeff. were derived(C006), 20100928
static const int NTIMER_AQ = 121; // number of times for which Aqua/MODIS polar. coeff. were derived(C006), 20101012
static const int NBANDLUTR = 3;    // number of bands in LUT, 9 bands (B08-B16)
static const int NCOEFF = 6; // number of polynomial coefficients (3rd order polynomical)

extern const int chindx[8]; //{ 1,2,3,4,5,6,10,11 }
extern const int bindx[3]; //{ 1,3,8 } DeepBlue targetted bands in MODIS
extern const string str_season[NUM_SEASONS];

struct dbLUT {
    short xxyear[NTIMER]; // year and DoY for which the pol. coeff. were derived
    short xxday[NTIMER];
    short xxyear_aq[NTIMER_AQ];
    short xxday_aq[NTIMER_AQ];
    float angle_of_incident[NAOI];
    float am12[NAOI];
    float am13[NSIDE][NAOI][NDET][NBANDLUT];
    float xxwave[NBANDLUTR]; // wavelengths for ocean bands
    float xxwave_aq[NBANDLUTR];
    double sectab[NTIMER]; // reference time for LUT (seconds)
    double xxrvs[NCOEFF][NDET][NSIDE][NBANDLUTR][NTIMER]; // pol. coeff. for ocn bands
    double xxam12[NCOEFF][NDET][NSIDE][NBANDLUTR][NTIMER];
    double xxam13[NCOEFF][NDET][NSIDE][NBANDLUTR][NTIMER];
    double xxrvs_aq[NCOEFF][NDET][NSIDE][NBANDLUTR][NTIMER_AQ];
    double sectab_aq[NTIMER_AQ];
};

static const int ENTRIES_412 = 121;
static const int ENTRIES_470 = 104;
static const int XCAL_COEF = 6;
static const int MODIS_DETECTORS = 10;
static const int MIRROR_SIDES = 2;
static const int NLATS = 180;
static const int NLONS = 360;
static const int NSEASONS = 4;
static const int NMONTHS = 12;
static const int NNDVI = 3;
static const int NTERMS = 4;

struct dbModisCorrectionsLUT {
    short YEAR_412[ENTRIES_412];
    short DAY_412[ENTRIES_412];
    short YEAR_470[ENTRIES_470];
    short DAY_470[ENTRIES_470];
    float GAIN_412[ENTRIES_412];
    float GAIN_470[ENTRIES_470];
    double XCAL_412_M11[XCAL_COEF][MODIS_DETECTORS][MIRROR_SIDES][ENTRIES_412-1];
    double XCAL_412_M12[XCAL_COEF][MODIS_DETECTORS][MIRROR_SIDES][ENTRIES_412-1];
    double XCAL_412_M13[XCAL_COEF][MODIS_DETECTORS][MIRROR_SIDES][ENTRIES_412-1];
    double XCAL_470_M11[XCAL_COEF][MODIS_DETECTORS][MIRROR_SIDES][ENTRIES_470-1];
    double XCAL_470_M12[XCAL_COEF][MODIS_DETECTORS][MIRROR_SIDES][ENTRIES_470-1];
    double XCAL_470_M13[XCAL_COEF][MODIS_DETECTORS][MIRROR_SIDES][ENTRIES_470-1];
};

struct dbModisSurfReflLUT {
    float SR412_ALL[NSEASONS][NLATS*10][NLONS*10];
    float SR412_FWD[NSEASONS][NLATS*10][NLONS*10];
    float SR470_ALL[NSEASONS][NLATS*10][NLONS*10];
    float SR470_FWD[NSEASONS][NLATS*10][NLONS*10];
    float SR650_ALL[NSEASONS][NLATS*10][NLONS*10];
    float SR650_FWD[NSEASONS][NLATS*10][NLONS*10];
    float SR865_ALL[NSEASONS][NLATS*10][NLONS*10];
};

struct dbModisSurfReflLimited {
    boost::multi_array<float,2> SR412_ALL_L;
    boost::multi_array<float,2> SR412_FWD_L;
    boost::multi_array<float,2> SR470_ALL_L;
    boost::multi_array<float,2> SR470_FWD_L;
    boost::multi_array<float,2> SR650_ALL_L;
    boost::multi_array<float,2> SR650_FWD_L;
    boost::multi_array<float,2> SR865_ALL_L;
};

struct dbModisSwirVsVisLUT {
    float latitude[NLATS*10][NLONS*10];
    float longitude[NLATS*10][NLONS*10];
    float coeffs_2130_to_412[NLATS*10][NLONS*10];
    float coeffs_2130_to_470[NLATS*10][NLONS*10];
    float min_2130_for_412[NLATS*10][NLONS*10];
    float max_2130_for_412[NLATS*10][NLONS*10];
    float min_2130_for_470[NLATS*10][NLONS*10];
    float max_2130_for_470[NLATS*10][NLONS*10];
    float data_num_total[NLATS*10][NLONS*10];
    float data_num_fitting[NLATS*10][NLONS*10];
    float stderr_412[NLATS*10][NLONS*10];
    float stderr_470[NLATS*10][NLONS*10];
};
// Surface Pressure LUT
static const int SP_720 = 720;
static const int SP_360 = 360;
static const int SP_SETS = 6;

struct dbSurfacePressureLUT {
    float PS[SP_720][SP_360];
    float SURFACE_PRESSURE[SP_SETS * SP_720][SP_SETS * SP_360];
    float SURFACE_ELEVATION[SP_SETS * SP_720][SP_SETS * SP_360];
};

// Nvalx LUT
static const int NSZAV = 10;
static const int NSZA = 12;
static const int NVZA = 46;
static const int NRAA = 30;
static const int NTAU = 10;
static const int SSA412 = 8;
static const int SSA470 = 4;
static const int SSA650 = 1;
static const int SR412 = 20;
static const int SR470 = 24;
static const int SR650 = 24;

// Nvalx LUT

struct dbNvalxLUT {
    float NVALX_412[SR412][SSA412][NTAU][NRAA][NVZA][NSZA];
    float NVALX_470[SR470][SSA470][NTAU][NRAA][NVZA][NSZA];
    float NVALX_650[SR650][NTAU][NRAA][NVZA][NSZA];
};

static const int NSTOKES  = 3;
static const int NRRAA    = 31;

struct dbRayleighLUT {
    float RAYL_412[NRRAA][NVZA][NSZA][NSTOKES];
    float RAYL_470[NRRAA][NVZA][NSZA][NSTOKES];
    float RAYL_650[NRRAA][NVZA][NSZA][NSTOKES];
};


// MODIS/SeaWifs LER Tables LUT

static const int MTABLE_20800 = 20800;
static const int MTABLE_260 = 260;
static const int MTABLE_160 = 160;
static const int MTABLE_2 = 2;

struct dbTablesLUT {
    float LOGI0[MTABLE_20800];
    float Z1I0[MTABLE_20800];
    float Z2I0[MTABLE_20800];
    float TI0[MTABLE_20800];
    float SB[MTABLE_260];
    float LOGI0R[MTABLE_160];
    float Z1I0R[MTABLE_160];
    float Z2I0R[MTABLE_160];
    float TI0R[MTABLE_160];
    float SBR[MTABLE_2];
};

struct dbVeg_21sfcLUT {
    float NVALX21_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float R0X21_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float SX21_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float TX21_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float NVALX672_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float R0X672_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float SX672_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float TX672_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float NVALX865_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float R0X865_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float SX865_SFC[NSEASONS][NRAA][NVZA][NSZAV];
    float TX865_SFC[NSEASONS][NRAA][NVZA][NSZAV];
};

struct dbViirsSurfReflLUT {
    float SR412_ALL[NSEASONS][NLATS*10][NLONS*10];
    float SR488_ALL[NSEASONS][NLATS*10][NLONS*10];
    float SR670_ALL[NSEASONS][NLATS*10][NLONS*10];
    float BRDF_650[NLATS*10][NLONS*10];
};

struct dbViirsSurfReflLimited {
    boost::multi_array<float,2> SR412_ALL_L;
    boost::multi_array<float,2> SR488_ALL_L;
    boost::multi_array<float,2> SR670_ALL_L;
    boost::multi_array<float,2> BRDF_650_L;
};

struct dbSurfCoeffLUT {
    float SC412_ALL[NSEASONS][NNDVI][NTERMS][NLATS*10][NLONS*10];
    float SC412_FWD[NSEASONS][NNDVI][NTERMS][NLATS*10][NLONS*10];
    float SC470_ALL[NSEASONS][NNDVI][NTERMS][NLATS*10][NLONS*10];
    float SC470_FWD[NSEASONS][NNDVI][NTERMS][NLATS*10][NLONS*10];
    float SC650_ALL[NSEASONS][NNDVI][NTERMS][NLATS*10][NLONS*10];
    float SC650_FWD[NSEASONS][NNDVI][NTERMS][NLATS*10][NLONS*10];
};

struct dbSurfCoeffLimited {
    boost::multi_array<float,4> SC412_ALL_L;
    boost::multi_array<float,4> SC412_FWD_L;
    boost::multi_array<float,4> SC470_ALL_L;
    boost::multi_array<float,4> SC470_FWD_L;
    boost::multi_array<float,4> SC650_ALL_L;
    boost::multi_array<float,4> SC650_FWD_L;
};

struct dbLandcoverLUT {
    int   VEGETATION[4*NLATS*10][NLONS*10];
    short IGBP[NLATS*10][NLONS*10];
    short REGION_INDEX[NLATS*10][NLONS*10];
    float DESERTS_FLAG[NLATS*10][NLONS*10];
};

struct dbGeozoneLUT {
    float GEOZONE_FLAG[NLATS*10][NLONS*10];
    float ELEVATION_STDV[NLATS*10][NLONS*10];
    float BACKGROUND_AOD[4][NLATS][NLONS];
};

static const int NSCOEF = 3;
static const int NSLATS = 3000;
static const int NSLONS = 6000;
struct dbViirsSwirVsVisLUT {
    float latitude[NSEASONS][NSLATS][NSLONS];
    float longitude[NSEASONS][NSLATS][NSLONS];
    float coeffs_2250_to_412[NSEASONS][NSCOEF][NSLATS][NSLONS];
    float coeffs_2250_to_488[NSEASONS][NSCOEF][NSLATS][NSLONS];
    float coeffs_2250_to_670[NSEASONS][NSCOEF][NSLATS][NSLONS];
    float min_2250_for_412[NSEASONS][NSLATS][NSLONS];
    float max_2250_for_412[NSEASONS][NSLATS][NSLONS];
    float min_2250_for_488[NSEASONS][NSLATS][NSLONS];
    float max_2250_for_488[NSEASONS][NSLATS][NSLONS];
    float min_2250_for_670[NSEASONS][NSLATS][NSLONS];
    float max_2250_for_670[NSEASONS][NSLATS][NSLONS];
    float data_num_total[NSEASONS][NSLATS][NSLONS];
    float data_num_fitting[NSEASONS][NSLATS][NSLONS];
    float stderr_412[NSEASONS][NSLATS][NSLONS];
    float stderr_488[NSEASONS][NSLATS][NSLONS];
    float stderr_670[NSEASONS][NSLATS][NSLONS];
};

static const int NVSZA = 22;
static const int NVVZA = 20;
static const int NVRAA = 21;
static const int NAOT1 = 14;
static const int NAOT2 = 7;
static const int NFMF1 = 7;
static const int NFMF2 = 5;
static const int NFMF3 = 8;
static const int NFMF4 = 4;
static const int NWS = 6;
static const int NCHL = 4;
static const int NDBOWL = 7;
static const int NDBMDL = 4;
static const int NDBLWL = 3;

struct dbOceanAerosolLUMA {
    const static size_t nsza = NVSZA;
    const static size_t nvza = NVVZA;
    const static size_t nraa = NVRAA;
    size_t naot;
    size_t nfmf;
    const static size_t nwspd = NWS;
    const static size_t nchl = NCHL;
    const static size_t nwave = NDBOWL;
    boost::multi_array<float,7> m03;
    boost::multi_array<float,7> m04;
    boost::multi_array<float,7> m05;
    boost::multi_array<float,7> m07;
    boost::multi_array<float,6> m08;
    boost::multi_array<float,6> m10;
    boost::multi_array<float,6> m11;
    boost::multi_array<float,1> sza;
    boost::multi_array<float,1> vza;
    boost::multi_array<float,1> raa;
    boost::multi_array<float,1> aot550;
    boost::multi_array<float,1> fmf;
    boost::multi_array<float,1> wspd;
    boost::multi_array<float,1> chl;
    boost::multi_array<float,1> wave;
    boost::multi_array<float,2> ae;
    boost::multi_array<float,3> aot;
    boost::multi_array<float,3> fine_aot;
    boost::multi_array<float,3> coarse_aot;
    boost::multi_array<float,3> ssa;
    boost::multi_array<float,3> fine_ssa;
    boost::multi_array<float,3> coarse_ssa;
    boost::multi_array<float,3> asy;
    boost::multi_array<float,3> fine_asy;
    boost::multi_array<float,3> coarse_asy;
};

struct dbOceanAerosolLUT {
    const static size_t nsza = NVSZA;
    const static size_t nvza = NVVZA;
    const static size_t nraa = NVRAA;
    size_t naot;
    size_t nfmf;
    const static size_t nwspd = NWS;
    const static size_t nchl = NCHL;
    const static size_t nwave = NDBOWL;
    float m03[NCHL][NWS][NFMF3][NAOT1][NVRAA][NVVZA][NVSZA];
    float m04[NCHL][NWS][NFMF3][NAOT1][NVRAA][NVVZA][NVSZA];
    float m05[NCHL][NWS][NFMF3][NAOT1][NVRAA][NVVZA][NVSZA];
    float m07[NCHL][NWS][NFMF3][NAOT1][NVRAA][NVVZA][NVSZA];
    float m08[NWS][NFMF3][NAOT1][NVRAA][NVVZA][NVSZA];
    float m10[NWS][NFMF3][NAOT1][NVRAA][NVVZA][NVSZA];
    float m11[NWS][NFMF3][NAOT1][NVRAA][NVVZA][NVSZA];
    float sza[NVSZA];
    float vza[NVVZA];
    float raa[NVRAA];
    float aot550[NAOT1];
    float fmf[NFMF3];
    float wspd[NWS];
    float chl[NCHL];
    float wave[NDBOWL];
    float ae[NFMF3][NAOT1];
    float aot[NDBOWL][NFMF3][NAOT1];
    float fine_aot[NDBOWL][NFMF3][NAOT1];
    float coarse_aot[NDBOWL][NFMF3][NAOT1];
    float ssa[NDBOWL][NFMF3][NAOT1];
    float fine_ssa[NDBOWL][NFMF3][NAOT1];
    float coarse_ssa[NDBOWL][NFMF3][NAOT1];
    float asy[NDBOWL][NFMF3][NAOT1];
    float fine_asy[NDBOWL][NFMF3][NAOT1];
    float coarse_asy[NDBOWL][NFMF3][NAOT1];
};


struct dbLandAerosolLUT {
    float SZA412_Nodes[NSZA];
    float VZA412_Nodes[NVZA];
    float RAA412_Nodes[NRAA];
    float AOT412_Nodes[NTAU];
    float SSA412_Nodes[SSA412];
    float SR412_Nodes[SR412];
    float nvalx412[SR412][SSA412][NTAU][NRAA][NVZA][NSZA];
    float SZA488_Nodes[NSZA];
    float VZA488_Nodes[NVZA];
    float RAA488_Nodes[NRAA];
    float AOT488_Nodes[NTAU];
    float SSA488_Nodes[SSA470];
    float SR488_Nodes[SR470];
    float nvalx488[SR470][SSA470][NTAU][NRAA][NVZA][NSZA];
    float SZA672_Nodes[NSZA];
    float VZA672_Nodes[NVZA];
    float RAA672_Nodes[NRAA];
    float AOT672_Nodes[NTAU];
    float SSA672_Nodes[SSA650];
    float SR672_Nodes[SR650];
    float nvalx672[SR650][SSA650][NTAU][NRAA][NVZA][NSZA];
};

static const int NUMX = 21601;
static const int NUMY = 10801;

struct dbBathymetryLUT {
    double x[NUMX];
    double y[NUMY];
    int    z[NUMY][NUMX];
};

struct dbChlLUT {
    float  longitude[NLATS*10][NLONS*10];
    float  latitude[NLATS*10][NLONS*10];
    double time[NMONTHS];
    float  log_chl[NLATS*10][NLONS*10][NMONTHS];
};


static const int NDAYS = 95;

struct dbXCalLUT {
//    int GEOZONE[NLATS*10][NLONS*10];
};

class DbLutNetcdf {

public:

    // file attributes
    string lut_title_;
    string lut_prod_name_;

    // global attributes:
    string sensor_;
    string processing_version_;
    string Conventions_;
    string institution_;
    string license_;
    string naming_authority_;
    string date_created_;
    string keywords_vocabulary_;
    string stdname_vocabulary_;
    string creator_name_;
    string creator_email_;
    string creator_url_;
    string project_;
    string publisher_name_;
    string publisher_url_;
    string publisher_email_;
    string processing_level_;
    string cdm_data_type_;
    int orbit_number_;
    string history_;
    string source_files_;
    string time_coverage_start_;
    string time_coverage_end_;
    int format_version_;
    int instrument_number_;

    /**
     *  Class constructor
     */

    DbLutNetcdf();

    /**
     *  Class destructor
     */

    ~DbLutNetcdf();

    /**
     *  Initialize L1A data
     */

    int initialize();

    /**
     *  Create dark target aerosol netCDF4 LUT
     */

    int create_db_nc4_lut();

    /**
     *  Read aerosol netCDF4 LUT
     */

    int read_tables_lut(dbTablesLUT*  lut);
    int read_surface_pressure_lut(dbSurfacePressureLUT*  lut);
    int read_nvalx_lut(dbNvalxLUT*  lut);
    int read_veg_21sfc_lut(dbVeg_21sfcLUT*  lut);
    int read_swir_lut(dbViirsSwirVsVisLUT*  lut);
    int read_modis_surf_refl_lut(dbModisSurfReflLimited* lut, int* start,
            int* edge, int& season, int &dateline);
    int read_viirs_surf_refl_lut(dbViirsSurfReflLimited* lut, int* start,
            int* edge, int& season, int &dateline);
    int read_surf_coeff_lut(dbSurfCoeffLimited* lut, int* start,
            int* edge, int&  season, int &dateline);
    int read_landcover_lut(dbLandcoverLUT*  lut);
    int read_ocean_aero_lut(dbOceanAerosolLUMA* lut, const string sType);
    int read_land_aero_lut(dbLandAerosolLUT* lut, const string sType);
    int read_bathymetry_lut(dbBathymetryLUT* lut);
    int read_chl_lut(dbChlLUT* lut);
    int read_geozone_lut(dbGeozoneLUT* lut);
    int read_rayleigh_lut(dbRayleighLUT*  lut);

    /**
     *  Set/Get History
     */

    void setHistory(std::string history) {
        history_ = history;
    }
    std::string getHistory() {
        return history_;
    }

protected:

    /**
     *  Read LUT files.
     */

    int read_tables_file(dbTablesLUT* lut);
    int read_surface_pressure_file(dbSurfacePressureLUT* lut);
    int read_nvalx_files(dbNvalxLUT* lut);
    int read_veg_21sfc_files(dbVeg_21sfcLUT* lut);
    int read_swir_file( dbViirsSwirVsVisLUT* lut );
    int read_modis_surf_refl_files(dbModisSurfReflLUT* lut);
    int read_viirs_surf_refl_files(dbViirsSurfReflLUT* lut);
    int read_surf_coeff_files(dbSurfCoeffLUT* lut);
    int read_landcover_files(dbLandcoverLUT* lut);
    int read_ocean_aero_file( dbOceanAerosolLUT* lut, const string strType );
    int read_land_aero_file( dbLandAerosolLUT* lut, const string strType );
    int read_bathymetry_files( dbBathymetryLUT* lut );
    int read_chl_files( dbChlLUT* lut );
    int read_geozone_files( dbGeozoneLUT* lut );
    int read_rayleigh_files(dbRayleighLUT* lut);

    /**
     *  Write LUT files to NetCDF4 lut file.
     */

    int write_tables_lut(NcFile* nc_output, dbTablesLUT* lut);
    int write_surface_pressure_lut(NcFile* nc_output,
            dbSurfacePressureLUT* mt_lut);
    int write_nvalx_lut(NcFile* nc_output, dbNvalxLUT* lut);
    int write_veg_21sfc_lut(NcFile* nc_output, dbVeg_21sfcLUT* lut);
    int write_swir_lut(NcFile* nc_output, dbViirsSwirVsVisLUT* lut);
    int write_modis_surf_refl_lut(NcFile* nc_output, dbModisSurfReflLUT* lut);
    int write_viirs_surf_refl_lut(NcFile* nc_output, dbViirsSurfReflLUT* lut);
    int write_surf_coeff_lut(NcFile* nc_output, dbSurfCoeffLUT* lut);
    int write_landcover_lut(NcFile* nc_output, dbLandcoverLUT* lut);
    int write_ocean_aero_lut(NcFile* nc_output, dbOceanAerosolLUT* lut,
            const string strType);
    int write_land_aero_lut(NcFile* nc_output, dbLandAerosolLUT* lut,
            const string strType);
    int write_bathymetry_lut(NcFile* nc_output, dbBathymetryLUT* lut);
    int write_chl_lut(NcFile* nc_output, dbChlLUT* lut);
    int write_geozone_lut(NcFile* nc_output, dbGeozoneLUT* lut);
    int write_rayleigh_lut(NcFile* nc_output, dbRayleighLUT* lut);

    /**
     * Write global attributes to file
     */

    int write_global_attributes(NcFile* nc_output);

    /**
     * Determine if platform is little endian
     */

    bool isPlatformLittleEndian();

    // Dimensions

    NcDim scalar_dim_;
    NcDim dim_seasons_;
    NcDim dim_months_;
    NcDim dim_20800_;
    NcDim dim_260_;
    NcDim dim_160_;
    NcDim dim_2_;
    NcDim dim_3_;
    NcDim dim_4_;
    NcDim dim_720_;
    NcDim dim_180_;
    NcDim dim_360_;
    NcDim dim_4320_;
    NcDim dim_2160_;
    NcDim dim_nszav_;
    NcDim dim_nsza_;
    NcDim dim_nvza_;
    NcDim dim_nraa_;
    NcDim dim_ntau_;
    NcDim dim_nssa_;
    NcDim dim_nsr_;
    NcDim dim_ssa412_;
    NcDim dim_ssa470_;
    NcDim dim_sr412_;
    NcDim dim_sr470_;
    NcDim dim_sr650_;
    NcDim dim_ndvi_;
    NcDim dim_terms_;
    NcDim dim_3000_;
    NcDim dim_6000_;
    NcDim dim_1800_;
    NcDim dim_3600_;
    NcDim dim_7200_;
    NcDim dim_nrraa_;
    NcDim dim_nstokes_;
    NcDim dim_xcoeff_;
    NcDim dim_det_;
    NcDim dim_ms_;
    NcDim dim_entries_412_;
    NcDim dim_entries_470_;
    NcDim dim_entries_412x_;
    NcDim dim_entries_470x_;
    NcDim dim_nfmf_;
    NcDim dim_nws_;
    NcDim dim_nchl_;
    NcDim dim_nwl_;
    NcDim dim_nx_;
    NcDim dim_ny_;
    NcDim dim_nz_;

    /**
     * Converts the endianness of the parameter by performing the appropriate
     * byte swapping.
     *
     */

    template<typename T>
    static void byteSwap(T& aValue);

};

//---------------------------------------------------------------------------
// Converts the endianness of the parameter by performing the appropriate
// byte swapping.
//---------------------------------------------------------------------------
template<typename T>
void DbLutNetcdf::byteSwap(T& aValue) {
    T tempValue = aValue;  // a temporary copy of the value

    // Pointers to the first byte of both variables
    unsigned char* aValuePtr = reinterpret_cast<unsigned char*>(&aValue);
    unsigned char* tempValuePtr = reinterpret_cast<unsigned char*>(&tempValue);

    // Swap the byte order
    for (unsigned int byte = 0; byte < sizeof(aValue); ++byte) {
        aValuePtr[byte] = tempValuePtr[(sizeof(aValue) - 1) - byte];
    }
}

#endif /* DtLutNetcdf_H_ */
