/******************************************************************************
 *  NAME: DtLutNetcdf.h
 *
 *  DESCRIPTION: Object class that generates a netCDF4 LUT.
 *
 *  Created on: April 25, 2017
 *      Author: Sam Anderson
 *
 *
 ******************************************************************************/


#ifndef DtLutNetcdf_H_
#define DtLutNetcdf_H_

#include <netcdf>
#include <DDProcess.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

static const int LON_BINS = 360;
static const int LAT_BINS = 181;
static const int DATA_BINS = 55;
static const int GRIB_ARRAY_SIZE = LON_BINS*LAT_BINS*DATA_BINS*4;
static const int GRIB_ROW_SIZE = LON_BINS*LAT_BINS*4;
static const int NUM_CASES_SMALL = 4;
static const int NUM_CASES_BIG = 5;
static const int NUM_MOMENTS = 4;
static const int NUM_LUTS = 4;
static const int NUM_LATS = 180;
static const int NUM_LONS = 360;
static const int NUMCASES = 4;
static const int NUMCASEB = 5;
static const int NWAV = 7;
static const int NAOT = 6;
static const int NTH0 = 11;
static const int NTHET = 16;
static const int NPHI = 16;
static const int NUM_SOLUTIONS = 2;
static const int NLUTWAV = 4;
static const int NLTHET0 = 11;
static const int NLTHE = 15;
static const int NLPHI = 16;
static const int NLTAU = 7;
static const int NLTABLE = 5;
static const int NLETA = 13;
static const int NLSIZE = 2;
static const int DTABLE = 5;
static const int O3_COEFS = 2;
static const int H2O_COEFS = 3;
static const int IW550 = 1;

//GRIB Data LUT
struct GDAS
{
	float data[DATA_BINS][LAT_BINS][LON_BINS];
};

struct dtGribLUT
{
	float ugrd[LAT_BINS][LON_BINS];
	float vgrd[LAT_BINS][LON_BINS];
	float pwat[LAT_BINS][LON_BINS];
	float ozone[LAT_BINS][LON_BINS];
};

// Gas Correction data from LUT
constexpr int NUM_DT_BANDS = 10;

struct dtGasCorrectionLUT
{
    int     MBAND[NUM_DT_BANDS];
    int     VBAND[NUM_DT_BANDS];
    float   WAVE[NUM_DT_BANDS];
    float   MOL[NUM_DT_BANDS];
    float   OPT_O3_CLIM[NUM_DT_BANDS];
    float   OPT_H2O_CLIM[NUM_DT_BANDS];
    float   OPT_CO2_CLIM[NUM_DT_BANDS];
    float   O3_COEF[NUM_DT_BANDS][2];
    float   H2O_COEF[NUM_DT_BANDS][3];
};

// Land aerosol LUT
struct dtLandAerosolLUT
{
	int	  AEROSOL_ALL[NUM_SEASONS][NUM_LATS][NUM_LONS];
	float PHI_NL[NLPHI];
	float THE_NL[NLTHE];
	float THET0_NL[NLTHET0];
	float MU0_NL[NLTHET0];
	float WAV_NL[NLUTWAV];
	float OPTH_NL0[NLTABLE][NLUTWAV][NLTAU];
	float MASSCOEF_NL0[NLTABLE][NLUTWAV][NLTAU]; // Mass Coefficient
	float EXTNORM_NL0[NLTABLE][NLUTWAV][NLTAU]; // Extinction
	float SSA_NL0[NLTABLE][NLUTWAV][NLTAU];
	float QEXT_NL0[NLTABLE][NLUTWAV][NLTAU];
	float BEXT_NL0[NLTABLE][NLUTWAV][NLTAU];
	float VEXT_NL0[NLTABLE][NLUTWAV][NLTAU];
	float MEXT_NL0[NLTABLE][NLUTWAV][NLTAU];
	float SBAR_NL0[NLTABLE][NLUTWAV][NLTAU][NLTHET0];
	float INT_NL0[NLTABLE][NLUTWAV][NLTAU][NLTHET0][NLTHE][NLPHI];
	float Fd_NL0[NLTABLE][NLUTWAV][NLTAU][NLTHET0]; // Downward Flux
	float T_NL0[NLTABLE][NLUTWAV][NLTAU][NLTHET0][NLTHE]; // Atmospheric transmission
	float OMEGA0[NLTABLE][NLUTWAV][NLTAU];
	float ROD[NLUTWAV];
	float GOD[NLUTWAV];
	float OPTH_NL1[NLSIZE][NLUTWAV][NLTAU];
	float MASSCOEF_NL1[NLSIZE][NLUTWAV][NLTAU];
	float EXTNORM_NL1[NLSIZE][NLUTWAV][NLTAU];
	float SSA_NL1[NLSIZE][NLUTWAV][NLTAU];
	float QEXT_NL1[NLSIZE][NLUTWAV][NLTAU];
	float BEXT_NL1[NLSIZE][NLUTWAV][NLTAU];
	float VEXT_NL1[NLSIZE][NLUTWAV][NLTAU];
	float MEXT_NL1[NLSIZE][NLUTWAV][NLTAU];
	float SBAR_NL1[NLSIZE][NLUTWAV][NLTAU][NLTHET0];
	float INT_NL1[NLSIZE][NLUTWAV][NLTAU][NLTHET0][NLTHE][NLPHI];
	float Fd_NL1[NLSIZE][NLUTWAV][NLTAU][NLTHET0];
	float T_NL1[NLSIZE][NLUTWAV][NLTAU][NLTHET0][NLTHE];
};

// Ocean aerosol LUTs
struct dtOceanAerosolLUT
{
    int   JPHI[NPHI];
    float PHC[NPHI];
    float THET[NTHET];
    float THET0[NTH0];
    float WAVE[NWAV];
    float REF_RAYALL[NWAV][NUM_LUTS][NTH0][NTHET][NPHI];
    float EXTSMALL[NWAV][NUM_CASES_SMALL][NUM_LUTS];
    float RGSS[NUM_CASES_SMALL];
    float SIGMAS[NUM_CASES_SMALL];
    float MOMENTSSMALL[NUM_CASES_SMALL][NUM_MOMENTS][NUM_LUTS];
    float CCNSMALL[NUM_CASES_SMALL][NUM_LUTS];
    float BACKSCTTSMALL[NWAV][NUM_CASES_SMALL][NUM_LUTS];
    float ASSYMSMALL[NWAV][NUM_CASES_SMALL][NUM_LUTS];
    float ALBEDOSMALL[NWAV][NUM_CASES_SMALL][NUM_LUTS];
    float ALBEDO_R_SMALL[NWAV][NUM_CASES_SMALL][NAOT][NUM_LUTS][NTH0];
    float ALBEDO_T_SMALL[NWAV][NUM_CASES_SMALL][NAOT][NUM_LUTS][NTH0];
    float ALBEDO_R_RAY[NWAV][NTH0];
    float ALBEDO_T_RAY[NWAV][NTH0];
    float AINTS[NWAV][NUM_CASES_SMALL][NAOT][NUM_LUTS][NTH0][NTHET][NPHI];
    float TAUAS[NWAV][NUM_CASES_SMALL][NAOT];
    float EFFRADSMALL[NUM_CASES_SMALL];
    float EXTBIG[NWAV][NUM_CASES_BIG][NUM_LUTS];
    float RGSB[NUM_CASES_BIG];
    float SIGMAB[NUM_CASES_BIG];
    float MOMENTSBIG[NUM_CASES_BIG][NUM_MOMENTS][NUM_LUTS];
    float BACKSCTTBIG[NWAV][NUM_CASES_BIG][NUM_LUTS];
    float ASSYMBIG[NWAV][NUM_CASES_BIG][NUM_LUTS];
    float ALBEDOBIG[NWAV][NUM_CASES_BIG][NUM_LUTS];
    float ALBEDO_R_BIG[NWAV][NUM_CASES_BIG][NAOT][NUM_LUTS][NTH0];
    float ALBEDO_T_BIG[NWAV][NUM_CASES_BIG][NAOT][NUM_LUTS][NTH0];
    float AINTB[NWAV][NUM_CASES_BIG][NAOT][NUM_LUTS][NTH0][NTHET][NPHI];
    float TAUAB[NWAV][NUM_CASES_BIG][NAOT];
    float EFFRADBIG[NUM_CASES_BIG];
};

// Water Vapor LUTs
static const int TRANSM_H2O_TABLES = 6;
static const int TRANSM_H2O_ROWS = 220;
static const int TRANSM_H2O_VALS = 6;
static const int REFL_CH2_ROWS = 158080;
static const int REFL_CH2_VALS = 7;
static const int WEIGHT_VALS = 3;

struct dtWaterVaporLUT
{
    float TRANSM_H20[TRANSM_H2O_TABLES][TRANSM_H2O_ROWS][TRANSM_H2O_VALS];
    float WEIGHTS[TRANSM_H2O_ROWS][WEIGHT_VALS];
    float REFL_CH2[REFL_CH2_ROWS][REFL_CH2_VALS];
    float RATIO_CH19_TO_CH2[REFL_CH2_ROWS][REFL_CH2_VALS];
};

class DtLutNetcdf
{

public:

	// file attributes
	string lut_title_ ;
	string lut_prod_name_ ;

	// global attributes:
	string sensor_ ;
	string processing_version_;
	string Conventions_ ;
	string institution_;
	string license_ ;
	string naming_authority_;
	string date_created_ ;
	string keywords_vocabulary_;
	string stdname_vocabulary_ ;
	string creator_name_;
	string creator_email_;
	string creator_url_;
	string project_;
	string publisher_name_;
	string publisher_url_;
	string publisher_email_;
	string processing_level_;
	string cdm_data_type_;
	int	   orbit_number_;
	string history_;
	string source_files_;
	string time_coverage_start_;
	string time_coverage_end_;
	int    format_version_;
	int    instrument_number_;


	/**
	*  Class constructor
	*/

	DtLutNetcdf();

	/**
	*  Class destructor
	*/

	~DtLutNetcdf();

	/**
	*  Initialize L1A data
	*/

	int initialize();

	/**
	*  Create dark target aerosol netCDF4 LUT
	*/

	int create_dt_nc4_lut();

	/**
	*  Read aerosol netCDF4 LUT
	*/

	int read_grib_lut( dtGribLUT &grib_lut );
	int read_gas_correction_lut( dtGasCorrectionLUT &gc_lut );
	int read_land_aerosol_lut( dtLandAerosolLUT  &la_lut );
	int read_ocean_aerosol_lut( dtOceanAerosolLUT  &lo_lut );
	int read_water_vapor_lut( dtWaterVaporLUT &wv_lut );

	void setHistory(std::string history) {history_ = history;}
	std::string getHistory() {return history_;}

protected:

	/**
	 *  Read binary GRIB file.
	 */

	int read_grib_bin( const string filepath, dtGribLUT* grib_lut );

	/**
	 *  Read HDF GRIB file.
	 */

	int read_grib_hdf( const string filepath, dtGribLUT* grib_lut );

	/**
	 *  Read HDF ozone file.
	 */

	int read_ozone( const string filepath, dtGribLUT* grib_lut );

	/**
	 *  Read gas correction file data from text file.
	 */

	int read_gas_correction_file( dtGasCorrectionLUT* gc_lut );

	/**
	 *  Read land aerosol LUT file data from text file.
	 */

	int read_land_aerosol_file( const string groupname, int wnum,
								dtLandAerosolLUT* la_lut );

	/**
	 *  Read land aerosol map LUT file data from text file.
	 */

	int read_land_aerosol_map( const string groupname,
								dtLandAerosolLUT* la_lut );

	/**
	 *  Read big ocean aerosol LUT data from text file.
	 */

	int read_ocean_big_aerosol_file(const string groupname, int wnum,
									dtOceanAerosolLUT* boa_lut );

    /**
     *  Read small ocean aerosol LUT data from text file.
     */

    int read_ocean_small_aerosol_file(const string groupname, int wnum,
                                        dtOceanAerosolLUT* soa_lut );

    /**
     *  Read transmission H20 LUT data from text file.
     */

    int read_transm_h2o_file(const string groupname, int num,
                                        dtWaterVaporLUT* wv_lut );

    /**
     *  Read ch2 reflectance data from text file.
     */

    int read_ch2_reflectance_file( const string groupname,
                                             dtWaterVaporLUT* la_lut );

    /**
     *  Read Reads ch19-to-ch2 ratio data from text file.
     */

    int read_ch19_to_ch2_ratio_file( const string groupname,
                                             dtWaterVaporLUT* la_lut );

    /**
     *  Read weight table data from text file.
     */

    int read_weight_table_file( const string groupname,
                                             dtWaterVaporLUT* la_lut );

    /**
	* Write global attributes to file
	*/

	int write_global_attributes( NcFile* nc_output );

	/**
	* Determine if platform is little endian
	*/

	bool isPlatformLittleEndian();

	/**
	 *  Write binary GRIB file.
	 */

	int write_grib_lut( NcFile* nc_output, dtGribLUT* grib_lut );

	/**
	 *  Write gas correction file data to netcdf LUT.
	 */

	int write_gas_correction_lut( NcFile* nc_output, dtGasCorrectionLUT* gc_lut );

	/**
	 *  Write land aerosol LUT file data to netcdf LUT.
	 */

	int write_land_aerosol_lut( NcFile* nc_output, dtLandAerosolLUT* la_lut  );

    /**
     *  Write ocean aerosol LUT data to netcdf LUT.
     */

    int write_ocean_aerosol_lut( NcFile* nc_output, dtOceanAerosolLUT* oa_lut );

    /**
     *  Write water vapor LUT data to netcdf LUT.
     */

    int write_water_vapor_lut( NcFile* nc_output, dtWaterVaporLUT* oa_lut );


	// Dimensions

	NcDim 	scalar_dim_;
	NcDim 	num_grib_data_bins_dim_;
	NcDim 	num_grib_lat_bins_dim_;
	NcDim 	num_grib_lon_bins_dim_;
	NcDim 	num_gc_dt_bands_dim_;
	NcDim 	num_gc_O3_coef_dim_;
	NcDim 	num_gc_H2O_coef_dim_;
	NcDim 	num_land_lats_dim_;
	NcDim 	num_land_lons_dim_;
	NcDim 	num_land_phi_dim_;
	NcDim 	num_land_the_dim_;
	NcDim 	num_land_thet0_dim_;
	NcDim 	num_land_tau_dim_;
	NcDim 	num_land_wav_dim_;
	NcDim 	num_land_table_dim_;
	NcDim 	num_land_size_dim_;
	NcDim 	num_land_season_dim_;
	NcDim 	num_ocean_phi_dim_;
	NcDim 	num_ocean_the_dim_;
	NcDim 	num_ocean_thet0_dim_;
	NcDim 	num_ocean_tau_dim_;
	NcDim 	num_ocean_wave_dim_;
	NcDim 	num_ocean_cases_dim_;
	NcDim 	num_ocean_caseb_dim_;
	NcDim 	num_ocean_wslut_dim_;
	NcDim 	num_ocean_moments_dim_;
	NcDim   num_h2o_tables_dim_;
	NcDim   num_h2o_rows_dim_;
	NcDim   num_h2o_vals_dim_;
	NcDim   num_weight_vals_dim_;
	NcDim   num_ch2_rows_dim_;
	NcDim   num_ch2_vals_dim_;

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
void DtLutNetcdf::byteSwap(T& aValue)
{
	T tempValue = aValue;  // a temporary copy of the value

	// Pointers to the first byte of both variables
	unsigned char* aValuePtr = reinterpret_cast<unsigned char*>(&aValue);
	unsigned char* tempValuePtr = reinterpret_cast<unsigned char*>(&tempValue);

	// Swap the byte order
	for (unsigned int byte = 0; byte < sizeof(aValue); ++byte)
	{
		aValuePtr[byte] = tempValuePtr[(sizeof(aValue) - 1) - byte];
	}
}


#endif /* DtLutNetcdf_H_ */
