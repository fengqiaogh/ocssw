/******************************************************************************
 *  NAME: DbLutNetcdf.cpp
 *
 *  DESCRIPTION: Object class that generates a netCDF4 LUT for NASA Deep Blue
 *  aerosols algorithm
 *
 *  Created on: April 25, 2017
 *      Author: Sam Anderson
 *
 *
 ******************************************************************************/


#include "deepblue/DbLutNetcdf.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <mfhdf.h>

#include <libgen.h>

#include <DDAlgorithm.h>
#include <DDOptions.h>

const int chindx[8] = { 1,2,3,4,5,6,10,11 };
const int bindx[3] = { 1,3,8 };// DeepBlue targetted bands in MODIS
const std::string str_season[NUM_SEASONS] = {"winter", "spring", "summer", "fall"};
/**************************************************************************
 * NAME: DbLutNetcdf()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DbLutNetcdf::DbLutNetcdf() {
}

/**************************************************************************
 * NAME: ~DbLutNetcdf()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DbLutNetcdf::~DbLutNetcdf() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Initializes data and object classes for granule
 *
 *************************************************************************/

int DbLutNetcdf::initialize() {
    return DTDB_SUCCESS;
}


/**************************************************************************
 * NAME: create_db_nc4_lut()
 *
 * DESCRIPTION: Create deep blue aerosol netCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::create_db_nc4_lut() {

	int status = DTDB_SUCCESS;
	int istatus = DTDB_SUCCESS;
	int num_good = 0;

	string path = get_option(NETCDF_LUT_PATH);
	string filepath = path + "/VIIRS_DEEPBLUE_LUT_" + "version_source" + ".nc";

	NcFile* nc_output;

	try {
		nc_output = new NcFile( filepath, NcFile::replace );
	}
	catch( NcException& e) {
		e.what();
		cerr << "DbLutNetcdf:: Failure creating netCDF4 LUT file: " + filepath + "." << endl;
		return DTDB_FAIL;
	}
	nc_output->putAtt( "title", "VIIRS DEEP BLUE LUTs" );

	write_global_attributes( nc_output );

// Read input files and create LUTs
    cerr << "DbLutNetcdf:: Begin moving LUTs to NetCDF4 LUT file"  << endl;

    dbOceanAerosolLUT* oa_fine = new dbOceanAerosolLUT;
    istatus = read_ocean_aero_file( oa_fine, INPUT_AERO_OCEAN_FINE );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading VIIRS Ocean Aerosol LUT file " << endl;
        status = istatus;
    }
    istatus = write_ocean_aero_lut( nc_output, oa_fine, LUT_OCEAN_AEROSOL_FINE );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing VIIRS Ocean Aerosol data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete oa_fine;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_OCEAN_AEROSOL_FINE  << endl;
    }
    dbOceanAerosolLUT* oa_dust = new dbOceanAerosolLUT;
    istatus = read_ocean_aero_file( oa_dust, INPUT_AERO_OCEAN_DUST );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading VIIRS Ocean Aerosol LUT file " << endl;
        status = istatus;
    }
    istatus = write_ocean_aero_lut( nc_output, oa_dust, LUT_OCEAN_AEROSOL_DUST );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing VIIRS Ocean Aerosol data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete oa_dust;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_OCEAN_AEROSOL_DUST  << endl;
    }
    dbOceanAerosolLUT* oa_mari = new dbOceanAerosolLUT;
    istatus = read_ocean_aero_file( oa_mari, INPUT_AERO_OCEAN_MARI );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading VIIRS Ocean Aerosol LUT file " << endl;
        status = istatus;
    }
    istatus = write_ocean_aero_lut( nc_output, oa_mari, LUT_OCEAN_AEROSOL_MARI );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing VIIRS Ocean Aerosol data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete oa_mari;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
       cerr << "DbLutNetcdf:: Created " + LUT_OCEAN_AEROSOL_MARI  << endl;
    }
    dbOceanAerosolLUT* oa_mix = new dbOceanAerosolLUT;
    istatus = read_ocean_aero_file( oa_mix, INPUT_AERO_OCEAN_MIX );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading VIIRS Ocean Aerosol LUT file " << endl;
        status = istatus;
    }
    istatus = write_ocean_aero_lut( nc_output, oa_mix, LUT_OCEAN_AEROSOL_MIX );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing VIIRS Ocean Aerosol data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete oa_mix;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_OCEAN_AEROSOL_MIX  << endl;
    }
    dbLandAerosolLUT* la_lut = new dbLandAerosolLUT;
    istatus = read_land_aero_file( la_lut, INPUT_AERO_LAND_FINE );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading Land Aerosol LUT file " << endl;
        status = istatus;
    }
    istatus = write_land_aero_lut( nc_output, la_lut,
            LUT_LAND_AEROSOL_FINE );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing Land Aerosol data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete la_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_LAND_AEROSOL_FINE  << endl;
    }
    la_lut = new dbLandAerosolLUT;
    istatus = read_land_aero_file( la_lut, INPUT_AERO_LAND_DUST );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading Land Aerosol LUT file " << endl;
        status = istatus;
    }
    istatus = write_land_aero_lut( nc_output, la_lut,
            LUT_LAND_AEROSOL_DUST );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing Land Aerosol data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete la_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_LAND_AEROSOL_DUST  << endl;
    }
    dbBathymetryLUT* b_lut = new dbBathymetryLUT;
    istatus = read_bathymetry_files( b_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading bathymetry LUT file " << endl;
        status = istatus;
    }
    istatus = write_bathymetry_lut( nc_output, b_lut  );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing bathymetry data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete b_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_BATHYMETRY  << endl;
    }
    dbChlLUT* c_lut = new dbChlLUT;
    istatus = read_chl_files( c_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading CHL LUT file " << endl;
        status = istatus;
    }
    istatus = write_chl_lut( nc_output, c_lut  );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing CHL data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete c_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_CHL  << endl;
    }
    dbTablesLUT* t_lut = new dbTablesLUT;
    istatus = read_tables_file( t_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading LER tables file " << endl;
        status = istatus;
    }
    istatus = write_tables_lut( nc_output, t_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing LER tables to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete t_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_LER_TABLES  << endl;
    }
    dbLandcoverLUT* lc_lut = new dbLandcoverLUT;
    istatus = read_landcover_files( lc_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading Landcover LUT file " << endl;
        status = istatus;
    }
    istatus = write_landcover_lut( nc_output, lc_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing Landcover data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete lc_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_LANDCOVER  << endl;
    }
    dbGeozoneLUT* gz_lut = new dbGeozoneLUT;
    istatus = read_geozone_files( gz_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading geozone LUT file " << endl;
        status = istatus;
    }
    istatus = write_geozone_lut( nc_output, gz_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing geozone data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete gz_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_GEOZONE  << endl;
    }
    dbModisSurfReflLUT* msr_lut = new dbModisSurfReflLUT;
    istatus = read_modis_surf_refl_files( msr_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading surface reflectance LUT file " << endl;
        status = istatus;
    }
    istatus = write_modis_surf_refl_lut( nc_output, msr_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing surface reflectance data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete msr_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
       cerr << "DbLutNetcdf:: Created " + LUT_MODIS_SURFACE_REFL  << endl;
    }
    dbViirsSurfReflLUT* vsr_lut = new dbViirsSurfReflLUT;
    istatus = read_viirs_surf_refl_files( vsr_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading surface reflectance LUT file " << endl;
        status = istatus;
    }
    istatus = write_viirs_surf_refl_lut( nc_output, vsr_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing surface reflectance data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete vsr_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
       cerr << "DbLutNetcdf:: Created " + LUT_VIIRS_SURFACE_REFL  << endl;
    }
    dbSurfCoeffLUT* sc_lut = new dbSurfCoeffLUT;
    istatus = read_surf_coeff_files( sc_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading surface coefficients LUT file " << endl;
        status = istatus;
    }
    istatus = write_surf_coeff_lut( nc_output, sc_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing surface coefficients data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete sc_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_SURFACE_COEFF  << endl;
    }
    dbSurfacePressureLUT* sp_lut = new dbSurfacePressureLUT;
    istatus = read_surface_pressure_file( sp_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading surface pressure LUT file " << endl;
        status = istatus;
    }
    istatus = write_surface_pressure_lut( nc_output, sp_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing surface pressure data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete sp_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_SURFACE_PRESSURE  << endl;
    }
    dbViirsSwirVsVisLUT* vsw_lut = new dbViirsSwirVsVisLUT;
    istatus = read_swir_file( vsw_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading VIIRS SWIR LUT file " << endl;
        status = istatus;
    }
    istatus = write_swir_lut( nc_output, vsw_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing VIIRS SWIR data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete vsw_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_SWIR  << endl;
    }

    dbVeg_21sfcLUT* v_lut = new dbVeg_21sfcLUT;
    istatus = read_veg_21sfc_files( v_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure reading veg_21sfc LUT files " << endl;
        status = istatus;
    }
    istatus = write_veg_21sfc_lut( nc_output, v_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DbLutNetcdf:: Failure writing veg_21sfc data to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete v_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DbLutNetcdf:: Created " + LUT_NVALX21  << endl;
    }
    delete nc_output;

    if (status == DTDB_SUCCESS) {
        return status;
    } else {
        return num_good++;
    }
}

/**************************************************************************
 * NAME: read_tables_file()
 *
 * DESCRIPTION: Read MODIS/SEAWIFS tables file.
 *
 *************************************************************************/

int DbLutNetcdf::read_tables_file( dbTablesLUT* mt_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_LER_TABLE);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_tables_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    bool isFileBigEndian = true;
    std::ifstream fin(filepath.c_str(), std::ios::in | std::ios::binary);
    if(fin.is_open()) {
        const int RECORD_DELIMITER_LENGTH = 4;
        int r1_length = (4*20800 + 260)*sizeof(float);
        int r2_length = (4*160 + 2)*sizeof(float);
        fin.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
        fin.read( (char*) mt_lut->LOGI0, r1_length);
        fin.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
        fin.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
        fin.read( (char*) mt_lut->LOGI0R, r2_length);
        bool good = fin.good();
        fin.close();
        if(!good) {
          cerr <<
          "DbLutNetcdf::read_tables_file() Error reading binary file "
                  << INPUT_LER_TABLE << endl;
          return DTDB_FAIL;
       }
    }
    else {
       cerr << "DbLutNetcdf::read_tables_file() Error opening binary file "
               << INPUT_LER_TABLE << endl;
       return DTDB_FAIL;
    }
    if ( isPlatformLittleEndian() && isFileBigEndian) {
        for ( int i=0; i<MTABLE_20800; i++) {
            byteSwap( mt_lut->LOGI0[i] );
            byteSwap( mt_lut->Z1I0[i] );
            byteSwap( mt_lut->Z2I0[i] );
            byteSwap( mt_lut->TI0[i] );
        }
        for ( int i=0; i<MTABLE_260; i++) {
            byteSwap( mt_lut->SB[i] );
        }
        for ( int i=0; i<MTABLE_160; i++) {
            byteSwap( mt_lut->LOGI0R[i] );
            byteSwap( mt_lut->Z1I0R[i] );
            byteSwap( mt_lut->Z2I0R[i] );
            byteSwap( mt_lut->TI0R[i] );
        }
        for ( int i=0; i<MTABLE_2; i++) {
            byteSwap( mt_lut->SBR[i] );
        }
    }

    return status;
}

/**************************************************************************
 * NAME: write_tables_lut()
 *
 * DESCRIPTION: Write modis tables LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_tables_lut( NcFile* nc_output,
        dbTablesLUT* mt_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_LER_TABLES );

    dim_20800_ = lut_grp.addDim( "Dim_20800", MTABLE_20800 );
    dim_260_ = lut_grp.addDim( "Dim_260", MTABLE_260 );
    dim_160_ = lut_grp.addDim( "Dim_160", MTABLE_160 );
    dim_2_ = lut_grp.addDim( "Dim_2", MTABLE_2 );

    NcVar var = lut_grp.addVar( "LOGI0", ncFloat, dim_20800_ );
    var.putVar( mt_lut->LOGI0 );

    var = lut_grp.addVar( "Z1I0", ncFloat, dim_20800_ );
    var.putVar( mt_lut->Z1I0 );

    var = lut_grp.addVar( "Z2I0", ncFloat, dim_20800_ );
    var.putVar( mt_lut->Z2I0 );

    var = lut_grp.addVar( "TI0", ncFloat, dim_20800_ );
    var.putVar( mt_lut->TI0 );

    var = lut_grp.addVar( "SB", ncFloat, dim_260_ );
    var.putVar( mt_lut->SB );

    var = lut_grp.addVar( "LOGI0R", ncFloat, dim_160_ );
    var.putVar( mt_lut->LOGI0R );

    var = lut_grp.addVar( "Z1I0R", ncFloat, dim_160_ );
    var.putVar( mt_lut->Z1I0R );

    var = lut_grp.addVar( "Z2I0R", ncFloat, dim_160_ );
    var.putVar( mt_lut->Z2I0R );

    var = lut_grp.addVar( "TI0R", ncFloat, dim_160_ );
    var.putVar( mt_lut->TI0R );

    var = lut_grp.addVar( "SBR", ncFloat, dim_2_ );
    var.putVar( mt_lut->SBR );

    return DTDB_SUCCESS;
}



/**************************************************************************
 * NAME: read_tables_lut()
 *
 * DESCRIPTION: Read MODIS Tables NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_tables_lut( dbTablesLUT* mt_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_LER_TABLES );
    NcVar var = lut_grp.getVar( "LOGI0" );
    var.getVar( mt_lut->LOGI0 );
    var = lut_grp.getVar( "Z1I0" );
    var.getVar( mt_lut->Z1I0 );
    var = lut_grp.getVar( "Z2I0" );
    var.getVar( mt_lut->Z2I0 );
    var = lut_grp.getVar( "TI0" );
    var.getVar( mt_lut->TI0 );
    var = lut_grp.getVar( "SB" );
    var.getVar( mt_lut->SB );
    var = lut_grp.getVar( "LOGI0R" );
    var.getVar( mt_lut->LOGI0R );
    var = lut_grp.getVar( "Z1I0R" );
    var.getVar( mt_lut->Z1I0R );
    var = lut_grp.getVar( "Z2I0R" );
    var.getVar( mt_lut->Z2I0R );
    var = lut_grp.getVar( "TI0R" );
    var.getVar( mt_lut->TI0R );
    var = lut_grp.getVar( "SBR" );
    var.getVar( mt_lut->SBR );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_surface_pressure_file()
 *
 * DESCRIPTION: Read surface pressure hdf4 file.
 *
 *************************************************************************/

int DbLutNetcdf::read_surface_pressure_file( dbSurfacePressureLUT* sp_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_SURFACE_PRESSURE);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_surface_pressure_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    int fileID;
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening surface pressure file: " + filepath << endl;
        return DTDB_FAIL;
    }
    int sds_index, sds_id, numtype, rank, nattrs;
    int start[2], edges[2], dims[2];
    start[0] = 0;
    start[1] = 0;
    edges[0] = SP_720;
    edges[1] = SP_360;
    string sds_name = "ps";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_surface_pressure_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &sp_lut->PS[0][0]);

    edges[0] = SP_SETS*SP_720;
    edges[1] = SP_SETS*SP_360;
    sds_name = "surface_pressure";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_surface_pressure_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &sp_lut->SURFACE_PRESSURE[0][0]);

    sds_name = "surface_elevation";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_surface_pressure_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &sp_lut->SURFACE_ELEVATION[0][0]);

    return status;
}

/**************************************************************************
 * NAME: write_surface_pressure_lut()
 *
 * DESCRIPTION: Write surface pressure LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_surface_pressure_lut( NcFile* nc_output,
        dbSurfacePressureLUT* sp_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_SURFACE_PRESSURE );

    dim_720_ = lut_grp.addDim( "Dim_720", SP_720 );
    dim_360_ = lut_grp.addDim( "Dim_360", SP_360 );
    dim_4320_ = lut_grp.addDim( "Dim_4320", SP_SETS*SP_720 );
    dim_2160_ = lut_grp.addDim( "Dim_2160", SP_SETS*SP_360 );

    vector<NcDim> small_dims;
    small_dims.push_back(dim_720_);
    small_dims.push_back(dim_360_);

    vector<NcDim> big_dims;
    big_dims.push_back(dim_4320_);
    big_dims.push_back(dim_2160_);

    NcVar var = lut_grp.addVar( "PS", ncFloat, small_dims );
    var.putVar( sp_lut->PS );

    var = lut_grp.addVar( "SURFACE_PRESSURE", ncFloat, big_dims );
    var.putVar( sp_lut->SURFACE_PRESSURE );

    var = lut_grp.addVar( "SURFACE_ELEVATION", ncFloat, big_dims );
    var.putVar( sp_lut->SURFACE_ELEVATION );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_surface_pressure_lut()
 *
 * DESCRIPTION: Read surface pressure NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_surface_pressure_lut( dbSurfacePressureLUT* sp_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_SURFACE_PRESSURE );

    NcVar var = lut_grp.getVar( "PS" );
    var.getVar( sp_lut->PS );

    var = lut_grp.getVar( "SURFACE_PRESSURE" );
    var.getVar( sp_lut->SURFACE_PRESSURE );

    var = lut_grp.getVar( "SURFACE_ELEVATION" );
    var.getVar( sp_lut->SURFACE_ELEVATION );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_nvalx_files()
 *
 * DESCRIPTION: Read nvalx files.
 *
 *************************************************************************/

int DbLutNetcdf::read_nvalx_files( dbNvalxLUT* nv_lut )
{
    int status = DTDB_SUCCESS;

    const int RECORD_DELIMITER_LENGTH = 4;

    string filepath = get_option(INPUT_NVALX_412);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_nvalx_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    bool isFileBigEndian = false;
    std::ifstream fin1(filepath.c_str(), std::ios::in | std::ios::binary);
    if(fin1.is_open()) {
        int length = SR412*SSA412*NTAU*NRAA*NVZA*NSZA*sizeof(float);
        fin1.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
        fin1.read( (char*) nv_lut->NVALX_412, length);
        bool good = fin1.good();
        fin1.close();
        if(!good) {
          cerr <<
          "DbLutNetcdf::read_nvalx_file() Error reading binary file "
                  << INPUT_LER_TABLE << endl;
          return DTDB_FAIL;
       }
    }
    else {
       cerr << "DbLutNetcdf::read_nvalx_file() Error opening binary file "
               << INPUT_NVALX_412 << endl;
       return DTDB_FAIL;
    }

    filepath = get_option(INPUT_NVALX_470);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_nvalx_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    isFileBigEndian = false;
    std::ifstream fin2(filepath.c_str(), std::ios::in | std::ios::binary);
    if(fin2.is_open()) {
        int length = SR470*SSA470*NTAU*NRAA*NVZA*NSZA*sizeof(float);
        fin2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
        fin2.read( (char*) nv_lut->NVALX_470, length);
        bool good = fin2.good();
        fin2.close();
        if(!good) {
          cerr <<
          "DbLutNetcdf::read_nvalx_file() Error reading binary file "
                  << INPUT_NVALX_470 << endl;
          return DTDB_FAIL;
       }
    }
    else {
       cerr << "DbLutNetcdf::read_nvalx_file() Error opening binary file "
               << INPUT_NVALX_470 << endl;
       return DTDB_FAIL;
    }

    filepath = get_option(INPUT_NVALX_650);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_nvalx_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    isFileBigEndian = false;
    std::ifstream fin3(filepath.c_str(), std::ios::in | std::ios::binary);
    if(fin3.is_open()) {
        int length = SR650*NTAU*NRAA*NVZA*NSZA*sizeof(float);
        fin3.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
        fin3.read( (char*) nv_lut->NVALX_650, length);
        bool good = fin3.good();
        fin3.close();
        if(!good) {
          cerr <<
          "DbLutNetcdf::read_nvalx_file() Error reading binary file "
                  << INPUT_NVALX_650 << endl;
          return DTDB_FAIL;
       }
    }
    else {
       cerr << "DbLutNetcdf::read_nvalx_file() Error opening binary file "
               << INPUT_NVALX_650 << endl;
       return DTDB_FAIL;
    }

    if ( isPlatformLittleEndian() && isFileBigEndian) {
    }

    return status;
}

/**************************************************************************
 * NAME: write_nvalx_lut()
 *
 * DESCRIPTION: Write nvalx LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_nvalx_lut( NcFile* nc_output,
        dbNvalxLUT* nv_lut  )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_NVALX );

    dim_nsza_ = lut_grp.addDim( "Dim_NSZA", NSZA );
    dim_nvza_ = lut_grp.addDim( "Dim_NVZA", NVZA );
    dim_nraa_ = lut_grp.addDim( "Dim_NRAA", NRAA );
    dim_ntau_ = lut_grp.addDim( "Dim_NTAU", NTAU );
    dim_ssa412_ = lut_grp.addDim( "Dim_SSA412", SSA412 );
    dim_ssa470_ = lut_grp.addDim( "Dim_SSA470", SSA470 );
    dim_sr412_ = lut_grp.addDim( "Dim_SR412", SR412 );
    dim_sr470_ = lut_grp.addDim( "Dim_SR470", SR470 );
    dim_sr650_ = lut_grp.addDim( "Dim_SR650", SR650 );

    vector<NcDim> nvalx412_dims;
    nvalx412_dims.push_back(dim_sr412_);
    nvalx412_dims.push_back(dim_ssa412_);
    nvalx412_dims.push_back(dim_ntau_);
    nvalx412_dims.push_back(dim_nraa_);
    nvalx412_dims.push_back(dim_nvza_);
    nvalx412_dims.push_back(dim_nsza_);

    vector<NcDim> nvalx470_dims;
    nvalx470_dims.push_back(dim_sr470_);
    nvalx470_dims.push_back(dim_ssa470_);
    nvalx470_dims.push_back(dim_ntau_);
    nvalx470_dims.push_back(dim_nraa_);
    nvalx470_dims.push_back(dim_nvza_);
    nvalx470_dims.push_back(dim_nsza_);

    vector<NcDim> nvalx650_dims;
    nvalx650_dims.push_back(dim_sr650_);
    nvalx650_dims.push_back(dim_ntau_);
    nvalx650_dims.push_back(dim_nraa_);
    nvalx650_dims.push_back(dim_nvza_);
    nvalx650_dims.push_back(dim_nsza_);

    NcVar var = lut_grp.addVar( "NVALX_412", ncFloat, nvalx412_dims );
    var.putVar( nv_lut->NVALX_412 );

    var = lut_grp.addVar( "NVALX_470", ncFloat, nvalx470_dims );
    var.putVar( nv_lut->NVALX_470 );

    var = lut_grp.addVar( "NVALX_650", ncFloat, nvalx650_dims );
    var.putVar( nv_lut->NVALX_650 );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_nvalx_lut()
 *
 * DESCRIPTION: Read nvalx NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_nvalx_lut(  dbNvalxLUT* nv_lut  )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_NVALX );

    NcVar var = lut_grp.getVar( "NVALX_412" );
    var.getVar( nv_lut->NVALX_412 );

    var = lut_grp.getVar( "NVALX_470" );
    var.getVar( nv_lut->NVALX_470 );

    var = lut_grp.getVar( "NVALX_650" );
    var.getVar( nv_lut->NVALX_650 );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_veg_21sfc_files()
 *
 * DESCRIPTION: Read veg_21sfc files.
 *
 *************************************************************************/

int DbLutNetcdf::read_veg_21sfc_files( dbVeg_21sfcLUT* nv_lut )
{
    int status = DTDB_SUCCESS;

// NVALX21
    string filepath = get_option(INPUT_VEG_21SFC);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_nvalx21_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    int fileID;
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening nvalx21 file: " + filepath << endl;
        return DTDB_FAIL;
    }
    int sds_index, sds_id, numtype, rank, nattrs;
    int start[4], edges[4], dims[4];
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    edges[0] = NSEASONS;
    edges[1] = NRAA;
    edges[2] = NVZA;
    edges[3] = NSZAV;

    string sds_name = "NVALX21_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->NVALX21_SFC[0]);

    sds_name = "R0X21_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->R0X21_SFC[0]);

    sds_name = "SX21_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->SX21_SFC[0]);

    sds_name = "TX21_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->TX21_SFC[0]);

    sds_name = "NVALX672_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->NVALX672_SFC[0]);

    sds_name = "R0X672_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->R0X672_SFC[0]);

    sds_name = "SX672_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->SX672_SFC[0]);

    sds_name = "TX672_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->TX672_SFC[0]);

    sds_name = "NVALX865_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->NVALX865_SFC[0]);

    sds_name = "R0X865_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->R0X865_SFC[0]);

    sds_name = "SX865_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->SX865_SFC[0]);

    sds_name = "TX865_SFC";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_nvalx21_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &nv_lut->TX865_SFC[0]);

    return status;
}

/**************************************************************************
 * NAME: write_veg_21sfc_lut()
 *
 * DESCRIPTION: Write veg_21sfc LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_veg_21sfc_lut( NcFile* nc_output,
        dbVeg_21sfcLUT* nv_lut  )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_NVALX21 );

    dim_nszav_ = lut_grp.addDim( "Dim_NSZA", NSZAV );
    dim_nvza_ = lut_grp.addDim( "Dim_NVZA", NVZA );
    dim_nraa_ = lut_grp.addDim( "Dim_NRAA", NRAA );
    dim_4_ = lut_grp.addDim( "Dim_SEASONS", 4 );

    vector<NcDim> nvalx21_dims;
    nvalx21_dims.push_back(dim_4_);
    nvalx21_dims.push_back(dim_nraa_);
    nvalx21_dims.push_back(dim_nvza_);
    nvalx21_dims.push_back(dim_nszav_);

    NcVar var = lut_grp.addVar( "NVALX21_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->NVALX21_SFC );

    var = lut_grp.addVar( "R0X21_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->R0X21_SFC );

    var = lut_grp.addVar( "SX21_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->SX21_SFC );

    var = lut_grp.addVar( "TX21_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->TX21_SFC );

    var = lut_grp.addVar( "NVALX672_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->NVALX672_SFC );

    var = lut_grp.addVar( "R0X672_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->R0X672_SFC );

    var = lut_grp.addVar( "SX672_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->SX672_SFC );

    var = lut_grp.addVar( "TX672_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->TX672_SFC );

    var = lut_grp.addVar( "NVALX865_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->NVALX865_SFC );

    var = lut_grp.addVar( "R0X865_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->R0X865_SFC );

    var = lut_grp.addVar( "SX865_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->SX865_SFC );

    var = lut_grp.addVar( "TX865_SFC", ncFloat, nvalx21_dims );
    var.putVar( nv_lut->TX865_SFC );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_veg_21sfc_lut()
 *
 * DESCRIPTION: Read veg_21sfc NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_veg_21sfc_lut( dbVeg_21sfcLUT* nv_lut  )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_NVALX21 );

    NcVar var = lut_grp.getVar( "NVALX21_SFC" );
    var.getVar( nv_lut->NVALX21_SFC );

    var = lut_grp.getVar( "R0X21_SFC" );
    var.getVar( nv_lut->R0X21_SFC );

    var = lut_grp.getVar( "SX21_SFC" );
    var.getVar( nv_lut->SX21_SFC );

    var = lut_grp.getVar( "TX21_SFC" );
    var.getVar( nv_lut->TX21_SFC );

    var = lut_grp.getVar( "NVALX672_SFC" );
    var.getVar( nv_lut->NVALX672_SFC );

    var = lut_grp.getVar( "R0X672_SFC" );
    var.getVar( nv_lut->R0X672_SFC );

    var = lut_grp.getVar( "SX672_SFC" );
    var.getVar( nv_lut->SX672_SFC );

    var = lut_grp.getVar( "TX672_SFC" );
    var.getVar( nv_lut->TX672_SFC );

    var = lut_grp.getVar( "NVALX865_SFC" );
    var.getVar( nv_lut->NVALX865_SFC );

    var = lut_grp.getVar( "R0X865_SFC" );
    var.getVar( nv_lut->R0X865_SFC );

    var = lut_grp.getVar( "SX865_SFC" );
    var.getVar( nv_lut->SX865_SFC );

    var = lut_grp.getVar( "TX865_SFC" );
    var.getVar( nv_lut->TX865_SFC );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_modis_surf_refl_file()
 *
 * DESCRIPTION: Read seasonal surface reflectance hdf4 file.
 *
 *************************************************************************/

int DbLutNetcdf::read_modis_surf_refl_files( dbModisSurfReflLUT* sr_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_MODIS_SURF_REFL);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_modis_surf_refl_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    SEASON iseason = SEASON::NEVER;
    size_t rpos = 0;
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        size_t pos = filepath.find(str_season[iS]);
        if (pos != string::npos) {
            iseason = (SEASON) iS;
            rpos = pos;
        }
    }
    if (iseason == SEASON::NEVER) {
        cerr << "DbLutNetcdf:: Could not identify season in file name: " + filepath << endl;
        return DTDB_FAIL;
    }
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        string rfilepath = filepath;
        rfilepath.replace(rpos, str_season[(int)iseason].size(), str_season[iS]);
        int fileID;
        try {
            fileID = SDstart(rfilepath.c_str(), DFACC_READ );
        }
        catch( std::exception& e) {
            e.what();
            cerr << "DbLutNetcdf:: Failure opening modis surface refl file: " + rfilepath << endl;
            return DTDB_FAIL;
        }
        int sds_index, sds_id, numtype, rank, nattrs;
        int start[2], edges[2], dims[2];
        start[0] = 0;
        start[1] = 0;
        edges[0] = NLATS*10;
        edges[1] = NLONS*10;

        string sds_name = "412_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_modis_surf_refl_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->SR412_ALL[iS][0][0]);
// 470
        sds_name = "470_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_modis_surf_refl_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->SR470_ALL[iS][0][0]);
// 650
        sds_name = "650_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_modis_surf_refl_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->SR650_ALL[iS][0][0]);

// 865
        sds_name = "865_all_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_modis_surf_refl_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->SR865_ALL[iS][0][0]);

        SDend(fileID);
    }

    return status;
}

/**************************************************************************
 * NAME: write_modis_surf_refl_lut()
 *
 * DESCRIPTION: Write surf_refl LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_modis_surf_refl_lut( NcFile* nc_output,
        dbModisSurfReflLUT* sr_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_MODIS_SURFACE_REFL );

    dim_seasons_ = lut_grp.addDim( "Dim_seasons", NUM_SEASONS );
    dim_1800_ = lut_grp.addDim( "Dim_1800", NLATS*10 );
    dim_3600_ = lut_grp.addDim( "Dim_3600", NLONS*10 );

    vector<NcDim> latlon_dims;
    latlon_dims.push_back(dim_1800_);
    latlon_dims.push_back(dim_3600_);

    vector<NcDim> slatlon_dims;
    slatlon_dims.push_back(dim_seasons_);
    slatlon_dims.push_back(dim_1800_);
    slatlon_dims.push_back(dim_3600_);

    NcVar var = lut_grp.addVar( "SR412_ALL", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR412_ALL );

    var = lut_grp.addVar( "SR470_ALL", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR470_ALL );

    var = lut_grp.addVar( "SR650_ALL", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR650_ALL );

    var = lut_grp.addVar( "SR412_FWD", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR412_FWD );

    var = lut_grp.addVar( "SR470_FWD", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR470_FWD );

    var = lut_grp.addVar( "SR650_FWD", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR650_FWD );

    var = lut_grp.addVar( "SR865_ALL", ncFloat, latlon_dims );
    var.putVar( sr_lut->SR865_ALL );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_modis_surf_refl_lut()
 *
 * DESCRIPTION: Read surf_refl NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_modis_surf_refl_lut( dbModisSurfReflLimited* sr_lut,
        int* start, int* edge, int &season, int &dateline )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_MODIS_SURFACE_REFL );

    sr_lut->SR412_ALL_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR470_ALL_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR650_ALL_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR412_FWD_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR470_FWD_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR650_FWD_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR865_ALL_L.resize(boost::extents[edge[1]][edge[0]]);

    vector<size_t> startp;
    vector<size_t> countp;
    startp.push_back(season);
    startp.push_back(start[1]);
    startp.push_back(start[0]);
    countp.push_back(1);
    countp.push_back(edge[1]);
    countp.push_back(edge[0]);

    if (dateline ==0) {
        NcVar var = lut_grp.getVar( "SR412_ALL" );
        var.getVar( startp, countp, &sr_lut->SR412_ALL_L[0][0] );
        var = lut_grp.getVar( "SR470_ALL" );
        var.getVar( startp, countp, &sr_lut->SR470_ALL_L[0][0] );
        var = lut_grp.getVar( "SR650_ALL" );
        var.getVar( startp, countp, &sr_lut->SR650_ALL_L[0][0] );
        var = lut_grp.getVar( "SR412_FWD" );
        var.getVar( startp, countp, &sr_lut->SR412_FWD_L[0][0] );
        var = lut_grp.getVar( "SR470_FWD" );
        var.getVar( startp, countp, &sr_lut->SR470_FWD_L[0][0] );
        var = lut_grp.getVar( "SR650_FWD" );
        var.getVar( startp, countp, &sr_lut->SR650_FWD_L[0][0] );
        var = lut_grp.getVar( "SR865_ALL" );
        var.getVar( startp, countp, &sr_lut->SR865_ALL_L[0][0] );
    } else {
        countp[2] = dateline;
        NcVar var = lut_grp.getVar( "SR412_ALL" );
        var.getVar( startp, countp, &sr_lut->SR412_ALL_L[0][0] );
        var = lut_grp.getVar( "SR470_ALL" );
        var.getVar( startp, countp, &sr_lut->SR470_ALL_L[0][0] );
        var = lut_grp.getVar( "SR650_ALL" );
        var.getVar( startp, countp, &sr_lut->SR650_ALL_L[0][0] );
        var = lut_grp.getVar( "SR412_FWD" );
        var.getVar( startp, countp, &sr_lut->SR412_FWD_L[0][0] );
        var = lut_grp.getVar( "SR470_FWD" );
        var.getVar( startp, countp, &sr_lut->SR470_FWD_L[0][0] );
        var = lut_grp.getVar( "SR650_FWD" );
        var.getVar( startp, countp, &sr_lut->SR650_FWD_L[0][0] );
        var = lut_grp.getVar( "SR865_ALL" );
        var.getVar( startp, countp, &sr_lut->SR865_ALL_L[0][0] );
        startp[2] = 0;
        countp[2] = edge[0] - dateline;
        var = lut_grp.getVar( "SR412_ALL" );
        var.getVar( startp, countp, &sr_lut->SR412_ALL_L[0][dateline] );
        var = lut_grp.getVar( "SR470_ALL" );
        var.getVar( startp, countp, &sr_lut->SR470_ALL_L[0][dateline] );
        var = lut_grp.getVar( "SR650_ALL" );
        var.getVar( startp, countp, &sr_lut->SR650_ALL_L[0][dateline] );
        var = lut_grp.getVar( "SR412_FWD" );
        var.getVar( startp, countp, &sr_lut->SR412_FWD_L[0][dateline] );
        var = lut_grp.getVar( "SR470_FWD" );
        var.getVar( startp, countp, &sr_lut->SR470_FWD_L[0][dateline] );
        var = lut_grp.getVar( "SR650_FWD" );
        var.getVar( startp, countp, &sr_lut->SR650_FWD_L[0][dateline] );
        var = lut_grp.getVar( "SR865_ALL" );
        var.getVar( startp, countp, &sr_lut->SR865_ALL_L[0][dateline] );
    }
    delete nc_input;

    return status;
}

/**************************************************************************
 * NAME: read_viirs_surf_refl_file()
 *
 * DESCRIPTION: Read seasonal surface reflectance hdf4 file.
 *
 *************************************************************************/

int DbLutNetcdf::read_viirs_surf_refl_files( dbViirsSurfReflLUT* sr_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_VIIRS_SURF_REFL);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_surf_refl_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    SEASON iseason = SEASON::NEVER;
    size_t rpos = 0;
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        size_t pos = filepath.find(str_season[iS]);
        if (pos != string::npos) {
            iseason = (SEASON) iS;
            rpos = pos;
        }
    }
    if (iseason == SEASON::NEVER) {
        cerr << "DbLutNetcdf:: Could not identify season in file name: " + filepath << endl;
        return DTDB_FAIL;
    }
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        string rfilepath = filepath;
        rfilepath.replace(rpos, str_season[(int) iseason].size(), str_season[iS]);
        int fileID;
        try {
            fileID = SDstart(rfilepath.c_str(), DFACC_READ );
        }
        catch( std::exception& e) {
            e.what();
            cerr << "DbLutNetcdf:: Failure opening surface refl file: " + rfilepath << endl;
            return DTDB_FAIL;
        }
        int sds_index, sds_id, numtype, rank, nattrs;
        int start[2], edges[2], dims[2];
        start[0] = 0;
        start[1] = 0;
        edges[0] = NLATS*10;
        edges[1] = NLONS*10;

        string sds_name = "412_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_refl_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->SR412_ALL[iS][0][0]);
// 470
        sds_name = "488_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_refl_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->SR488_ALL[iS][0][0]);
// 650
        sds_name = "670_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_refl_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->SR670_ALL[iS][0][0]);

        string filepath = get_option(INPUT_BRDF);
        if(filepath.empty()) {
            cerr << "DbLutNetcdf::read_landcover_file() Invalid path." << endl;
            return DTDB_FAIL;
        }
        try {
            fileID = SDstart(filepath.c_str(), DFACC_READ );
        }
        catch( std::exception& e) {
            e.what();
            cerr << "DbLutNetcdf:: Failure opening brdf file: " + filepath << endl;
            return DTDB_FAIL;
        }
        start[0] = 0;
        start[1] = 0;
        edges[0] = NLATS*10;
        edges[1] = NLONS*10;
        sds_name = "brdf_base_650";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_landcover_files() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sr_lut->BRDF_650[0][0]);

        SDend(fileID);
    }

    return status;
}

/**************************************************************************
 * NAME: write_viirs_surf_refl_lut()
 *
 * DESCRIPTION: Write surf_refl LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_viirs_surf_refl_lut( NcFile* nc_output,
        dbViirsSurfReflLUT* sr_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_VIIRS_SURFACE_REFL );

    dim_seasons_ = lut_grp.addDim( "Dim_seasons", NUM_SEASONS );
    dim_1800_ = lut_grp.addDim( "Dim_1800", NLATS*10 );
    dim_3600_ = lut_grp.addDim( "Dim_3600", NLONS*10 );

    vector<NcDim> latlon_dims;
    latlon_dims.push_back(dim_1800_);
    latlon_dims.push_back(dim_3600_);

    vector<NcDim> slatlon_dims;
    slatlon_dims.push_back(dim_seasons_);
    slatlon_dims.push_back(dim_1800_);
    slatlon_dims.push_back(dim_3600_);

    NcVar var = lut_grp.addVar( "SR412_ALL", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR412_ALL );

    var = lut_grp.addVar( "SR488_ALL", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR488_ALL );

    var = lut_grp.addVar( "SR670_ALL", ncFloat, slatlon_dims );
    var.putVar( sr_lut->SR670_ALL );

    var = lut_grp.addVar( "BRDF_650", ncFloat, latlon_dims );
    var.putVar( sr_lut->BRDF_650 );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_viirs_surf_refl_lut()
 *
 * DESCRIPTION: Read surf_refl NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_viirs_surf_refl_lut( dbViirsSurfReflLimited* sr_lut,
        int* start, int* edge, int &season, int &dateline )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_VIIRS_SURFACE_REFL );

    sr_lut->SR412_ALL_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR488_ALL_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->SR670_ALL_L.resize(boost::extents[edge[1]][edge[0]]);
    sr_lut->BRDF_650_L.resize(boost::extents[edge[1]][edge[0]]);

    vector<size_t> startp;
    vector<size_t> countp;
    startp.push_back(season);
    startp.push_back(start[1]);
    startp.push_back(start[0]);
    countp.push_back(1);
    countp.push_back(edge[1]);
    countp.push_back(edge[0]);

    if (dateline ==0) {
        NcVar var = lut_grp.getVar( "SR412_ALL" );
        var.getVar( startp, countp, &sr_lut->SR412_ALL_L[0][0] );
        var = lut_grp.getVar( "SR488_ALL" );
        var.getVar( startp, countp, &sr_lut->SR488_ALL_L[0][0] );
        var = lut_grp.getVar( "SR670_ALL" );
        var.getVar( startp, countp, &sr_lut->SR670_ALL_L[0][0] );
        var = lut_grp.getVar( "BRDF_650" );
        var.getVar( startp, countp, &sr_lut->BRDF_650_L[0][0] );
    } else {
        countp[2] = dateline;
        NcVar var = lut_grp.getVar( "SR412_ALL" );
        var.getVar( startp, countp, &sr_lut->SR412_ALL_L[0][0] );
        var = lut_grp.getVar( "SR488_ALL" );
        var.getVar( startp, countp, &sr_lut->SR488_ALL_L[0][0] );
        var = lut_grp.getVar( "SR670_ALL" );
        var.getVar( startp, countp, &sr_lut->SR670_ALL_L[0][0] );
        var = lut_grp.getVar( "BRDF_650" );
        var.getVar( startp, countp, &sr_lut->BRDF_650_L[0][0] );
        startp[2] = 0;
        countp[2] = edge[0] - dateline;
        var = lut_grp.getVar( "SR412_ALL" );
        var.getVar( startp, countp, &sr_lut->SR412_ALL_L[0][dateline] );
        var = lut_grp.getVar( "SR488_ALL" );
        var.getVar( startp, countp, &sr_lut->SR488_ALL_L[0][dateline] );
        var = lut_grp.getVar( "SR670_ALL" );
        var.getVar( startp, countp, &sr_lut->SR670_ALL_L[0][dateline] );
        var = lut_grp.getVar( "BRDF_650" );
        var.getVar( startp, countp, &sr_lut->BRDF_650_L[0][dateline] );
    }
    delete nc_input;

    return status;
}

/**************************************************************************
 * NAME: read_surf_coeff_files()
 *
 * DESCRIPTION: Read seasonal surface coefficients hdf4 files.
 *
 *************************************************************************/

int DbLutNetcdf::read_surf_coeff_files( dbSurfCoeffLUT* sc_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_SURF_COEFF);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_surf_coeff_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    SEASON iseason = SEASON::NEVER;
    size_t rpos = 0;
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        size_t pos = filepath.find(str_season[iS]);
        if (pos != string::npos) {
            iseason = (SEASON) iS;
            rpos = pos;
        }
    }
    if (iseason == SEASON::NEVER) {
        cerr << "DbLutNetcdf:: Could not identify season in file name: " + filepath << endl;
        return DTDB_FAIL;
    }
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        string rfilepath = filepath;
        rfilepath.replace(rpos, str_season[(int) iseason].size(), str_season[iS]);
        int fileID;
        try {
            fileID = SDstart(rfilepath.c_str(), DFACC_READ );
        }
        catch( std::exception& e) {
            e.what();
            cerr << "DbLutNetcdf:: Failure opening surface coeff file: " + rfilepath << endl;
            return DTDB_FAIL;
        }
        int sds_index, sds_id, numtype, rank, nattrs;
        int start[4], edges[4], dims[4];
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        start[3] = 0;
        edges[0] = NNDVI;
        edges[1] = NTERMS;
        edges[2] = NLATS*10;
        edges[3] = NLONS*10;

        string sds_name = "412_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_coeff_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sc_lut->SC412_ALL[iS][0][0][0][0]);

        sds_name = "412_fwd";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_coeff_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sc_lut->SC412_FWD[iS][0][0][0][0]);

// 470
        sds_name = "470_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_coeff_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sc_lut->SC470_ALL[iS][0][0][0][0]);

        sds_name = "470_fwd";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_coeff_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sc_lut->SC470_FWD[iS][0][0][0][0]);

// 650
        sds_name = "650_all";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_coeff_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sc_lut->SC650_ALL[iS][0][0][0][0]);

        sds_name = "650_fwd";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_surf_coeff_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &sc_lut->SC650_FWD[iS][0][0][0][0]);

        SDend(fileID);
    }

    return status;
}

/**************************************************************************
 * NAME: write_surf_coeff_lut()
 *
 * DESCRIPTION: Write surf_coeff LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_surf_coeff_lut( NcFile* nc_output,
        dbSurfCoeffLUT* sc_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_SURFACE_COEFF );

    dim_seasons_ = lut_grp.addDim( "Dim_seasons", NUM_SEASONS );
    dim_ndvi_ = lut_grp.addDim( "Dim_ndvi", NNDVI );
    dim_terms_ = lut_grp.addDim( "Dim_terms", NTERMS );
    dim_1800_ = lut_grp.addDim( "Dim_1800", NLATS*10 );
    dim_3600_ = lut_grp.addDim( "Dim_3600", NLONS*10 );

    vector<NcDim> latlon_dims;
    latlon_dims.push_back(dim_seasons_);
    latlon_dims.push_back(dim_ndvi_);
    latlon_dims.push_back(dim_terms_);
    latlon_dims.push_back(dim_1800_);
    latlon_dims.push_back(dim_3600_);

    NcVar var = lut_grp.addVar( "SC412_ALL", ncFloat, latlon_dims );
    var.putVar( sc_lut->SC412_ALL );

    var = lut_grp.addVar( "SC412_FWD", ncFloat, latlon_dims );
    var.putVar( sc_lut->SC412_FWD );

    var = lut_grp.addVar( "SC470_ALL", ncFloat, latlon_dims );
    var.putVar( sc_lut->SC470_ALL );

    var = lut_grp.addVar( "SC470_FWD", ncFloat, latlon_dims );
    var.putVar( sc_lut->SC470_FWD );

    var = lut_grp.addVar( "SC650_ALL", ncFloat, latlon_dims );
    var.putVar( sc_lut->SC650_ALL );

    var = lut_grp.addVar( "SC650_FWD", ncFloat, latlon_dims );
    var.putVar( sc_lut->SC650_FWD );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_surf_coeff_lut()
 *
 * DESCRIPTION: Read surface coefficients NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_surf_coeff_lut( dbSurfCoeffLimited* sc_lut, int* start,
        int* edge, int &season, int &dateline )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_SURFACE_COEFF );

    sc_lut->SC412_ALL_L.resize(boost::extents[3][4][edge[1]][edge[0]]);
    sc_lut->SC412_FWD_L.resize(boost::extents[3][4][edge[1]][edge[0]]);
    sc_lut->SC470_ALL_L.resize(boost::extents[3][4][edge[1]][edge[0]]);
    sc_lut->SC470_FWD_L.resize(boost::extents[3][4][edge[1]][edge[0]]);
    sc_lut->SC650_ALL_L.resize(boost::extents[3][4][edge[1]][edge[0]]);
    sc_lut->SC650_FWD_L.resize(boost::extents[3][4][edge[1]][edge[0]]);

    vector<size_t> startp;
    vector<size_t> countp;
    startp.push_back(season);
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(start[1]);
    startp.push_back(start[0]);
    countp.push_back(1);
    countp.push_back(3);
    countp.push_back(4);
    countp.push_back(edge[1]);
    countp.push_back(edge[0]);

    if (dateline ==0) {
        NcVar var = lut_grp.getVar( "SC412_ALL" );
        var.getVar( startp, countp, &sc_lut->SC412_ALL_L[0][0][0][0] );
        var = lut_grp.getVar( "SC412_FWD" );
        var.getVar( startp, countp, &sc_lut->SC412_FWD_L[0][0][0][0] );
        var = lut_grp.getVar( "SC470_ALL" );
        var.getVar( startp, countp, &sc_lut->SC470_ALL_L[0][0][0][0] );
        var = lut_grp.getVar( "SC470_FWD" );
        var.getVar( startp, countp, &sc_lut->SC470_FWD_L[0][0][0][0] );
        var = lut_grp.getVar( "SC650_ALL" );
        var.getVar( startp, countp, &sc_lut->SC650_ALL_L[0][0][0][0] );
        var = lut_grp.getVar( "SC650_FWD" );
        var.getVar( startp, countp, &sc_lut->SC650_FWD_L[0][0][0][0] );
    } else {
        countp[4] = dateline;
        NcVar var = lut_grp.getVar( "SC412_ALL" );
        var.getVar( startp, countp, &sc_lut->SC412_ALL_L[0][0][0][0] );
        var = lut_grp.getVar( "SC412_FWD" );
        var.getVar( startp, countp, &sc_lut->SC412_FWD_L[0][0][0][0] );
        var = lut_grp.getVar( "SC470_ALL" );
        var.getVar( startp, countp, &sc_lut->SC470_ALL_L[0][0][0][0] );
        var = lut_grp.getVar( "SC470_FWD" );
        var.getVar( startp, countp, &sc_lut->SC470_FWD_L[0][0][0][0] );
        var = lut_grp.getVar( "SC650_ALL" );
        var.getVar( startp, countp, &sc_lut->SC650_ALL_L[0][0][0][0] );
        var = lut_grp.getVar( "SC650_FWD" );
        var.getVar( startp, countp, &sc_lut->SC650_FWD_L[0][0][0][0] );
        startp[4] = 0;
        countp[4] = edge[0] - dateline;
        var = lut_grp.getVar( "SC412_ALL" );
        var.getVar( startp, countp, &sc_lut->SC412_ALL_L[0][0][0][dateline] );
        var = lut_grp.getVar( "SC412_FWD" );
        var.getVar( startp, countp, &sc_lut->SC412_FWD_L[0][0][0][dateline] );
        var = lut_grp.getVar( "SC470_ALL" );
        var.getVar( startp, countp, &sc_lut->SC470_ALL_L[0][0][0][dateline] );
        var = lut_grp.getVar( "SC470_FWD" );
        var.getVar( startp, countp, &sc_lut->SC470_FWD_L[0][0][0][dateline] );
        var = lut_grp.getVar( "SC650_ALL" );
        var.getVar( startp, countp, &sc_lut->SC650_ALL_L[0][0][0][dateline] );
        var = lut_grp.getVar( "SC650_FWD" );
        var.getVar( startp, countp, &sc_lut->SC650_FWD_L[0][0][0][dateline] );
    }
        delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_rayleigh_files()
 *
 * DESCRIPTION: Read rayleigh files.
 *
 *************************************************************************/

int DbLutNetcdf::read_rayleigh_files( dbRayleighLUT* rl_lut )
{
    int status = DTDB_SUCCESS;

    string filepath1 = get_option(INPUT_RAYL_412);
    string filepath2 = get_option(INPUT_RAYL_470);
    string filepath3 = get_option(INPUT_RAYL_650);
    if(filepath1.empty() || filepath2.empty() || filepath3.empty()) {
        cerr << "DbLutNetcdf::read_rayleigh_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    bool isFileBigEndian = false;
    string   line1, line2, line3;
    int      nlines = 8556;
    int      nrows = 5;
    float    data1[nlines*nrows];
    float    data2[nlines*nrows];
    float    data3[nlines*nrows];
    ifstream fin1(filepath1.c_str());
    ifstream fin2(filepath2.c_str());
    ifstream fin3(filepath3.c_str());
    if(fin1.is_open() && fin2.is_open() && fin3.is_open()) {
        for (int iL=0; iL<nlines; iL++) {
            getline(fin1, line1);
            getline(fin2, line2);
            getline(fin3, line3);
            stringstream ss1(line1);
            stringstream ss2(line2);
            stringstream ss3(line3);
            for (int iR=0; iR<nrows; iR++) {
                ss1 >> data1[iL*nrows + iR];
                ss2 >> data2[iL*nrows + iR];
                ss3 >> data3[iL*nrows + iR];
            }
        }
        int iL = 0;
        for (int iR=0; iR<NRRAA; iR++) {
            for (int iT=0; iT<NVZA; iT++) {
                for (int iZ=0; iZ<NSZA; iZ++) {
                    for (int iS=0; iS<NSTOKES; iS++) {
                        rl_lut->RAYL_412[iR][iT][iZ][iS] = data1[iL];
                        rl_lut->RAYL_470[iR][iT][iZ][iS] = data2[iL];
                        rl_lut->RAYL_650[iR][iT][iZ][iS] = data3[iL];
                        iL++;
                   }
               }
           }
       }
    } else {
       cerr << "DbLutNetcdf::read_rayleigh_file() Error opening file "
               << endl;
       return DTDB_FAIL;
    }

    if ( isPlatformLittleEndian() && isFileBigEndian) {
        for( int iRAA=0; iRAA<NRAA; iRAA++) {
            for( int iTHE=0; iTHE<NVZA; iTHE++) {
                for( int iSZA=0; iSZA<NSZA; iSZA++) {
                    for( int iSTOKES=0; iSTOKES<NSTOKES; iSTOKES++) {
                        byteSwap(rl_lut->RAYL_412[iRAA][iTHE][iSZA][iSTOKES]);
                        byteSwap(rl_lut->RAYL_470[iRAA][iTHE][iSZA][iSTOKES]);
                        byteSwap(rl_lut->RAYL_650[iRAA][iTHE][iSZA][iSTOKES]);
                    }
                }
            }
        }
    }

    return status;
}

/**************************************************************************
 * NAME: write_rayleigh_lut()
 *
 * DESCRIPTION: Write rayleigh LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_rayleigh_lut( NcFile* nc_output,
        dbRayleighLUT* rl_lut  )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_RAYLEIGH );

    dim_nraa_ = lut_grp.addDim( "Dim_NRAA", NRRAA );
    dim_nvza_ = lut_grp.addDim( "Dim_NVZA", NVZA );
    dim_nsza_ = lut_grp.addDim( "Dim_NSZA", NSZA );
    dim_nstokes_ = lut_grp.addDim( "Dim_NSTOKES", NSTOKES );

    vector<NcDim> rayl_dims;
    rayl_dims.push_back(dim_nrraa_);
    rayl_dims.push_back(dim_nvza_);
    rayl_dims.push_back(dim_nsza_);
    rayl_dims.push_back(dim_nstokes_);

    NcVar var = lut_grp.addVar( "RAYL_412", ncFloat, rayl_dims );
    var.putVar( rl_lut->RAYL_412 );

    var = lut_grp.addVar( "RAYL_470", ncFloat, rayl_dims );
    var.putVar( rl_lut->RAYL_470 );

    var = lut_grp.addVar( "RAYL_650", ncFloat, rayl_dims );
    var.putVar( rl_lut->RAYL_650 );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_rayleigh_lut()
 *
 * DESCRIPTION: Read rayleigh NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_rayleigh_lut(  dbRayleighLUT* rl_lut  )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_RAYLEIGH );

    NcVar var = lut_grp.getVar( "RAYL_412" );
    var.getVar( rl_lut->RAYL_412 );

    var = lut_grp.getVar( "NVALX_470" );
    var.getVar( rl_lut->RAYL_470 );

    var = lut_grp.getVar( "NVALX_650" );
    var.getVar( rl_lut->RAYL_650 );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_landcover_file()
 *
 * DESCRIPTION: Read seasonal landcover hdf4 file.
 *
 *************************************************************************/

int DbLutNetcdf::read_landcover_files( dbLandcoverLUT* lc_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_LANDCOVER);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_landcover_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    int fileID;
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening landcover file: " + filepath << endl;
        return DTDB_FAIL;
    }
    int sds_index, sds_id, numtype, rank, nattrs;
    int start[2], edges[2], dims[2];
    start[0] = 0;
    start[1] = 0;
    edges[0] = 4*NLATS*10;
    edges[1] = NLONS*10;
    string sds_name = "Land_Vegetation_Type";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_landcover_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lc_lut->VEGETATION[0][0]);

    filepath = get_option(INPUT_VEG_LANDCOVER);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_landcover_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening global IGBP file: " + filepath << endl;
        return DTDB_FAIL;
    }
    start[0] = 0;
    start[1] = 0;
    edges[0] = NLATS*10;
    edges[1] = NLONS*10;
    sds_name = "IGBP_Land_Cover";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_landcover_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lc_lut->IGBP[0][0]);

    sds_name = "Region_Index";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_landcover_files() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lc_lut->REGION_INDEX[0][0]);

    filepath = get_option(INPUT_SEASONAL_DESERTS);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_landcover_files() Invalid path." << endl;
        return DTDB_FAIL;
    }
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening seasonal deserts file: " + filepath << endl;
        return DTDB_FAIL;
    }
    start[0] = 0;
    start[1] = 0;
    edges[0] = NLATS*10;
    edges[1] = NLONS*10;
    sds_name = "seasonal_desert_flag";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_deserts_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lc_lut->DESERTS_FLAG[0][0]);

    return status;
}

/**************************************************************************
 * NAME: write_landcover_lut()
 *
 * DESCRIPTION: Write landcover LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_landcover_lut( NcFile* nc_output,
        dbLandcoverLUT* lc_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_LANDCOVER );

    dim_1800_ = lut_grp.addDim( "Dim_1800", NLATS*10 );
    dim_3600_ = lut_grp.addDim( "Dim_3600", NLONS*10 );
    dim_7200_ = lut_grp.addDim( "Dim_7200", 4*NLATS*10 );

    vector<NcDim> vegie_dims;
    vegie_dims.push_back(dim_7200_);
    vegie_dims.push_back(dim_3600_);

    vector<NcDim> latlon_dims;
    latlon_dims.push_back(dim_1800_);
    latlon_dims.push_back(dim_3600_);

    NcVar var = lut_grp.addVar( "VEGETATION", ncInt, vegie_dims );
    var.putVar( lc_lut->VEGETATION );

    var = lut_grp.addVar( "IGBP", ncShort, latlon_dims );
    var.putVar( lc_lut->IGBP );

    var = lut_grp.addVar( "REGION_INDEX", ncShort, latlon_dims );
    var.putVar( lc_lut->REGION_INDEX );

    var = lut_grp.addVar( "DESERTS_FLAG", ncFloat, latlon_dims );
    var.putVar( lc_lut->DESERTS_FLAG );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_landcover_lut()
 *
 * DESCRIPTION: Read landcover NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_landcover_lut( dbLandcoverLUT* lc_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_LANDCOVER );

    NcVar var = lut_grp.getVar( "VEGETATION" );
    var.getVar( lc_lut->VEGETATION );

    var = lut_grp.getVar( "IGBP" );
    var.getVar( lc_lut->IGBP );

    var = lut_grp.getVar( "REGION_INDEX" );
    var.getVar( lc_lut->REGION_INDEX );

    var = lut_grp.getVar( "DESERTS_FLAG" );
    var.getVar( lc_lut->DESERTS_FLAG );

   delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_swir_file()
 *
 * DESCRIPTION: Read swir vs. vis hdf4 file.
 *
 *************************************************************************/

int DbLutNetcdf::read_swir_file( dbViirsSwirVsVisLUT* vsw_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_SWIR);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_swir_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    SEASON iseason = SEASON::NEVER;
    size_t rpos = 0;
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        size_t pos = filepath.find(str_season[iS]);
        if (pos != string::npos) {
            iseason = (SEASON) iS;
            rpos = pos;
        }
    }
    if (iseason == SEASON::NEVER) {
        cerr << "DbLutNetcdf:: Could not identify season in file name: " + filepath << endl;
        return DTDB_FAIL;
    }
    for (int iS=0; iS<NUM_SEASONS; iS++) {
        string rfilepath = filepath;
        rfilepath.replace(rpos, str_season[(int) iseason].size(), str_season[iS]);
        int fileID;
        try {
            fileID = SDstart(rfilepath.c_str(), DFACC_READ );
        }
        catch( std::exception& e) {
            e.what();
            cerr << "DbLutNetcdf:: Failure opening swir file: " + rfilepath << endl;
            return DTDB_FAIL;
        }
        int sds_index, sds_id, numtype, rank, nattrs;
        int start[3], edges[3], dims[3];
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        edges[0] = NSLATS;
        edges[1] = NSLONS;
        string sds_name = "latitude";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->latitude[iS][0][0]);
        sds_name = "longitude";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->longitude[iS][0][0]);
        sds_name = "coeffs_2250_to_412";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        edges[0] = NSCOEF;
        edges[1] = NSLATS;
        edges[2] = NSLONS;
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->coeffs_2250_to_412[iS][0][0][0]);
        sds_name = "coeffs_2250_to_488";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->coeffs_2250_to_488[iS][0][0][0]);
        sds_name = "coeffs_2250_to_670";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->coeffs_2250_to_670[iS][0][0][0]);
        sds_name = "min_2250_for_412";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        edges[0] = NSLATS;
        edges[1] = NSLONS;
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->min_2250_for_412[iS][0][0]);
        sds_name = "max_2250_for_412";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->max_2250_for_412[iS][0][0]);
        sds_name = "min_2250_for_488";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->min_2250_for_488[iS][0][0]);
        sds_name = "max_2250_for_488";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->max_2250_for_488[iS][0][0]);
        sds_name = "min_2250_for_670";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->min_2250_for_670[iS][0][0]);
        sds_name = "max_2250_for_670";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->max_2250_for_670[iS][0][0]);
        sds_name = "data_num_total";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->data_num_total[iS][0][0]);
        sds_name = "data_num_fitting";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->data_num_fitting[iS][0][0]);
        sds_name = "stderr_412";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->stderr_412[iS][0][0]);
        sds_name = "stderr_488";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->stderr_488[iS][0][0]);
        sds_name = "stderr_670";
        sds_index = SDnametoindex(fileID, sds_name.c_str());
        if (sds_index < 0) {
            cerr << "DbLutNetcdf::read_swir_file() " <<
                    "SDnametoindex() failure for "<< sds_name << endl;
            SDend(fileID);
            return DTDB_FAIL;
        }
        sds_id = SDselect(fileID, sds_index);
        SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
        SDreaddata(sds_id, start, NULL, edges, &vsw_lut->stderr_670[iS][0][0]);
    }

    return status;
}

/**************************************************************************
 * NAME: write_swir_lut()
 *
 * DESCRIPTION: Write swir LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_swir_lut( NcFile* nc_output,
        dbViirsSwirVsVisLUT* vsw_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_SWIR );

    dim_seasons_ = lut_grp.addDim( "Dim_seasons", NUM_SEASONS );
    dim_3_ = lut_grp.addDim( "Dim_3", NSCOEF );
    dim_3000_ = lut_grp.addDim( "Dim_3000", NSLATS );
    dim_6000_ = lut_grp.addDim( "Dim_6000", NSLONS );

    vector<NcDim> coef_dims;
    coef_dims.push_back(dim_seasons_);
    coef_dims.push_back(dim_3_);
    coef_dims.push_back(dim_3000_);
    coef_dims.push_back(dim_6000_);
    vector<NcDim> latlon_dims;
    latlon_dims.push_back(dim_seasons_);
    latlon_dims.push_back(dim_3000_);
    latlon_dims.push_back(dim_6000_);

    NcVar var = lut_grp.addVar( "latitude", ncFloat, latlon_dims );
    var.putVar( vsw_lut->latitude );

    var = lut_grp.addVar( "longitude", ncFloat, latlon_dims );
    var.putVar( vsw_lut->longitude );

    var = lut_grp.addVar( "coeffs_2250_to_412", ncFloat, coef_dims );
    var.putVar( vsw_lut->coeffs_2250_to_412 );

    var = lut_grp.addVar( "coeffs_2250_to_488", ncFloat, coef_dims );
    var.putVar( vsw_lut->coeffs_2250_to_488 );

    var = lut_grp.addVar( "coeffs_2250_to_670", ncFloat, coef_dims );
    var.putVar( vsw_lut->coeffs_2250_to_670 );

    var = lut_grp.addVar( "min_2250_for_412", ncFloat, latlon_dims );
    var.putVar( vsw_lut->min_2250_for_412 );

    var = lut_grp.addVar( "max_2250_for_412", ncFloat, latlon_dims );
    var.putVar( vsw_lut->max_2250_for_412 );

    var = lut_grp.addVar( "min_2250_for_488", ncFloat, latlon_dims );
    var.putVar( vsw_lut->min_2250_for_488 );

    var = lut_grp.addVar( "max_2250_for_488", ncFloat, latlon_dims );
    var.putVar( vsw_lut->max_2250_for_488 );

    var = lut_grp.addVar( "min_2250_for_670", ncFloat, latlon_dims );
    var.putVar( vsw_lut->min_2250_for_670 );

    var = lut_grp.addVar( "max_2250_for_670", ncFloat, latlon_dims );
    var.putVar( vsw_lut->max_2250_for_670 );

    var = lut_grp.addVar( "data_num_total", ncFloat, latlon_dims );
    var.putVar( vsw_lut->data_num_total );

    var = lut_grp.addVar( "data_num_fitting", ncFloat, latlon_dims );
    var.putVar( vsw_lut->data_num_fitting );

    var = lut_grp.addVar( "stderr_412", ncFloat, latlon_dims );
    var.putVar( vsw_lut->stderr_412 );

    var = lut_grp.addVar( "stderr_488", ncFloat, latlon_dims );
    var.putVar( vsw_lut->stderr_488 );

    var = lut_grp.addVar( "stderr_670", ncFloat, latlon_dims );
    var.putVar( vsw_lut->stderr_670 );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_swir_lut()
 *
 * DESCRIPTION: Read swir NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_swir_lut( dbViirsSwirVsVisLUT* vsw_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_SWIR );
    NcVar var = lut_grp.getVar( "latitude" );
    var.getVar( vsw_lut->latitude );
    var = lut_grp.getVar( "longitude" );
    var.getVar( vsw_lut->longitude );
    var = lut_grp.getVar( "coeffs_2250_to_412" );
    var.getVar( vsw_lut->coeffs_2250_to_412 );
    var = lut_grp.getVar( "coeffs_2250_to_488" );
    var.getVar( vsw_lut->coeffs_2250_to_488 );
    var = lut_grp.getVar( "coeffs_2250_to_670" );
    var.getVar( vsw_lut->coeffs_2250_to_670 );
    var = lut_grp.getVar( "min_2250_for_412" );
    var.getVar( vsw_lut->min_2250_for_412 );
    var = lut_grp.getVar( "max_2250_for_412" );
    var.getVar( vsw_lut->max_2250_for_412 );
    var = lut_grp.getVar( "min_2250_for_488" );
    var.getVar( vsw_lut->min_2250_for_488 );
    var = lut_grp.getVar( "max_2250_for_488" );
    var.getVar( vsw_lut->max_2250_for_488 );
    var = lut_grp.getVar( "min_2250_for_670" );
    var.getVar( vsw_lut->min_2250_for_670 );
    var = lut_grp.getVar( "max_2250_for_670" );
    var.getVar( vsw_lut->max_2250_for_670 );
    var = lut_grp.getVar( "data_num_total" );
    var.getVar( vsw_lut->data_num_total );
    var = lut_grp.getVar( "data_num_fitting" );
    var.getVar( vsw_lut->data_num_fitting );
    var = lut_grp.getVar( "stderr_412" );
    var.getVar( vsw_lut->stderr_412 );
    var = lut_grp.getVar( "stderr_488" );
    var.getVar( vsw_lut->stderr_488 );
    var = lut_grp.getVar( "stderr_670" );
    var.getVar( vsw_lut->stderr_670 );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_ocean_aero_file()
 *
 * DESCRIPTION: Read all viirs fine mode aerosol files.
 *
 *************************************************************************/

int DbLutNetcdf::read_ocean_aero_file( dbOceanAerosolLUT* lut,
        const string sType )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(sType);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_ocean_aero_file() Invalid path." << endl;
        return DTDB_FAIL;
    }
    int fileID;
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening ocean aerosol file: " + filepath << endl;
        return DTDB_FAIL;
    }

    if (sType == INPUT_AERO_OCEAN_FINE) {
        lut->nfmf = NFMF1;
        lut->naot = NAOT1;
    } else if (sType == INPUT_AERO_OCEAN_DUST) {
        lut->nfmf = NFMF2;
        lut->naot = NAOT1;
    } else if (sType == INPUT_AERO_OCEAN_MARI) {
        lut->nfmf = NFMF3;
        lut->naot = NAOT2;
    }else if (sType == INPUT_AERO_OCEAN_MIX) {
        lut->nfmf = NFMF4;
        lut->naot = NAOT1;
    } else {
        return DTDB_FAIL;
    };
    int sds_index, sds_id, numtype, rank, nattrs;
    int dims[7];
    int start[7] = {0,0,0,0,0,0,0};
    int edges[7] = {(int)NCHL,(int)NWS,(int)lut->nfmf,(int)lut->naot,(int)NVRAA,(int)NVVZA,(int)NVSZA};

    string sds_name = "IoverF_m03";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->m03[0][0][0][0][0][0][0]);
    sds_name = "IoverF_m04";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->m04[0][0][0][0][0][0][0]);
    sds_name = "IoverF_m05";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->m05[0][0][0][0][0][0][0]);
    sds_name = "IoverF_m07";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->m07[0][0][0][0][0][0][0]);
    sds_name = "IoverF_m08";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    edges[0] = NWS;
    edges[1] = lut->nfmf;
    edges[2] = lut->naot;
    edges[3] = NVRAA;
    edges[4] = NVVZA;
    edges[5] = NVSZA;
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->m08[0][0][0][0][0][0]);
    sds_name = "IoverF_m10";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->m10[0][0][0][0][0][0]);
    sds_name = "IoverF_m11";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->m11[0][0][0][0][0][0]);
    edges[0] = NVSZA;
    sds_name = "Solar_Zenith_Angle";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->sza[0]);
    edges[0] = NVVZA;
    sds_name = "View_Zenith_Angle";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->vza[0]);
    edges[0] = NVRAA;
    sds_name = "Relative_Azimuth_Angle";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->raa[0]);
    edges[0] = lut->naot;
    sds_name = "Aerosol_Optical_Depth_550";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->aot550[0]);
    edges[0] = lut->nfmf;
    sds_name = "Fine_Mode_Fraction_550";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->fmf[0]);
    edges[0] = NWS;
    sds_name = "Wind_Speed";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->wspd[0]);
    edges[0] = NCHL;
    sds_name = "Chl_Conc";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->chl[0]);
    edges[0] = NDBOWL;
    sds_name = "Band_Central_Wavelength";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->wave[0]);

    edges[0] = lut->nfmf;
    edges[1] = lut->naot;
    sds_name = "Angstrom_Exponent";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->ae[0][0]);
    edges[0] = NDBOWL;
    edges[1] = lut->nfmf;
    edges[2] = lut->naot;
    sds_name = "Spectral_Total_AOD";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->aot[0][0][0]);
    sds_name = "Spectral_Fine_AOD";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->fine_aot[0][0][0]);
    sds_name = "Spectral_Coarse_AOD";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->coarse_aot[0][0][0]);
    sds_name = "Spectral_Total_SSA";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->ssa[0][0][0]);
    sds_name = "Spectral_Fine_SSA";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->fine_ssa[0][0][0]);
    sds_name = "Spectral_Coarse_SSA";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->coarse_ssa[0][0][0]);
    sds_name = "Spectral_Total_ASY";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->asy[0][0][0]);
    sds_name = "Spectral_Fine_ASY";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->fine_asy[0][0][0]);
    sds_name = "Spectral_Coarse_ASY";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_ocean_aero_fine_file() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &lut->coarse_asy[0][0][0]);

    return status;
}

/**************************************************************************
 * NAME: write_ocean_aero_lut()
 *
 * DESCRIPTION: Write VIIRS Aerosol LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_ocean_aero_lut( NcFile* nc_output,
        dbOceanAerosolLUT* lut, const string sType )
{

    NcGroup lut_grp = nc_output->addGroup( sType );
    dim_nsza_ = lut_grp.addDim( "Dim_Solar_Zenith_Angle", NVSZA );
    dim_nvza_ = lut_grp.addDim( "Dim_View_Zenith_Angle", NVVZA );
    dim_nraa_ = lut_grp.addDim( "Dim_Relative_Azimuth_Angle", NVRAA );
    dim_ntau_ = lut_grp.addDim( "Dim_Aerosol_Optical_Depth_550", lut->naot );
    dim_nfmf_ = lut_grp.addDim( "Dim_Fine_Mode_Fraction_550", lut->nfmf );
    dim_nws_ = lut_grp.addDim( "Dim_Wind_Speed", NWS );
    dim_nchl_ = lut_grp.addDim( "Dim_Chl_Conc", NCHL );
    dim_nwl_ = lut_grp.addDim( "Dim_Band_Central_Wavelength", NDBOWL );
    vector<NcDim> ioverf_dims;
    ioverf_dims = {dim_nchl_,dim_nws_,dim_nfmf_,dim_ntau_,dim_nraa_,dim_nvza_,dim_nsza_};
    vector<size_t> start, edges;
    start = {0,0,0,0,0,0,0};
    edges = {NCHL,NWS,lut->nfmf,lut->naot,NVRAA,NVVZA,NVSZA};
    NcVar var = lut_grp.addVar( "IoverF_m03", ncFloat, ioverf_dims );
    var.putVar( start, edges, &lut->m03[0][0][0][0][0][0][0] );
    var = lut_grp.addVar( "IoverF_m04", ncFloat, ioverf_dims );
    var.putVar( start, edges, &lut->m04[0][0][0][0][0][0][0] );
    var = lut_grp.addVar( "IoverF_m05", ncFloat, ioverf_dims );
    var.putVar( start, edges, &lut->m05[0][0][0][0][0][0][0] );
    var = lut_grp.addVar( "IoverF_m07", ncFloat, ioverf_dims );
    var.putVar( start, edges, &lut->m07[0][0][0][0][0][0][0] );
    ioverf_dims = {dim_nws_,dim_nfmf_,dim_ntau_,dim_nraa_,dim_nvza_,dim_nsza_};
    edges = {NWS,lut->nfmf,lut->naot,NVRAA,NVVZA,NVSZA};
    var = lut_grp.addVar( "IoverF_m08", ncFloat, ioverf_dims );
    var.putVar( start, edges, &lut->m08[0][0][0][0][0][0] );
    var = lut_grp.addVar( "IoverF_m10", ncFloat, ioverf_dims );
    var.putVar( start, edges, &lut->m10[0][0][0][0][0][0] );
    var = lut_grp.addVar( "IoverF_m11", ncFloat, ioverf_dims );
    var.putVar( start, edges, &lut->m11[0][0][0][0][0][0] );
    var = lut_grp.addVar( "Solar_Zenith_Angle", ncFloat, dim_nsza_ );
    var.putVar( &lut->sza[0] );
    var = lut_grp.addVar( "View_Zenith_Angle", ncFloat, dim_nvza_ );
    var.putVar( &lut->vza[0] );
    var = lut_grp.addVar( "Relative_Azimuth_Angle", ncFloat, dim_nraa_ );
    var.putVar( &lut->raa[0] );
    var = lut_grp.addVar( "Aerosol_Optical_Depth_550", ncFloat, dim_ntau_ );
    edges = {lut->naot};
    var.putVar( start, edges, &lut->aot550[0] );
    var = lut_grp.addVar( "Fine_Mode_Fraction_550", ncFloat, dim_nfmf_ );
    edges = {lut->nfmf};
    var.putVar( start, edges, &lut->fmf[0] );
    var = lut_grp.addVar( "Wind_Speed", ncFloat, dim_nws_ );
    var.putVar( &lut->wspd[0] );
    var = lut_grp.addVar( "Chl_Conc", ncFloat, dim_nchl_ );
    var.putVar( &lut->chl[0] );
    var = lut_grp.addVar( "Band_Central_Wavelength", ncFloat, dim_nwl_ );
    var.putVar( &lut->wave[0] );
    vector<NcDim> ae_dims;
    ae_dims = {dim_nfmf_,dim_ntau_};
    edges = {lut->nfmf,lut->naot};
    var = lut_grp.addVar( "Angstrom_Exponent", ncFloat, ae_dims );
    var.putVar( start, edges, &lut->ae[0][0] );
    vector<NcDim> spectral_dims;
    spectral_dims = {dim_nwl_,dim_nfmf_,dim_ntau_};
    edges = {lut->nwave,lut->nfmf,lut->naot};
    var = lut_grp.addVar( "Spectral_Total_AOD", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->aot[0][0][0] );
    var = lut_grp.addVar( "Spectral_Fine_AOD", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->fine_aot[0][0][0] );
    var = lut_grp.addVar( "Spectral_Coarse_AOD", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->coarse_aot[0][0][0] );
    var = lut_grp.addVar( "Spectral_Total_SSA", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->ssa[0][0][0] );
    var = lut_grp.addVar( "Spectral_Fine_SSA", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->fine_ssa[0][0][0] );
    var = lut_grp.addVar( "Spectral_Coarse_SSA", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->coarse_ssa[0][0][0] );
    var = lut_grp.addVar( "Spectral_Total_ASY", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->asy[0][0][0] );
    var = lut_grp.addVar( "Spectral_Fine_ASY", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->fine_asy[0][0][0] );
    var = lut_grp.addVar( "Spectral_Coarse_ASY", ncFloat, spectral_dims );
    var.putVar( start, edges, &lut->coarse_asy[0][0][0] );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_ocean_aero_lut()
 *
 * DESCRIPTION: Read VIIRS aerosol NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_ocean_aero_lut( dbOceanAerosolLUMA* lut,
        const string sType )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    if (sType == LUT_OCEAN_AEROSOL_FINE) {
        lut->nfmf = NFMF1;
        lut->naot = NAOT1;
    } else if (sType == LUT_OCEAN_AEROSOL_DUST) {
        lut->nfmf = NFMF2;
        lut->naot = NAOT1;
    } else if (sType == LUT_OCEAN_AEROSOL_MARI) {
        lut->nfmf = NFMF3;
        lut->naot = NAOT2;
    }else if (sType == LUT_OCEAN_AEROSOL_MIX) {
        lut->nfmf = NFMF4;
        lut->naot = NAOT1;
    } else {
        return DTDB_FAIL;
    };

    vector<size_t> start = {0,0,0,0,0,0,0};
    vector<size_t> edges = {NCHL,NWS,lut->nfmf,lut->naot,NVRAA,NVVZA,NVSZA};
    NcGroup lut_grp = nc_input->getGroup( sType );
    lut->m03.resize(boost::extents[NCHL][NWS][lut->nfmf][lut->naot][NVRAA][NVVZA][NVSZA]);
    NcVar var = lut_grp.getVar( "IoverF_m03" );
    var.getVar( start, edges, &lut->m03[0][0][0][0][0][0][0] );
    lut->m04.resize(boost::extents[NCHL][NWS][lut->nfmf][lut->naot][NVRAA][NVVZA][NVSZA]);
    var = lut_grp.getVar( "IoverF_m04" );
    var.getVar( start, edges, &lut->m04[0][0][0][0][0][0][0] );
    lut->m05.resize(boost::extents[NCHL][NWS][lut->nfmf][lut->naot][NVRAA][NVVZA][NVSZA]);
    var = lut_grp.getVar( "IoverF_m05" );
    var.getVar( start, edges, &lut->m05[0][0][0][0][0][0][0] );
    lut->m07.resize(boost::extents[NCHL][NWS][lut->nfmf][lut->naot][NVRAA][NVVZA][NVSZA]);
    var = lut_grp.getVar( "IoverF_m07" );
    var.getVar( start, edges, &lut->m07[0][0][0][0][0][0][0] );

    start = {0,0,0,0,0,0,0};
    edges = {NWS,lut->nfmf,lut->naot,NVRAA,NVVZA,NVSZA};
    lut->m08.resize(boost::extents[NWS][lut->nfmf][lut->naot][NVRAA][NVVZA][NVSZA]);
    var = lut_grp.getVar( "IoverF_m08" );
    var.getVar( start, edges, &lut->m08[0][0][0][0][0][0] );
    lut->m10.resize(boost::extents[NWS][lut->nfmf][lut->naot][NVRAA][NVVZA][NVSZA]);
    var = lut_grp.getVar( "IoverF_m10" );
    var.getVar( start, edges, &lut->m10[0][0][0][0][0][0] );
    lut->m11.resize(boost::extents[NWS][lut->nfmf][lut->naot][NVRAA][NVVZA][NVSZA]);
    var = lut_grp.getVar( "IoverF_m11" );
    var.getVar( start, edges, &lut->m11[0][0][0][0][0][0] );

    lut->sza.resize(boost::extents[NVSZA]);
    var = lut_grp.getVar( "Solar_Zenith_Angle" );
    var.getVar( &lut->sza[0] );
    lut->vza.resize(boost::extents[NVVZA]);
    var = lut_grp.getVar( "View_Zenith_Angle" );
    var.getVar( &lut->vza[0] );
    lut->raa.resize(boost::extents[NVRAA]);
    var = lut_grp.getVar( "Relative_Azimuth_Angle" );
    var.getVar( &lut->raa[0] );
    lut->aot550.resize(boost::extents[lut->naot]);
    var = lut_grp.getVar( "Aerosol_Optical_Depth_550" );
    var.getVar( &lut->aot550[0] );
    lut->fmf.resize(boost::extents[lut->nfmf]);
    var = lut_grp.getVar( "Fine_Mode_Fraction_550" );
    var.getVar( &lut->fmf[0] );
    lut->wspd.resize(boost::extents[NWS]);
    var = lut_grp.getVar( "Wind_Speed" );
    var.getVar( &lut->wspd[0] );
    lut->chl.resize(boost::extents[NCHL]);
    var = lut_grp.getVar( "Chl_Conc" );
    var.getVar( &lut->chl[0] );
    lut->wave.resize(boost::extents[NDBOWL]);
    var = lut_grp.getVar( "Band_Central_Wavelength" );
    var.getVar( &lut->wave[0] );

    start = {0,0};
    edges = {lut->nfmf,lut->naot};
    lut->ae.resize(boost::extents[lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Angstrom_Exponent" );
    var.getVar( start, edges, &lut->ae[0][0] );
    start = {0,0,0};
    edges = {lut->nwave,lut->nfmf,lut->naot};
    lut->aot.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Total_AOD" );
    var.getVar( start, edges, &lut->aot[0][0][0] );
    lut->fine_aot.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Fine_AOD" );
    var.getVar( start, edges, &lut->fine_aot[0][0][0] );
    lut->coarse_aot.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Coarse_AOD" );
    var.getVar( start, edges, &lut->coarse_aot[0][0][0] );
    lut->ssa.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Total_SSA" );
    var.getVar( start, edges, &lut->ssa[0][0][0] );
    lut->fine_ssa.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Fine_SSA" );
    var.getVar( start, edges, &lut->fine_ssa[0][0][0] );
    lut->coarse_ssa.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Coarse_SSA" );
    var.getVar( start, edges, &lut->coarse_ssa[0][0][0] );
    lut->asy.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Total_ASY" );
    var.getVar( start, edges, &lut->asy[0][0][0] );
    lut->fine_asy.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Fine_ASY" );
    var.getVar( start, edges, &lut->fine_asy[0][0][0] );
    lut->coarse_asy.resize(boost::extents[lut->nwave][lut->nfmf][lut->naot]);
    var = lut_grp.getVar( "Spectral_Coarse_ASY" );
    var.getVar( start, edges, &lut->coarse_asy[0][0][0] );

    delete nc_input;

   return status;
}


/**************************************************************************
 * NAME: read_land_aero_file()
 *
 * DESCRIPTION: Read all land aerosol files.
 *
 *************************************************************************/

int DbLutNetcdf::read_land_aero_file( dbLandAerosolLUT* va_lut,
        const string aero_input)
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(aero_input);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_aero_land_files() Invalid path." << endl;
        return DTDB_FAIL;
    }
    int fileID;
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening land aerosol file: " + filepath << endl;
        return DTDB_FAIL;
    }
    int sds_index, sds_id, numtype, rank, nattrs;
    int start[6], edges[6], dims[6];
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    start[4] = 0;
    start[5] = 0;
    edges[0] = NSZA;
    string sds_name = "SZA412_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SZA412_Nodes[0]);
    edges[0] = NVZA;
    sds_name = "VZA412_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->VZA412_Nodes[0]);
    edges[0] = NTAU;
    sds_name = "AOT412_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->AOT412_Nodes[0]);
    edges[0] = NRAA;
    sds_name = "RAA412_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->RAA412_Nodes[0]);
    edges[0] = SSA412;
    sds_name = "SSA412_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SSA412_Nodes[0]);
    edges[0] = SR412;
    sds_name = "SR412_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SR412_Nodes[0]);
    edges[0] = SR412;
    edges[1] = SSA412;
    edges[2] = NTAU;
    edges[3] = NRAA;
    edges[4] = NVZA;
    edges[5] = NSZA;
    sds_name = "nvalx412";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_land_aero_files() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->nvalx412[0][0][0][0][0][0]);

    edges[0] = NSZA;
    sds_name = "SZA488_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SZA488_Nodes[0]);
    edges[0] = NVZA;
    sds_name = "VZA488_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->VZA488_Nodes[0]);
    edges[0] = NTAU;
    sds_name = "AOT488_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->AOT488_Nodes[0]);
    edges[0] = NRAA;
    sds_name = "RAA488_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->RAA488_Nodes[0]);
    edges[0] = SSA470;
    sds_name = "SSA488_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SSA488_Nodes[0]);
    edges[0] = SR470;
    sds_name = "SR488_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SR488_Nodes[0]);
    edges[0] = SR470;
    edges[1] = SSA470;
    edges[2] = NTAU;
    edges[3] = NRAA;
    edges[4] = NVZA;
    edges[5] = NSZA;
    sds_name = "nvalx488";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_land_aero_files() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->nvalx488[0][0][0][0][0][0]);

    edges[0] = NSZA;
    sds_name = "SZA672_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SZA672_Nodes[0]);
    edges[0] = NVZA;
    sds_name = "VZA672_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->VZA672_Nodes[0]);
    edges[0] = NTAU;
    sds_name = "AOT672_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->AOT672_Nodes[0]);
    edges[0] = NRAA;
    sds_name = "RAA672_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->RAA672_Nodes[0]);
    edges[0] = SSA650;
    sds_name = "SSA672_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SSA672_Nodes[0]);
    edges[0] = SR650;
    sds_name = "SR672_Nodes";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->SR672_Nodes[0]);
    edges[0] = SR650;
    edges[1] = SSA650;
    edges[2] = NTAU;
    edges[3] = NRAA;
    edges[4] = NVZA;
    edges[5] = NSZA;
    sds_name = "nvalx672";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    if (sds_index < 0) {
        cerr << "DbLutNetcdf::read_land_aero_files() " <<
                "SDnametoindex() failure for "<< sds_name << endl;
        SDend(fileID);
        return DTDB_FAIL;
    }
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &va_lut->nvalx672[0][0][0][0][0][0]);

    SDend(fileID);
    return status;
}

/**************************************************************************
 * NAME: write_land_aero_lut()
 *
 * DESCRIPTION: Write Land aerosol LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_land_aero_lut( NcFile* nc_output,
        dbLandAerosolLUT* va_lut, const string aero_lut )
{
    NcGroup lut_grp = nc_output->addGroup( aero_lut );

    dim_nsza_ = lut_grp.addDim( "Dim_Solar_Zenith_Angle", NSZA );
    dim_nvza_ = lut_grp.addDim( "Dim_View_Zenith_Angle", NVZA );
    dim_nraa_ = lut_grp.addDim( "Dim_Relative_Azimuth_Angle", NRAA );
    dim_ntau_ = lut_grp.addDim( "Dim_Aerosol_Optical_Depth", NTAU );
    dim_nssa_ = lut_grp.addDim( "Dim_SSA_412", SSA412 );
    dim_nsr_ = lut_grp.addDim( "Dim_SR_412", SR412 );

    vector<NcDim> nvalx_dims;
    nvalx_dims.push_back(dim_nsr_);
    nvalx_dims.push_back(dim_nssa_);
    nvalx_dims.push_back(dim_ntau_);
    nvalx_dims.push_back(dim_nraa_);
    nvalx_dims.push_back(dim_nvza_);
    nvalx_dims.push_back(dim_nsza_);
    NcVar var = lut_grp.addVar( "SZA412_Nodes", ncFloat, dim_nsza_ );
    var.putVar( va_lut->SZA412_Nodes );
    var = lut_grp.addVar( "VZA412_Nodes", ncFloat, dim_nvza_ );
    var.putVar( va_lut->VZA412_Nodes );
    var = lut_grp.addVar( "RAA412_Nodes", ncFloat, dim_nraa_ );
    var.putVar( va_lut->RAA412_Nodes );
    var = lut_grp.addVar( "AOT412_Nodes", ncFloat, dim_ntau_ );
    var.putVar( va_lut->AOT412_Nodes );
    var = lut_grp.addVar( "SSA412_Nodes", ncFloat, dim_nssa_ );
    var.putVar( va_lut->SSA412_Nodes );
    var = lut_grp.addVar( "SR412_Nodes", ncFloat, dim_nsr_ );
    var.putVar( va_lut->SR412_Nodes );
    var = lut_grp.addVar( "NVALX412", ncFloat, nvalx_dims );
    var.putVar( va_lut->nvalx412 );

    dim_nssa_ = lut_grp.addDim( "Dim_SSA_488", SSA470 );
    dim_nsr_ = lut_grp.addDim( "Dim_SR_488", SR470 );
    nvalx_dims.clear();
    nvalx_dims.push_back(dim_nsr_);
    nvalx_dims.push_back(dim_nssa_);
    nvalx_dims.push_back(dim_ntau_);
    nvalx_dims.push_back(dim_nraa_);
    nvalx_dims.push_back(dim_nvza_);
    nvalx_dims.push_back(dim_nsza_);
    var = lut_grp.addVar( "SZA488_Nodes", ncFloat, dim_nsza_ );
    var.putVar( va_lut->SZA488_Nodes );
    var = lut_grp.addVar( "VZA488_Nodes", ncFloat, dim_nvza_ );
    var.putVar( va_lut->VZA488_Nodes );
    var = lut_grp.addVar( "RAA488_Nodes", ncFloat, dim_nraa_ );
    var.putVar( va_lut->RAA488_Nodes );
    var = lut_grp.addVar( "AOT488_Nodes", ncFloat, dim_ntau_ );
    var.putVar( va_lut->AOT488_Nodes );
    var = lut_grp.addVar( "SSA488_Nodes", ncFloat, dim_nssa_ );
    var.putVar( va_lut->SSA488_Nodes );
    var = lut_grp.addVar( "SR488_Nodes", ncFloat, dim_nsr_ );
    var.putVar( va_lut->SR488_Nodes );
    var = lut_grp.addVar( "NVALX488", ncFloat, nvalx_dims );
    var.putVar( va_lut->nvalx488 );

    dim_nssa_ = lut_grp.addDim( "Dim_SSA_672", SSA650 );
    dim_nsr_ = lut_grp.addDim( "Dim_SR_672", SR650 );
    nvalx_dims.clear();
    nvalx_dims.push_back(dim_nsr_);
    nvalx_dims.push_back(dim_nssa_);
    nvalx_dims.push_back(dim_ntau_);
    nvalx_dims.push_back(dim_nraa_);
    nvalx_dims.push_back(dim_nvza_);
    nvalx_dims.push_back(dim_nsza_);
    var = lut_grp.addVar( "SZA672_Nodes", ncFloat, dim_nsza_ );
    var.putVar( va_lut->SZA672_Nodes );
    var = lut_grp.addVar( "VZA672_Nodes", ncFloat, dim_nvza_ );
    var.putVar( va_lut->VZA672_Nodes );
    var = lut_grp.addVar( "RAA672_Nodes", ncFloat, dim_nraa_ );
    var.putVar( va_lut->RAA672_Nodes );
    var = lut_grp.addVar( "AOT672_Nodes", ncFloat, dim_ntau_ );
    var.putVar( va_lut->AOT672_Nodes );
    var = lut_grp.addVar( "SSA672_Nodes", ncFloat, dim_nssa_ );
    var.putVar( va_lut->SSA672_Nodes );
    var = lut_grp.addVar( "SR672_Nodes", ncFloat, dim_nsr_ );
    var.putVar( va_lut->SR672_Nodes );
    var = lut_grp.addVar( "NVALX672", ncFloat, nvalx_dims );
    var.putVar( va_lut->nvalx672 );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_land_aero_lut()
 *
 * DESCRIPTION: Read land aerosol NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_land_aero_lut( dbLandAerosolLUT* va_lut,
        const string aero_lut)
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( aero_lut );

    NcVar var = lut_grp.getVar( "SZA412_Nodes" );
    var.getVar( va_lut->SZA412_Nodes );
    var = lut_grp.getVar( "VZA412_Nodes" );
    var.getVar( va_lut->VZA412_Nodes );
    var = lut_grp.getVar( "RAA412_Nodes" );
    var.getVar( va_lut->RAA412_Nodes );
    var = lut_grp.getVar( "AOT412_Nodes" );
    var.getVar( va_lut->AOT412_Nodes );
    var = lut_grp.getVar( "SSA412_Nodes" );
    var.getVar( va_lut->SSA412_Nodes );
    var = lut_grp.getVar( "SR412_Nodes" );
    var.getVar( va_lut->SR412_Nodes );
    var = lut_grp.getVar( "NVALX412" );
    var.getVar( va_lut->nvalx412 );

    var = lut_grp.getVar( "SZA488_Nodes" );
    var.getVar( va_lut->SZA488_Nodes );
    var = lut_grp.getVar( "VZA488_Nodes" );
    var.getVar( va_lut->VZA488_Nodes );
    var = lut_grp.getVar( "RAA488_Nodes" );
    var.getVar( va_lut->RAA488_Nodes );
    var = lut_grp.getVar( "AOT488_Nodes" );
    var.getVar( va_lut->AOT488_Nodes );
    var = lut_grp.getVar( "SSA488_Nodes" );
    var.getVar( va_lut->SSA488_Nodes );
    var = lut_grp.getVar( "SR488_Nodes" );
    var.getVar( va_lut->SR488_Nodes );
    var = lut_grp.getVar( "NVALX488" );
    var.getVar( va_lut->nvalx488 );

    var = lut_grp.getVar( "SZA672_Nodes" );
    var.getVar( va_lut->SZA672_Nodes );
    var = lut_grp.getVar( "VZA672_Nodes" );
    var.getVar( va_lut->VZA672_Nodes );
    var = lut_grp.getVar( "RAA672_Nodes" );
    var.getVar( va_lut->RAA672_Nodes );
    var = lut_grp.getVar( "AOT672_Nodes" );
    var.getVar( va_lut->AOT672_Nodes );
    var = lut_grp.getVar( "SSA672_Nodes" );
    var.getVar( va_lut->SSA672_Nodes );
    var = lut_grp.getVar( "SR672_Nodes" );
    var.getVar( va_lut->SR672_Nodes );
    var = lut_grp.getVar( "NVALX672" );
    var.getVar( va_lut->nvalx672 );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_bathymetry_files()
 *
 * DESCRIPTION: Read all bathymetry files.
 *
 *************************************************************************/

int DbLutNetcdf::read_bathymetry_files( dbBathymetryLUT* b_lut)
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_BATHYMETRY);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_bathymetry_files() Invalid path." << endl;
        return DTDB_FAIL;
    }
    int fileID;
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening bathymetry file: " + filepath << endl;
        return DTDB_FAIL;
    }
    int sds_index, sds_id, numtype, rank, nattrs;
    int start[2], edges[2], dims[2];
    start[0] = 0;
    start[1] = 0;
    edges[0] = NUMY;
    edges[1] = NUMX;
    string sds_name = "z";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &b_lut->z[0][0]);
    edges[0] = NUMX;
    sds_name = "x";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &b_lut->x[0]);
    edges[0] = NUMY;
    sds_name = "y";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &b_lut->y[0]);

    SDend(fileID);
    return status;
}

/**************************************************************************
 * NAME: write_bathymetry_lut()
 *
 * DESCRIPTION: Write bathymetry LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_bathymetry_lut( NcFile* nc_output,
        dbBathymetryLUT* b_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_BATHYMETRY );

    dim_nx_ = lut_grp.addDim( "Dim_Bathymetry_X", NUMX );
    dim_ny_ = lut_grp.addDim( "Dim_Bathymetry_Y", NUMY );

    vector<NcDim> bath_dims;
    bath_dims.push_back(dim_ny_);
    bath_dims.push_back(dim_nx_);
    NcVar var = lut_grp.addVar( "X", ncDouble, dim_nx_ );
    var.putVar( b_lut->x );
    var = lut_grp.addVar( "Y", ncDouble, dim_ny_ );
    var.putVar( b_lut->y );
    var = lut_grp.addVar( "Z", ncInt, bath_dims );
    var.putVar( b_lut->z );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_bathymetry_lut()
 *
 * DESCRIPTION: Read bathymetry NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_bathymetry_lut( dbBathymetryLUT* b_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_BATHYMETRY );

    NcVar var = lut_grp.getVar( "X" );
    var.getVar( b_lut->x );
    var = lut_grp.getVar( "Y" );
    var.getVar( b_lut->y );
    var = lut_grp.getVar( "Z" );
    var.getVar( b_lut->z );
    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_chl_files()
 *
 * DESCRIPTION: Read all chl files.
 *
 *************************************************************************/

int DbLutNetcdf::read_chl_files( dbChlLUT* c_lut)
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_CHL);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_chl_files() Invalid path." << endl;
        return DTDB_FAIL;
    }
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening chl file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcVar var = nc_input->getVar( "latitude" );
    var.getVar( c_lut->latitude );
    var = nc_input->getVar( "longitude" );
    var.getVar( c_lut->longitude );
    var = nc_input->getVar( "time" );
    var.getVar( c_lut->time );
    var = nc_input->getVar( "log_chl" );
    var.getVar( c_lut->log_chl );

    delete nc_input;
    return status;
}

/**************************************************************************
 * NAME: write_chl_lut()
 *
 * DESCRIPTION: Write chl LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_chl_lut( NcFile* nc_output,
        dbChlLUT* c_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_CHL );

    dim_1800_ = lut_grp.addDim( "Dim_Latitude", NLATS*10 );
    dim_3600_ = lut_grp.addDim( "Dim_Longitude", NLONS*10 );
    dim_months_ = lut_grp.addDim( "Dim_Time", NMONTHS );
    vector<NcDim> map_dims;
    map_dims.push_back(dim_1800_);
    map_dims.push_back(dim_3600_);
    vector<NcDim> chl_dims;
    chl_dims.push_back(dim_1800_);
    chl_dims.push_back(dim_3600_);
    chl_dims.push_back(dim_months_);
    NcVar var = lut_grp.addVar( "LATITUDE", ncFloat, map_dims );
    var.putVar( c_lut->latitude );
    var = lut_grp.addVar( "LONGITUDE", ncFloat, map_dims );
    var.putVar( c_lut->latitude );
    var = lut_grp.addVar( "TIME", ncDouble, dim_months_ );
    var.putVar( c_lut->time );
    var = lut_grp.addVar( "LOG_CHL", ncFloat, chl_dims );
    var.putVar( c_lut->log_chl );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_chl_lut()
 *
 * DESCRIPTION: Read chl NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_chl_lut( dbChlLUT* c_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_CHL );

    NcVar var = lut_grp.getVar( "LATITUDE" );
    var.getVar( c_lut->latitude);
    var = lut_grp.getVar( "LONGITUDE" );
    var.getVar( c_lut->longitude );
    var = lut_grp.getVar( "TIME" );
    var.getVar( c_lut->time );
    var = lut_grp.getVar( "LOG_CHL" );
    var.getVar( c_lut->log_chl );

    delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_geozone_files()
 *
 * DESCRIPTION: Read all geozone files.
 *
 *************************************************************************/

int DbLutNetcdf::read_geozone_files( dbGeozoneLUT* g_lut)
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_GEOZONE);
    if(filepath.empty()) {
        cerr << "DbLutNetcdf::read_geozone_files() Invalid path." << endl;
        return DTDB_FAIL;
    }
    int fileID;
    try {
        fileID = SDstart(filepath.c_str(), DFACC_READ );
    }
    catch( std::exception& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening geozone file: " + filepath << endl;
        return DTDB_FAIL;
    }
    int sds_index, sds_id, numtype, rank, nattrs;
    int start[3], edges[3], dims[3];
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    edges[0] = NLATS*10;
    edges[1] = NLONS*10;
    string sds_name = "geographical_zone_flag";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &g_lut->GEOZONE_FLAG[0][0]);
    sds_name = "surface_elevation_stddev";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &g_lut->ELEVATION_STDV[0][0]);
    edges[0] = NSEASONS;
    edges[1] = NLATS;
    edges[2] = NLONS;
    sds_name = "background_aod";
    sds_index = SDnametoindex(fileID, sds_name.c_str());
    sds_id = SDselect(fileID, sds_index);
    SDgetinfo(sds_id, (char*) sds_name.c_str(), &rank, dims, &numtype, &nattrs);
    SDreaddata(sds_id, start, NULL, edges, &g_lut->BACKGROUND_AOD[0]);

    SDend(fileID);
    return status;
}

/**************************************************************************
 * NAME: write_geozone_lut()
 *
 * DESCRIPTION: Write geozone LUT to NetCDF4 file.
 *
 *************************************************************************/

int DbLutNetcdf::write_geozone_lut( NcFile* nc_output,
        dbGeozoneLUT* g_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_GEOZONE );

    dim_1800_ = lut_grp.addDim( "Dim_Latx10", NLATS*10 );
    dim_3600_ = lut_grp.addDim( "Dim_Lonx10", NLONS*10 );
    dim_180_ = lut_grp.addDim( "Dim_Lat", NLATS );
    dim_360_ = lut_grp.addDim( "Dim_Lon", NLONS );
    dim_seasons_ = lut_grp.addDim( "Dim_Seasons", NSEASONS );
    vector<NcDim> map_dims;
    map_dims.push_back(dim_1800_);
    map_dims.push_back(dim_3600_);
    vector<NcDim> aod_dims;
    aod_dims.push_back(dim_seasons_);
    aod_dims.push_back(dim_180_);
    aod_dims.push_back(dim_360_);
    NcVar var = lut_grp.addVar( "GEOZONE_FLAG", ncFloat, map_dims );
    var.putVar( g_lut->GEOZONE_FLAG );
    var = lut_grp.addVar( "ELEVATION_STDV", ncFloat, map_dims );
    var.putVar( g_lut->ELEVATION_STDV );
    var = lut_grp.addVar( "BACKGROUND_AOD", ncFloat, aod_dims );
    var.putVar( g_lut->BACKGROUND_AOD );

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_geozone_lut()
 *
 * DESCRIPTION: Read geozone NetCDF4 LUT.
 *
 *************************************************************************/

int DbLutNetcdf::read_geozone_lut( dbGeozoneLUT* g_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DB_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DbLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_GEOZONE );

    NcVar var = lut_grp.getVar( "GEOZONE_FLAG" );
    var.getVar( g_lut->GEOZONE_FLAG);
    var = lut_grp.getVar( "ELEVATION_STDV" );
    var.getVar( g_lut->ELEVATION_STDV );
    var = lut_grp.getVar( "BACKGROUND_AOD" );
    var.getVar( g_lut->BACKGROUND_AOD );

    delete nc_input;

   return status;
}

/**************************************************************************
 * write_global_attributes()
 *
 * Write global attributes to specified netCDF file ID
 *
 **************************************************************************/

int DbLutNetcdf::write_global_attributes( NcFile* nc_output )
{
	nc_output->putAtt( "processing_version", processing_version_);
	nc_output->putAtt( "Conventions", Conventions_);
	nc_output->putAtt( "institution", institution_);
	nc_output->putAtt( "license", license_);
	nc_output->putAtt( "naming_authority", naming_authority_);
	nc_output->putAtt( "date_created", date_created_);
    nc_output->putAtt( "ProductionTime", date_created_);
	nc_output->putAtt( "keywords_vocabulary", keywords_vocabulary_);
	nc_output->putAtt( "stdname_vocabulary", stdname_vocabulary_);
	nc_output->putAtt( "creator_name", creator_name_);
	nc_output->putAtt( "creator_email", creator_email_);
	nc_output->putAtt( "creator_url", creator_url_);
	nc_output->putAtt( "project", project_);
	nc_output->putAtt( "publisher_name", publisher_name_);
	nc_output->putAtt( "publisher_url", publisher_url_);
	nc_output->putAtt( "publisher_email", publisher_email_);
	nc_output->putAtt( "processing_level", processing_level_);
	nc_output->putAtt( "cdm_data_type", cdm_data_type_);
	nc_output->putAtt( "orbit_number", ncInt, orbit_number_);
	nc_output->putAtt( "history", history_);
	nc_output->putAtt( "source", source_files_);
	nc_output->putAtt( "time_coverage_start", time_coverage_start_);
	nc_output->putAtt( "time_coverage_end", time_coverage_end_);
	string pge_name_ = basename((char*)get_option("PGE_Name").c_str());
	if (!pge_name_.empty()) {
		nc_output->putAtt("PGE_Name",pge_name_);
	}
	string versionid_ =  basename((char*)get_option("VersionID").c_str());
	if (!versionid_.empty()) {
		nc_output->putAtt("VersionId",versionid_);
	}

	nc_output->putAtt( "format_version", ncInt, format_version_);
	nc_output->putAtt( "instrument_number", ncInt, instrument_number_);

	return DTDB_SUCCESS;
}


/**************************************************************************
 * NAME: isPlatformLittleEndian()
 *
 * DESCRIPTION: Determine if target platform is little endian.
 * Return true if platform is little endian.
 *
 *************************************************************************/

bool DbLutNetcdf::isPlatformLittleEndian()
{
    unsigned short checkValue = 0xAABB;
    unsigned char* bytePtr = reinterpret_cast<unsigned char*>(&checkValue);

    if (bytePtr[0] == 0xAA) {  // big-endian
        return false;
    } else {  // little-endian
        return true;
    }
}



