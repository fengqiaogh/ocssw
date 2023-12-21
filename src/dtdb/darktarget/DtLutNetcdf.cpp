/******************************************************************************
 *  NAME: DtLutNetcdf.cpp
 *
 *  DESCRIPTION: Object class that generates a netCDF4 LUT for NASA Dark Target
 *  aerosols algorithm
 *
 *  Created on: April 25, 2017
 *      Author: Sam Anderson
 *
 *
 ******************************************************************************/

#include "darktarget/DtLutNetcdf.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <mfhdf.h>
#include <libgen.h>

#include <DDProcess.h>
#include <DDOptions.h>

/**************************************************************************
 * NAME: DtLutNetcdf()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtLutNetcdf::DtLutNetcdf() {
}

/**************************************************************************
 * NAME: ~DtLutNetcdf()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DtLutNetcdf::~DtLutNetcdf() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Initializes data and object classes for granule
 *
 *************************************************************************/

int DtLutNetcdf::initialize() {
    return DTDB_SUCCESS;
}


/**************************************************************************
 * NAME: create_dt_nc4_lut()
 *
 * DESCRIPTION: Create dark target aerosol netCDF4 LUT.
 *
 *************************************************************************/

int DtLutNetcdf::create_dt_nc4_lut() {

	int status = DTDB_SUCCESS;
	int istatus = DTDB_SUCCESS;
    int num_good = 0;

	string path = get_option(NETCDF_LUT_PATH);
	string filepath = path + "/VIIRS_DARKTARGET_LUT_" + "version_source" + ".nc";

	NcFile* nc_output;

	try {
		nc_output = new NcFile( filepath, NcFile::replace );
	}
	catch( NcException& e) {
		e.what();
		cerr << "DtLutNetcdf:: Failure creating granule netCDF4 LUT file: " + filepath + "." << endl;
		return DTDB_FAIL;
	}

	nc_output->putAtt( "title", "VIIRS DARKTARGET LUTs" );

	write_global_attributes( nc_output );

// Read input files and create LUTs

	dtGasCorrectionLUT* gc_lut = new dtGasCorrectionLUT;
	istatus = read_gas_correction_file( gc_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading Gas Correction file " << endl;
		status = istatus;
	}
	istatus = write_gas_correction_lut( nc_output, gc_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure writing Gas Correction file to netCDF4 LUT file " << endl;
		status = istatus;
	}
	delete gc_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DtLutNetcdf:: Created " + LUT_GAS_CORRECTION  << endl;
    }

	dtLandAerosolLUT* la_lut = new dtLandAerosolLUT;
	istatus = read_land_aerosol_file( INPUT_LAND_W0466, 0, la_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading land aerosol file " << endl;
		status = istatus;
	}
	istatus = read_land_aerosol_file( INPUT_LAND_W0554, 1, la_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading land aerosol file " << endl;
		status = istatus;
	}
	istatus = read_land_aerosol_file( INPUT_LAND_W0645, 2, la_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading land aerosol file " << endl;
		status = istatus;
	}
	istatus = read_land_aerosol_file( INPUT_LAND_W2113, 3, la_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading land aerosol file " << endl;
		status = istatus;
	}
	istatus = read_land_aerosol_map( INPUT_LAND_MAP, la_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading land aerosol map file " << endl;
		status = istatus;
	}
	istatus = write_land_aerosol_lut( nc_output, la_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure writing land aerosol file to netCDF4 LUT file " << endl;
		status = istatus;
	}
	delete la_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DtLutNetcdf:: Created " + LUT_LAND_AEROSOL  << endl;
    }

	dtOceanAerosolLUT* oa_lut = new dtOceanAerosolLUT;
	istatus = read_ocean_small_aerosol_file( INPUT_OCEAN_SMALL1, 0, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean small aerosol file " << endl;
		status = istatus;
	}
	istatus = read_ocean_small_aerosol_file( INPUT_OCEAN_SMALL2, 1, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean small aerosol file " << endl;
		status = istatus;
	}
	istatus = read_ocean_small_aerosol_file( INPUT_OCEAN_SMALL3, 2, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean small aerosol file " << endl;
		status = istatus;
	}
	istatus = read_ocean_small_aerosol_file( INPUT_OCEAN_SMALL4, 3, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean small aerosol file " << endl;
		status = istatus;
	}

	istatus = read_ocean_big_aerosol_file( INPUT_OCEAN_BIG1, 0, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean big aerosol file " << endl;
		status = istatus;
	}
	istatus = read_ocean_big_aerosol_file( INPUT_OCEAN_BIG2, 1, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean big aerosol file " << endl;
		status = istatus;
	}
	istatus = read_ocean_big_aerosol_file( INPUT_OCEAN_BIG3, 2, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean big aerosol file " << endl;
		status = istatus;
	}
	istatus = read_ocean_big_aerosol_file( INPUT_OCEAN_BIG4, 3, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure reading ocean big aerosol file " << endl;
		status = istatus;
	}
	istatus = write_ocean_aerosol_lut( nc_output, oa_lut );
	if ( istatus != DTDB_SUCCESS ) {
		cerr << "DtLutNetcdf:: Failure writing ocean big aerosol file to netCDF4 LUT file " << endl;
		status = istatus;
	}
	delete oa_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DtLutNetcdf:: Created " + LUT_OCEAN_AEROSOL  << endl;
    }

    dtWaterVaporLUT* wv_lut = new dtWaterVaporLUT;
    istatus = read_transm_h2o_file( INPUT_TRANSM_H2O_1, 1, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading Transmission H20 file " << endl;
        status = istatus;
    }

    istatus = read_transm_h2o_file( INPUT_TRANSM_H2O_2, 2, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading Transmission H20 file " << endl;
        status = istatus;
    }

    istatus = read_transm_h2o_file( INPUT_TRANSM_H2O_3, 3, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading Transmission H20 file " << endl;
        status = istatus;
    }

    istatus = read_transm_h2o_file( INPUT_TRANSM_H2O_4, 4, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading Transmission H20 file " << endl;
        status = istatus;
    }

    istatus = read_transm_h2o_file( INPUT_TRANSM_H2O_5, 5, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading Transmission H20 file " << endl;
        status = istatus;
    }

    istatus = read_transm_h2o_file( INPUT_TRANSM_H2O_6, 6, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading Transmission H20 file " << endl;
        status = istatus;
    }
/*
    istatus = read_weight_table_file( INPUT_WEIGHT_TABLE, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading weight table file " << endl;
        status = istatus;
    }

    istatus = read_ch2_reflectance_file( INPUT_REFL_CH2, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading CH2 reflectance file " << endl;
        status = istatus;
    }

    istatus = read_ch19_to_ch2_ratio_file( INPUT_RATIO_CH19_TO_CH2, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure reading ch19_to_ch2_ratio file " << endl;
        status = istatus;
    }
*/
    istatus = write_water_vapor_lut( nc_output, wv_lut );
    if ( istatus != DTDB_SUCCESS ) {
        cerr << "DtLutNetcdf:: Failure writing Water Vapor files to netCDF4 LUT file " << endl;
        status = istatus;
    }
    delete wv_lut;
    if (istatus == DTDB_SUCCESS) {
        num_good++;
        cerr << "DtLutNetcdf:: Created " + LUT_WATER_VAPOR  << endl;
    }

	delete nc_output;

	return status;
}

/**************************************************************************
 * NAME: read_grib_lut()
 *
 * DESCRIPTION: Read darktarget GRIB LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_grib_lut( dtGribLUT& gc_lut )
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_GRIB);
    if(! filepath.empty()) {
        if (!Hishdf(filepath.c_str())) {
            status = read_grib_bin( filepath, &gc_lut );
            if (status != DTDB_SUCCESS) {
                cerr << "DtLutNetcdf::read_grib_lut(): " <<
                        "Failure reading GRIB binary ancillary file." << endl;
                return status;
            }
        }
        else {
            status = read_grib_hdf( filepath, &gc_lut );
            if (status != DTDB_SUCCESS) {
                cerr << "DtLutNetcdf::read_grib_lut(): " <<
                        "Failure reading GRIB hdf4 ancillary file." << endl;
                return status;
            }
            string filepath = get_option(INPUT_OZONE);
            if(! filepath.empty()) {
                status = read_ozone( filepath, &gc_lut );
                if (status != DTDB_SUCCESS) {
                    cerr << "DtLutNetcdf::read_ozone(): " <<
                            "Failure reading hdf4 ozone ancillary file." << endl;
                    return status;
                }
            }
            else {
                cerr << "DtLutNetcdf::read_grib_lut(): " <<
                    "Failure obtaining GRIB ozone input file path." << endl;
                return DTDB_FAIL;
            }
        }
    }
    else {
        cerr << "DtLutNetcdf::read_grib_lut(): " <<
            "Failure obtaining GRIB input file path." << endl;
        return DTDB_FAIL;
    }

    return status;
}

/**************************************************************************
 * NAME: read_grib_bin()
 *
 * DESCRIPTION: Read GRIB binary file.
 *
 *************************************************************************/

int DtLutNetcdf::read_grib_bin( const string filepath, dtGribLUT* grib_lut )
{
	int status = DTDB_SUCCESS;

	GDAS* gdas = new GDAS;
	bool isFileBigEndian = false;
	std::ifstream fin(filepath.c_str(), std::ios::in | std::ios::binary);
	if(fin.is_open()) {
		fin.read( (char*) gdas->data, GRIB_ARRAY_SIZE);
	   bool good = fin.good();
	   fin.close();
	   if(!good) {
		  cerr << "DtLutNetcdf::read_grib_bin() Error reading LUT data " << INPUT_GRIB << endl;
		  return DTDB_FAIL;
	   }
	}
	else {
	   cerr << "DtLutNetcdf::read_grib_bin() Error opening LUT file " << INPUT_GRIB << endl;
	   return DTDB_FAIL;
	}
	if ( isPlatformLittleEndian() && isFileBigEndian) {
		for ( int i=0; i<DATA_BINS; i++) {
			for ( int j=0; j<LAT_BINS; j++) {
				for ( int k=0; k<LON_BINS; k++) {
					byteSwap( gdas->data[i][j][k] );
				}
			}
		}
	}
	memcpy(grib_lut->pwat, gdas->data[50], GRIB_ROW_SIZE);
	memcpy(grib_lut->ugrd, gdas->data[51], GRIB_ROW_SIZE);
	memcpy(grib_lut->vgrd, gdas->data[52], GRIB_ROW_SIZE);
	memcpy(grib_lut->ozone, gdas->data[54], GRIB_ROW_SIZE);

	delete gdas;
    return status;
}

/**************************************************************************
 * NAME: read_grib_hdf()
 *
 * DESCRIPTION: Read GRIB HDF file.
 *
 *************************************************************************/

int DtLutNetcdf::read_grib_hdf( const string filepath, dtGribLUT* grib_lut )
{
	int status = DTDB_SUCCESS;

    int sdfid = SDstart(filepath.c_str(), DFACC_RDONLY);
    int vfid = Hopen(filepath.c_str(), DFACC_RDONLY, 0);
    Vstart(vfid);
    int vgroup_ref = Vfind(vfid, "Geophysical Data");
    int vgid = Vattach(vfid, vgroup_ref, "r");
	int nrefs = Vntagrefs(vgid);
    int tag[32];
    int ref[32];
	Vgettagrefs(vgid, tag, ref, nrefs);
	for (int iRef=0; iRef<nrefs; iRef++) {
		if ((tag[iRef]!=DFTAG_NDG) && (tag[iRef]!=DFTAG_SD)) {
			continue;
		}
		int index = SDreftoindex(sdfid, ref[iRef]);
		int sdsid = SDselect(sdfid, index);
		int rank = 0;
		int numtype = 0;
		int numatt = 0;
		int dims[2] = {0,0};
	    int start[2] = {0,0};
		char name[32];
		SDgetinfo(sdsid, name, &rank, dims, &numtype, &numatt);
		string sname = string(name);
		if (sname == "z_wind") {
			status = SDreaddata(sdsid, start, (int*)NULL, dims, (void*) grib_lut->ugrd);
			if (status != DTDB_SUCCESS) {
				cerr << "DtLutNetcdf::read_grib_hdf(): " <<
						"Failure reading hdf4 z_wind ancillary data." << endl;
				return status;
			}
		}
		if (sname == "m_wind") {
			status = SDreaddata(sdsid, start, (int*)NULL, dims, (void*) grib_lut->vgrd);
			if (status != DTDB_SUCCESS) {
				cerr << "DtLutNetcdf::read_grib_hdf(): " <<
						"Failure reading hdf4 m_wind ancillary data." << endl;
				return status;
			}
		}
		if (sname == "p_water") {
			status = SDreaddata(sdsid, start, (int*)NULL, dims, (void*) grib_lut->pwat);
			if (status != DTDB_SUCCESS) {
				cerr << "DtLutNetcdf::read_grib_hdf(): " <<
						"Failure reading hdf4 p_water ancillary data." << endl;
				return status;
			}
		}
	}
	Vdetach(vgid);
	Vend (vfid);
	SDend(sdfid);
	Hclose (vfid);

	return status;
}

/**************************************************************************
 * NAME: read_ozone()
 *
 * DESCRIPTION: Read ozone HDF file.
 *
 *************************************************************************/

int DtLutNetcdf::read_ozone( const string filepath, dtGribLUT* grib_lut )
{
	int status = DTDB_SUCCESS;

    int sdfid = SDstart(filepath.c_str(), DFACC_RDONLY);
    int vfid = Hopen(filepath.c_str(), DFACC_RDONLY, 0);
    status = Vstart(vfid);
    int vgroup_ref = Vfind(vfid, "Geophysical Data");
    int vgid = Vattach(vfid, vgroup_ref, "r");
	int nrefs = Vntagrefs(vgid);
    int tag[32];
    int ref[32];
	Vgettagrefs(vgid, tag, ref, nrefs);
	for (int iRef=0; iRef<nrefs; iRef++) {
		if ((tag[iRef]!=DFTAG_NDG) && (tag[iRef]!=DFTAG_SD)) {
			continue;
		}
		int index = SDreftoindex(sdfid, ref[iRef]);
		int sdsid = SDselect(sdfid, index);
		int rank = 0;
		int numtype = 0;
		int numatt = 0;
		int dims[2] = {0,0};
	    int start[2] = {0,0};
		char name[32];
		SDgetinfo(sdsid, name, &rank, dims, &numtype, &numatt);
		string sname = string(name);
		if (sname == "ozone") {
			status = SDreaddata(sdsid, start, (int*)NULL, dims, (void*) grib_lut->ozone);
			if (status != DTDB_SUCCESS) {
				cerr << "DtLutNetcdf::read_ozone(): " <<
						"Failure reading hdf4 ozone ancillary data." << endl;
				return status;
			}
		}
	}
	Vdetach(vgid);
	Vend (vfid);
	SDend(sdfid);
	Hclose (vfid);

	return status;
}

/**************************************************************************
 * NAME: read_gas_correction_file()
 *
 * DESCRIPTION: Read darktarget gas correction LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_gas_correction_file( dtGasCorrectionLUT* gc_lut )
{
	int status = DTDB_SUCCESS;

	string filePath = get_option(INPUT_GAS_CORRECTION);
	if(filePath.empty()) {
		cerr << "DtGranule::read_gas_correction_file() Invalid path for gas correction file." << endl;
		return DTDB_FAIL;
	}
   ifstream fin(filePath.c_str());
   if(!fin) {
       std::cout << "DtGranule::read_gas_correction_file() Error opening gas correction file." << endl;
       return DTDB_FAIL;
   }
   string   line;
   getline(fin, line);
   int i = 0;
   while (getline(fin, line)) {
       stringstream ss(line);
       ss >> gc_lut->MBAND[i] >> gc_lut->VBAND[i] >> gc_lut->WAVE[i] >> gc_lut->MOL[i] >>
	   gc_lut->OPT_O3_CLIM[i] >> gc_lut->OPT_H2O_CLIM[i] >>
	   gc_lut->OPT_CO2_CLIM[i] >> gc_lut->O3_COEF[i][0] >> gc_lut->O3_COEF[i][1] >>
	   gc_lut->H2O_COEF[i][0] >> gc_lut->H2O_COEF[i][1] >>  gc_lut->H2O_COEF[i][2];
	   i++;
   }

   return status;
}

/**************************************************************************
 * NAME: read_land_aerosol_file()
 *
 * DESCRIPTION: Read Land Aerosol LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_land_aerosol_file( const string groupname, int wnum,
										  dtLandAerosolLUT* la_lut )
{
	int status = DTDB_SUCCESS;

	string filePath = get_option(groupname);

	if(filePath.empty()) {
		cerr << "DtLutNetcdf::read_aerosol_file() Invalid path for aerosol LUT file." << endl;
		return DTDB_FAIL;
	}
	ifstream fin(filePath.c_str());
	if(!fin) {
       std::cout << "DtLutNetcdf::read_aerosol_file() Error opening aerosol LUT file." << endl;
       return DTDB_FAIL;
	}
   string   line;
   string   dummy;
   for (int iTab=0; iTab<NLTABLE; iTab++) {
//  Reads observation zenith angles(the),and observation azimuth angles
//  from look-up tables
		getline(fin, line);
		stringstream ss(line);
		ss >> dummy;
		for (int i=0; i<NLTHE; i++) {
		  ss >> la_lut->THE_NL[i];
		}
		getline(fin, line);
		ss.clear();
		ss.str(line);
		ss >> dummy;
		for (int i=0; i<NLPHI; i++) {
		  ss >> la_lut->PHI_NL[i];
		}
//  Reads wavelength(wav),optical thickness(opth),solar zenith angle(thet0),
//  reflectance ofatmosphere(sbarw) from look-up tables.
		for (int iTau=0; iTau<NLTAU; iTau++) {
			getline(fin, line);
			getline(fin, line);
			ss.clear();
			ss.str(line);
			ss >> dummy >> la_lut->SSA_NL0[iTab][wnum][iTau];
			ss >> dummy >> la_lut->QEXT_NL0[iTab][wnum][iTau];
			ss >> dummy >> la_lut->BEXT_NL0[iTab][wnum][iTau];
			ss >> dummy >> la_lut->VEXT_NL0[iTab][wnum][iTau];
			ss >> dummy >> la_lut->MEXT_NL0[iTab][wnum][iTau];
			ss >> dummy >> la_lut->MASSCOEF_NL0[iTab][wnum][iTau];
			getline(fin, line);
			ss.clear();
			ss.str(line);
			ss >> dummy >> la_lut->WAV_NL[wnum];
			ss >> dummy >> la_lut->OPTH_NL0[iTab][wnum][iTau];
			ss >> dummy >> la_lut->ROD[wnum];
			ss >> dummy >> la_lut->GOD[wnum];
			for (int iThet0=0; iThet0<NLTHET0; iThet0++) {
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> la_lut->THET0_NL[iThet0];
				ss >> dummy >> la_lut->MU0_NL[iThet0];
				ss >> dummy >> la_lut->SBAR_NL0[iTab][wnum][iTau][iThet0];
				ss >> dummy >> la_lut->Fd_NL0[iTab][wnum][iTau][iThet0];
//  Reads transmission as a function of observation zenith angle,
//  optical thickness
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy;
				for (int iThe=0; iThe<NLTHE; iThe++) {
					ss >> la_lut->T_NL0[iTab][wnum][iTau][iThet0][iThe];
				}
				getline(fin, line);
				for (int iThe=0; iThe<NLTHE; iThe++) {
					getline(fin, line);
					ss.clear();
					ss.str(line);
//  Reads atmospheric radiance (int) as a function of solar zenith angle,
//  optical thickness, height, observation zenith angle and azimuth angle
//  from look-up table.
					for (int iPhi=0; iPhi<NLPHI; iPhi++) {
						ss >> la_lut->INT_NL0[iTab][wnum][iTau][iThet0][iThe][iPhi];
					}
				}
			}
		}
//  Set extinction parameters for "AOD = 0.0" case
		la_lut->QEXT_NL0[iTab][wnum][0] = la_lut->QEXT_NL0[iTab][wnum][1];
		la_lut->BEXT_NL0[iTab][wnum][0] = la_lut->BEXT_NL0[iTab][wnum][1];
		la_lut->VEXT_NL0[iTab][wnum][0] = la_lut->VEXT_NL0[iTab][wnum][1];
		la_lut->MEXT_NL0[iTab][wnum][0] = la_lut->MEXT_NL0[iTab][wnum][1];
		la_lut->MASSCOEF_NL0[iTab][wnum][0] = la_lut->MASSCOEF_NL0[iTab][wnum][1];
   	}
	for (int iTab=0; iTab<NLTABLE; iTab++) {
		for (int iTau=0; iTau<NLTAU; iTau++) {
			for (int iWav=0; iWav<NLUTWAV; iWav++) {
				la_lut->EXTNORM_NL0[iTab][iWav][iTau] =
					 la_lut->QEXT_NL0[iTab][iWav][iTau] /
					 la_lut->QEXT_NL0[iTab][IW550][iTau];
			}
		}
	}

   return status;
}

/**************************************************************************
 * NAME: read_land_aerosol_map()
 *
 * DESCRIPTION: Subroutine AEROSOL_MAP reads
 *               the aerosol map (1 degree resolution)
 *               to determine expected aerosol type at
 *               a given location and season
 *
 *************************************************************************/

int DtLutNetcdf::read_land_aerosol_map( const string groupname,
										 dtLandAerosolLUT* la_lut )
{
	int status = DTDB_SUCCESS;

	string filePath = get_option(groupname);
	if(filePath.empty()) {
	cerr << "DtLutNetcdf::read_aerosol_map() Invalid path for aerosol LUT file." << endl;
	return DTDB_FAIL;
	}
	ifstream fin(filePath.c_str());
	if(!fin) {
	   std::cout << "DtLutNetcdf::read_aerosol_map() Error opening aerosol LUT file." << endl;
	   return DTDB_FAIL;
	}

	string   line;

	for (int iS=0; iS<NUM_SEASONS; iS++) {
	   getline(fin, line);
	   for (int iLat=0; iLat<NUM_LATS; iLat++) {
		   for (int iLon=0; iLon<NUM_LONS-1; iLon++)  {
			   getline(fin, line, ',');
			   stringstream ss(line);
			   ss >> la_lut->AEROSOL_ALL[iS][iLat][iLon];
		   }
           getline(fin, line);
           stringstream ss(line);
           ss >> la_lut->AEROSOL_ALL[iS][iLat][NUM_LONS-1];
	   }
	}

	return status;
}

/**************************************************************************
 * NAME: read_smalll_aerosol_file()
 *
 * DESCRIPTION: Read Small Ocean Aerosol LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_ocean_small_aerosol_file(const string groupname, int iLut,
												dtOceanAerosolLUT* oa_lut )
{
	int status = DTDB_SUCCESS;

	string filePath = get_option(groupname);
	if(filePath.empty()) {
	cerr << "DtProcessOcean::read_small_aerosol_file() Invalid path for aerosol LUT file." << endl;
	return DTDB_FAIL;
	}
    ifstream fin(filePath.c_str());
    if(!fin) {
       std::cout << "DtProcessOcean::read_small_aerosol_file() Error opening aerosol LUT file." << endl;
       return DTDB_FAIL;
    }
    string   line;
    string 	 aline;
    string   dummy;
    float	 TR, QSCT;
    for (int iWav=0; iWav<NWAV; iWav++) {
		getline(fin, line);
		getline(fin, line);
		stringstream ss(line);
		ss >> dummy >> oa_lut->WAVE[iWav];
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		ss.clear();
		ss.str(line);
		ss >> dummy >> TR;
		ss >> dummy >> oa_lut->TAUAS[iWav][0][0];
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		for (int il=0; il<2; il++) {
			getline(fin, line);
			aline += line;
		}
		ss.clear();
		ss.str(aline);
		ss >> dummy;
		for (int iTho=0; iTho<NTH0; iTho++) {
		   ss >> oa_lut->THET0[iTho];
		}
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		getline(fin, line);
		aline.clear();
		for (int il=0; il<2; il++) {
			getline(fin, line);
			aline += line;
		}
		ss.clear();
		ss.str(aline);
		ss >> dummy;
		for (int iTho=0; iTho<NTH0; iTho++) {
		   ss >> oa_lut->ALBEDO_R_RAY[iWav][iTho];
		}
		aline.clear();
		for (int il=0; il<2; il++) {
			getline(fin, line);
			aline += line;
		}
		ss.clear();
		ss.str(aline);
		ss >> dummy;
		for (int iTho=0; iTho<NTH0; iTho++) {
		   ss >> oa_lut->ALBEDO_T_RAY[iWav][iTho];
		}
		getline(fin, line);
		for (int iTho=0; iTho<NTH0; iTho++) {
			getline(fin, line);
			aline.clear();
			for (int il=0; il<3; il++) {
				getline(fin, line);
				aline += line;
			}
			ss.clear();
			ss.str(aline);
			ss >> dummy >> dummy;
			for (int iTh=0; iTh<NTHET; iTh++) {
			   ss >> oa_lut->THET[iTh];
			}
			for (int iPhi=0; iPhi<NPHI; iPhi++) {
				aline.clear();
				for (int il=0; il<3; il++) {
					getline(fin, line);
					aline += line;
				}
				ss.clear();
				ss.str(aline);
				ss >> oa_lut->JPHI[iPhi];
				for (int iTh=0; iTh<NTHET; iTh++) {
				   ss >> oa_lut->REF_RAYALL[iWav][iLut][iTho][iTh][iPhi];
				}
			}
		}
    }
	for (int iCase=0; iCase<NUMCASES; iCase++) {
		for (int iWav=0; iWav<NWAV; iWav++) {
			for (int iTau=1; iTau<NAOT; iTau++) {
				getline(fin, line);
				getline(fin, line);
				stringstream ss(line);
				ss >> dummy >> oa_lut->WAVE[iWav];
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> oa_lut->RGSS[iCase] >> dummy >> oa_lut->SIGMAS[iCase];
				getline(fin, line);
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> oa_lut->EFFRADSMALL[iCase];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> dummy >> oa_lut->MOMENTSSMALL[iCase][0][iLut] >>
				      dummy >> dummy >> dummy >> oa_lut->MOMENTSSMALL[iCase][1][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> dummy >> oa_lut->MOMENTSSMALL[iCase][2][iLut] >>
				      dummy >> dummy >> dummy >> oa_lut->MOMENTSSMALL[iCase][3][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> dummy >> oa_lut->ALBEDOSMALL[iWav][iCase][iLut] >>
				      dummy >> dummy >> dummy >> oa_lut->ASSYMSMALL[iWav][iCase][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> oa_lut->CCNSMALL[iCase][iLut] >>
				      dummy >> dummy >> dummy >> oa_lut->BACKSCTTSMALL[iWav][iCase][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> QSCT >>
				dummy >> dummy >> oa_lut->EXTSMALL[iWav][iCase][iLut];
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> TR ;
				ss >> dummy >> oa_lut->TAUAS[iWav][iCase][iTau];
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				aline.clear();
				for (int il=0; il<2; il++) {
					getline(fin, line);
					aline += line;
				}
				ss.clear();
				ss.str(aline);
				ss >> dummy;
				for (int iTho=0; iTho<NTH0; iTho++) {
				   ss >> oa_lut->THET0[iTho];
				}
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				aline.clear();
				for (int il=0; il<2; il++) {
					getline(fin, line);
					aline += line;
				}
				ss.clear();
				ss.str(aline);
				ss >> dummy;
				for (int iTho=0; iTho<NTH0; iTho++) {
				   ss >> oa_lut->ALBEDO_R_SMALL[iWav][iCase][iTau][iLut][iTho];
				}
				aline.clear();
				for (int il=0; il<2; il++) {
					getline(fin, line);
					aline += line;
				}
				ss.clear();
				ss.str(aline);
				ss >> dummy;
				for (int iTho=0; iTho<NTH0; iTho++) {
				   ss >> oa_lut->ALBEDO_T_SMALL[iWav][iCase][iTau][iLut][iTho];
				}
				getline(fin, line);
				for (int iTho=0; iTho<NTH0; iTho++) {
					getline(fin, line);
					aline.clear();
					for (int il=0; il<3; il++) {
						getline(fin, line);
						aline += line;
					}
					ss.clear();
					ss.str(aline);
					ss >> dummy >> dummy;
					for (int iTh=0; iTh<NTHET; iTh++) {
					   ss >> oa_lut->THET[iTh];
					}
					for (int iPhi=0; iPhi<NPHI; iPhi++) {
						aline.clear();
						for (int il=0; il<3; il++) {
							getline(fin, line);
							aline += line;
						}
						ss.clear();
						ss.str(aline);
						ss >> oa_lut->JPHI[iPhi];
						for (int iTh=0; iTh<NTHET; iTh++) {
						   ss >> oa_lut->AINTS[iWav][iCase][iTau][iLut][iTho][iTh][iPhi];
						}
					}
				}
			}
// Fill the array for albedo and transmission for all cases tau=0
			oa_lut->TAUAS[iWav][iCase][0]=oa_lut->TAUAS[iWav][0][0];
			for (int iTho=0; iTho<NTH0; iTho++) {
				oa_lut->ALBEDO_R_SMALL[iWav][iCase][0][iLut][iTho] = oa_lut->ALBEDO_R_RAY[iWav][iTho];
				oa_lut->ALBEDO_T_SMALL[iWav][iCase][0][iLut][iTho] = oa_lut->ALBEDO_T_RAY[iWav][iTho];
			}
		}
	}
	for (int iPhi=0; iPhi<NPHI; iPhi++) {
		oa_lut->PHC[iPhi] = (float) oa_lut->JPHI[iPhi];
	}

   return status;
}

/**************************************************************************
 * NAME: read_big_aerosol_file()
 *
 * DESCRIPTION: Read Big Ocean Aerosol LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_ocean_big_aerosol_file(const string groupname, int iLut,
											  dtOceanAerosolLUT* oa_lut )
{
	int status = DTDB_SUCCESS;

	string filePath = get_option(groupname);
	if(filePath.empty()) {
		cerr << "DtProcessOcean::read_ocean_big_aerosol_file() Invalid path for aerosol LUT file." << endl;
		return DTDB_FAIL;
	}
	ifstream fin(filePath.c_str());
	if(!fin) {
	   std::cout << "DtProcessOcean::big_ocean_small_aerosol_file() Error opening aerosol LUT file." << endl;
	   return DTDB_FAIL;
	}
	string   line;
	string 	 aline;
	string   dummy;
    float	 TR, QSCT;
	for (int iCase=0; iCase<NUMCASEB; iCase++) {
		for (int iWav=0; iWav<NWAV; iWav++) {
			for (int iTau=1; iTau<NAOT; iTau++) {
				getline(fin, line);
				getline(fin, line);
				stringstream ss(line);
				ss >> dummy >> oa_lut->WAVE[iWav];
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> oa_lut->RGSB[iCase]
				   >> dummy >> oa_lut->SIGMAB[iCase];
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> oa_lut->EFFRADBIG[iCase];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> dummy >> oa_lut->MOMENTSBIG[iCase][0][iLut] >>
				      dummy >> dummy >> dummy >> oa_lut->MOMENTSBIG[iCase][1][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> dummy >> oa_lut->MOMENTSBIG[iCase][2][iLut] >>
				      dummy >> dummy >> dummy >> oa_lut->MOMENTSBIG[iCase][3][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> dummy >> oa_lut->ALBEDOBIG[iWav][iCase][iLut] >>
				      dummy >> dummy >> dummy >> oa_lut->ASSYMBIG[iWav][iCase][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> dummy >> dummy >>
					  dummy >> oa_lut->BACKSCTTBIG[iWav][iCase][iLut];
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> dummy >> QSCT >>
				      dummy >> dummy >> oa_lut->EXTBIG[iWav][iCase][iLut];
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				ss.clear();
				ss.str(line);
				ss >> dummy >> TR ;
				ss >> dummy >> oa_lut->TAUAB[iWav][iCase][iTau];
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				aline.clear();
				for (int il=0; il<2; il++) {
					getline(fin, line);
					aline += line;
				}
				ss.clear();
				ss.str(aline);
				ss >> dummy;
				for (int iTho=0; iTho<NTH0; iTho++) {
				   ss >> oa_lut->THET0[iTho];
				}
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);
				aline.clear();
				for (int il=0; il<2; il++) {
					getline(fin, line);
					aline += line;
				}
				ss.clear();
				ss.str(aline);
				ss >> dummy;
				for (int iTho=0; iTho<NTH0; iTho++) {
				   ss >> oa_lut->ALBEDO_R_BIG[iWav][iCase][iTau][iLut][iTho];
				}
				aline.clear();
				for (int il=0; il<2; il++) {
					getline(fin, line);
					aline += line;
				}
				ss.clear();
				ss.str(aline);
				ss >> dummy;
				for (int iTho=0; iTho<NTH0; iTho++) {
				   ss >> oa_lut->ALBEDO_T_BIG[iWav][iCase][iTau][iLut][iTho];
				}
				getline(fin, line);
				for (int iTho=0; iTho<NTH0; iTho++) {
					getline(fin, line);
					aline.clear();
					for (int il=0; il<3; il++) {
						getline(fin, line);
						aline += line;
					}
					ss.clear();
					ss.str(aline);
					ss >> dummy >> dummy;
					for (int iTh=0; iTh<NTHET; iTh++) {
					   ss >> oa_lut->THET[iTh];
					}
					for (int iPhi=0; iPhi<NPHI; iPhi++) {
						aline.clear();
						for (int il=0; il<3; il++) {
							getline(fin, line);
							aline += line;
						}
						ss.clear();
						ss.str(aline);
						ss >> oa_lut->JPHI[iPhi];
						for (int iTh=0; iTh<NTHET; iTh++) {
						   ss >> oa_lut->AINTB[iWav][iCase][iTau][iLut][iTho][iTh][iPhi];
						}
					}
				}
			}
		}
	}

   return status;
}


/**************************************************************************
 * NAME: read_transm_h2o_file()
 *
 * DESCRIPTION: Reads TRANSM_H2O ascii file
 *
 *************************************************************************/

int DtLutNetcdf::read_transm_h2o_file( const string groupname, const int num,
                                         dtWaterVaporLUT* wv_lut )
{
    int status = DTDB_SUCCESS;

    string filePath = get_option(groupname);
    if(filePath.empty()) {
    cerr << "DtLutNetcdf::read_transm_h2o_file() Invalid path for LUT file." << endl;
    return DTDB_FAIL;
    }
    ifstream fin(filePath.c_str());
    if(!fin) {
       std::cout << "DtLutNetcdf::read_transm_h2o_file() Error opening LUT file." << endl;
       return DTDB_FAIL;
    }
    string   line;
    for (int iR=0; iR<TRANSM_H2O_ROWS; iR++) {
       getline(fin, line);
       stringstream ss(line);
       for (int iV=0; iV<TRANSM_H2O_VALS; iV++) {
           ss >> wv_lut->TRANSM_H20[num-1][iR][iV];
       }
    }

    return status;
}

/**************************************************************************
 * NAME: read_ch19_to_ch2_ratio_file()
 *
 * DESCRIPTION: Reads ch19-to-ch2 ratio ascii file
 *
 *************************************************************************/

int DtLutNetcdf::read_ch19_to_ch2_ratio_file( const string groupname,
                                         dtWaterVaporLUT* wv_lut )
{
    int status = DTDB_SUCCESS;

    string filePath = get_option(groupname);
    if(filePath.empty()) {
        std::cout <<
    "DtLutNetcdf::read_ch19_to_ch2_ratio_file() Invalid path for LUT file." << endl;
    return DTDB_FAIL;
    }
    ifstream fin(filePath.c_str());
    if(!fin) {
       std::cout <<
       "DtLutNetcdf::read_ch19_to_ch2_ratio_file() Error opening LUT file." << endl;
       return DTDB_FAIL;
    }
    string   line;
    for (int iR=0; iR<REFL_CH2_ROWS; iR++) {
       getline(fin, line);
       stringstream ss(line);
       for (int iV=0; iV<REFL_CH2_VALS; iV++) {
           ss >> wv_lut->RATIO_CH19_TO_CH2[iR][iV];
       }
    }

    return status;
}

/**************************************************************************
 * NAME: read_ch2_reflectance_file()
 *
 * DESCRIPTION: Reads ch2 reflectance ascii file
 *
 *************************************************************************/

int DtLutNetcdf::read_ch2_reflectance_file( const string groupname,
                                         dtWaterVaporLUT* wv_lut )
{
    int status = DTDB_SUCCESS;

    string filePath = get_option(groupname);
    if(filePath.empty()) {
        std::cout <<
    "DtLutNetcdf::read_ch2_file() Invalid path for LUT file." << endl;
    return DTDB_FAIL;
    }
    ifstream fin(filePath.c_str());
    if(!fin) {
       std::cout <<
       "DtLutNetcdf::read_ch2_file() Error opening LUT file." << endl;
       return DTDB_FAIL;
    }
    string   line;
    for (int iR=0; iR<REFL_CH2_ROWS; iR++) {
       getline(fin, line);
       stringstream ss(line);
       for (int iV=0; iV<REFL_CH2_VALS; iV++) {
           ss >> wv_lut->REFL_CH2[iR][iV];
       }
    }

    return status;
}

/**************************************************************************
 * NAME: read_weight_table_file()
 *
 * DESCRIPTION: Reads weight table ascii file
 *
 *************************************************************************/

int DtLutNetcdf::read_weight_table_file( const string groupname,
                                         dtWaterVaporLUT* wv_lut )
{
    int status = DTDB_SUCCESS;

    string filePath = get_option(groupname);
    if(filePath.empty()) {
        std::cout <<
    "DtLutNetcdf::read_weight_table_file() Invalid path for LUT file." << endl;
    return DTDB_FAIL;
    }
    ifstream fin(filePath.c_str());
    if(!fin) {
       std::cout <<
       "DtLutNetcdf::read_weight_table_file() Error opening LUT file." << endl;
       return DTDB_FAIL;
    }
    string   line;
    getline(fin, line);
    for (int iR=0; iR<TRANSM_H2O_ROWS; iR++) {
       getline(fin, line);
       stringstream ss(line);
       for (int iV=0; iV<WEIGHT_VALS; iV++) {
           ss >> wv_lut->WEIGHTS[iR][iV];
       }
    }

    return status;
}

/**************************************************************************
 * NAME: write_gas_correction_lut()
 *
 * DESCRIPTION: Write gas correction LUT to NetCDF4 file.
 *
 *************************************************************************/

int DtLutNetcdf::write_gas_correction_lut( NcFile* nc_output,
											dtGasCorrectionLUT* gc_lut )
{
	NcGroup lut_grp = nc_output->addGroup( LUT_GAS_CORRECTION );

	num_gc_dt_bands_dim_ = lut_grp.addDim( "Dim_Bands", NUM_DT_BANDS );
	num_gc_O3_coef_dim_ = lut_grp.addDim( "Dim_O3_Coefficients", O3_COEFS );
	num_gc_H2O_coef_dim_ = lut_grp.addDim( "Dim_H2O_Coefficients", H2O_COEFS );

	vector<NcDim> O3_dims;
	O3_dims.push_back(num_gc_dt_bands_dim_);
	O3_dims.push_back(num_gc_O3_coef_dim_);

	vector<NcDim> H2O_dims;
	H2O_dims.push_back(num_gc_dt_bands_dim_);
	H2O_dims.push_back(num_gc_H2O_coef_dim_);

	NcVar var = lut_grp.addVar( "MBAND", ncInt, num_gc_dt_bands_dim_ );
	var.putVar( gc_lut->MBAND );

	var = lut_grp.addVar( "VBAND", ncInt, num_gc_dt_bands_dim_ );
	var.putVar( gc_lut->VBAND );

	var = lut_grp.addVar( "WAVE", ncFloat, num_gc_dt_bands_dim_ );
	var.putVar( gc_lut->WAVE );

	var = lut_grp.addVar( "MOL", ncFloat, num_gc_dt_bands_dim_ );
	var.putVar( gc_lut->MOL );

	var = lut_grp.addVar( "OPT_O3_CLIM", ncFloat, num_gc_dt_bands_dim_ );
	var.putVar( gc_lut->OPT_O3_CLIM );

	var = lut_grp.addVar( "OPT_H2O_CLIM", ncFloat, num_gc_dt_bands_dim_ );
	var.putVar( gc_lut->OPT_H2O_CLIM );

	var = lut_grp.addVar( "OPT_CO2_CLIM", ncFloat, num_gc_dt_bands_dim_ );
	var.putVar( gc_lut->OPT_CO2_CLIM );

	var = lut_grp.addVar( "O3_COEF", ncFloat, O3_dims );
	var.putVar( gc_lut->O3_COEF );

	var = lut_grp.addVar( "H2O_COEF", ncFloat, H2O_dims );
	var.putVar( gc_lut->H2O_COEF );

	return DTDB_SUCCESS;
}


/**************************************************************************
 * NAME: write_land_aerosol_lut()
 *
 * DESCRIPTION: Write land aerosol LUT to NetCDF4 file.
 *
 *************************************************************************/

int DtLutNetcdf::write_land_aerosol_lut( NcFile* nc_output,
										  dtLandAerosolLUT* la_lut )
{
	NcGroup lut_grp = nc_output->addGroup( LUT_LAND_AEROSOL );

	num_land_lats_dim_ = lut_grp.addDim( "Dim_Latitude", NUM_LATS );
	num_land_lons_dim_ = lut_grp.addDim( "Dim_Longitude", NUM_LONS );
	num_land_phi_dim_ = lut_grp.addDim( "Dim_Phi", NLPHI );
	num_land_the_dim_ = lut_grp.addDim( "Dim_Theta", NLTHE );
	num_land_thet0_dim_ = lut_grp.addDim( "Dim_SZA", NLTHET0 );
	num_land_tau_dim_ = lut_grp.addDim( "Dim_Tau", NLTAU );
	num_land_wav_dim_ = lut_grp.addDim( "Dim_Bands", NLUTWAV );
	num_land_table_dim_ = lut_grp.addDim( "Dim_Tables", NLTABLE );
	num_land_size_dim_ = lut_grp.addDim( "Dim_Sizes", NLSIZE );
	num_land_season_dim_ = lut_grp.addDim( "Dim_Seasons", NUM_SEASONS );

	vector<NcDim> map_dims;
	map_dims.push_back(num_land_season_dim_);
	map_dims.push_back(num_land_lats_dim_);
	map_dims.push_back(num_land_lons_dim_);

	vector<NcDim> tab3_dims;
	tab3_dims.push_back(num_land_table_dim_);
	tab3_dims.push_back(num_land_wav_dim_);
	tab3_dims.push_back(num_land_tau_dim_);

	vector<NcDim> tab4_dims;
	tab4_dims.push_back(num_land_table_dim_);
	tab4_dims.push_back(num_land_wav_dim_);
	tab4_dims.push_back(num_land_tau_dim_);
	tab4_dims.push_back(num_land_thet0_dim_);

	vector<NcDim> tab5_dims;
	tab5_dims.push_back(num_land_table_dim_);
	tab5_dims.push_back(num_land_wav_dim_);
	tab5_dims.push_back(num_land_tau_dim_);
	tab5_dims.push_back(num_land_thet0_dim_);
	tab5_dims.push_back(num_land_the_dim_);

	vector<NcDim> tab6_dims;
	tab6_dims.push_back(num_land_table_dim_);
	tab6_dims.push_back(num_land_wav_dim_);
	tab6_dims.push_back(num_land_tau_dim_);
	tab6_dims.push_back(num_land_thet0_dim_);
	tab6_dims.push_back(num_land_the_dim_);
	tab6_dims.push_back(num_land_phi_dim_);

	vector<NcDim> siz3_dims;
	siz3_dims.push_back(num_land_size_dim_);
	siz3_dims.push_back(num_land_wav_dim_);
	siz3_dims.push_back(num_land_tau_dim_);

	vector<NcDim> siz4_dims;
	siz4_dims.push_back(num_land_size_dim_);
	siz4_dims.push_back(num_land_wav_dim_);
	siz4_dims.push_back(num_land_tau_dim_);
	siz4_dims.push_back(num_land_thet0_dim_);

	vector<NcDim> siz5_dims;
	siz5_dims.push_back(num_land_size_dim_);
	siz5_dims.push_back(num_land_wav_dim_);
	siz5_dims.push_back(num_land_tau_dim_);
	siz5_dims.push_back(num_land_thet0_dim_);
	siz5_dims.push_back(num_land_the_dim_);

	vector<NcDim> siz6_dims;
	siz6_dims.push_back(num_land_size_dim_);
	siz6_dims.push_back(num_land_wav_dim_);
	siz6_dims.push_back(num_land_tau_dim_);
	siz6_dims.push_back(num_land_thet0_dim_);
	siz6_dims.push_back(num_land_the_dim_);
	siz6_dims.push_back(num_land_phi_dim_);

	NcVar var = lut_grp.addVar( "AEROSOL_ALL", ncInt, map_dims );
	var.putVar( la_lut->AEROSOL_ALL );

	var = lut_grp.addVar( "PHI_NL", ncFloat, num_land_phi_dim_ );
	var.putVar( la_lut->PHI_NL );

	var = lut_grp.addVar( "THE_NL", ncFloat, num_land_the_dim_ );
	var.putVar( la_lut->THE_NL );

	var = lut_grp.addVar( "THET0_NL", ncFloat, num_land_thet0_dim_ );
	var.putVar( la_lut->THET0_NL );

	var = lut_grp.addVar( "MU0_NL", ncFloat, num_land_thet0_dim_ );
	var.putVar( la_lut->MU0_NL );

	var = lut_grp.addVar( "WAV_NL", ncFloat, num_land_wav_dim_ );
	var.putVar( la_lut->WAV_NL );

	var = lut_grp.addVar( "OPTH_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->OPTH_NL0 );

	var = lut_grp.addVar( "MASSCOEF_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->MASSCOEF_NL0 );

	var = lut_grp.addVar( "EXTNORM_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->EXTNORM_NL0 );

	var = lut_grp.addVar( "SSA_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->SSA_NL0 );

	var = lut_grp.addVar( "QEXT_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->QEXT_NL0 );

	var = lut_grp.addVar( "BEXT_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->BEXT_NL0 );

	var = lut_grp.addVar( "VEXT_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->VEXT_NL0 );

	var = lut_grp.addVar( "MEXT_NL0", ncFloat, tab3_dims );
	var.putVar( la_lut->MEXT_NL0 );

	var = lut_grp.addVar( "SBAR_NL0", ncFloat, tab4_dims );
	var.putVar( la_lut->SBAR_NL0 );

	var = lut_grp.addVar( "INT_NL0", ncFloat, tab6_dims );
	var.putVar( la_lut->INT_NL0 );

	var = lut_grp.addVar( "Fd_NL0", ncFloat, tab4_dims );
	var.putVar( la_lut->Fd_NL0 );

	var = lut_grp.addVar( "T_NL0", ncFloat, tab5_dims );
	var.putVar( la_lut->T_NL0 );

	var = lut_grp.addVar( "OMEGA0", ncFloat, tab3_dims );
	var.putVar( la_lut->OMEGA0 );

	var = lut_grp.addVar( "ROD", ncFloat, num_land_wav_dim_ );
	var.putVar( la_lut->ROD );

	var = lut_grp.addVar( "GOD", ncFloat, num_land_wav_dim_ );
	var.putVar( la_lut->GOD );

	return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: write_ocean_aerosol_lut()
 *
 * DESCRIPTION: Write ocean aerosol LUT to NetCDF4 file.
 *
 *************************************************************************/

int DtLutNetcdf::write_ocean_aerosol_lut( NcFile* nc_output,
											dtOceanAerosolLUT* oa_lut )
{
	NcGroup lut_grp = nc_output->addGroup( LUT_OCEAN_AEROSOL );

	num_ocean_phi_dim_ = lut_grp.addDim( "Dim_Phi", NPHI );
	num_ocean_the_dim_ = lut_grp.addDim( "Dim_Theta", NTHET );
	num_ocean_thet0_dim_ = lut_grp.addDim( "Dim_SZA", NTH0 );
	num_ocean_tau_dim_ = lut_grp.addDim( "Dim_Tau", NAOT );
	num_ocean_wave_dim_ = lut_grp.addDim( "Dim_Bands", NWAV );
	num_ocean_cases_dim_ = lut_grp.addDim( "Dim_Small_Cases", NUM_CASES_SMALL );
	num_ocean_caseb_dim_ = lut_grp.addDim( "Dim_Large_Cases", NUM_CASES_BIG );
	num_ocean_wslut_dim_ = lut_grp.addDim( "Dim_Wind_LUT_bins", NUM_LUTS );
	num_ocean_moments_dim_ = lut_grp.addDim( "Dim_Moments", 4 );

	vector<NcDim> moments_dims;
    moments_dims.push_back(num_ocean_cases_dim_);
	moments_dims.push_back(num_ocean_moments_dim_);
	moments_dims.push_back(num_ocean_wslut_dim_);

	vector<NcDim> tab2s_dims;
	tab2s_dims.push_back(num_ocean_cases_dim_);
	tab2s_dims.push_back(num_ocean_wslut_dim_);

	vector<NcDim> tab3s_dims;
	tab3s_dims.push_back(num_ocean_wave_dim_);
	tab3s_dims.push_back(num_ocean_cases_dim_);
	tab3s_dims.push_back(num_ocean_wslut_dim_);

	vector<NcDim> tau3s_dims;
    tau3s_dims.push_back(num_ocean_wave_dim_);
    tau3s_dims.push_back(num_ocean_cases_dim_);
    tau3s_dims.push_back(num_ocean_tau_dim_);

	vector<NcDim> tab5s_dims;
    tab5s_dims.push_back(num_ocean_wave_dim_);
	tab5s_dims.push_back(num_ocean_cases_dim_);
	tab5s_dims.push_back(num_ocean_tau_dim_);
	tab5s_dims.push_back(num_ocean_wslut_dim_);
    tab5s_dims.push_back(num_ocean_thet0_dim_);

	vector<NcDim> tab7s_dims;
    tab7s_dims.push_back(num_ocean_wave_dim_);
	tab7s_dims.push_back(num_ocean_cases_dim_);
	tab7s_dims.push_back(num_ocean_tau_dim_);
    tab7s_dims.push_back(num_ocean_wslut_dim_);
	tab7s_dims.push_back(num_ocean_thet0_dim_);
	tab7s_dims.push_back(num_ocean_the_dim_);
	tab7s_dims.push_back(num_ocean_phi_dim_);

	NcVar var = lut_grp.addVar( "JPHI", ncInt, num_ocean_phi_dim_ );
	var.putVar( oa_lut->JPHI );

	var = lut_grp.addVar( "PHC", ncFloat, num_ocean_phi_dim_ );
	var.putVar( oa_lut->PHC );

	var = lut_grp.addVar( "THET", ncFloat, num_ocean_the_dim_ );
	var.putVar( oa_lut->THET );

	var = lut_grp.addVar( "THET0", ncFloat, num_ocean_thet0_dim_ );
	var.putVar( oa_lut->THET0 );

	var = lut_grp.addVar( "WAVE", ncFloat, num_ocean_wave_dim_ );
	var.putVar( oa_lut->WAVE );

	var = lut_grp.addVar( "EXTSMALL", ncFloat, tab3s_dims );
	var.putVar( oa_lut->EXTSMALL );

	var = lut_grp.addVar( "RGSS", ncFloat, num_ocean_cases_dim_ );
	var.putVar( oa_lut->RGSS );

	var = lut_grp.addVar( "SIGMAS", ncFloat, num_ocean_cases_dim_ );
	var.putVar( oa_lut->SIGMAS );

	var = lut_grp.addVar( "MOMENTSSMALL", ncFloat, moments_dims );
	var.putVar( oa_lut->MOMENTSSMALL );

	var = lut_grp.addVar( "CCNSMALL", ncFloat, tab2s_dims );
	var.putVar( oa_lut->CCNSMALL );

	var = lut_grp.addVar( "BACKSCTTSMALL", ncFloat, tab3s_dims );
	var.putVar( oa_lut->BACKSCTTSMALL );

	var = lut_grp.addVar( "ASSYMSMALL", ncFloat, tab3s_dims );
	var.putVar( oa_lut->ASSYMSMALL );

	var = lut_grp.addVar( "ALBEDOSMALL", ncFloat, tab3s_dims );
	var.putVar( oa_lut->ALBEDOSMALL );

	var = lut_grp.addVar( "ALBEDO_R_SMALL", ncFloat, tab5s_dims );
	var.putVar( oa_lut->ALBEDO_R_SMALL );

	var = lut_grp.addVar( "ALBEDO_T_SMALL", ncFloat, tab5s_dims );
	var.putVar( oa_lut->ALBEDO_T_SMALL );

	var = lut_grp.addVar( "AINTS", ncFloat, tab7s_dims );
	var.putVar( oa_lut->AINTS );

	var = lut_grp.addVar( "TAUAS", ncFloat, tau3s_dims );
	var.putVar( oa_lut->TAUAS );

	var = lut_grp.addVar( "EFFRADSMALL", ncFloat, num_ocean_cases_dim_ );
	var.putVar( oa_lut->EFFRADSMALL );

	vector<NcDim> momentb_dims;
    momentb_dims.push_back(num_ocean_caseb_dim_);
	momentb_dims.push_back(num_ocean_moments_dim_);
	momentb_dims.push_back(num_ocean_wslut_dim_);

	vector<NcDim> tab2b_dims;
	tab2b_dims.push_back(num_ocean_caseb_dim_);
	tab2b_dims.push_back(num_ocean_wslut_dim_);

	vector<NcDim> tab3b_dims;
	tab3b_dims.push_back(num_ocean_wave_dim_);
	tab3b_dims.push_back(num_ocean_caseb_dim_);
	tab3b_dims.push_back(num_ocean_wslut_dim_);

	vector<NcDim> tau3b_dims;
    tau3b_dims.push_back(num_ocean_wave_dim_);
    tau3b_dims.push_back(num_ocean_caseb_dim_);
    tau3b_dims.push_back(num_ocean_tau_dim_);

	vector<NcDim> tab5b_dims;
    tab5b_dims.push_back(num_ocean_wave_dim_);
	tab5b_dims.push_back(num_ocean_caseb_dim_);
	tab5b_dims.push_back(num_ocean_tau_dim_);
	tab5b_dims.push_back(num_ocean_wslut_dim_);
    tab5b_dims.push_back(num_ocean_thet0_dim_);

	vector<NcDim> tab7b_dims;
    tab7b_dims.push_back(num_ocean_wave_dim_);
	tab7b_dims.push_back(num_ocean_caseb_dim_);
	tab7b_dims.push_back(num_ocean_tau_dim_);
    tab7b_dims.push_back(num_ocean_wslut_dim_);
	tab7b_dims.push_back(num_ocean_thet0_dim_);
	tab7b_dims.push_back(num_ocean_the_dim_);
	tab7b_dims.push_back(num_ocean_phi_dim_);

	vector<NcDim> alb2_dims;
    alb2_dims.push_back(num_ocean_wave_dim_);
	alb2_dims.push_back(num_ocean_thet0_dim_);

	vector<NcDim> ray5_dims;
	ray5_dims.push_back(num_ocean_wave_dim_);
    ray5_dims.push_back(num_ocean_wslut_dim_);
	ray5_dims.push_back(num_ocean_thet0_dim_);
	ray5_dims.push_back(num_ocean_the_dim_);
	ray5_dims.push_back(num_ocean_phi_dim_);

	var = lut_grp.addVar( "EXTBIG", ncFloat, tab3b_dims );
	var.putVar( oa_lut->EXTBIG );

	var = lut_grp.addVar( "RGSB", ncFloat, num_ocean_caseb_dim_ );
	var.putVar( oa_lut->RGSB );

	var = lut_grp.addVar( "SIGMAB", ncFloat, num_ocean_caseb_dim_ );
	var.putVar( oa_lut->SIGMAB );

	var = lut_grp.addVar( "MOMENTSBIG", ncFloat, momentb_dims );
	var.putVar( oa_lut->MOMENTSBIG );

	var = lut_grp.addVar( "BACKSCTTBIG", ncFloat, tab3b_dims );
	var.putVar( oa_lut->BACKSCTTBIG );

	var = lut_grp.addVar( "ASSYMBIG", ncFloat, tab3b_dims );
	var.putVar( oa_lut->ASSYMBIG );

	var = lut_grp.addVar( "ALBEDOBIG", ncFloat, tab3b_dims );
	var.putVar( oa_lut->ALBEDOBIG );

	var = lut_grp.addVar( "ALBEDO_R_BIG", ncFloat, tab5b_dims );
	var.putVar( oa_lut->ALBEDO_R_BIG );

	var = lut_grp.addVar( "ALBEDO_T_BIG", ncFloat, tab5b_dims );
	var.putVar( oa_lut->ALBEDO_T_BIG );

	var = lut_grp.addVar( "AINTB", ncFloat, tab7b_dims );
	var.putVar( oa_lut->AINTB );

	var = lut_grp.addVar( "TAUAB", ncFloat, tau3b_dims );
	var.putVar( oa_lut->TAUAB );

	var = lut_grp.addVar( "EFFRADBIG", ncFloat, num_ocean_caseb_dim_ );
	var.putVar( oa_lut->EFFRADBIG );

	var = lut_grp.addVar( "ALBEDO_R_RAY", ncFloat, alb2_dims );
	var.putVar( oa_lut->ALBEDO_R_RAY );

	var = lut_grp.addVar( "ALBEDO_T_RAY", ncFloat, alb2_dims );
	var.putVar( oa_lut->ALBEDO_T_RAY );

	var = lut_grp.addVar( "REF_RAYALL", ncFloat, ray5_dims );
	var.putVar( oa_lut->REF_RAYALL );

	return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: write_water_vapor_lut()
 *
 * DESCRIPTION: Write water vapor LUT to NetCDF4 file.
 *
 *************************************************************************/

int DtLutNetcdf::write_water_vapor_lut( NcFile* nc_output,
                                    dtWaterVaporLUT* wv_lut )
{
    NcGroup lut_grp = nc_output->addGroup( LUT_WATER_VAPOR );

    num_h2o_tables_dim_ = lut_grp.addDim( "Dim_Transm_H2O_tables", TRANSM_H2O_TABLES );
    num_h2o_rows_dim_ = lut_grp.addDim( "Dim_H2O_Rows", TRANSM_H2O_ROWS );
    num_h2o_vals_dim_ = lut_grp.addDim( "Dim_H2O_Vals", TRANSM_H2O_VALS );
//    num_weight_vals_dim_ = lut_grp.addDim( "Dim_Weight_Vals", WEIGHT_VALS );
//    num_ch2_rows_dim_ = lut_grp.addDim( "Dim_Refl_Ch2_Rows", REFL_CH2_ROWS );
//    num_ch2_vals_dim_ = lut_grp.addDim( "Dim_Refl_Ch2_Vals", REFL_CH2_VALS );

    vector<NcDim> H2O_dims;
    H2O_dims.push_back(num_h2o_tables_dim_);
    H2O_dims.push_back(num_h2o_rows_dim_);
    H2O_dims.push_back(num_h2o_vals_dim_);

//    vector<NcDim> weight_dims;
//    weight_dims.push_back(num_h2o_rows_dim_);
//    weight_dims.push_back(num_weight_vals_dim_);

//    vector<NcDim> ch2_dims;
//    ch2_dims.push_back(num_ch2_rows_dim_);
//    ch2_dims.push_back(num_ch2_vals_dim_);

    NcVar var = lut_grp.addVar( "TRANSM_H20", ncFloat, H2O_dims );
    var.putVar( wv_lut->TRANSM_H20 );
/*
    var = lut_grp.addVar( "WEIGHTS", ncFloat, weight_dims );
    var.putVar( wv_lut->WEIGHTS );

    var = lut_grp.addVar( "REFL_CH2", ncFloat, ch2_dims );
    var.putVar( wv_lut->REFL_CH2 );

    var = lut_grp.addVar( "RATIO_CH19_TO_CH2", ncFloat, ch2_dims );
    var.putVar( wv_lut->RATIO_CH19_TO_CH2 );
*/
    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: read_gas_correction_lut()
 *
 * DESCRIPTION: Read darktarget gas correction NetCDF4 LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_gas_correction_lut( dtGasCorrectionLUT &gc_lut )
{
	int status = DTDB_SUCCESS;

	std::string filepath = get_option(INPUT_NC4_LUT);
	if (filepath.empty()) {
		filepath = get_option(INPUT_DT_NC4_LUT);
	}
	NcFile* nc_input;
	try {
		nc_input = new NcFile(filepath, NcFile::read );
	}
	catch( NcException& e) {
		e.what();
		cerr << "DtLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
		return DTDB_FAIL;
	}
	NcGroup lut_grp = nc_input->getGroup( LUT_GAS_CORRECTION );

	NcVar var = lut_grp.getVar ( "MBAND" );
	var.getVar( gc_lut.MBAND );

	var = lut_grp.getVar ( "VBAND" );
	var.getVar( gc_lut.VBAND );

	var = lut_grp.getVar ( "WAVE" );
	var.getVar( gc_lut.WAVE );

	var = lut_grp.getVar ( "MOL" );
	var.getVar( gc_lut.MOL );

	var = lut_grp.getVar ( "OPT_O3_CLIM" );
	var.getVar( gc_lut.OPT_O3_CLIM );

	var = lut_grp.getVar ( "OPT_H2O_CLIM" );
	var.getVar( gc_lut.OPT_H2O_CLIM );

	var = lut_grp.getVar ( "OPT_CO2_CLIM" );
	var.getVar( gc_lut.OPT_CO2_CLIM );

	var = lut_grp.getVar ( "O3_COEF" );
	var.getVar( gc_lut.O3_COEF );

	var = lut_grp.getVar ( "H2O_COEF" );
	var.getVar( gc_lut.H2O_COEF );

	delete nc_input;

	return status;
}

/**************************************************************************
 * NAME: read_land_aerosol_lut()
 *
 * DESCRIPTION: Read Land Aerosol NetCDF4 LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_land_aerosol_lut( dtLandAerosolLUT  &la_lut )
{
	int status = DTDB_SUCCESS;

	std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DT_NC4_LUT);
	}
	NcFile* nc_input;
	try {
		nc_input = new NcFile(filepath, NcFile::read );
	}
	catch( NcException& e) {
		e.what();
		cerr << "DtLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
		return DTDB_FAIL;
	}
	NcGroup lut_grp = nc_input->getGroup( LUT_LAND_AEROSOL );

	NcVar var = lut_grp.getVar( "AEROSOL_ALL");
	var.getVar( la_lut.AEROSOL_ALL );

	var = lut_grp.getVar( "PHI_NL" );
	var.getVar( la_lut.PHI_NL );

	var = lut_grp.getVar( "THE_NL" );
	var.getVar( la_lut.THE_NL );

	var = lut_grp.getVar( "THET0_NL" );
	var.getVar( la_lut.THET0_NL );

	var = lut_grp.getVar( "MU0_NL" );
	var.getVar( la_lut.MU0_NL );

	var = lut_grp.getVar( "WAV_NL" );
	var.getVar( la_lut.WAV_NL );

	var = lut_grp.getVar( "OPTH_NL0" );
	var.getVar( la_lut.OPTH_NL0 );

	var = lut_grp.getVar( "MASSCOEF_NL0" );
	var.getVar( la_lut.MASSCOEF_NL0 );

	var = lut_grp.getVar( "EXTNORM_NL0" );
	var.getVar( la_lut.EXTNORM_NL0 );

	var = lut_grp.getVar( "SSA_NL0" );
	var.getVar( la_lut.SSA_NL0 );

	var = lut_grp.getVar( "QEXT_NL0" );
	var.getVar( la_lut.QEXT_NL0 );

	var = lut_grp.getVar( "BEXT_NL0" );
	var.getVar( la_lut.BEXT_NL0 );

	var = lut_grp.getVar( "VEXT_NL0" );
	var.getVar( la_lut.VEXT_NL0 );

	var = lut_grp.getVar( "MEXT_NL0" );
	var.getVar( la_lut.MEXT_NL0 );

	var = lut_grp.getVar( "SBAR_NL0" );
	var.getVar( la_lut.SBAR_NL0 );

	var = lut_grp.getVar( "INT_NL0" );
	var.getVar( la_lut.INT_NL0 );

	var = lut_grp.getVar( "Fd_NL0" );
	var.getVar( la_lut.Fd_NL0 );

	var = lut_grp.getVar( "T_NL0" );
	var.getVar( la_lut.T_NL0 );

	var = lut_grp.getVar( "OMEGA0" );
	var.getVar( la_lut.OMEGA0 );

	var = lut_grp.getVar( "ROD" );
	var.getVar( la_lut.ROD );

	var = lut_grp.getVar( "GOD" );
	var.getVar( la_lut.GOD );

	delete nc_input;

   return status;
}


/**************************************************************************
 * NAME: read_ocean_aerosol_lut()
 *
 * DESCRIPTION: Read Ocean Aerosol NetCDF4 LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_ocean_aerosol_lut( dtOceanAerosolLUT  &oa_lut )
{
	int status = DTDB_SUCCESS;

	std::string filepath = get_option( INPUT_NC4_LUT );
	if (filepath.empty()) {
		filepath = get_option(INPUT_DT_NC4_LUT);
	}
	NcFile* nc_input;
	try {
		nc_input = new NcFile(filepath, NcFile::read );
	}
	catch( NcException& e) {
		e.what();
		cerr << "DtLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
		return DTDB_FAIL;
	}
	NcGroup lut_grp = nc_input->getGroup( LUT_OCEAN_AEROSOL );

	NcVar var = lut_grp.getVar ( "JPHI" );
	var.getVar( oa_lut.JPHI );

	var = lut_grp.getVar ( "PHC" );
	var.getVar( oa_lut.PHC );

	var = lut_grp.getVar ( "THET" );
	var.getVar( oa_lut.THET );

	var = lut_grp.getVar ( "THET0" );
	var.getVar( oa_lut.THET0 );

	var = lut_grp.getVar ( "WAVE" );
	var.getVar( oa_lut.WAVE );

	var = lut_grp.getVar ( "EXTSMALL" );
	var.getVar( oa_lut.EXTSMALL );

	var = lut_grp.getVar ( "RGSS" );
	var.getVar( oa_lut.RGSS );

	var = lut_grp.getVar( "SIGMAS" );
	var.getVar( oa_lut.SIGMAS );

	var = lut_grp.getVar( "MOMENTSSMALL" );
	var.getVar( oa_lut.MOMENTSSMALL );

	var = lut_grp.getVar( "CCNSMALL" );
	var.getVar( oa_lut.CCNSMALL );

	var = lut_grp.getVar( "BACKSCTTSMALL" );
	var.getVar( oa_lut.BACKSCTTSMALL );

	var = lut_grp.getVar( "ASSYMSMALL" );
	var.getVar( oa_lut.ASSYMSMALL );

	var = lut_grp.getVar( "ALBEDOSMALL" );
	var.getVar( oa_lut.ALBEDOSMALL );

	var = lut_grp.getVar( "ALBEDO_R_SMALL" );
	var.getVar( oa_lut.ALBEDO_R_SMALL );

	var = lut_grp.getVar( "ALBEDO_T_SMALL" );
	var.getVar( oa_lut.ALBEDO_T_SMALL );

	var = lut_grp.getVar( "AINTS" );
	var.getVar( oa_lut.AINTS );

	var = lut_grp.getVar( "TAUAS" );
	var.getVar( oa_lut.TAUAS );

	var = lut_grp.getVar( "EFFRADSMALL" );
	var.getVar( oa_lut.EFFRADSMALL );

	var = lut_grp.getVar( "EXTBIG" );
	var.getVar( oa_lut.EXTBIG );

	var = lut_grp.getVar( "RGSB" );
	var.getVar( oa_lut.RGSB );

	var = lut_grp.getVar( "SIGMAB" );
	var.getVar( oa_lut.SIGMAB );

	var = lut_grp.getVar( "MOMENTSBIG" );
	var.getVar( oa_lut.MOMENTSBIG );

	var = lut_grp.getVar( "BACKSCTTBIG" );
	var.getVar( oa_lut.BACKSCTTBIG );

	var = lut_grp.getVar( "ASSYMBIG" );
	var.getVar( oa_lut.ASSYMBIG );

	var = lut_grp.getVar( "ALBEDOBIG" );
	var.getVar( oa_lut.ALBEDOBIG );

	var = lut_grp.getVar( "ALBEDO_R_BIG" );
	var.getVar( oa_lut.ALBEDO_R_BIG );

	var = lut_grp.getVar( "ALBEDO_T_BIG" );
	var.getVar( oa_lut.ALBEDO_T_BIG );

	var = lut_grp.getVar( "AINTB" );
	var.getVar( oa_lut.AINTB );

	var = lut_grp.getVar( "TAUAB" );
	var.getVar( oa_lut.TAUAB );

	var = lut_grp.getVar( "EFFRADBIG" );
	var.getVar( oa_lut.EFFRADBIG );

	var = lut_grp.getVar( "ALBEDO_R_RAY" );
	var.getVar( oa_lut.ALBEDO_R_RAY );

	var = lut_grp.getVar( "ALBEDO_T_RAY" );
	var.getVar( oa_lut.ALBEDO_T_RAY );

	var = lut_grp.getVar( "REF_RAYALL" );
	var.getVar( oa_lut.REF_RAYALL );

	delete nc_input;

   return status;
}

/**************************************************************************
 * NAME: read_water_vapor_lut()
 *
 * DESCRIPTION: Read water vapor NetCDF4 LUT.
 *
 *************************************************************************/

int DtLutNetcdf::read_water_vapor_lut( dtWaterVaporLUT &wv_lut )
{
    int status = DTDB_SUCCESS;

    std::string filepath = get_option(INPUT_NC4_LUT);
	if (filepath.empty()) {
		filepath = get_option(INPUT_DT_NC4_LUT);
	}
    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "DtLutNetcdf:: Failure opening netcdf LUT file: " + filepath << endl;
        return DTDB_FAIL;
    }
    NcGroup lut_grp = nc_input->getGroup( LUT_WATER_VAPOR );

    NcVar var = lut_grp.getVar ( "TRANSM_H20" );
    var.getVar( wv_lut.TRANSM_H20 );
/*
    var = lut_grp.getVar ( "WEIGHTS" );
    var.getVar( wv_lut.WEIGHTS );

    var = lut_grp.getVar ( "REFL_CH2" );
    var.getVar( wv_lut.REFL_CH2 );

    var = lut_grp.getVar ( "RATIO_CH19_TO_CH2" );
    var.getVar( wv_lut.RATIO_CH19_TO_CH2 );
*/
    delete nc_input;

    return status;
}


/**************************************************************************
 * write_global_attributes()
 *
 * Write global attributes to specified netCDF file ID
 *
 **************************************************************************/

int DtLutNetcdf::write_global_attributes( NcFile* nc_output )
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

bool DtLutNetcdf::isPlatformLittleEndian()
{
	unsigned short checkValue = 0xAABB;
	unsigned char* bytePtr = reinterpret_cast<unsigned char*>(&checkValue);

	if (bytePtr[0] == 0xAA)  // big-endian
	{
		return false;
	}
	else  // little-endian
	{
		return true;
	}
}


