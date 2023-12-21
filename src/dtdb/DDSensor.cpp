/**************************************************************************
 *
 * NAME: DDSensor.cpp
 *
 * DESCRIPTION: Object class that reads in L1B radiance and geolocation data
 * for for processing by selected algorithms.
 *
 *  Created on: August 24, 2020
 *      Author: Sam Anderson
 *
 **************************************************************************/

#include <netcdf>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <netcdf>

#include <timeutils.h>
#include <DDataset.hpp>
#include <DDProcess.h>
#include <DDOptions.h>
#include <DDSensor.h>
#include "resam_viirs/RsViirs.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

string nodestr [3] = {"Ascending", "Descending", "Unknown"};
string daynightstr[4] = {"Day", "Night", "Mixed", "Unknown"};

map<string,float> F0 = {
	    {"rhot_410", 167.403},
	    {"rhot_445", 194.417},
	    {"rhot_490", 195.211},
	    {"rhot_550", 187.053},
	    {"rhot_670", 151.462},
	    {"rhot_865", 94.948},
	    {"rhot_1240", 44.830},
	    {"rhot_1380", 36.468},
	    {"rhot_1610", 24.439},
	    {"rhot_2250", 7.686}
};

/**************************************************************************
 * NAME: DDSensor()
 * DESCRIPTION: Class Constructor
 *************************************************************************/
DDSensor::DDSensor()
{
	brad_ = false;
}

/**************************************************************************
 * NAME: ~DDSensor()
 * DESCRIPTION: Class Destructor
 *************************************************************************/
DDSensor::~DDSensor()
{
}

/**************************************************************************
 * NAME: create()
 * DESCRIPTION: Create l2 template data.
 *************************************************************************/

map<string,ddata*> DDSensor::create(vector<size_t> start, vector<size_t> count)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    string filepath_l1b = get_option(INPUT_L1B);
    string filepath_geo = get_option(INPUT_GEO);
    filepath_geo = (filepath_geo=="") ? filepath_l1b : filepath_geo;
    map<string,ddata*> gmap = read_geo( filepath_geo, start, count );
    if (static_cast<ddval<int>*>(gmap["status"])->val != DTDB_SUCCESS) {
        cerr << "DDSensor:: Read GEO failure. " << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete gmap["status"];
    gmap.erase("status");
    omap.insert(gmap.begin(), gmap.end());

    map<string,ddata*> lmap = create_l2( gmap );
    if (static_cast<ddval<int>*>(lmap["status"])->val != DTDB_SUCCESS) {
        cerr << "DDSensor:: Create L2 failure. " << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete lmap["status"];
    lmap.erase("status");
    omap.insert(lmap.begin(), lmap.end());
    return omap;
}

/**************************************************************************
 * NAME: read()
 * DESCRIPTION: Read selection of l1 data.
 *************************************************************************/

map<string,ddata*> DDSensor::read(vector<size_t> start, vector<size_t> count)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    string filepath_l1b = get_option(INPUT_L1B);
    string filepath_geo = get_option(INPUT_GEO);
    filepath_geo = (filepath_geo=="") ? filepath_l1b : filepath_geo;

	map<string,ddata*> lmap = read_l1b( filepath_l1b, start, count );
    if (static_cast<ddval<int>*>(lmap["status"])->val != DTDB_SUCCESS) {
        cerr << "DDSensor::Error: Read L1B Failure. " << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete lmap["status"];
    lmap.erase("status");
    omap.insert(lmap.begin(), lmap.end());

    map<string,ddata*> gmap = read_geo( filepath_geo, start, count );
    if (static_cast<ddval<int>*>(gmap["status"])->val != DTDB_SUCCESS) {
        cerr << "DDSensor::Error: Read GEO Failure. " << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete gmap["status"];
    gmap.erase("status");
    omap.insert(gmap.begin(), gmap.end());

    // Apply solar zenith angle corrections for VIIRS
    if (instrument_ == SENSOR::VIIRS) {
        for ( size_t iy=0; iy<count[0]; iy++) {
            for ( size_t ix=start[1]; ix<count[1]; ix++) {
                double sza = static_cast<ddma<float,2>*>(gmap["solar_zenith"])->pts[iy][ix];
                double cossza = cos(sza*DEGtoRAD);
                for (auto &it : DDProcess::rhot_band_names) {
                    string name = (string) it.first;
                    float rfl = static_cast<ddma<float,2>*>(lmap[name])->pts[iy][ix];
                    if (rfl > 0) {
                         static_cast<ddma<float,2>*>(lmap[name])->pts[iy][ix] = rfl/cossza;
                    }
                }
            }
        }
    }
    return omap;
}

/**************************************************************************
 * NAME: read_land_mask()
 *
 * DESCRIPTION: Read global Land/Sea mask LUT.
 *
 *************************************************************************/
map<string,ddata*>  DDSensor::read_landmask( map<string,ddata*> gmap,
						vector<size_t> start, vector<size_t> count )
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    LandMaskLUT* lm_lut = new LandMaskLUT;
    string filepath = get_option(INPUT_LANDMASK);
    if(! filepath.empty()) {
        NcFile* nc_input;
        try {
            nc_input = new NcFile(filepath, NcFile::read );
        }
        catch( NcException& e) {
            e.what();
            cerr << "DDSensor:: Failure opening Land Mask LUT file: "
                    + filepath << endl;
            pstat->val = DTDB_FAIL;
            return omap;
        }
        NcVar var = nc_input->getVar ( "lat" );
        var.getVar( lm_lut->lat );
        var = nc_input->getVar ( "lon" );
        var.getVar( lm_lut->lon );
        var = nc_input->getVar ( "watermask" );
        var.getVar( lm_lut->watermask );
        delete nc_input;
    }
    else {
        cerr << "DDSensor:: Failure, bad Land Mask file path." << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }

    ddma<unsigned char,2>* dlw = new ddma<unsigned char,2>(dtype::UBYTE, start, count);

    for ( size_t iy=start[0]; iy<count[0]; iy++) {
        for ( size_t ix=start[1]; ix<count[1]; ix++) {
            double lat = static_cast<ddma<float,2>*>(gmap["latitude"])->pts[iy][ix];
            double lon = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[iy][ix];
            int iLat = (90.0+lat)*(LANDMASK_LAT-1)/180.0;
            int iLon = (180.0+lon)*(LANDMASK_LON-1)/360.0;
            iLat = max (iLat, 0);
            iLat = min (iLat, LANDMASK_LAT-1);
            iLon = max (iLon, 0);
            iLon = min (iLon, LANDMASK_LON-1);
            unsigned char c = (unsigned char) lm_lut->watermask[iLat][iLon];
            dlw->pts[iy][ix] = (c < 0.5) ? 1 : 0;
        }
    }
    delete lm_lut;

    omap.insert({"land_water", dlw});
    return omap;
}

/**************************************************************************
 * NAME: random_seed()
 *
 * DESCRIPTION: Obtain a random seed for
 *
 *************************************************************************/
unsigned long int DDSensor::random_seed() {
    /* Seed generator for gsl. */
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (tv.tv_sec + tv.tv_usec);
}

/**************************************************************************
 * NAME: make_noise()
 *
 * DESCRIPTION: Noise generator adapted from OCSSW loadl1.c
 *
 *************************************************************************/
float DDSensor::make_noise(float sigma) {
    unsigned long randSeed = random_seed();
    float noise;
    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, randSeed);
    noise = gsl_ran_gaussian(rng, sigma);
    gsl_rng_free(rng);
    return (noise);
}

/**************************************************************************
 * NAME: noise_model()
 *
 * DESCRIPTION: Virtual instrument-dependent noise model.
 *
 *************************************************************************/
float DDSensor::noise_model(float lt, int iw, float snr_mult)
{ return 0.0; }

/**************************************************************************
 * NAME: read_l1b()
 *
 * DESCRIPTION: Virtual method.
 *
 *************************************************************************/
map<string,ddata*> DDSensor::read_l1b(const string & filepath,
		vector<size_t> start, vector<size_t> count)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});
	return omap;
}

/**************************************************************************
 * NAME: read_l1b_attributes()
 *
 * DESCRIPTION: Virtual method.
 *
 *************************************************************************/
map<string,ddata*> DDSensor::read_l1b_attributes(NcFile* nc_input)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});
	return omap;
}

/**************************************************************************
 * NAME: create_l2()
 *
 * DESCRIPTION: Populate input map for generic granule data.
 *
 *************************************************************************/
map<string, ddata*> DDSensor::create_l2(map<string, ddata*> gmap)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    int spixl = get_option_int(INPUT_SPIX);
    int epixl = get_option_int(INPUT_EPIX);
    int sline = get_option_int(INPUT_SLINE);
    int eline = get_option_int(INPUT_ELINE);

    time_t rawtime;
    struct tm * timeinfo;
    char timebuffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(timebuffer, 80, "%Y-%m-%dT%H:%M:%SZ", timeinfo);
    string date_created = string(timebuffer);

    string filepath_l1b = get_option(INPUT_L1B);
    NcFile* nc_input = new NcFile(filepath_l1b, NcFile::read );
    map<string,ddata*> amap = read_l1b_attributes(nc_input);
    delete nc_input;
    if (static_cast<ddval<int>*>(amap["status"])->val != DTDB_SUCCESS) {
        cerr << "DDSensor::Error: Read L1B attributes failure. " << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete amap["status"];
    amap.erase("status");
    size_t numx = static_cast<ddval<int>*>(amap["num_pixels"])->val;
    size_t numy = static_cast<ddval<int>*>(amap["num_lines"])->val;
    vector<size_t> start = {0,0};
    vector<size_t> count = {numy,numx};
    size_t cpix = numx / 2;
    size_t npix = numx;
    size_t epix = (epixl < 1) ? npix-1 : epixl-1;
    size_t spix = (spixl < 1) ? 0 : spixl-1;
    size_t nscan = numy;
    if (eline < 1) {
        nscan = numy - max(sline - 1, 0);
    } else {
        nscan = (eline - 1) - max(sline - 1, 0);
    }
    size_t cscan = nscan / 2;
    size_t line_cnt = 0;
    vector<size_t> strt = {0};
    vector<size_t> cntr = {nscan};
    vector<size_t> cntc = {numx};

    map<string,ddata*> lmap = read_landmask(gmap, start, count);
    if (static_cast<ddval<int>*>(lmap["status"])->val != DTDB_SUCCESS) {
        cerr << "DDSensor::Error: Read landmask failure. " << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete lmap["status"];
    lmap.erase("status");
    omap.insert(lmap.begin(), lmap.end());
    omap.insert(gmap.begin(), gmap.end());
    omap.insert(amap.begin(), amap.end());
    omap.insert({"mside", amap["mside"]});
    ddma<int,1>* dyear = new ddma<int,1>(dtype::INT, strt, cntr);
    omap.insert({"year", dyear});
    ddma<int,1>* dday = new ddma<int,1>(dtype::INT, strt, cntr);
    omap.insert({"day", dday});
    ddma<int,1>* dmsec = new ddma<int,1>(dtype::INT, strt, cntr);
    omap.insert({"msec", dmsec});
    ddma<int,1>* drows = new ddma<int,1>(dtype::INT, strt, cntr);
    omap.insert({"cntl_pt_rows", drows});
    ddma<int,1>* dcols = new ddma<int,1>(dtype::INT, strt, cntc);
    omap.insert({"cntl_pt_cols", dcols});
    ddma<float,1>* dslon = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"slon", dslon});
    ddma<float,1>* dclon = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"clon", dclon});
    ddma<float,1>* delon = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"elon", delon});
    ddma<float,1>* dslat = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"slat", dslat});
    ddma<float,1>* dclat = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"clat", dclat});
    ddma<float,1>* delat = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"elat", delat});
    ddma<float,1>* dcsol_z = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"csol_z", dcsol_z});
    ddma<float,1>* dtilt = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"tilt", dtilt});

    size_t pix_cnt = 0;
    for ( size_t ix=0; ix<numx; ix++) {
        dcols->pts[ix] = pix_cnt++;
    }
    for ( size_t iy=0; iy<nscan; iy++) {
        double utime = (double) static_cast<ddma<double,1>*>(amap["scan_time"])->pts[iy];
        int16_t year = 0;
        int16_t day = 0;
        double secs = 0.0;
        unix2yds(utime, &year, &day, &secs);
        dyear->pts[iy] = (int32_t) year;
        dday->pts[iy] = (int32_t) day;
        dmsec->pts[iy] = (int32_t) (secs * 1.e3);
        dslon->pts[iy] = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[iy][spix];
        dclon->pts[iy] = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[iy][cpix];
        delon->pts[iy] = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[iy][epix];
        dslat->pts[iy] = static_cast<ddma<float,2>*>(gmap["latitude"])->pts[iy][spix];
        dclat->pts[iy] = static_cast<ddma<float,2>*>(gmap["latitude"])->pts[iy][cpix];
        delat->pts[iy] = static_cast<ddma<float,2>*>(gmap["latitude"])->pts[iy][epix];
        dcsol_z->pts[iy] = static_cast<ddma<float,2>*>(gmap["solar_zenith"])->pts[iy][cpix];
        drows->pts[iy] = line_cnt++;
    }

    double esd = esdist_(&dyear->pts[cscan], &dday->pts[cscan], &dmsec->pts[cscan]);
    double esdist_correction = pow(1.0 / esd, 2);
    float start_center_lon = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[0][cpix];
    float start_center_lat = static_cast<ddma<float,2>*>(gmap["latitude"])->pts[0][cpix];
    float end_center_lon = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[nscan-1][cpix];
    float end_center_lat = static_cast<ddma<float,2>*>(gmap["latitude"])->pts[nscan-1][cpix];

    auto mm = minmax_element(static_cast<ddma<float,2>*>(gmap["latitude"])->pts.data(),
    		static_cast<ddma<float,2>*>(gmap["latitude"])->pts.data() +
			static_cast<ddma<float,2>*>(gmap["latitude"])->pts.num_elements());
    float geospatial_lat_max = *mm.second;
    float geospatial_lat_min = *mm.first;
    mm = minmax_element(static_cast<ddma<float,2>*>(gmap["longitude"])->pts.data(),
        		static_cast<ddma<float,2>*>(gmap["longitude"])->pts.data() +
    			static_cast<ddma<float,2>*>(gmap["longitude"])->pts.num_elements());
    float geospatial_lon_max = *mm.second;
    float geospatial_lon_min = *mm.first;

    float lastLat = start_center_lat;
    int daynight = UNKNOWNSCENE;
    float northern_lat = -90.0;
    float southern_lat = +90.0;
    float eastern_lon = +180.0;
    float western_lon = -180.0;
    for (size_t is=0; is<nscan; is++) {
        for (size_t ip = spix; ip <= epix; ip++) {
            northern_lat = max(northern_lat, (static_cast<ddma<float,2>*>(gmap["latitude"]))->pts[is][ip]);
            southern_lat = min(southern_lat, (static_cast<ddma<float,2>*>(gmap["latitude"]))->pts[is][ip]);
            if (daynight != DAYANDNIGHT) {
                if (static_cast<ddma<float,2>*>(gmap["solar_zenith"])->pts[is][ip] > SOLZNIGHT) {
                    daynight = (daynight == DAYSCENE) ? DAYANDNIGHT : NIGHTSCENE;
                } else {
                    daynight = (daynight == NIGHTSCENE) ? DAYANDNIGHT : DAYSCENE;
                }
            }
        }
        float wl = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[is][spix];
        float el = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[is][epix];
        if (static_cast<ddma<float,2>*>(gmap["latitude"])->pts[is][cpix] < lastLat) {
            wl = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[is][epix];
            el = static_cast<ddma<float,2>*>(gmap["longitude"])->pts[is][spix];
        }
        western_lon = (wl > western_lon) ? wl : western_lon;
        eastern_lon = (el < eastern_lon) ? el : eastern_lon;
        lastLat = static_cast<ddma<float,2>*>(gmap["latitude"])->pts[is][cpix];
    }

    string start_direction = (static_cast<ddma<float,2>*>(gmap["latitude"])->pts[1][cpix] >
    		static_cast<ddma<float,2>*>(gmap["latitude"])->pts[0][cpix]) ?
            nodestr[ASCENDING] : nodestr[DSCENDING];
    string end_direction = (static_cast<ddma<float,2>*>(gmap["latitude"])->pts[nscan-1][cpix] >
    		static_cast<ddma<float,2>*>(gmap["latitude"])->pts[nscan-2][cpix]) ?
            nodestr[ASCENDING] : nodestr[DSCENDING];
    string day_night_flag = daynightstr[daynight];

    // fill map
    ddstr* dstr = new ddstr(date_created.c_str());
    omap.insert({"date_created", dstr});
    dstr = new ddstr(start_direction.c_str());
    omap.insert({"start_direction", dstr});
    dstr = new ddstr(end_direction.c_str());
    omap.insert({"end_direction", dstr});
    dstr = new ddstr(day_night_flag.c_str());
    omap.insert({"day_night_flag", dstr});
    ddval<double>* ddoub = new ddval<double>(dtype::DOUBLE, esdist_correction);
    omap.insert({"earth_sun_distance_correction", ddoub});
    ddval<float>* dflt = new ddval<float>(dtype::FLOAT, start_center_lon);
    omap.insert({"start_center_longitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, start_center_lat);
    omap.insert({"start_center_latitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, end_center_lon);
    omap.insert({"end_center_longitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, end_center_lat);
    omap.insert({"end_center_latitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, northern_lat);
    omap.insert({"northernmost_latitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, southern_lat);
    omap.insert({"southernmost_latitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, eastern_lon);
    omap.insert({"easternmost_longitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, western_lon);
    omap.insert({"westernmost_longitude", dflt});
    dflt = new ddval<float>(dtype::FLOAT, geospatial_lat_max);
    omap.insert({"geospatial_lat_max", dflt});
    dflt = new ddval<float>(dtype::FLOAT, geospatial_lat_min);
    omap.insert({"geospatial_lat_min", dflt});
    dflt = new ddval<float>(dtype::FLOAT, geospatial_lon_max);
    omap.insert({"geospatial_lon_max", dflt});
    dflt = new ddval<float>(dtype::FLOAT, geospatial_lon_min);
    omap.insert({"geospatial_lon_min", dflt});

    ddma<float,1>* gringlat = new ddma<float,1>(dtype::FLOAT, {0}, {8});
    omap.insert({"GRingPointLatitude", gringlat});
    gringlat->pts[0] = dslat->pts[0];
    gringlat->pts[1] = dclat->pts[0];
    gringlat->pts[2] = delat->pts[0];
    gringlat->pts[3] = delat->pts[cscan];
    gringlat->pts[4] = delat->pts[nscan-1];
    gringlat->pts[5] = dclat->pts[nscan-1];
    gringlat->pts[6] = dslat->pts[nscan-1];
    gringlat->pts[7] = dslat->pts[cscan];
    ddma<float,1>* gringlon = new ddma<float,1>(dtype::FLOAT, {0}, {8});
    omap.insert({"GRingPointLongitude", gringlon});
    gringlon->pts[0] = dslon->pts[0];
    gringlon->pts[1] = dclon->pts[0];
    gringlon->pts[2] = delon->pts[0];
    gringlon->pts[3] = delon->pts[cscan];
    gringlon->pts[4] = delon->pts[nscan-1];
    gringlon->pts[5] = dclon->pts[nscan-1];
    gringlon->pts[6] = dslon->pts[nscan-1];
    gringlon->pts[7] = dslon->pts[cscan];
    ddma<int,1>* gringseq = new ddma<int,1>(dtype::INT, {0}, {8});
    omap.insert({"GRingPointSequenceNo", gringseq});
    gringseq->pts[0] = 0;
    gringseq->pts[1] = 1;
    gringseq->pts[2] = 2;
    gringseq->pts[3] = 3;
    gringseq->pts[4] = 4;
    gringseq->pts[5] = 5;
    gringseq->pts[6] = 6;
    gringseq->pts[7] = 7;

    strt = {0};
    cntr = {NTWL};
    ddma<float,1>* ddat = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"vcal_gain", ddat});
    ddat = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"vcal_offset", ddat});
    ddat = new ddma<float,1>(dtype::FLOAT, strt, cntr);
    omap.insert({"F0", ddat});

    fsol_ = esdist_correction;
	for (auto &it : DDProcess::rhot_band_names) {
		string name = (string) it.first;
		ddat->pts[(size_t)it.second] = F0[name]*fsol_;
	}

    return omap;
}

/**************************************************************************
 * NAME: POCI::read_geo()
 *
 * DESCRIPTION: Virtual function to read geolocation data.
 *
 *************************************************************************/
map<string,ddata*> DDSensor::read_geo(const string & filepath,
		vector<size_t> start, vector<size_t> count)
{
	map<string,ddata*> omap;
	return omap;
}

/**************************************************************************
 * NAME: POCI()
 * DESCRIPTION: Class Constructor
 *************************************************************************/
POCI::POCI()
{
}

/**************************************************************************
 * NAME: ~POCI()
 * DESCRIPTION: Class Destructor
 *************************************************************************/
POCI::~POCI()
{
}

/**************************************************************************
 * NAME: noise_model()
 *
 * DESCRIPTION: Virtual instrument-dependent noise model. Noise coefficients
 * entered in order: {C0,C1} such that C0 + C1 * lt
 *
 *************************************************************************/
float POCI::noise_model(float lt, int iw, float snr_mult)
{
    float sigma, noise, snr, scaled_lt;
    float coef[9][2] = {
        {/*412:C0,C1*/2.1525676E-4, 1.283166E-5},
        {/*488:*/2.4496944E-4, 9.037615E-6},
        {/*550*/2.3188611E-4, 9.68949E-6},
        {/*670*/2.011162E-5, 8.404901E-6},
        {/*865*/3.2728847E-5, 2.2284337E-6},
        {/*1240*/0.07877732, 0.00049940},
        {/*1380*/0.0, 0.0},
        {/*1610*/0.26743281, 0.01044864},
        {/*2250*/0.00628912, 0.00021160}
    };
    scaled_lt = lt;
    noise = coef[iw][0] + coef[iw][1] * scaled_lt;
    snr = scaled_lt * snr_mult / noise; //noise model based on
    sigma = 1 / snr;

    return sigma;
}

/**************************************************************************
 * NAME: POCI::read_l1b()
 *
 * DESCRIPTION: Read granule-level PACE L1B data and compute derived status
 * values.
 *
 *************************************************************************/

map<string,ddata*> POCI::read_l1b(const string & filepath,
		vector<size_t> start, vector<size_t> count)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

	instrument_ = SENSOR::POCI;

    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "POCI:: Failure opening PACE L1B input file: " + filepath << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    try {
        NcGroup nc_group = nc_input->getGroup("observation_data");

        vector<size_t> countp = {1,count[0],count[1]};

        vector<size_t> startp = {bindex(410),start[0],start[1]};
        ddma<float,2>* dda = new ddma<float,2>(nc_group,
        		"rhot_blue", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_410", dda});

        startp = {bindex(445),start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_blue", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_445", dda});

        startp = {bindex(490),start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_blue", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_490", dda});

        startp = {bindex(550),start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_blue", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_550", dda});

        startp = {rindex(670),start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_red", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_670", dda});

        startp = {rindex(865),start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_red", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_865", dda});

        startp = {PS1250,start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_SWIR", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_1240", dda});

        startp = {PS1378,start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_SWIR", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_1380", dda});

        startp = {PS1615,start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_SWIR", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_1610", dda});

        startp = {PS2260,start[0],start[1]};
        dda = new ddma<float,2>(nc_group,
        		"rhot_SWIR", dtype::FLOAT, startp, countp);
        omap.insert({"rhot_2250", dda});

// Convert to VIIRS units
        for ( size_t iy=0; iy<count[0]; iy++) {
            for ( size_t ix=start[1]; ix<count[1]; ix++) {
                for (auto &it : DDProcess::rhot_band_names) {
                    const string name = (string) it.first;
                    float rfl = static_cast<ddma<float,2>*>(omap[name])->pts[iy][ix];
                    if (rfl > 0) {
                        if (brad_) {
                            static_cast<ddma<float,2>*>(omap[name])->pts[iy][ix] = rfl;
                        } else {
                            static_cast<ddma<float,2>*>(omap[name])->pts[iy][ix] = M_PI*rfl;
                        }
                    }
                }
            }
        }
    }
    catch( NcException& e) {
        e.what();
        cerr << "POCI::Failure reading PACE L1B input data." << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete nc_input;

    return omap;
}

/**************************************************************************
 * POCI::read_l1b_attributes()
 *
 * Read granule global attributes contained in input file
 *
 *************************************************************************/
map<string,ddata*> POCI::read_l1b_attributes(NcFile* nc_input)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    try {
    	multimap <string, NcGroupAtt> attributes;

		attributes = nc_input->getAtts(NcGroup::Current);
		string tcs,tce,tstr;
		(attributes.find("time_coverage_start")->second).getValues(tcs);
		ddstr* dstr = new ddstr(tcs);
		omap.insert({"time_coverage_start", dstr});
		(attributes.find("time_coverage_end")->second).getValues(tce);
		dstr = new ddstr(tce);
		omap.insert({"time_coverage_end", dstr});
		(attributes.find("instrument")->second).getValues(tstr);
		dstr = new ddstr(tstr);
		omap.insert({"instrument", dstr});
		tstr = "PACE";
		dstr = new ddstr(tstr);
		omap.insert({"platform", dstr});

		int normLt = -1;
		(attributes.find("normalizedLt")->second).getValues(&normLt);
		brad_ = (normLt==0) ? true : false;

	//  2015-08-18T17:30:00.000Z

		ddval<int>* dval = new ddval<int>(dtype::INT, atoi(tcs.substr(0,4).c_str()));
		omap.insert({"start_year", dval});
		int start_month = atoi(tcs.substr(5,2).c_str());
		dval = new ddval<int>(dtype::INT, start_month);
		omap.insert({"start_month", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(8,2).c_str()));
		omap.insert({"start_day", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(11,2).c_str()));
		omap.insert({"start_hour", dval});
		dval = new ddval<int>(dtype::INT, atoi(tce.substr(11,2).c_str()));
		omap.insert({"end_hour", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(14,2).c_str()));
		omap.insert({"start_minute", dval});
		dval = new ddval<int>(dtype::INT, atoi(tce.substr(14,2).c_str()));
		omap.insert({"end_minute", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(17,2).c_str()));
		omap.insert({"start_second", dval});
		dval = new ddval<int>(dtype::INT, atoi(tce.substr(17,2).c_str()));
		omap.insert({"end_second", dval});
		int season = (start_month == 12) ? 0 : start_month/3;
		dval = new ddval<int>(dtype::INT, season);
		omap.insert({"season", dval});

		int num;
		(attributes.find("orbit_number")->second).getValues(&num);
		dval = new ddval<int>(dtype::INT, num);
		omap.insert({"orbit_number", dval});

		dval = new ddval<int>(dtype::INT, nc_input->getDim("ccd_pixels").getSize());
		omap.insert({"num_pixels", dval});
		dval = new ddval<int>(dtype::INT, nc_input->getDim("number_of_scans").getSize());
		omap.insert({"num_lines", dval});

		NcGroup nc_group = nc_input->getGroup("scan_line_attributes");
		vector<size_t> starts = {0};
		vector<size_t> counts = {(size_t)dval->val};
		ddma<unsigned char,1>* ddh = new ddma<unsigned char,1>(nc_group,
				"HAM_side", dtype::UBYTE, starts, counts);
		omap.insert({"mside", ddh});
		ddma<double,1>* dde = new ddma<double,1>(nc_group,
				"time", dtype::DOUBLE, starts, counts);
		omap.insert({"scan_time", dde});

		NcVar var = nc_group.getVar("time");
		std::map<std::string,NcVarAtt> tatts = var.getAtts();
		(tatts.find("units")->second).getValues(tstr);
		size_t pos = tstr.find("20");
		tstr = tstr.substr (pos);
		double lepoch = isodate2unix(tstr.c_str());
	    for ( size_t iy=0; iy<dde->count[0]; iy++) {
	        dde->pts[iy] += lepoch;
	    }
		vector<size_t> startw = {0};
		vector<size_t> countw = {NTWL};
		ddma<int,1>* ddw = new ddma<int,1>(dtype::INT, startw, countw);
		omap.insert({"wavelength", ddw});
		ddw->pts[(size_t)rhot_band::W410] = 410;
		ddw->pts[(size_t)rhot_band::W445] = 445;
		ddw->pts[(size_t)rhot_band::W490] = 490;
		ddw->pts[(size_t)rhot_band::W550] = 550;
		ddw->pts[(size_t)rhot_band::W670] = 670;
		ddw->pts[(size_t)rhot_band::W865] = 865;
		ddw->pts[(size_t)rhot_band::W1240] = 1240;
		ddw->pts[(size_t)rhot_band::W1380] = 1378;
		ddw->pts[(size_t)rhot_band::W1610] = 1615;
		ddw->pts[(size_t)rhot_band::W2250] = 2260;
	}
	catch( NcException& e) {
		e.what();
		cerr << "POCI::Failure reading PACE L1B attributes." << endl;
		pstat->val = DTDB_FAIL;
		return omap;
	}

	return omap;
}

/**************************************************************************
 * POCI::read_geo()
 *
 * Read geolocation input file
 *
 *************************************************************************/
map<string,ddata*> POCI::read_geo(const string & filepath,
		vector<size_t> start, vector<size_t> count)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "POCI:: Failure opening PACE L1B input file: " + filepath << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }

    ddma<float,2>* ddsen;
    ddma<float,2>* ddsol;
    try {
        NcGroup nc_group = nc_input->getGroup("geolocation_data");

        ddma<float,2>* dda = new ddma<float,2>(nc_group,
        		"latitude", dtype::FLOAT, start, count);
        omap.insert({"latitude", dda});

        dda = new ddma<float,2>(nc_group,
        		"longitude", dtype::FLOAT, start, count);
        omap.insert({"longitude", dda});

        dda = new ddma<float,2>(nc_group,
        		"height", dtype::FLOAT, start, count);
        omap.insert({"elevation", dda});

        ddsen = new ddma<float,2>(nc_group,
        		"sensor_azimuth", dtype::FLOAT, start, count);
        omap.insert({"sensor_azimuth", ddsen});

        dda = new ddma<float,2>(nc_group,
        		"sensor_zenith", dtype::FLOAT, start, count);
        omap.insert({"sensor_zenith", dda});

        ddsol = new ddma<float,2>(nc_group,
        		"solar_azimuth", dtype::FLOAT, start, count);
        omap.insert({"solar_azimuth", ddsol});

        dda = new ddma<float,2>(nc_group,
        		"solar_zenith", dtype::FLOAT, start, count);
        omap.insert({"solar_zenith", dda});
    }
    catch( NcException& e) {
        e.what();
        cerr << "POCI::Failure reading PACE L1B input data." << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete nc_input;

    ddma<float,2>* ddra = new ddma<float,2>(dtype::FLOAT, start, count);
    for ( size_t iy=0; iy<count[0]; iy++) {
        for ( size_t ix=start[1]; ix<count[1]; ix++) {
			float raa = ddsen->pts[iy][ix] - ddsol->pts[iy][ix]  - 180.0;
			if (raa > 180.0) raa = raa - 360.0;
			if (raa < -180.0) raa = raa + 360.0;
			if (raa < 0.0) raa = -raa;
			ddra->pts[iy][ix] = raa;
        }
    }
    omap.insert({"relative_azimuth", ddra});
    return omap;
}

/**************************************************************************
 * NAME: VIIRS()
 * DESCRIPTION: Class Constructor
 *************************************************************************/
VIIRS::VIIRS()
{
}

/**************************************************************************
 * NAME: ~VIIRS()
 * DESCRIPTION: Class Destructor
 *************************************************************************/
VIIRS::~VIIRS()
{
}

/**************************************************************************
 * NAME: noise_model()
 *
 * DESCRIPTION: Virtual instrument-dependent noise model. Noise coefficients
 * entered in order: {C0,C1} such that C0 + C1 * lt
 *
 *************************************************************************/
float VIIRS::noise_model(float lt, int iw, float snr_mult)
{
    float sigma, noise, snr, scaled_lt;
    float coef[9][2] = {
        {/*412:C0,C1*/0.05499859, 0.00008340},
        {/*488:*/0.01927545, 0.00009450},
        {/*550*/0.08769538, 0.00007000},
        {/*670*/0.00496291, 0.00014050},
        {/*865*/0.00312263, 0.00018600},
        {/*1240*/0.07877732, 0.00049940},
        {/*1380*/0.0, 0.0},
        {/*1610*/0.26743281, 0.01044864},
        {/*2250*/0.00628912, 0.00021160}
    };
    scaled_lt = lt;
    noise = coef[iw][0] + coef[iw][1] * scaled_lt;
    snr = scaled_lt * snr_mult / noise; //noise model based on
    sigma = 1 / snr;

    return sigma;
}

/**************************************************************************
 * NAME: VIIRS::read_l1b()
 *
 * DESCRIPTION: Read granule-level VIIRS L1B data and compute derived status
 * values.
 *
 *************************************************************************/
map<string,ddata*> VIIRS::read_l1b(const string & filepath,
		vector<size_t> start, vector<size_t> count)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    instrument_ = SENSOR::VIIRS;

    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "VIIRS:: Failure opening VIIRS L1B file: " + filepath << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    try {
        NcGroup nc_group = nc_input->getGroup("observation_data");

        ddma<float,2>* dda = new ddma<float,2>(nc_group, "M01", dtype::FLOAT, start, count);
        omap.insert({"rhot_410", dda});

        dda = new ddma<float,2>(nc_group, "M02", dtype::FLOAT, start, count);
        omap.insert({"rhot_445", dda});

        dda = new ddma<float,2>(nc_group, "M03", dtype::FLOAT, start, count);
        omap.insert({"rhot_490", dda});

        dda = new ddma<float,2>(nc_group, "M04", dtype::FLOAT, start, count);
        omap.insert({"rhot_550", dda});

        dda = new ddma<float,2>(nc_group, "M05", dtype::FLOAT, start, count);
        omap.insert({"rhot_670", dda});

        dda = new ddma<float,2>(nc_group, "M07", dtype::FLOAT, start, count);
        omap.insert({"rhot_865", dda});

        dda = new ddma<float,2>(nc_group, "M08", dtype::FLOAT, start, count);
        omap.insert({"rhot_1240", dda});

        dda = new ddma<float,2>(nc_group, "M09", dtype::FLOAT, start, count);
        omap.insert({"rhot_1380", dda});

        dda = new ddma<float,2>(nc_group, "M10", dtype::FLOAT, start, count);
        omap.insert({"rhot_1610", dda});

        dda = new ddma<float,2>(nc_group, "M11", dtype::FLOAT, start, count);
        omap.insert({"rhot_2250", dda});
    }
    catch( NcException& e) {
        e.what();
        cerr << "VIIRS::Failure reading VIIRS L1B data." << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete nc_input;
/*
    RsViirs* rv = new RsViirs( count[0], count[1]);
    rv->generate_sort_index();
	for (auto &it : omap) {
		string name = (string) it.first;
		rv->resort( (static_cast<ddma<float,2>*>(omap[name]))->pts );
		rv->fill_the_fills( (static_cast<ddma<float,2>*>(omap[name]))->pts );
	}
	delete rv;
*/
    return omap;
}

/**************************************************************************
 * VIIRS::read_l1b_attributes()
 *
 * Read granule global attributes contained in input file
 *
 *************************************************************************/
map<string,ddata*> VIIRS::read_l1b_attributes(NcFile* nc_input)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    try {
		multimap <string, NcGroupAtt> attributes;

		attributes = nc_input->getAtts(NcGroup::Current);
		string tcs,tce,tstr;
		(attributes.find("time_coverage_start")->second).getValues(tcs);
		ddstr* dstr = new ddstr(tcs);
		omap.insert({"time_coverage_start", dstr});
		(attributes.find("time_coverage_end")->second).getValues(tce);
		dstr = new ddstr(tce);
		omap.insert({"time_coverage_end", dstr});
		(attributes.find("instrument")->second).getValues(tstr);
		dstr = new ddstr(tstr);
		omap.insert({"instrument", dstr});
		(attributes.find("platform")->second).getValues(tstr);
		dstr = new ddstr(tstr);
		omap.insert({"platform", dstr});

	//  2015-08-18T17:30:00.000Z

		ddval<int>* dval = new ddval<int>(dtype::INT, atoi(tcs.substr(0,4).c_str()));
		omap.insert({"start_year", dval});
		int start_month = atoi(tcs.substr(5,2).c_str());
		dval = new ddval<int>(dtype::INT, start_month);
		omap.insert({"start_month", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(8,2).c_str()));
		omap.insert({"start_day", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(11,2).c_str()));
		omap.insert({"start_hour", dval});
		dval = new ddval<int>(dtype::INT, atoi(tce.substr(11,2).c_str()));
		omap.insert({"end_hour", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(14,2).c_str()));
		omap.insert({"start_minute", dval});
		dval = new ddval<int>(dtype::INT, atoi(tce.substr(14,2).c_str()));
		omap.insert({"end_minute", dval});
		dval = new ddval<int>(dtype::INT, atoi(tcs.substr(17,2).c_str()));
		omap.insert({"start_second", dval});
		dval = new ddval<int>(dtype::INT, atoi(tce.substr(17,2).c_str()));
		omap.insert({"end_second", dval});
		int season = (start_month == 12) ? 0 : start_month/3;
		dval = new ddval<int>(dtype::INT, season);
		omap.insert({"season", dval});

		int num;
		(attributes.find("orbit_number")->second).getValues(&num);
		dval = new ddval<int>(dtype::INT, num);
		omap.insert({"orbit_number", dval});

		dval = new ddval<int>(dtype::INT, nc_input->getDim("number_of_pixels").getSize());
		omap.insert({"num_pixels", dval});
		dval = new ddval<int>(dtype::INT, nc_input->getDim("number_of_lines").getSize());
		omap.insert({"num_lines", dval});

		NcGroup nc_group = nc_input->getGroup("scan_line_attributes");
		vector<size_t> starts = {0};
		vector<size_t> countl = {(size_t)dval->val};
		vector<size_t> counts = {(size_t)dval->val/VDETECTORS};
		ddma<unsigned char,1>* ddh = new ddma<unsigned char,1>(dtype::UBYTE, starts, countl);
		omap.insert({"mside", ddh});
		ddma<double,1>* dde = new ddma<double,1>(dtype::DOUBLE, starts, countl);
		omap.insert({"scan_time", dde});
		double scantimes[VIIRSSCANS];
		nc_group.getVar("ev_mid_time").getVar( &scantimes[0] );
		for (size_t i=0; i<countl[0]; i++) {
			ddh->pts[i] = 255;
			dde->pts[i] = tai58_to_unix(scantimes[i/VIIRSSCANS]);
		}

		vector<size_t> startw = {0};
		vector<size_t> countw = {NTWL};
		ddma<int,1>* ddw = new ddma<int,1>(dtype::INT, startw, countw);
		omap.insert({"wavelength", ddw});
		ddw->pts[(size_t)rhot_band::W410] = 410;
		ddw->pts[(size_t)rhot_band::W445] = 445;
		ddw->pts[(size_t)rhot_band::W490] = 490;
		ddw->pts[(size_t)rhot_band::W550] = 550;
		ddw->pts[(size_t)rhot_band::W670] = 670;
		ddw->pts[(size_t)rhot_band::W865] = 865;
		ddw->pts[(size_t)rhot_band::W1240] = 1240;
		ddw->pts[(size_t)rhot_band::W1380] = 1380;
		ddw->pts[(size_t)rhot_band::W1610] = 1610;
		ddw->pts[(size_t)rhot_band::W2250] = 2250;
	}
	catch( NcException& e) {
		e.what();
		cerr << "VIIRS::Failure reading VIIRS L1B attributes." << endl;
		pstat->val = DTDB_FAIL;
		return omap;
	}

	return omap;
}

/**************************************************************************
 * NAME: VIIRS::read_geo()
 *
 * DESCRIPTION: Read VIIRS geolocation data and compute derived status
 * values.
 *
 *************************************************************************/

map<string,ddata*> VIIRS::read_geo(const string & filepath,
		vector<size_t> start, vector<size_t> count)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

    NcFile* nc_input;
    try {
        nc_input = new NcFile(filepath, NcFile::read );
    }
    catch( NcException& e) {
        e.what();
        cerr << "VIIRS:: Failure opening VIIRS geo file: " + filepath << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }

    try {
        NcGroup nc_group = nc_input->getGroup("geolocation_data");

        ddma<float,2>* dda = new ddma<float,2>(nc_group,
        		"latitude", dtype::FLOAT, start, count);
        omap.insert({"latitude", dda});

        dda = new ddma<float,2>(nc_group,
        		"longitude", dtype::FLOAT, start, count);
        omap.insert({"longitude", dda});

        dda = new ddma<float,2>(nc_group,
        		"height", dtype::FLOAT, start, count);
        omap.insert({"elevation", dda});

        dda = new ddma<float,2>(nc_group,
        		"sensor_azimuth", dtype::FLOAT, start, count);
        omap.insert({"sensor_azimuth", dda});

        dda = new ddma<float,2>(nc_group,
        		"sensor_zenith", dtype::FLOAT, start, count);
        omap.insert({"sensor_zenith", dda});

        dda = new ddma<float,2>(nc_group,
        		"solar_azimuth", dtype::FLOAT, start, count);
        omap.insert({"solar_azimuth", dda});

        dda = new ddma<float,2>(nc_group,
        		"solar_zenith", dtype::FLOAT, start, count);
        omap.insert({"solar_zenith", dda});
    }
    catch( NcException& e) {
        e.what();
        cerr << "VIIRS::Failure reading VIIRS geo data." << endl;
        pstat->val = DTDB_FAIL;
        return omap;
    }
    delete nc_input;

    ddma<float,2>* ddra = new ddma<float,2>(dtype::FLOAT, start, count);
    for ( size_t iy=0; iy<count[0]; iy++) {
        for ( size_t ix=start[1]; ix<count[1]; ix++) {
        	float sen = static_cast<ddma<float,2>*>(omap["sensor_azimuth"])->pts[iy][ix];
        	float sol = static_cast<ddma<float,2>*>(omap["solar_azimuth"])->pts[iy][ix];
			float raa = sen - sol  - 180.0;
			if (raa > 180.0) raa = raa - 360.0;
			if (raa < -180.0) raa = raa + 360.0;
			if (raa < 0.0) raa = -raa;
			ddra->pts[iy][ix] = raa;
        }
    }
    omap.insert({"relative_azimuth", ddra});

    return omap;
}

