/*******************************************************************************
 *
 * NAME: DDAncillary.cpp
 *
 * DESCRIPTION: Reads ancillary meteorological data
 *
 *  Created on: April 17, 2019
 *      Author: Sam Anderson
 *
 *******************************************************************************/

#include <cmath>
#include <DDAncillary.h>
#include <DDOptions.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

/**************************************************************************
 * NAME: Ancillary()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DDAncillary::DDAncillary() {
    lut_ = nullptr;
}

/**************************************************************************
 * NAME: ~DDAncillary()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DDAncillary::~DDAncillary() {
    if (lut_ != nullptr) {
        delete lut_;
        lut_ = nullptr;
    }
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Read and time-interpolate meteorology lut for granule
 *
 *************************************************************************/

int DDAncillary::initialize(int hr, int min)
{
    int status = DTDB_SUCCESS;

    string filepath = get_option(INPUT_MET);
    if ( !filepath.empty()) {
    } else {
        string filepath_1 = get_option(INPUT_MET1);
        string filepath_2 = get_option(INPUT_MET2);

        if (filepath_1.empty() || filepath_2.empty()) {
            cerr << "DDAncillary:: Ancillary files missing: "
                            + filepath_1 << endl;
            return DTDB_FAIL;
        }

        met_lut* lut_1 = new met_lut;
        memset(lut_1, 0, sizeof(met_lut));
        status = read_merra2_file( filepath_1, lut_1);
        if (status != DTDB_SUCCESS) {
            cerr << "DDAncillary:: failure reading file 1: "
                            + filepath_1 << endl;
            return status;
        }
        met_lut* lut_2 = new met_lut;
        memset(lut_2, 0, sizeof(met_lut));
        status = read_merra2_file( filepath_2, lut_2);
        if (status != DTDB_SUCCESS) {
            cerr << "DDAncillary:: failure reading file 2: "
                            + filepath_2 << endl;
            return status;
        }

        lut_ = new met_lut;
        memset(lut_, 0, sizeof(met_lut));

        for (int iy=0; iy<NMLATS; iy++) {
            for (int ix=0; ix<NMLONS; ix++) {
                lut_->ps[iy][ix] = lut_1->ps[iy][ix] +
                        (lut_2->ps[iy][ix]-lut_1->ps[iy][ix])*min/60.0;
                lut_->oz[iy][ix] = lut_1->oz[iy][ix] +
                        (lut_2->oz[iy][ix]-lut_1->oz[iy][ix])*min/60.0;
                lut_->pw[iy][ix] = lut_1->pw[iy][ix] +
                        (lut_2->pw[iy][ix]-lut_1->pw[iy][ix])*min/60.0;
                lut_->u_ws[iy][ix] = lut_1->u_ws[iy][ix] +
                        (lut_2->u_ws[iy][ix]-lut_1->u_ws[iy][ix])*min/60.0;
                lut_->v_ws[iy][ix] = lut_1->v_ws[iy][ix] +
                        (lut_2->v_ws[iy][ix]-lut_1->v_ws[iy][ix])*min/60.0;
            }
        }

        delete lut_1;
        delete lut_2;
    }

    return status;
}

/**************************************************************************
 * NAME: read_merra2_file()
 *
 * DESCRIPTION: Open and read a MERRA2 file
 *
 *************************************************************************/

int DDAncillary::read_merra2_file(const string filepath, met_lut* in)
{
    int status = DTDB_SUCCESS;

    try {
        NcFile* fio = new NcFile(filepath, NcFile::read);
        NcVar var = fio->getVar("PS");
        var.getVar(&in->ps[0][0]);
        var = fio->getVar("TO3");
        var.getVar(&in->oz[0][0]);
        var = fio->getVar("TQV");
        var.getVar(&in->pw[0][0]);
        var = fio->getVar("U10M");
        var.getVar(&in->u_ws[0][0]);
        var = fio->getVar("V10M");
        var.getVar(&in->v_ws[0][0]);
        delete fio;
    } catch (NcException& e) {
        e.what();
        cerr << "DDAncillary:: Failure opening ancillary file: "
                        + filepath << endl;
        return DTDB_FAIL;
    }

    return status;
}

/**************************************************************************
 * NAME: get_meteorology()
 *
 * DESCRIPTION: Get ancillary data from look-up table.
 * Initialize granule-level output data attributes
 *
 *************************************************************************/

map<string,ddata*> DDAncillary::read( map<string,ddata*> imap )
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});
	metio io;

	ddma<float,2>* plon = static_cast<ddma<float,2>*>(imap["longitude"]);
	ddma<float,2>* plat = static_cast<ddma<float,2>*>(imap["latitude"]);
	ddma<float,2>* pws = new ddma<float,2>(dtype::FLOAT, plon->start, plon->count);
	omap.insert({"windspeed", pws});
	ddma<float,2>* pwd = new ddma<float,2>(dtype::FLOAT, plon->start, plon->count);
	omap.insert({"windangle", pwd});
	ddma<float,2>* ppw = new ddma<float,2>(dtype::FLOAT, plon->start, plon->count);
	omap.insert({"water_vapor", ppw});
	ddma<float,2>* poz = new ddma<float,2>(dtype::FLOAT, plon->start, plon->count);
	omap.insert({"ozone", poz});
	ddma<float,2>* pps = new ddma<float,2>(dtype::FLOAT, plon->start, plon->count);
	omap.insert({"pressure", pps});
	ddma<unsigned char,2>* puc = new ddma<unsigned char,2>(dtype::UBYTE, plon->start, plon->count);
	omap.insert({"cloud_mask", puc});

    string filepath = get_option(INPUT_MET);
    if ( !filepath.empty()) {
        NcFile* fio = new NcFile(filepath, NcFile::read);
        NcVar var = fio->getVar("PS");
        var.getVar(plon->start, plon->count, (float*)pps->ptr);
        var = fio->getVar("TO3");
        var.getVar(plon->start, plon->count, (float*)poz->ptr);
        var = fio->getVar("TQV");
        var.getVar(plon->start, plon->count, (float*)ppw->ptr);

    	ddma<float,2>* pu_ws = new ddma<float,2>(dtype::FLOAT, plon->start, plon->count);
    	ddma<float,2>* pv_ws = new ddma<float,2>(dtype::FLOAT, plon->start, plon->count);
        var = fio->getVar("U10M");
        var.getVar(plon->start, plon->count, (float*)pu_ws->ptr);
        var = fio->getVar("V10M");
        var.getVar(plon->start, plon->count, (float*)pv_ws->ptr);

		for ( size_t iy=0; iy<plon->count[0]; iy++) {
			for ( size_t ix=0; ix<plon->count[1]; ix++) {
				// wind speed
				pws->pts[iy][ix] = sqrt(pu_ws->pts[iy][ix]*pu_ws->pts[iy][ix] +
						pv_ws->pts[iy][ix]*pv_ws->pts[iy][ix]);
				// wind direction
				float wd = atan2(-1.0*pu_ws->pts[iy][ix],-1.0*pv_ws->pts[iy][ix]) *
						180.0/M_PI;
				pwd->pts[iy][ix] = (wd < 0) ? wd + 360.0 : wd;
			}
		}
        delete fio;
        delete pu_ws;
        delete pv_ws;

    } else {
		for ( size_t iy=0; iy<plon->count[0]; iy++) {
			for ( size_t ix=0; ix<plon->count[1]; ix++) {
				get_met(plat->pts[iy][ix], plon->pts[iy][ix], &io);
				pws->pts[iy][ix] = io.ws;
				pwd->pts[iy][ix] = io.wd;
				ppw->pts[iy][ix] = io.pw;
				poz->pts[iy][ix] = io.oz;
				pps->pts[iy][ix] = io.ps;
			}
		}
    }
    filepath = get_option(INPUT_CLDMASK);
    if ( !filepath.empty()) {
        NcFile* fio = new NcFile(filepath, NcFile::read);
        NcVar var = fio->getVar("ADJ_MASK");
        var.getVar(plon->start, plon->count, (unsigned char*)puc->ptr);
        delete fio;
    } else {
		for ( size_t iy=0; iy<plon->count[0]; iy++) {
			for ( size_t ix=0; ix<plon->count[1]; ix++) {
				puc->pts[iy][ix] = DFILL_UBYTE;
			}
		}
    }
    return omap;
}

int DDAncillary::get_met( const float lat, const float lon, metio* io )
{
    int status = DTDB_SUCCESS;

    if(lat > 90.0 || lat < -90.0) {
        cerr << "DDAncillary:: lat out of range, must be 90 to -90:" <<  lat << endl;
        return DTDB_FAIL;
    }
    if (lon < -180.0 || lon > 180.0) {
        cerr << "DDAncillary:: lon out of range, must be [-180,180):" <<  lon << endl;
        return DTDB_FAIL;
    }

    int   i1 = 0;
    int   i2 = 0;
    int   j1 = 0;
    int   j2 = 0;

    status = get_interp_indexes(lat, lon, i1, i2, j1, j2);
    if (status != DTDB_SUCCESS) {
        cerr << "DDAncillary:: Failed to get indexes for interpolation: " << status << endl;
        return status;
    }

    float f[4] = {0,0,0,0};
    float x1 = index2lon(i1);
    float x2 = index2lon(i2);
    float y1 = index2lat(j1);
    float y2 = index2lat(j2);

// perform 2D bilinear interpolation (http://en.wikipedia.org/wiki/Bilinear_interpolation)
// wind speed u-direction

    f[0] = lut_->u_ws[j1][i1];
    f[1] = lut_->u_ws[j1][i2];
    f[2] = lut_->u_ws[j2][i2];
    f[3] = lut_->u_ws[j2][i1];

    float u = f[0]*(x2-lon)*(y2-lat) + f[1]*(lon-x1)*(y2-lat) +
          f[2]*(lon-x1)*(lat-y1) + f[3]*(x2-lon)*(lat-y1);
    u /= ((x2-x1)*(y2-y1));

// wind speed v-direction

    f[0] = lut_->v_ws[j1][i1];
    f[1] = lut_->v_ws[j1][i2];
    f[2] = lut_->v_ws[j2][i2];
    f[3] = lut_->v_ws[j2][i1];

    float v = f[0]*(x2-lon)*(y2-lat) + f[1]*(lon-x1)*(y2-lat) +
          f[2]*(lon-x1)*(lat-y1) + f[3]*(x2-lon)*(lat-y1);
    v /= ((x2-x1)*(y2-y1));

// wind speed
    io->ws = sqrt(u*u + v*v);

// wind direction
    float wd = atan2(-1.0*u,-1.0*v) * 180.0/M_PI;
    io->wd = (wd < 0) ? wd + 360.0 : wd;

// precipitable water vapor
    f[0] = lut_->pw[j1][i1];
    f[1] = lut_->pw[j1][i2];
    f[2] = lut_->pw[j2][i2];
    f[3] = lut_->pw[j2][i1];

    float pwat = f[0]*(x2-lon)*(y2-lat) + f[1]*(lon-x1)*(y2-lat) +
         f[2]*(lon-x1)*(lat-y1) + f[3]*(x2-lon)*(lat-y1);

    io->pw = pwat / ((x2-x1)*(y2-y1));
    if (io->pw <= 0.0) {
        io->pw = 1.0e-5;
    }

// ozone
    f[0] = lut_->oz[j1][i1];
    f[1] = lut_->oz[j1][i2];
    f[2] = lut_->oz[j2][i2];
    f[3] = lut_->oz[j2][i1];

    float oz = f[0]*(x2-lon)*(y2-lat) + f[1]*(lon-x1)*(y2-lat) +
        f[2]*(lon-x1)*(lat-y1) + f[3]*(x2-lon)*(lat-y1);

    io->oz = oz / ((x2-x1)*(y2-y1));

// surface pressure
	f[0] = lut_->ps[j1][i1];
	f[1] = lut_->ps[j1][i2];
	f[2] = lut_->ps[j2][i2];
	f[3] = lut_->ps[j2][i1];

	float ps = f[0]*(x2-lon)*(y2-lat) + f[1]*(lon-x1)*(y2-lat) +
		f[2]*(lon-x1)*(lat-y1) + f[3]*(x2-lon)*(lat-y1);

	io->ps = ps / ((x2-x1)*(y2-y1));


    return status;
}

/**************************************************************************
 * NAME: lat2index()
 *
 * DESCRIPTION: Converts latitudes into indexes into the 2D
 * lat/lon data arrays. Assumes a 1x1 degree resolution. Index returned is
 * for the lower-left corner of the grid box containing the location specified.
 *
 *************************************************************************/

int DDAncillary::lat2index(const float lat) {

    if (lat < -90.0 || lat > 90.0) {
        cerr <<  "DDAncillary:: Invalid latitude: " << lat << endl;
        return -1;
    }
    int index = floor((lat + 90.0) / 0.5);

    return index;
}

/**************************************************************************
 * NAME: index2lat()
 *
 * DESCRIPTION: Converts indexes into latitude.Assumes a 1x1 degree resolution.
 *
 *************************************************************************/

float DDAncillary::index2lat(const int index) {

    if (index < 0 || index > NMLATS) {
        cerr <<  "DDAncillary:: Index is out of bounds: " << index << endl;
        return -1;
    }
    float lat = -90.0 + 0.5*(index);

    return lat;
}

/**************************************************************************
 * NAME: lon2index()
 *
 * DESCRIPTION: Converts longitudes into indexes into the 2D
 * lat/lon data arrays. Assumes a 1x1 degree resolution. Index returned is
 * for the lower-left corner of the grid box containing the location specified.
 *
 *************************************************************************/

int DDAncillary::lon2index(const float lon) {

    if (lon < -180.0 || lon >= 180.0) {
        cerr <<  "DDAncillary:: Invalid longitude: " << lon << endl;
        return -1;
    }
    int index = floor((lon + 180.0)/0.625);

    return index;
}

/**************************************************************************
 * NAME: index2lon()
 *
 * DESCRIPTION: Converts indexes into longitude. Assumes a 1x1 degree resolution.
 *
 *************************************************************************/

float DDAncillary::index2lon(const int index) {

    if (index < 0 || index > NMLONS) {
        cerr <<  "DDAncillary:: Index is out of bounds of windspeed array: " << index << endl;
        return -1;
    }
    float lon = -180.0 + 0.625*index;

    return lon;
}

/**************************************************************************
 * NAME: get_interp_indexes()
 *
 * DESCRIPTION: Based on values of lat and lon, returns indexes of surrounding
 * points for 2D linear interpolation.
 * i1 and i2 are the indexes to the west and east of the given longitude.
 * j1 and j2 are the indexes to the south and north of the given latitude.
 * Corner cases involving dateline and north pole (exactly 90.0N)
 * should be handled correctly.
 *
 *************************************************************************/

int DDAncillary::get_interp_indexes(const float lat, const float lon,
                       int& i1, int& i2, int& j1, int& j2)
{

// specifically exclude locations where lat==90.0 exactly. This would break
//  the interpolation below.
    if(lat >= 90.0 || lat < -90.0) {
        cerr << "DDAncillary:: lat out of range, must be [-90,90]:" << lat << endl;
        return DTDB_FAIL;
    }
    if (lon < -180.0 || lon >= 180.0) {
        cerr << "DDAncillary:: lon out of range, must be [-180,180]:" << lon << endl;
        return DTDB_FAIL;
    }

// get indexes into the windspeed table for lower left point.
    i1 = lon2index(lon);
    if (i1 == -1) {
        cerr <<  "DDAncillary:: Failed to convert longitude to index: " << i1 << endl;
        return DTDB_FAIL;
    }

    j1 = lat2index(lat);
    if (j1 == -1) {
        cerr <<  "DDAncillary:: Failed to convert latitude to index: " << j1 << endl;
        return DTDB_FAIL;
    }

// simply calculate the next indexes. Watch for instances where we're
// crossing the dateline (i2 == NMLONS) and circle back to index 0.
    i2 = i1 + 1;
    if (i2 == NMLONS) i2 = 0;
    j2 = j1 + 1;

    return DTDB_SUCCESS;
}

