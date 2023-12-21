/*******************************************************************************
 *
 * NAME: DDAncillary.h
 *
 * DESCRIPTION: Object class that supports reading and interpolation of
 * ancillary data extracted from MERA2 files.
 *
 *  Created on: April 11, 2019
 *      Author: Sam Anderson, NASA Ocean Color
 *
 *******************************************************************************/

#ifndef DDAncillary_H_
#define DDAncillary_H_

#include <DDProcess.h>

using namespace std;

const int NMLONS = 576;
const int NMLATS = 361;

struct metio {
    float ps;
    float oz;
    float pw;
    float ws;
    float wd;
};

struct met_lut {
    float ps[NMLATS][NMLONS];
    float oz[NMLATS][NMLONS];
    float pw[NMLATS][NMLONS];
    float u_ws[NMLATS][NMLONS];
    float v_ws[NMLATS][NMLONS];
};

class DDProcess;

class DDAncillary
{
public:

//  Class constructor
    DDAncillary ();

//  Class destructor
	~DDAncillary ();

//  Initialize lut
	int initialize( int hr, int min );

//  Retrieve meteorology
	map<string,ddata*> read( map<string,ddata*> imap );

protected:

// Meteorology LUT

    met_lut* lut_;

// Member functions

    int   read_merra2_file(const string filename, met_lut* lut);

    int   get_met(const float lat, const float lon, metio* io );

    int   get_interp_indexes(const float lat, const float lon,
                             int& i1, int& i2, int& j1, int& j2);

    float index2lat(const int index);
    float index2lon(const int index);
    int   lat2index(const float lat);
    int   lon2index(const float lon);
};

#endif /* DDAncillary_H_ */
