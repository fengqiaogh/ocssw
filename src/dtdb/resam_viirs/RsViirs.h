/*******************************************************************************
 *
 * NAME: RsViirs.h
 *
 * DESCRIPTION: Object class that supports the resampling of VIIRS GEO and L1B
 * files in order to eliminate unsightly bowtie fill regions in L1B.
 *
 *  Created on: February 17, 2019
 *      Author: Sam Anderson, NASA Ocean Color
 *
 *******************************************************************************/

#ifndef RsViirs_H_
#define RsViirs_H_

#include <string>
#include <math.h>

using namespace std;

static constexpr int 	 RS_SUCCESS = 0;
static constexpr int 	 RS_FAIL = 1;

enum RESAM_ENUM {
    M01,
    M02,
    M03,
    M04,
    M05,
    M06,
    M07,
    M08,
    M09,
    M10,
    M11,
    M12,
    M13,
    M14,
    M15,
    M16,
    NUM_MOD_BANDS
};

// Bands 1-11 reflectance, and bands 12, 14-16 brightnesss temperature
const unsigned short NA_UINT16_FILL = 65535;
const unsigned short MISS_UINT16_FILL = 65534;
const unsigned short ONBOARD_PT_UINT16_FILL = 65533;
const unsigned short ONGROUND_PT_UINT16_FILL = 65532;
const unsigned short ERR_UINT16_FILL = 65531;
const unsigned short VDNE_UINT16_FILL = 65529;
const unsigned short SOUB_UINT16_FILL = 65528;

// Band 13 brightnesss temperature
const float NA_FLOAT32_FILL = -999.9;
const float ONBOARD_PT_FLOAT32_FILL = -999.7;
const float ONGROUND_PT_FLOAT32_FILL = -999.6;
const float VDNE_FLOAT32_FILL = -999.3;

const unsigned short DELETION_ZONE_INT = ONBOARD_PT_UINT16_FILL;
const float DELETION_ZONE_FLOAT = ONBOARD_PT_FLOAT32_FILL;

const int NMBANDS = 16;
const int NSCANS = 203;
const int NPIXELS = 3200;
const int NDETECTORS = 16;
const int NBREAKS = 11;

const float D2R = M_PI/180.0;
const float R2D = 180.0/M_PI;

struct fio {
    float in[NSCANS*NDETECTORS][NPIXELS];
    float out[NSCANS*NDETECTORS][NPIXELS];
};
struct sio {
    short in[NSCANS*NDETECTORS][NPIXELS];
    short out[NSCANS*NDETECTORS][NPIXELS];
};
struct usio {
    unsigned short in[NSCANS*NDETECTORS][NPIXELS];
    unsigned short out[NSCANS*NDETECTORS][NPIXELS];
};
struct uio {
    unsigned int in[NSCANS*NDETECTORS][NPIXELS];
    unsigned int out[NSCANS*NDETECTORS][NPIXELS];
};


class RsViirs
{
public:

//  Class constructor
    RsViirs ();
    RsViirs ( int lines, int pixels );

//  Class destructor
	~RsViirs ();

//  Apply algorithm to data, produce product
	int process();

// Strings
    string config_file_;
    string history_;
    string source_files_;

// Instrument-specific dimensions
    int lines_;
    int pixels_;
    short  si_[NSCANS*NDETECTORS][NPIXELS];

	int generate_sort_index ();
    int resort (sio* io);
    int resort (usio* io);
    int resort (uio* io);
    int resort (fio* io);
    int fill_the_fills (usio* io);
    int fill_the_fills (uio* io);
    int fill_the_fills (fio* io);

};

#endif /* RsViirs_H_ */
