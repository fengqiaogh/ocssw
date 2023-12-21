/**************************************************************************
 *
 * NAME: DDSensor.h
 *
 * DESCRIPTION: Declaration for the DDSensor object class.
 *
 *  Created on: May 10, 2017
 *      Author: Sam Anderson
 *
 ***************************************************************************/

#ifndef INCLUDE_INPUT_H_
#define INCLUDE_INPUT_H_

#include <vector>
#include <map>
#include <DDOptions.h>

using namespace std;
using namespace netCDF;

extern "C" {
    double esdist_(int32_t *year, int32_t *day, int32_t *msec);
}

static const int BTLUTSIZE = 65536;
static const int VALIDMAX = 65527;
static const int SOLZNIGHT = 90.0;
static const int SOLZNIGHTA = 80.0;
static const int DAYSCENE = 0;
static const int NIGHTSCENE = 1;
static const int DAYANDNIGHT = 2;
static const int UNKNOWNSCENE = 3;
static const int ASCENDING = 0;
static const int DSCENDING = 1;
static const int UNKNOWNNODE = 2;

class DDProcess;

struct btlut_vector
{
    float btlut[BTLUTSIZE];
};

// VIIRS input structures
struct vfio {
    float in[VLINES][VELEMS];
    float out[VLINES][VELEMS];
};
struct vsio {
    short in[VLINES][VELEMS];
    short out[VLINES][VELEMS];
};
struct vusio {
    unsigned short in[VLINES][VELEMS];
    unsigned short out[VLINES][VELEMS];
};

// PACE input structures
struct pfio {
    boost::multi_array<float,2> in;
    boost::multi_array<float,2> out;
};
struct psio {
    boost::multi_array<short,2> in;
    boost::multi_array<short,2> out;
};
struct pusio {
    boost::multi_array<unsigned short,2> in;
    boost::multi_array<unsigned short,2> out;
};

#define LANDMASK_LAT 43201
#define LANDMASK_LON 86401

struct LandMaskLUT
{
    double lat[LANDMASK_LAT];
    double lon[LANDMASK_LON];
    char   watermask[LANDMASK_LAT][LANDMASK_LON];
};

const static int MAXMODIS250RSB = 2;
const static int MAXMODIS500RSB = 7;
const static int MAXMODIS1KMRSB = 19;

enum MODIS_BANDS {
    MW645, MW859, MW470, MW550, MW124, MW164,
    MW213, MW412, MW443, MW488, MW531,
    MW551, MW667L, MW667H, MW678L,
    MW678H, MW748, MW869, MW905, MW936, MW940,
    MW3750, MW3959H, MW3959L, MW4050,
    MW4465, MW4515, MW1375, MW6715,
    MW7325, MW8550, MW9730, MW11030,
    MW12020, MW13335, MW13635, MW13935, MW14235,
    MAXMODISALL
};

const static float WLS  = 2.5;
const static float WLB  = 310.0;
const static float WLR  = 600.0;

enum PACE_SWIR {
    PS940, PS1040, PS1250, PD1250, PS1378, PS1615,
    PD1615, PS2130, PS2260, MAXSWIRBAND
};

/*
 * BASE CLASS
 */

class DDSensor
{
public:
/**
 *  Class constructor
 */
    DDSensor();
/**
 *  Class destructor
 */
    virtual ~DDSensor ();

/**
 *  Create l2 required data and return on map
 */
    map<string, ddata*> create(vector<size_t> start, vector<size_t> count);

/**
 *  Read input and return on map
 */
    map<string, ddata*> read(vector<size_t> start, vector<size_t> count);

/**
 *  Read Global Land Mask LUT
 */
    map<string, ddata*>  read_landmask( map<string, ddata*> gmap,
    		vector<size_t> start, vector<size_t> count );

protected:
    friend class POCI;
    friend class VIIRS;

    SENSOR  instrument_;
    bool    brad_;
    float   fsol_;

/**
 *  Noise generation utilities
 */
    unsigned long int random_seed();
    float make_noise(float sigma);

/**
 *  Store datasets into STL maps
 */
     map<string, ddata*> create_l2( map<string, ddata*> );

/**
 *  Instrument-dependent noise model
 */
    virtual float noise_model(float rfl, int iw, float snr_mult);

/**
 *  Read L1B file
 */
    virtual map<string,ddata*> read_l1b(const string& filepath,
    		vector<size_t> start, vector<size_t> count);

/**
 *  Read L1B attributes
 */
    virtual map<string,ddata*> read_l1b_attributes(NcFile* nc_input);

 /**
 *  Read geolocation file
 */
    virtual map<string,ddata*> read_geo(const string& filepath,
    		vector<size_t> start, vector<size_t> count);
};

/*
 * PACE
 */

class POCI : public DDSensor
{
public:

/**
 *  Class constructor
 */
	POCI();
/**
 *  Class destructor
 */
    ~POCI ();

protected:

/**
 *  Instrument-dependent noise model
 */
    float noise_model(float rfl, int iw, float snr_mult);

/**
 *  Read POCI L1B file
 */
    map<string,ddata*> read_l1b(const string& filepath,
    		vector<size_t> start, vector<size_t> count);

/**
 *  Read POCI L1B attributes
 */
    map<string,ddata*> read_l1b_attributes(NcFile* nc_input);

/**
 *  Read POCI geolocation file
 */
	map<string,ddata*> read_geo(const string& filepath,
			vector<size_t> start, vector<size_t> count);

	size_t bindex(float wl) {return (size_t)((wl-WLB)/WLS);}
	size_t rindex(float wl) {return (size_t)((wl-WLR)/WLS);}

};

/*
 * VIIRS
 */

class VIIRS : public DDSensor
{
public:
/**
 *  Class constructor
 */
    VIIRS();

/**
 *  Class destructor
 */
    ~VIIRS ();

protected:

/**
 *  Instrument-dependent noise model
 */
    float noise_model(float rfl, int iw, float snr_mult);

/**
 *  Read VIIRS L1B file
 */
    map<string,ddata*> read_l1b(const string& filepath,
    		vector<size_t> start, vector<size_t> count);

/**
 *  Read VIIRS L1B attributes
 */
    map<string,ddata*> read_l1b_attributes(NcFile* nc_input);

/**
 *  Read VIIRS geolocation file
 */
    map<string,ddata*> read_geo(const string& filepath,
    		vector<size_t> start, vector<size_t> count);

};

#endif /* INCLUDE_INPUT_H_ */
