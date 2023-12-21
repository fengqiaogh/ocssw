/*******************************************************************************
 *
 * NAME: DDProcess.h
 *
 * DESCRIPTION: Declaration of the DDProcess object class.
 *
 *  Created on: August 25, 2020
 *      Author: Sam Anderson, NASA Ocean Color
 *
 *******************************************************************************/

#ifndef DDProcess_H_
#define DDProcess_H_

#include <pugixml.hpp>
#include <math.h>
#include <map>
#include <vector>
#include <netcdf>
#include <DDataset.hpp>

using namespace std;

// Enumerations

enum class rhot_band
{
    W410,           // m01
    W445,           // m02
    W490,           // m03
    W550,           // m04
    W670,           // m05
    W865,           // m07
    W1240,          // m08
    W1380,          // m09
    W1610,          // m10
    W2250,          // m11
    W8550,          // m14
    W11000,         // m15
    W12000,         // m16
    MAXRHOTBAND
};

enum class aot_band
{
    W410,           // m01
    W490,           // m03
    W550,           // m04
    W670,           // m05
    W865,           // m07
    W1240,          // m08
    W1610,          // m10
    W2250,          // m11
    WMAXAOTBAND
};

enum class srf_band
{
    W410,           // m01
    W490,           // m03
    W670,           // m05
    W2250,          // m11
    WMAXSRFBAND
};


enum class SEASON
{
    WINTER,
    SPRING,
    SUMMER,
    FALL,
    NEVER
};

enum class SENSOR
{
    VIIRS,
    POCI,
    NOTHING
};

enum class LW
{
    OCEAN,
    LAND,
    COAST
};

enum class ALGO
{
    DARKTARGET,
    DEEPBLUE,
    NOTHING
};

static constexpr int P_LEVELS = 26;
static constexpr int M_LEVELS = 21;
static constexpr int WIND_LUT_ENTRIES = 4;

static constexpr int   NOWL = 7;
static constexpr int   NLWL = 3;
static constexpr int   NTWL = 10;

static constexpr int 	 DTDB_SUCCESS = 0;
static constexpr int 	 DTDB_FAIL = 1;
static constexpr int 	 NUM_SEASONS = 4;
static constexpr double  DEGtoRAD =  M_PI/180.0L;
static constexpr double  RADtoDEG =  180.0L/M_PI;

// Compression
static constexpr  bool bShuffleFilter = true;
static constexpr  bool bDeflateFilter = true;
static constexpr  int deflateLevel = 5;

class DDSensor;
class DDAncillary;
class DDAlgorithm;
class DDProcess
{
public:

	string config_file_;
	string history_;
	string source_files_;

// Instrument-specific data

    string 		title_;
    bool 		bday_;
    bool 		bmaskcloud_;
    bool 		bmaskglint_;
    bool 		bmasksolz_;
    bool 		bmasksenz_;
    bool 		bgascorrect_;
	SENSOR		instrument_;
	size_t 		lines_;
	size_t 		pixels_;
	size_t 		lprw_;
    ALGO     	alg_;
    string		algStr_;

	/**
	 *  Class constructor / destructor
	 */
	DDProcess();
	virtual ~DDProcess ();

	/**
	 *  Initialize Input data
	 */
	int initialize();

	/**
	 *  Apply algorithm to data, produce product
	 */
	virtual int process();

	/**
	 *  Wavelength band names processed
	 */
	static const map<string, rhot_band> rhot_band_names;
	static const map<string, aot_band> aot_band_names;
	static const map<string, srf_band> srf_band_names;


protected:

    map<int,NcDim> dmap_;

    /**
     *  Objects
     */
    DDSensor*  		psensor_;
    DDAncillary*  	pancillary_;
    DDAlgorithm*  	pl_;
    DDAlgorithm*  	po_;

	vector<string> out_products;
	static const vector<string> out_groups;
	static const map<string, dtype> input_names;
	static const map<string, dtype> dtdb_names;
	static const map<string, string> odps2dtdb;

    /**
     *  Gather granule data essential for initialization.
     */
    int query_granule();

	/**
	*  Initialize output data
	*/
    int create_nc4( map<string, ddata*> imap );

	/**
	*  Write attributes to nc4 dataset from product xml node
	*/
    int write_nc4_attributes(string name, pugi::xml_node products,
    		NcVar var, ddata* dsp);

   /**
     *  Write a single line to the output file.
     */
    int write_nc4( map<string, ddata*> omap );

    /**
     *  Write decimated datasets to output file.
     */

	int write_decimated(vector<string> names, size_t& nwin);
};

#endif /* DDProcess_H_ */
