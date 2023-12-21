/*****************************************************************************
 *
 * NAME:  TmConstants.h
 *
 * DESCRIPTION: Defines constants
 *
 *****************************************************************************/

#ifndef TmConstants_h
#define TmConstants_h

#include <math.h>
#include <boost/multi_array.hpp>

#define BOOST_DISABLE_ASSERTS
typedef boost::multi_array<double, 1> double_1darray;
typedef boost::multi_array<float, 1> float_1darray;
typedef boost::multi_array<short, 1> short_1darray;
typedef boost::multi_array<double, 2> double_2darray;
typedef boost::multi_array<float, 2> float_2darray;
typedef boost::multi_array<short, 2> short_2darray;
typedef boost::multi_array<double, 3> double_3darray;
typedef boost::multi_array<float, 3> float_3darray;
typedef boost::multi_array<short, 3> short_3darray;
typedef boost::multi_array<double, 4> double_4darray;
typedef boost::multi_array<float, 4> float_4darray;
typedef boost::multi_array<short, 4> short_4darray;
typedef boost::multi_array<double, 5> double_5darray;
typedef boost::multi_array<float, 5> float_5darray;
typedef boost::multi_array<short, 5> short_5darray;
typedef boost::multi_array<double, 6> double_6darray;
typedef boost::multi_array<float, 6> float_6darray;
typedef boost::multi_array<short, 6> short_6darray;
typedef boost::multi_array<double, 7> double_7darray;
typedef boost::multi_array<float, 7> float_7darray;
typedef boost::multi_array<short, 7> short_7darray;
typedef boost::multi_array<double, 8> double_8darray;
typedef boost::multi_array<float, 8> float_8darray;
typedef boost::multi_array<short, 8> short_8darray;

static const float floatfill = -999.9;
static const short shortfill = -999;
static const short filltest = -990;

#define NUM_SEASONS 4

const int TM_SUCCESS = 0;
const int TM_FAIL = 1;

// Compression

const bool bShuffleFilter = true;
const bool bDeflateFilter = true;
const int  deflateLevel = 5;

// Enumerations

enum SEASON_ENUM
{
    WINTER,
    SPRING,
    SUMMER,
    FALL,
    NEVER
};

enum PROJECT_ENUM
{
    MODIS,
    VIIRS,
    PACE,
    NOTHING
};

// numeric constants
static const long long MSECPERDAY  = 86400000000ll; //24*60*60=86400 million
static const double TAI93_TAI58_SEC  = 1104537627.0;
static const double PlancksConstant = 6.6260755e-34;
static const double SpeedOfLight =  2.9979246e+8;
static const double BoltzConstant = 1.380658e-23;
static const double Wav2 = 11.0;


const double  PIO2 =     M_PI_2;              // pi/2
const double  PIO4 =     M_PI_4;              // pi/4
const double  TREPIO2 =  3.0L*M_PI/2.0L;      // 3*pi/2
const double  TWOPI =    2.0L*M_PI;           // 2*pi

// used to convert degrees to radians
const double  DEG2RAD =  M_PI/180.0L;         // (pi/180)

// used to convert radians to degrees
const double  RAD2DEG =  180.0L/M_PI;         // (180/pi)

// used to convert degrees to arcsec
const double  DEG2ARCSEC =  3600.0L;         // seconds in a degree

// Universal Gas Constant - used to compute virtual temp.  J/(kg*K)
const float  DRYGAS = 287.05;

// Average Atmospheric Pressure at Sea Level
const double  AVG_PRESS_SEALVL = 1013.25e0;

// Gravity (m/s^2)
const double  GRAVITY = 9.80665;

// Radius of the Area Weighting Reference Sphere (m)
// This is NOT the radius of the earth and should not be used for any
// calculations other than by Area Weighting to determine if a pixel
// is near the pole.
const double  EARTH_RADIUS_METERS = 6371007.181;

// The following constant definitions are the WGS84 Earth Ellipsoid
// constants.  The main reference for the WGS84 ellipsoid is
// NIMA Physical Geodesy web page: 164.214.2.59/GandG/wgs-84/egm96.htm
// See the updated page with the link to the NIMA TR8350.2 document:
// http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350_2.html
// The Flattening Factor(f) is computed from f = 1/298.257223563.
const double  EQUAT_RAD = 6.37813700000000e+6; // Equatoral rad., meters WGS84
const double  POLAR_RAD = 6.35675231424518e+6; // Polar radius, meters WGS84
const double  ECCEN_SQ = 6.69437999014132e-3;  // Eccentricity Squared WGS84
const double  FLATFAC = 3.35281066474748071e-3;   // Flattening Factor WGS84

// Earth Gravitational Parameter mu (G*M) in m^3 per s^2 from WGS84.
// The central term in the Earth's gravitational field (GM) is known with
// much greater accuracy than either 'G', the universal gravitational
// constant, or 'M', the mass of the Earth.  The refined value accounting
// for the mass of the atmosphere is (3986004.418 +/- 0.008) e+8 m^3/s^2.
// GPS OCS applications take advantage of the improved value, however,
// Section 3.2.3.2 and the ICD-GPS-200 recommends the original WGS 84 GM
// value to avoid introduction of error to the GPS user.
const double  EARTH_GRAV_mu = 3.986005000e+14;

// USAF Orbit Analyst Manuals, circa 1978.  1 - eccen_sqr
const double  DETIC2CENTRIC = 9.93305620009859e-1;
const double  CENTRIC2DETIC = 1.00673949674228e+0;

// The following constant definitions are for time conversions

const double  TAI2IET      = 1.0e+06;         // Conversion factor
const double  MIN_IN_HOUR  = 60.0e+0;         // Number of minutes in an hour
const double  SEC_IN_HOUR  = 3600.0e+0;       // Number of seconds in an hour
const double  MJD_CONV_FAC = 2.4000005e+6;    // Factor to convert AJD to MJD
const double  SEC_IN_DAY   = 8.64e+04;        // Number of seconds in a day
const double  UJD58        = 2.43620450e+06;  // Jan 1 1958  UJD format
const double  JAN012030    = 2.272147232e+09; // Jan 1 2030  TAI format
const double  TJD_CONV_FAC = 32.184e+0;       // Factor to convert TAI to TJD
const double  DEG_IN_HOUR  = 15.0e+0;         // Number of degrees in an hour

// The following constant definitions are for polarstereographic dataset
const double  MINUS30      = -0.523598775598299e0;  // -30 degrees in radians
const double  PLUS30       =  0.523598775598299e0;  // 30 degrees in radians

// The following constant definitions are for nwp ancillary granulation
// declare constant for calculation of water vapor mixing ratio (r)
const double GAS   = 621.97;  //-- ratio of the molecular weight
//-- of water vapor to dry air
//-- in units of grams/kilogram

// declare just a few of the more popular of the twenty SI prefixes
const double MICRO =   0.000001;  //-- scale by 1/1000000th
const double MILLI =   0.001;     //-- scale by 1/1000th
const double CENTI =   0.01;      //-- scale by 1/100th
const double DECI  =   0.1;       //-- scale by 1/10th
const double DEKA  =   10.0;      //-- scale by 10x
const double HECTO =   100.0;     //-- scale by 100x

// Kelvin/Celsius conversion factor
const double TCOEFF=273.15;

// Constant used to generate surface reflectance
// multiplier to convert pascal to atmospheres (1 atm = 101325 pascal)
const float PRESS_CONV = 1.0 / 101325.0;

// Standard Atmosphere Surface Pressure
const double STDPSL=1013.0;

// Moist air adiabatic lapse rate is 6.5 K/Km (equivalent to 6.5 C/Km)
// Converted value would be .0065 C/m
const double MOIST_AIR_LAPSE_RATE = 6.5/1000;

#endif
