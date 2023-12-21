/*******************************************************************************
 *
 * NAME: DbAlgorithm.cpp
 *
 *-----------------------------------------------------------------------
 * *F90
 *
 * *DESCRIPTION:
 *
 *   DeepBlue.f90 is the main driver for the MOD_PR04_DB process.
 *   MOD_PR04_DB reads in MOD02, MOD03, and MOD06 data, processes
 *   the data using the Deep Blue aerosol algorithm, and writes
 *   the Deep Blue data fields into the MOD04 product.
 *
 * *INPUT PARAMETERS:  none
 *
 * *OUTPUT PARAMETERS:  none
 *
 *
  * *REFERENCES AND CREDITS
 *
*   Modelled on code originally developed by Mark Gray and
*   modified by J. Wei.
 *
 *
 * Created on: November 3, 2017
 *     Author: Sam Anderson, Db
 *
 * Modified:
 *
 *******************************************************************************/

#include "deepblue/DbAlgorithm.h"

#include <math.h>
#include <new>          // nothrow
#include <algorithm>    // std::sort
#include <iostream>     // std::cout
#include <functional>   // std::bind
#include <vector>

#include <DDProcess.h>

using namespace std;

const float DbAlgorithm::xzlog[10] = { 0.0000000, 0.00977964, 0.0395086,
        0.0904221, 0.164818, 0.266515, 0.401776, 0.581261, 0.824689, 1.17436 };
//-- these numbers correspond to satellite zenith angle
//-- node points of 0.,16.,30.,40.,48.,54.,58.,60. degrees
const float DbAlgorithm::xlog[8] = { 0.000000, 0.0395086, 0.143841, 0.266515,
        0.401776, 0.531394, 0.635031, 0.693147 };

const float DbAlgorithm::htab[8] = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852};
const float DbAlgorithm::ttab[8] = {288.15, 216.65, 216.65, 228.65, 270.65,
        270.65, 214.65, 186.946};
const float DbAlgorithm::ptab[8] = {1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3,
        1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6};
const float DbAlgorithm::gtab[8] = {-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0};


/**************************************************************************
 * NAME: DbAlgorithm()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DbAlgorithm::DbAlgorithm() {
}

/**************************************************************************
*NAME: ~DbAlgorithm()
 *
*DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DbAlgorithm::~DbAlgorithm() {
}

/**************************************************************************
*NAME: initialize()
 *
*DESCRIPTION: Virtual function initializes data and object classes for
*process operations.
 *
 *************************************************************************/

int DbAlgorithm::initialize( map<string, ddata*> imap ) {
    int status = DTDB_SUCCESS;

    lines_ = static_cast<ddval<int>*>(imap["num_lines"])->val;
    pixels_ = static_cast<ddval<int>*>(imap["num_pixels"])->val;
    month_ = static_cast<ddval<int>*>(imap["start_month"])->val;
	bgascorrect_ = static_cast<ddval<bool>*>(imap["bgascorrect"])->val;
	bmaskglint_ = static_cast<ddval<bool>*>(imap["bmaskglint"])->val;
	bmaskcloud_ = static_cast<ddval<bool>*>(imap["bmaskcloud"])->val;
	bmasksolz_ = static_cast<ddval<bool>*>(imap["bmasksolz"])->val;
	bmasksenz_ = static_cast<ddval<bool>*>(imap["bmasksenz"])->val;
	threshsolz_ = static_cast<ddval<float>*>(imap["threshsolz"])->val;
	threshsenz_ = static_cast<ddval<float>*>(imap["threshsenz"])->val;
	threshglint_ = static_cast<ddval<float>*>(imap["threshglint"])->val;
	btest_ = false;

    status = initialize_LUT_data( imap );
    if (status != DTDB_SUCCESS) {
        std::cerr << "DtProcess:: LUT initialization failure" << std::endl;
    }

    return status;
}

/**************************************************************************
*NAME: initialize_LUT_data()
 *
*DESCRIPTION: Read deep blue  LUTs for land.
 *
 *************************************************************************/

int DbAlgorithm::initialize_LUT_data( map<string, ddata*> imap ) {
    int status = DTDB_SUCCESS;

    size_t num_lines = (size_t) static_cast<ddval<int>*>(imap["num_lines"])->val;
    size_t num_pixels = (size_t) static_cast<ddval<int>*>(imap["num_pixels"])->val;
    ddma<float,2>* plat = static_cast<ddma<float,2>*>(imap["latitude"]);
    ddma<float,2>* plon = static_cast<ddma<float,2>*>(imap["longitude"]);
    dateline_ = 0;
    float eastedge = DFILL_FLOAT;
    float westedge = DFILL_FLOAT;
    float minlat = 90.0;
    float maxlat = -90.0;
    float minlon = 180.0;
    float maxlon = -180.0;
    for (size_t i = 0; i < num_lines; i++) {
        for (size_t j = 0; j < num_pixels; j++) {
            float lat = plat->pts[i][j];
            float lon = plon->pts[i][j];
            if (lat < DFILL_TEST || lon < DFILL_TEST)
                continue;
            minlat = (lat < minlat) ? lat : minlat;
            maxlat = (lat > maxlat) ? lat : maxlat;
            minlon = (lon < minlon) ? lon : minlon;
            maxlon = (lon > maxlon) ? lon : maxlon;
        }
    }
    if (minlon < -175.0 && maxlon > 175.0) {
        eastedge = 180.0;
        westedge = -180.0;
        for (size_t i = 0; i < num_lines; i++) {
            for (size_t j = 0; j < num_pixels; j++) {
                float lon = plon->pts[i][j];
                if (lon <= DFILL_TEST)
                    continue;
                if (lon > 0.0 && lon < eastedge)
                    eastedge = lon;
                if (lon < 0.0 && lon > westedge)
                    westedge = lon;
            }
        }
        ler_start_[0] = 10*(180 + floor(eastedge) - 1);
        if (ler_start_[0] < 0) {
            ler_start_[0] = 0;
        }
        dateline_ = 3600 - ler_start_[0];
        ler_edge_[0] = 10*(180 + (floor(westedge) + 2)) + dateline_;
    } else {
        ler_start_[0] = 10*(180 + (floor(minlon) - 1));
        if (ler_start_[0] < 0) {
            ler_start_[0] = 0;
        }
        ler_edge_[0] = 10*(180 + (floor(maxlon) + 2)) - ler_start_[0];
        if (ler_edge_[0] + ler_start_[0] > 3600) {
            ler_edge_[0] = 3600 - ler_start_[0];
        }
    }
    ler_start_[1] = 10*(90 + (floor(minlat) - 1));
    if (ler_start_[1] < 0) {
        ler_start_[1] = 0;
    }
    ler_edge_[1] = 10*(90 + (floor(maxlat) + 2)) - ler_start_[1];
    if (ler_edge_[1] + ler_start_[1] > 1800) {
        ler_edge_[1] = 1800 - ler_start_[1];
    }

    return status;
}

/**************************************************************************
*NAME: process()
 *
*DESCRIPTION: Virtual function executes process algorithm.
 *
 *************************************************************************/

map<string, ddata*> DbAlgorithm::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap)
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});
	return omap;
}

/**************************************************************************
 * NAME: locate()
 *
 * DESCRIPTION: Perform a binary search of array xx for j such
 * that x lies between y[j) and y[j+1).
 * See Numerical Recipes in Fortran, Second Edition, p.110
 * Returns  values:
 *      0 if y[1] > x and size(xx) if y[size(xx)] < x
 *     -1, 1 if x < y[1] or x > y[size(xx)] respectively.
 *
 *************************************************************************/

int DbAlgorithm::locate( int size, float y[], float x, int& status)
{
// -- start binary search.
    int jl = 0;
    int ju = size;
    while (ju-jl > 1) {
        int jm = (ju + jl) / 2;
        if (x >= y[jm]) {
          jl = jm;
        } else {
          ju = jm;
        }
    }
// -- check endpoint equality, otherwise use jl from above.
    int j = 0;
    if (x == y[0]) {
      j = 0;
    } else if (x == y[size-1]) {
      j = size-2;
    } else {
      j = jl;
    }
//   -- set status, indicate success or appropriate failure condition.
    status = 0;
    if (j >= size-1) status = 1;
    if (j < 0) status = -1;

    return j;
}

/**************************************************************************
 * NAME: compute_gas_correction()
 *
 * DESCRIPTION: Calculate gas correction factors for each band for VIIRS.
 * SZA and VZA are in degrees. ozone is in atm/cm, water vapor should be in
 * cm, and surface pressure should be in atms.
 * good status = 0, bad status = -1.
 * Coefficients from Dr. Sayer's IDL script, viirs_new_gas_correctoin.pro.
 * Derived from gas correction expressions from paper by Dr. Robert Levy:
 * "The Collection 6 MODIS aerosol products over land and ocean"
 * http://www.atmos-meas-tech.net/6/2989/2013/amt-6-2989-20 *.html
 *
 * Water vapor reduced by half below based on Tanre paper:
 * D. Tanre, B. N. Holben and Y. J. Kaufman, "Atmospheric correction against
 * algorithm for NOAA-AVHRR products: theory and application," in IEEE
 * Transactions on Geoscience and Remote Sensing, vol. 30, no. 2, pp. 231-248, Mar 1992.
 * doi: 10.1109/36.134074

 *************************************************************************/

int DbAlgorithm::compute_gas_correction()  {
    int status = DTDB_SUCCESS;

    float xvza = senz_;
    float xsza = solz_;
    float xwv = pwv_/2.0;
    float xoz = oz_;
    float xsp = ps_/1013.25;  // convert to standard atmosphere
//    float sigma = 0.0;
//    float theta = 0.0;
//    compute_pressure(height_, sigma, xsp, theta);
    float amf_coeffs[4];
    float amf_oz;          // air mass factors for ozone
    float amf_wv;          // water vapor, and constant
    float amf_cs;          // species.
    float trans_oz;        // ozone transmittance
    float trans_wv;        // water vapor transmittance
    float trans_cs;        // constant species transmittance
    float mu;

    float oz_coeffs[NTWL][2] = {{-3.09E-07, 4.71E-07},
                                        {-5.72E-05, 2.99E-06,},
                                        {-1.25E-04, 1.98E-05},
                                        {-4.75E-05, 9.08E-05},
                                        {-4.79E-05, 4.37E-05},
                                        {4.18E-07,  2.24E-06},
                                        {1.19E-07,  5.17E-26},
                                        {0.0,       0.0},
                                        {1.19E-07,  1.03E-25},
                                        {-2.61E-08, 3.28E-09}};

    float wv_coeffs[NTWL][3] = {{-9.61E+00, 9.16E-01, -2.01E-02},
                                        {-8.52E+00, 9.90E-01, -1.49E-03},
                                        {-9.65E+00, 9.87E-01, 1.80E-04},
                                        {-7.50E+00, 9.84E-01, -3.87E-03},
                                        {-7.69E+00, 9.95E-01, -1.10E-02},
                                        {-6.05E+00, 9.65E-01, -1.53E-02},
                                        {-5.16E+00, 9.59E-01, -2.67E-02},
                                        {0.0,       0.0,      0.0},
                                        {-6.43E+00, 1.02E+00, -3.60E-03},
                                        {-5.85E+00, 1.28E+00, -5.04E-03}};

    float cs_coeffs[NTWL] = {2.47E-04, 3.77E-04, 1.84E-03, 8.34E-04,
            1.44E-03, 2.45E-05, 1.19E-02, 0.0, 2.13E-02, 5.32E-02};
// -- calculate tranmittance correction factor for VZA
// ---------------------------------------------------
    for (int ib=0; ib<NTWL; ib++) {
        if (ib == 7) {      // DB uses M09, but this table does not
            gasc_[ib] = 1;
            continue;
        }
        mu = cos(xvza*DEGtoRAD);
// -- calculate ozone air mass factor and transmittance correction factor
        amf_coeffs[0] = 268.45;
        amf_coeffs[1] = 0.5;
        amf_coeffs[2] = 115.42;
        amf_coeffs[3] = -3.2922;
        amf_oz = mu + amf_coeffs[0] * (pow(xvza, amf_coeffs[1]) *
                    pow((amf_coeffs[2]-xvza),amf_coeffs[3]));
        amf_oz = 1.0 / amf_oz;
        trans_oz = exp(oz_coeffs[ib][0] + oz_coeffs[ib][1]*amf_oz*xoz);
// -- calculate water vapor air mass factor
        amf_coeffs[0] = 0.0311;
        amf_coeffs[1] = 0.1;
        amf_coeffs[2] = 92.471;
        amf_coeffs[3] = -1.3814;
        amf_wv = mu + amf_coeffs[0] * (pow(xvza,amf_coeffs[1]) *
                    pow((amf_coeffs[2]-xvza),amf_coeffs[3]));
        amf_wv = 1.0 / amf_wv;
        trans_wv = exp(exp(wv_coeffs[ib][0] + wv_coeffs[ib][1]*log(amf_wv*xwv) +
                        wv_coeffs[ib][2]*pow(log(amf_wv*xwv),2)));
// -- calculate constant species air mass factor
        amf_coeffs[0] = 0.4567;
        amf_coeffs[1] = 0.07;
        amf_coeffs[2] = 96.484;
        amf_coeffs[3] = -1.697;
        amf_cs = mu + amf_coeffs[0] * (pow(xvza,amf_coeffs[1]) *
                    pow((amf_coeffs[2]-xvza),amf_coeffs[3]));
        amf_cs = 1.0 / amf_cs;
        trans_cs = exp(cs_coeffs[ib] * amf_cs * xsp);
        gasc_[ib] = trans_oz * trans_wv * trans_cs;
// -- calculate tranmittance correction factor for SZA
// ---------------------------------------------------
// -- calculate ozone air mass factor and transmittance correction factor
        mu = cos(xsza*DEGtoRAD);
        amf_coeffs[0] = 268.45;
        amf_coeffs[1] = 0.5;
        amf_coeffs[2] = 115.42;
        amf_coeffs[3] = -3.2922;
        amf_oz = mu + amf_coeffs[0] * (pow(xsza,amf_coeffs[1]) *
                pow((amf_coeffs[2]-xsza),amf_coeffs[3]));
        amf_oz = 1.0 / amf_oz;
        trans_oz = exp(oz_coeffs[ib][0]+oz_coeffs[ib][1]*amf_oz*xoz);
//       -- calculate water vapor air mass factor
        amf_coeffs[0] = 0.0311;
        amf_coeffs[1] = 0.1;
        amf_coeffs[2] = 92.471;
        amf_coeffs[3] = -1.3814;
        amf_wv = mu + amf_coeffs[0] * (pow(xsza,amf_coeffs[1]) *
                   pow((amf_coeffs[2]-xsza),amf_coeffs[3]));
        amf_wv = 1.0 / amf_wv;
        trans_wv = exp(exp(wv_coeffs[ib][0] + wv_coeffs[ib][1]*log(amf_wv*xwv) +
                       wv_coeffs[ib][2]*pow(log(amf_wv*xwv),2)));
//       -- calculate constant species air mass factor
        amf_coeffs[0] = 0.4567;
        amf_coeffs[1] = 0.07;
        amf_coeffs[2] = 96.484;
        amf_coeffs[3] = -1.697;
        amf_cs = mu + amf_coeffs[0] * (pow(xsza,amf_coeffs[1]) *
                          pow((amf_coeffs[2]-xsza),amf_coeffs[3]));
        amf_cs = 1.0 / amf_cs;
        trans_cs = exp(cs_coeffs[ib] * amf_cs * xsp);
        gasc_[ib] = gasc_[ib] * trans_oz * trans_wv * trans_cs;
     }

  return status;
}

/**************************************************************************
* NAME: compute_pressure()
*
* DESCRIPTION: Compute the properties of the 1976 standard atmosphere to 86 km.
* AUTHOR - Ralph Carmichael, Public Domain Aeronautical Software
* NOTE - If alt > 86, the values returned will not be correct, but they will
*   not be too far removed from the correct values for density.
*   The reference document does not use the terms pressure and temperature
*   above 86 km.
*
*   alt(in) = geometric altitude, km.
*   sigma(out) = density/sea-level standard density
*   ps(out) = pressure/sea-level standard pressure
*   theta(out) = temperature/sea-level standard temperature
*
 *************************************************************************/
int DbAlgorithm::compute_pressure (float height, float& sigma,
                                     float& ps, float& theta)
{
    int status = DTDB_SUCCESS;

    const float REARTH = 6369.0; // radius of the Earth (km)
    const float GMR = 34.163195; // hydrostatic constant
    const int   NTAB=8 ;         // number of entries in the defining tables

    float h = 0.0;      // geopotential altitude (km)
    float tgrad = 0.0;  // temperature gradient and base temp of this layer
    float tbase = 0.0;  // temperature gradient and base temp of this layer
    float tlocal = 0.0; // local temperature
    float deltah = 0.0; // height above base of this layer
//============================================================================
// ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
//============================================================================
    float htab[NTAB] = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852};
    float ttab[NTAB] = {288.15, 216.65, 216.65, 228.65, 270.65,
                      270.65, 214.65, 186.946};
    float ptab[NTAB] = {1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3,
                      1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6};
    float gtab[NTAB] = {-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0};

    float alt = height/1000;
    h=alt*REARTH/(alt+REARTH);  // convert geometric to geopotential altitude

    int ip = 0;
    do {
        ip++;
    } while (h > htab[ip] && ip < NTAB);
    ip--;

// temperature ratio
    tgrad=gtab[ip];
    tbase=ttab[ip];
    deltah=h-htab[ip];
    tlocal=tbase+tgrad*deltah;
    theta=tlocal/ttab[0];
// pressure ratio
    if (tgrad == 0.0) {
        ps=ptab[ip]*exp(-GMR*deltah/tbase);
    } else {
        ps=ptab[ip]*pow((tbase/tlocal),(GMR/tgrad));
    }
// density ratio
    sigma=ps/theta;

    return status;
}

/**************************************************************************
 * NAME: compute_glint_refl()
 *
 * DESCRIPTION: Compute the glint reflectance
 *
 *************************************************************************/

int DbAlgorithm::compute_glint_refl(float& glint_angle)
{
    int status = DTDB_SUCCESS;

    glint_angle = 0.0;
    if((solz_> 0.0) && (senz_> 0.0) && (raa_> 0.0)) {
        glint_angle = cos(solz_*DEGtoRAD)*cos(senz_*DEGtoRAD)
        + sin(solz_*DEGtoRAD)*sin(senz_*DEGtoRAD)*cos(raa_*DEGtoRAD);
        glint_angle = acos(glint_angle)*RADtoDEG;
    }

    double cc = DEGtoRAD;
    float nr = 1.341;

    double zx = (sin(senz_*cc)*sin(raa_*cc))/
            (cos(solz_*cc)+cos(senz_*cc));
    double zy = (sin(solz_*cc)-sin(senz_*cc)*cos(raa_*cc))/
            (cos(solz_*cc)+cos(senz_*cc));
    double sigx = sqrt(0.003+0.00192*ws_);
    double sigy = sqrt(0.00316*ws_);
    double zeta = zx / sigx;
    double eta  = zy / sigy;
    double p = (1.0/(2.0*M_PI*sigx*sigy))*exp(-0.5*(zeta*zeta + eta*eta));
    double costwoomega = cos(senz_*cc)*cos(solz_*cc) -
                         sin(senz_*cc)*sin(solz_*cc)*cos(raa_*cc);
    double cosbeta = (cos(solz_*cc)+cos(senz_*cc))/
            (sqrt(2.0 + 2.0*costwoomega));
    double w = 0.5 * acos(costwoomega);
    double wp = asin(sin(w)/nr);
    double a1 = sin(w - wp);
    double b1 = sin(w+wp);
    double c1 = tan(w-wp);
    double d1 = tan(w+wp);
    double R = 0.5*((a1*a1)/(b1*b1)+(c1*c1)/(d1*d1));
    glint_refl_ = p*R/(4*cos(senz_*cc)*pow(cosbeta,4.0));

    return status;
}

/**************************************************************************
 * NAME: compute_scatter_angle()
 *
 * DESCRIPTION: Compute land scatter angle.
 *
 *************************************************************************/

int DbAlgorithm::compute_scatter_angle(float& scatter_angle)
{
    int status = DTDB_SUCCESS;

    scatter_angle = -cos(solz_*DEGtoRAD)*cos(senz_*DEGtoRAD)
        +sin(solz_*DEGtoRAD)*sin(senz_*DEGtoRAD)*cos(raa_*DEGtoRAD);
    scatter_angle = acos(scatter_angle)*RADtoDEG;

    return status;
}

