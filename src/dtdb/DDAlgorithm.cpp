/*******************************************************************************
 *
 * NAME: DDAlgorithm.cpp
 *
 * DESCRIPTION: Base class for all algorithms.
 *
 *  Created on: October 19, 2020
 *      Author: Sam Anderson
 *
 *******************************************************************************/

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/timer/timer.hpp>
#include <algorithm>
#include <netcdf>

#include <DDataset.hpp>
#include <DDOptions.h>
#include <DDAlgorithm.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

/**************************************************************************
 * NAME: DDAlgorithm()
 *
 * DESCRIPTION: Constructor
 *
 *************************************************************************/

DDAlgorithm::DDAlgorithm()
{
	lines_ = 0;
	pixels_ = 0;
	for (size_t i=0; i< NTWL; i++) {
		rfl_[i] = 0.0;
		gasc_[i] = 0.0;
		for (size_t j=0; j<3; j++) {
			for (size_t k=0; k<3; k++) {
				rfla_[i][j][k] = 0.0;
			}
		}
	}
    lat_ = 0.0;
    lon_ = 0.0;
    solz_ = 0.0;
    senz_ = 0.0;
    raa_ = 0.0;
    height_ = 0.0;
    ws_ = 0.0;
    pwv_ = 0.0;
    oz_ = 0.0;
    ps_ = 0.0;
	month_ = 0;
	bgascorrect_ = false;
	bmaskglint_ = false;
	bmaskcloud_ = false;
	bmasksolz_ = false;
	bmasksenz_ = false;
	btest_ = false ;
	cloud_mask_ = DFILL_UBYTE;
	l2_flags_ = 0;
    qual_flag_ = 0;
    aerosol_type_ = 0;
    error_flag_ = 0;
    scatter_ang_ = 0.0;
    glint_ang_ = 0.0;
    sse_ = 0.0;
    fmf_ = 0.0;
    aot_550_ = 0.0;
    ae1_ = 0.0;
    ae2_ = 0.0;
    ndv_ = 0.0;
    chlor_ = 0.0;
	for (size_t i=0; i< NOWL+1; i++) {
	    aot_[i] = 0.0;
	}
	for (size_t i=0; i< NLWL+1; i++) {
	    ssa_[i] = 0.0;
	    sr_[i] = 0.0;
	}
	threshsolz_ = 90.0;
	threshsenz_ = 90.0;
	threshglint_ = 100.0;
};

DDAlgorithm::~DDAlgorithm()
{
};

/**************************************************************************
 * NAME: read_inputs()
 *
 * DESCRIPTION: Read inputs from map
 *
 *************************************************************************/

int DDAlgorithm::get_inputs( vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap )
{
	int status = DTDB_SUCCESS;
	size_t sy = start[0];
	size_t sx = start[1];
	lon_ = static_cast<ddma<float,2>*>(imap["longitude"])->pts[sy][sx];
	lon_ = (lon_ > 180.0) ? lon_ - 360.0 : lon_;
	lon_ = (lon_ < -180.0) ? lon_ + 360.0 : lon_;
	lat_ = static_cast<ddma<float,2>*>(imap["latitude"])->pts[sy][sx];
	lat_ = (lat_ > 90.0) ? lat_ - 180.0 : lat_;
	lat_ = (lat_ < -90.0) ? lat_ + 180.0 : lat_;
	solz_ = static_cast<ddma<float,2>*>(imap["solar_zenith"])->pts[sy][sx];
	senz_ = static_cast<ddma<float,2>*>(imap["sensor_zenith"])->pts[sy][sx];
	height_ = static_cast<ddma<float,2>*>(imap["elevation"])->pts[sy][sx]/1000.0;

	float sola = static_cast<ddma<float,2>*>(imap["solar_azimuth"])->pts[sy][sx];
	sola = (sola > 180.0) ? sola - 360.0 : sola;
	sola = (sola < -180.0) ? sola + 360.0 : sola;
	float sena = static_cast<ddma<float,2>*>(imap["sensor_azimuth"])->pts[sy][sx];
	sena = (sena > 180.0) ? sena - 360.0 : sena;
	sena = (sena < -180.0) ? sena + 360.0 : sena;
	raa_ = sena - sola - 180.0;
	if (raa_ < 0.0) raa_ = -raa_;
	raa_ = (raa_ > 180.0) ? raa_ - 360.0 : raa_;
	raa_ = (raa_ < -180.0) ? raa_ + 360.0 : raa_;
	if((solz_ > 0.0) && (senz_ > 0.0) && (raa_ > 0.0)) {
		scatter_ang_ = -cos (solz_*DEGtoRAD)*cos(senz_*DEGtoRAD)
						+ sin(solz_*DEGtoRAD)*sin(senz_*DEGtoRAD)
						* cos (raa_*DEGtoRAD);
		scatter_ang_ = acos (scatter_ang_)*RADtoDEG;
	}
	ws_ = static_cast<ddma<float,2>*>(imap["windspeed"])->pts[sy][sx];
	pwv_ = static_cast<ddma<float,2>*>(imap["water_vapor"])->pts[sy][sx];
	oz_ = static_cast<ddma<float,2>*>(imap["ozone"])->pts[sy][sx];
	ps_ = static_cast<ddma<float,2>*>(imap["pressure"])->pts[sy][sx]/100.0; // in units hPa
	cloud_mask_ = static_cast<ddma<unsigned char,2>*>(imap["cloud_mask"])->pts[sy][sx];

    size_t ilmin = (sy==0) ? 1 : 0;
    size_t ilmax = (sy==count[0]-1) ? 1 : 2;
    size_t ipmin = (sx==0) ? 1 : 0;
    size_t ipmax = (sx==count[1]-1) ? 1 : 2;
	for (auto &it : DDProcess::rhot_band_names) {
		string name = (string) it.first;
		rfl_[(size_t) it.second] = DFILL_FLOAT;
        for (size_t il=0; il<3; il++) {
            for (size_t ip=0; ip<3; ip++) {
            	rfla_[(size_t)it.second][il][ip] = DFILL_FLOAT;
            }
        }
		rfl_[(size_t) it.second] = (static_cast<ddma<float,2>*>(imap[name]))->pts[sy][sx];
        for (size_t il=ilmin; il<=ilmax; il++) {
            for (size_t ip=ipmin; ip<=ipmax; ip++) {
            	rfla_[(size_t)it.second][il][ip] =
            			(static_cast<ddma<float,2>*>(imap[name]))->pts[sy-1+il][sx-1+ip];
            }
        }
	}
	l2_flags_ = 0;

	return status;
}

/**************************************************************************
 * NAME: set_outputs()
 *
 * DESCRIPTION: Set outputs
 *
 *************************************************************************/

map<string, ddata*> DDAlgorithm::set_outputs()
{
	map<string, ddata*> omap;
	int status = DTDB_SUCCESS;
	ddval<int>* pstat = new ddval<int>(dtype::INT, status);
	omap.insert({"status", pstat});

	omap.insert({"status", pstat});
	ddval<unsigned char>* puc = new ddval<unsigned char>(dtype::UBYTE,
			(unsigned char) cloud_mask_);
	omap.insert({"cloud_mask", puc});
	ddval<short>* psh = new ddval<short>(dtype::SHORT, aerosol_type_);
	omap.insert({"aerosol_type", psh});
	ddval<float>* pfl = new ddval<float>(dtype::FLOAT, ae1_);
	omap.insert({"angstrom", pfl});
//	pfl = new ddval<float>(dtype::FLOAT, aot_550_);
//	if (!bcloudmask_ && (pfl->val > 0)) {
//		pfl->val = log10(pfl->val);
//	}
//	omap.insert({"aot_550", pfl});
	pfl = new ddval<float>(dtype::FLOAT, ndv_);
	omap.insert({"ndvi", pfl});
	pfl = new ddval<float>(dtype::FLOAT, sse_);
	omap.insert({"residual_error", pfl});
	short qual = 0;
	if (pfl->val > 0) {
		qual = (short) ceil(-log10(pfl->val));
	}
	qual = (qual>3) ? 3 : (qual<0) ? 0 : qual;
	psh = new ddval<short>(dtype::SHORT, qual);
	omap.insert({"quality", psh});

	pfl = new ddval<float>(dtype::FLOAT, fmf_);
	omap.insert({"fmf_550", pfl});
	pfl = new ddval<float>(dtype::FLOAT, chlor_);
	omap.insert({"chlorophyll", pfl});
	pfl = new ddval<float>(dtype::FLOAT, glint_ang_);
	omap.insert({"glint_angle", pfl});
	pfl = new ddval<float>(dtype::FLOAT, scatter_ang_);
	omap.insert({"scattang", pfl});

	for (auto &it : DDProcess::rhot_band_names) {
		string name = (string) it.first;
		pfl = new ddval<float>(dtype::FLOAT, rfl_[(size_t)it.second]);
		omap.insert({name, pfl});
	}
	for (auto &it : DDProcess::aot_band_names) {
		string name = (string) it.first;
		pfl = new ddval<float>(dtype::FLOAT, aot_[(size_t)it.second]);
		if (!bmaskcloud_ && (pfl->val > 0)) {
			pfl->val = log10(pfl->val);
		}
		omap.insert({name, pfl});
	}
	for (auto &it : DDProcess::srf_band_names) {
		string name = (string) it.first;
		pfl = new ddval<float>(dtype::FLOAT, sr_[(size_t)it.second]);
		omap.insert({name, pfl});
		pfl = new ddval<float>(dtype::FLOAT, ssa_[(size_t)it.second]);
		string oname = "ssa" + name.substr(7);
		omap.insert({oname, pfl});
	}

// for special aot_380
// Interpolation of ocean spectral optical depth to wavelength 380
// This is a placeholder for a planned UV model
	vector<float> tba = {490.0,550.0,670.0,865.0,1240.0,1610.0,2250.0};
	vector<float> yba(begin(aot_)+1, end(aot_));
	using boost::math::interpolators::barycentric_rational;
	barycentric_rational<float> interp(move(tba), move(yba));
	float uv = interp(380.0);
	pfl = new ddval<float>(dtype::FLOAT, uv);
	omap.insert({"aot_380", pfl});

//	if (qual == 0 && qual_flag_!= DFILL_SHORT) l2_flags_ = l2_flags_ | (unsigned int) flags::PRODFAIL;
//	if (qual == 1 && qual_flag_!= DFILL_SHORT) l2_flags_ = l2_flags_ | (unsigned int) flags::PRODWARN;

	ddval<int>* pflags = new ddval<int>(dtype::INT, l2_flags_);
	omap.insert({"l2_flags", pflags});

	return omap;
}

map<string, ddata*> DDAlgorithm::set_fills()
{
	l2_flags_ = l2_flags_ | (unsigned int) flags::PRODFAIL;
    qual_flag_ = DFILL_SHORT;
    aerosol_type_ = DFILL_SHORT;
    scatter_ang_ = DFILL_FLOAT;
    glint_ang_ = DFILL_FLOAT;
    sse_ = DFILL_FLOAT;
    fmf_ = DFILL_FLOAT;
    aot_550_ = DFILL_FLOAT;
    ae1_ = DFILL_FLOAT;
    ae2_ = DFILL_FLOAT;
    ndv_ = DFILL_FLOAT;
    chlor_ = DFILL_FLOAT;
	for (size_t i=0; i< NOWL+1; i++) {
	    aot_[i] = DFILL_FLOAT;
	}
	for (size_t i=0; i< NLWL+1; i++) {
	    ssa_[i] = DFILL_FLOAT;
	    sr_[i] = DFILL_FLOAT;
	}
	return set_outputs();
};



