/*******************************************************************************
 *
 * NAME: DtAlgorithm.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute cloud masks for a given DtDDProcess object class.
 *
 *  Created on: November 3, 2016
 *      Author: Sam Anderson, DT
 *
 *  Modified:
 *
 *******************************************************************************/

#include "darktarget/DtAlgorithm.h"

#include <new>    		// nothrow
#include <algorithm>    // std::sort
#include <iostream>     // std::cout
#include <functional>   // std::bind
#include <vector>

using namespace std;

const float DtAlgorithm::wind_[WIND_LUT_ENTRIES] = {
		2.0, 6.0, 10.0, 14.0 };

const float DtAlgorithm::pressure_[P_LEVELS] = {
		1000.0, 975.0, 950.0, 925.0, 900.0, 850.0, 800.0,
		750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0,
		350.0, 300.0, 250.0, 200.0, 150.0, 100.0,  70.0,  50.0,
		30.0,  20.0,  10.0 };

/**************************************************************************
 * NAME: DtAlgorithm()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

DtAlgorithm::DtAlgorithm()
{
    cloud_fraction_= 0;
    scatter_angle_= 0;
    glint_angle_= 0;
    glint_refl_= 0;
    ndvi_= 0;
	SZA_0_= 0;
    SZA_1_= 0;
    THE_0_= 0;
    THE_1_= 0;
    PHI_0_= 0;
    PHI_1_= 0;
    cmask_ = 1;
    season_ = 0;
}

/**************************************************************************
 * NAME: ~DtAlgorithm()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

DtAlgorithm::~DtAlgorithm()
{
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * process operations.
 *
 *************************************************************************/

int DtAlgorithm::initialize( map<string, ddata*> imap )
{
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

	season_= 0;
	if (month_==12 || month_==1 || month_==2 ) {
		season_=0;
	} else if (month_==3 || month_==4 || month_==5 ) {
		season_=1;
	} else if (month_==6 || month_==7 || month_==8 ) {
		season_=2;
	} else if (month_==9 || month_==10 || month_==11 ) {
		season_=3;
	}

    DtLutNetcdf* lutgen = new DtLutNetcdf();
    status = lutgen->read_gas_correction_lut( gc_lut_);
    delete lutgen;

	return status;
}


/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Virtual function executes process algorithm.
 *
 *************************************************************************/

map<string, ddata*> DtAlgorithm::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap)
{
	map<string, ddata*> omap;
	return omap;
}

/**************************************************************************
 * NAME: compute_gascorrection()
 *
 * DESCRIPTION: Compute gas correction
 *
 *************************************************************************/

int DtAlgorithm::compute_gas_correction()
{
    int status = DTDB_SUCCESS;

    float   RTrans_H2O[NUM_DT_BANDS];
    float   RTrans_O3[NUM_DT_BANDS];
    float   RTrans_CO2[NUM_DT_BANDS];
    for (int iWave = 0; iWave < NUM_DT_BANDS; iWave++) {
        RTrans_H2O[iWave] = 1.0;
        RTrans_O3[iWave] = 1.0;
        RTrans_CO2[iWave] = 1.0;
    }
//  Calculate geometric factor for 2-way transmission
    float degrad = acos(-1.0)/180.0;
    float G_factor = -1.0;
    float G_factor_H2O = -1.0;
    float G_factor_O3 = -1.0;
    float G_factor_CO2 = -1.0;
    if ((senz_ > 0.0) && (solz_ > 0.0)) {
        float c_VZ = cos(degrad*senz_);
        float c_SZ = cos(degrad*solz_);
//   Calculate spherical geometry g_factor
        float hh = 9.0;  //  9 km atmos scale height
        float r = 6371.0/hh;
        G_factor = sqrt(pow(r*c_VZ,2.0) + 2.0*r + 1.0) - r*c_VZ
                    + sqrt(pow(r*c_SZ,2.0) + 2.0*r + 1) - r*c_SZ;
//   Calculate Kasten and Young g_factors (  1. / cosz + a1*z**a2 * (a3-z)**a4 )
        G_factor_H2O =
        (1.0/(c_VZ + 0.0311*pow(senz_,0.1) * pow((92.471-senz_),-1.3814)))
        + (1.0/(c_SZ + 0.0311*pow(solz_,0.1) * pow((92.471-solz_),-1.3814)));
        G_factor_O3 =
        (1.0/(c_VZ + 268.45*pow(senz_,0.5) * pow((115.42-senz_),-3.2922)))
        + (1.0/(c_SZ + 268.45*pow(solz_,0.5) * pow((115.42-solz_),-3.2922)));
        G_factor_CO2 =
        (1.0/(c_VZ + 0.4567*pow(senz_,0.07) * pow((96.4836-senz_),-1.6970)))
        + (1.0/(c_SZ + 0.4567*pow(solz_,0.07) * pow((96.4836-solz_),-1.6970)));
    }
    float total_H20 = 0.0;
    if (pwv_ > 0.0){
       total_H20 = pwv_/10.0;
    }
    else {
       total_H20 = 0.0;
    }
//  If NCEP water is available compute Water transmission else use OptH20 from clim.
    if((total_H20 > 0.0) && (G_factor > 0.0)) {
        float logcon = log(total_H20*G_factor_H2O);
        float logcon2 = logcon*logcon;
        for (int iWave = 0; iWave < NUM_DT_BANDS; iWave++) {
            float exponent = gc_lut_.H2O_COEF[iWave][0] +
                    gc_lut_.H2O_COEF[iWave][1]*logcon +
                    gc_lut_.H2O_COEF[iWave][2]*logcon2;
            float exp1 = exp(exponent);
            RTrans_H2O[iWave]= exp(exp1);
        }
    }
    else {
        for (int iWave = 0; iWave < NUM_DT_BANDS; iWave++) {
            RTrans_H2O[iWave] = exp(gc_lut_.OPT_H2O_CLIM[iWave]*G_factor_H2O);
        }
    }
//  If NCEP Ozone is available compute Ozone transmission else use OptOzone from clim.
    float total_O3 = oz_/1000.0;
    if((total_O3 > 0.0) && (G_factor > 0.0)) {
        for (int iWave = 0; iWave < NUM_DT_BANDS; iWave++) {
            float exponent = total_O3*G_factor_O3;
            RTrans_O3[iWave] = exp(gc_lut_.O3_COEF[iWave][0] +
                    gc_lut_.O3_COEF[iWave][1]*exponent);
        }
    }
    else {
        for (int iWave = 0; iWave < NUM_DT_BANDS; iWave++) {
            RTrans_O3[iWave] = exp(gc_lut_.OPT_O3_CLIM[iWave]*G_factor_O3);
        }
    }
//  compute rest of gases from cli.
    for ( int iWave = 0; iWave < NUM_DT_BANDS; iWave++) {
        RTrans_CO2[iWave] = exp(gc_lut_.OPT_CO2_CLIM[iWave]*G_factor_CO2);
    }
//  compute total transmission
    for ( int iWave=0; iWave<NUM_DT_BANDS; iWave++) {
        gasc_[iWave] = RTrans_H2O[iWave]*RTrans_O3[iWave]*RTrans_CO2[iWave];
    }

    return status;
}

/**************************************************************************
 * NAME: mean_std()
 *
 * DESCRIPTION: Subroutine returns the mean and standard deviation of
 * input data array
 *
 *************************************************************************/

int DtAlgorithm::mean_std( int size, float* data, float& mean, float& sdev )
{
	int status = DTDB_SUCCESS;

    float s = 0.0;
    for ( int j=0; j<size; j++ ) {
    	s += data[j];
    }
    mean = s/size;
    float var = 0;
    for ( int j=0; j<size; j++ )  {
    	s = data[j] - mean;
    	float p = s*s;
    	var = var+p;
    }
    if (var > 0)  {
    	var = var/(size-1);
    	sdev = sqrt(var);
    }
    else {
    	sdev = 0;
    }

	return status;
}

/**************************************************************************
 * NAME: mean_std_weighted()
 *
 * DESCRIPTION: Subroutine returns the weighted mean and standard deviation of
 * input data array subject to a weighting mask
 *
 *************************************************************************/

int DtAlgorithm::mean_std_weighted( int size, float* data, float& mean,
		float& sdev, float* weight )
{
	int status = DTDB_SUCCESS;

	float  umean = 0.0;
	float* weighted_data;
	float  weight_sum = 0;
	weighted_data = new (nothrow) float[size];
	if (weighted_data == NULL) {
		cout << "DtAlgorithm::mean_std_weight: Memory could not be allocated";
		status = DTDB_FAIL;
	}
	else
	{
		for (int i=0; i < size; i++) {
			weighted_data[i] = data[i]*weight[i];
			weight_sum += weight[i];
		}

		mean_std( size, weighted_data, umean, sdev );

		if ( weight_sum > 0.0 ) {
			mean = umean * ((float)size) / weight_sum;
		}
	}
	delete[] weighted_data;

	return status;
}

/**************************************************************************
 * NAME: interp_extrap()
 *
 * DESCRIPTION: Perform linear interpolation or extrapolation.
 *
 *************************************************************************/

int DtAlgorithm::interp_extrap( int num, float xin,
                              float x[], float y[], float& yout)
{
	int status = DTDB_FAIL;

    if (xin <= x[0]) {
        // reverse extrapolation
        yout = y[0]+(xin-x[0])*(y[1]-y[0])/(x[1]-x[0]);
        status = DTDB_SUCCESS;
    } else if (xin >= x[num-1]) {
        // forward extrapolation
        yout = y[num-2]+(xin-x[num-2])*(y[num-1]-y[num-2])/(x[num-1]-x[num-2]);
        status = DTDB_SUCCESS;
    } else {
        // interpolation
        for (int i=0; i < num-1; i++) {
            if ((xin >= x[i]) && (xin <= x[i+1])) {
                yout=y[i]+(xin-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]);
                status = DTDB_SUCCESS;
                break;
            }
        }
    }

	return status;
}

/**************************************************************************
 * NAME: sort_index()
 *
 * DESCRIPTION: This subroutine finds the sorted index of an array.
 *
 *************************************************************************/

bool compare_pair (pair<float,int> x, pair<float,int> y)
{
	return (x.first < y.first);
}

int DtAlgorithm::sort_index( int numPts, float array[], int index[])
{
	int status = DTDB_SUCCESS;

	typedef pair<float,int> p;
	vector< pair<float,int> > index_pair;

	for (int i = 0; i < numPts; i++) {
		index_pair.push_back(make_pair(array[i], i));
	}

	stable_sort(index_pair.begin(), index_pair.end(), compare_pair);

	for (int i = 0; i < numPts; i++) {
		index[i] = index_pair[i].second;
	}

	return status;
}

/**************************************************************************
 * NAME: sort_inplace()
 *
 * DESCRIPTION: This subroutine performs sort in place.
 *
 *************************************************************************/

int DtAlgorithm::sort_inplace(int numPts, float array1[], float array2[])
{
	int status = DTDB_SUCCESS;

	typedef pair<float,float> p;
	vector< pair<float,float> > index_pair;

	for (int i = 0; i < numPts; i++) {
		index_pair.push_back(make_pair(array1[i], array2[i]));
	}

	stable_sort(index_pair.begin(), index_pair.end(), compare_pair);

	for (int i = 0; i < numPts; i++) {
		array1[i] = index_pair[i].first;
		array2[i] = index_pair[i].second;
	}

	return status;
}

/**************************************************************************
 * NAME: fit_line()
 *
 * DESCRIPTION: Linear fit data in x and y arrays to produce straight line
 * coefficients Y = A + BX.
 *
 *************************************************************************/

int DtAlgorithm::fit_line( float x[], float y[], float sig[], int ndata,
			  float& A, float& B )
{
	int status = DTDB_SUCCESS;

    float sx = 0.0;
    float sy = 0.0;
    float ST2 = 0.0;
    float WT = 0.0;
    float MWT = 0.0;
    float ss = 0.0;

    if (MWT != 0) {
        for (int i=0; i<ndata; i++) {
            WT = 1.0/pow(sig[i],2.0);
            ss += WT;
            sx += x[i]*WT;
            sy += y[i]*WT;
        }
    }
    else {
        for (int i=0; i<ndata; i++) {
            sx += x[i];
            sy += y[i];
        }
        ss = (float) ndata;
    }
    float sxOSS=sx/ss;
    if (MWT != 0) {
        for (int i=0; i<ndata; i++) {
            float T = (x[i]-sxOSS)/sig[i];
            ST2 += T*T;
            B += T*y[i]/sig[i];
        }
    }
    else {
        for (int i=0; i<ndata; i++) {
            float T = x[i]-sxOSS;
            ST2 += T*T;
            B += T*y[i];
        }
    }
    B = B/ST2;
    A = (sy-sx*B)/ss;
//    float SIGA = sqrt((1.0 + sx*sx/(ss*ST2))/ss);
//    float SIGB = sqrt(1.0/ST2);

    return status;
}

/**************************************************************************
 * NAME: set_byte()
 *
 * DESCRIPTION: Set bits in a quality control flag.
 **************************************************************************/

int DtAlgorithm::set_byte(short val, short bitPos, short& target)
{
	int status = DTDB_SUCCESS;

	short temp = 0;
	temp = val << bitPos;
	target |= temp;

    return status;
}

/**************************************************************************
 * NAME: compute_glint_angle()
 *
 * DESCRIPTION: Compute the glint angle
 *
 *************************************************************************/

int DtAlgorithm::compute_glint_refl()
{
    int status = DTDB_SUCCESS;

    glint_angle_ = 0.0;
    if((solz_> 0.0) && (senz_> 0.0) && (raa_> 0.0)) {
        glint_angle_ = cos(solz_*DEGtoRAD)*cos(senz_*DEGtoRAD)
        + sin(solz_*DEGtoRAD)*sin(senz_*DEGtoRAD)*cos(raa_*DEGtoRAD);
        glint_angle_ = acos(glint_angle_)*RADtoDEG;
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

int DtAlgorithm::compute_scatter_angle(float& scatter_angle)
{
    int status = DTDB_SUCCESS;

    scatter_angle = -cos(solz_*DEGtoRAD)*cos(senz_*DEGtoRAD)
        +sin(solz_*DEGtoRAD)*sin(senz_*DEGtoRAD)*cos(raa_*DEGtoRAD);
    scatter_angle = acos(scatter_angle)*RADtoDEG;

    return status;
}



