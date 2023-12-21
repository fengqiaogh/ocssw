/*
 * TmProcess.h
 *
 *  Created on: April, 2018
 *      Author: ssander2
 *

 * INPUT PARAMETERS:
 *
 *     RAT = 1 - particle size is specified in terms of the
 *               equal-volume-sphere radius
 *     RAT.NE.1 - particle size is specified in terms of the
 *               equal-surface-area-sphere radius
 *     NDISTR specifies the distribution of equivalent-sphere radii
 *     NDISTR = 1 - modified gamma distribution
 *          [Eq. (40) of Ref. 7]
 *              AXI=alpha
 *              B=r_c
 *              GAM=gamma
 *     NDISTR = 2 - log-normal distribution
 *          [Eq. 41) of Ref. 7]
 *              AXI=r_g
 *              B=[ln(sigma_g)]**2
 *     NDISTR = 3 - power law distribution
 *          [Eq. (42) of Ref. 7]
 *               AXI=r_eff (effective radius)
 *               B=v_eff (effective variance)
 *               Parameters R1 and R2 (see below) are calculated
 *               automatically for given AXI and B
 *     NDISTR = 4 - gamma distribution
 *          [Eq. (39) of Ref. 7]
 *               AXI=a
 *               B=b
 *     NDISTR = 5 - modified power law distribution
 *        [Eq. (24) in M. I. Mishchenko et al.,
 *        Bidirectional reflectance of flat,
 *        optically thick particulate laters: an efficient radiative
 *        transfer solution and applications to snow and soil surfaces,
 *        J. Quant. Spectrosc. Radiat. Transfer, Vol. 63, 409-432 (1999)].
 *               B=alpha
 *
 *     The code computes NPNAX size distributions of the same type
 *     and with the same values of B and GAM in one run.
 *     The parameter AXI varies from AXMAX to AXMAX/NPNAX in steps of
 *     AXMAX/NPNAX.  To compute a single size distribution, use
 *     NPNAX=1 and AXMAX equal to AXI of this size distribution.
 *
 *     R1 and R2 - minimum and maximum equivalent-sphere radii
 *          in the size distribution. They are calculated automatically
 *          for the power law distribution with given AXI and B
 *          but must be specified for other distributions
 *          after the lines
 *
 *            DO 600 IAX=1,NPNAX
 *               AXI=AXMAX-DAX*DFLOAT(IAX-1)
 *
 *          in the main program.
 *          For the modified power law distribution (NDISTR=5), the
 *          minimum radius is 0, R2 is the maximum radius,
 *          and R1 is the intermediate radius at which the
 *          n(r)=const dependence is replaced by the power law
 *          dependence.
 *
 *     NKMAX.LE.988 is such that NKMAX+2 is the
 *          number of Gaussian quadrature points used in
 *          integrating over the size distribution for particles with
 *          AXI=AXMAX.  For particles with AXI=AXMAX-AXMAX/NPNAX,
 *          AXMAX-2*AXMAX/NPNAX, etc. the number of Gaussian points
 *          linearly decreases.
 *          For the modified power law distribution, the number
 *          of integration points on the interval [0,R1] is also
 *          equal to NKMAX.
 *
 *     LAM - wavelength of light
 *     MRR and MRI - real and imaginary parts of the refractive
 *                 index (MRI.GE.0)
 *     EPS and NP - specify the shape of the particles.
 *            For spheroids NP=-1 and EPS is the ratio of the
 *                horizontal to rotational axes.  EPS is larger than
 *                1 for oblate spheroids and smaller than 1 for
 *                prolate spheroids.
 *            For cylinders NP=-2 and EPS is the ratio of the
 *                diameter to the length.
 *            For Chebyshev particles NP must be positive and
 *                is the degree of the Chebyshev polynomial, while
 *                EPS is the deformation parameter
 *                [Eq. (33) of Ref. 7].
 *     DDELT - accuracy of the computations
 *     NPNA - number of equidistant scattering angles (from 0
 *            to 180 deg) for which the scattering matrix is
 *            calculated.
 *     NDGS - parameter controlling the number of division points
 *            in computing integrals over the particle surface.
 *            For compact particles, the recommended value is 2.
 *            For highly aspherical particles larger values (3, 4,...)
 *            may be necessary to obtain convergence.
 *            The code does not check convergence over this parameter.
 *            Therefore, control comparisons of results obtained with
 *            different NDGS-values are recommended.
 *
 */

#ifndef INCLUDE_TmProcess_H_
#define INCLUDE_TmProcess_H_

#include <netcdf>
#include <TmConstants.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

static const double EQUAL_VOLUME = 1.0;
static const double EQUAL_AREA = 0.5;
static const int RADIUS_MAXIMUM = 2.0;

static const int SPHEROID = -1;
static const int CYLINDER = -2;
static const int CHEBYSHEV = 1;

static const int MODIFIED_GAMMA = 1;
static const int LOGNORMAL = 2;
static const int POWERLAW = 3;
static const int GAMMA = 4;
static const int MODIFIED_POWERLAW = 5;

struct tmatrix_in {
    double rat; // equal-volume-sphere radius (1) or equal-surface-area-sphere radius
    double* lam; // wavelength of light
    double eps_max; // maximum ratio of the horizontal to rotational axes
    double eps_min; // minimum ratio of the horizontal to rotational axes
    double eps_num; // number of ratios of the horizontal to rotational axes
    double mrr_max; // maximum real part of the refractive index
    double mrr_min; // minimum real part of the refractive index
    double mrr_num; // number of real part of the refractive index
    double mri_max; // maximum imaginary part of the refractive index
    double mri_min; // number of imaginary part of the refractive index
    double mri_num; // imaginary part of the refractive index
    double b_max;   // effective radius distribution parameter
    double b_min;   // effective radius distribution parameter
    double b_num;   // effective radius distribution parameter
    double ddelt; // accuracy of the computations
    double axmax; // AXI varies from AXMAX to AXMAX/NPNAX in steps of AXMAX/NPNAX
    double gam;   // effective variance distribution parameter
    int np;     // spheroid (-1), cylinder (-2), Chebyshev(degree)
    int ndgs;   // number of division points (default 2)
    int npnax;  // number of size distributions of the same type
    int ndistr; // specifies the distribution of equivalent-sphere radii
    int nkmax;  // number of Gaussian quadrature points used in integrating over the size distribution
    int npna;   // number of scattering angles
    int ncoeff; // number of expansion coefficients
};

static const int ALPHA_COEFF = 4;
static const int BETA_COEFF = 2;
static const int STOKES = 4;
static const int STOKES6 = 6;
static const int NPL = 801;
static const double fillvalue = -999.9;

struct tmatrix_out {
    double_1darray wl;
    double_1darray mrr;
    double_1darray mri;
    double_1darray eps;
    double_1darray a;
    double_1darray b;
    double_1darray angles;
    double_6darray reff;
    double_6darray veff;
    double_6darray cext;
    double_6darray cscat;
    double_6darray walb;
    double_6darray asymm;
    double_7darray alpha;
    double_7darray beta;
    double_8darray f;
};

static const int NUM_PACE_BANDS = 127;
static const int NUM_DTDB_BANDS = 12;

extern const double pace_wavelengths[NUM_PACE_BANDS];
extern const double dtdb_wavelengths[NUM_DTDB_BANDS];

class TmProcess
{
public:
    string  output_filepath_;

    // output file attributes

    string title_ ;
    string prod_name_ ;

    // output file global attributes

    string processing_version_;
    string Conventions_ ;
    string institution_;
    string license_ ;
    string naming_authority_;
    string date_created_ ;
    string source_files_;
    string versionid_;

    string shape_;
    string distribution_;
    string rad_type_;

/**
 *  Class constructor
 */
	TmProcess();

/**
 *  Class destructor
 */
	virtual ~TmProcess ();

/**
 *  Initialize Input data
 */
	virtual int initialize();

/**
 *  Compute scattering matrix array
 */
    virtual int compute_pace();
    virtual int compute_dtdb();

/**
 *  Write scattering matrix LUT to NetCDF4 file.
 */

    int create_output();

    int write_wavelength( const int wl );

/**
 *  Write global attributes to file..
 */

    int write_global_attributes( NcFile* nc_output );


/**
 *  Set/Get History
 */
    string history_;

    void setHistory(std::string history) {
        history_ = history;
    }
    std::string getHistory() {
        return history_;
    }


/**
 *  Public variables
 */

protected:

    int num_bands_;

    tmatrix_in*  tin_;
    tmatrix_out*  tout_;

};

#endif /* INCLUDE_TmProcess_H_ */
