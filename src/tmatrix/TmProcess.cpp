/*******************************************************************************
 *
 * NAME: TmProcess.cpp
 *
 * DESCRIPTION: .
 *
 *  Created on: November 3, 2016
 *      Author: Sam Anderson
 *
 *  Modified:
 *
 *******************************************************************************/

#include <new>    		// nothrow
#include <algorithm>    // sort
#include <iostream>     // cout
#include <functional>   // bind
#include <vector>

#include <libgen.h>

#include <TmProcess.h>
#include <TmParamsReader.h>

using namespace std;

extern "C" {
int calcrand_(double* RAT, double* LAM, double* MRR, double* MRI, double* EPS,
        int* NP, double* DDELT, int* NDGS, int* NPNAX, double* AXMAX, double* B,
        double* GAM, int* NDISTR, int* NKMAX, int* NPNA, int* NCOEFF,
        double* REFF, double* VEFF, double* CEXTIN, double* CSCAT, double* WALB,
        double* ASYMM, double* ALPHA, double* BETA, double* F, double* SZ);
}

const double pace_wavelengths[NUM_PACE_BANDS] = {
    310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0,
    350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0,
    390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0,
    430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0,
    470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0,
    510.0, 515.0, 520.0, 525.0, 530.0, 535.0, 540.0, 545.0,
    550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0,
    590.0, 595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0,
    630.0, 635.0, 640.0, 645.0, 650.0, 655.0, 660.0, 665.0,
    670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0, 705.0,
    710.0, 715.0, 720.0, 725.0, 730.0, 735.0, 740.0, 745.0,
    750.0, 755.0, 760.0, 765.0, 770.0, 775.0, 780.0, 785.0,
    790.0, 795.0, 800.0, 805.0, 810.0, 815.0, 820.0, 825.0,
    830.0, 835.0, 840.0, 845.0, 850.0, 855.0, 860.0, 865.0,
    870.0, 875.0, 880.0, 885.0, 890.0, 895.0, 940.0, 1040.0,
    1250.0, 1250.0, 1378.0, 1615.0, 1615.0, 2130.0, 2260.0
 };

const double dtdb_wavelengths[NUM_DTDB_BANDS] = {
    410.0, 490.0, 550.0, 670.0, 865.0, 1040.0,
    1250.0, 1378.0, 1615.0, 1615.0, 2130.0, 2260.0
 };

/**************************************************************************
 * NAME: TmProcess()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

TmProcess::TmProcess() {
}

/**************************************************************************
 * NAME: ~TmProcess()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

TmProcess::~TmProcess() {
    delete tin_;
    delete tout_;
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * process operations.
 *
 *************************************************************************/

int TmProcess::initialize() {
    int status = TM_SUCCESS;

    tin_ = new (tmatrix_in);

    string str = tm_get_option(INPUT_RAT);
    rad_type_ = str;
    if (str.empty() || (str == "EQUAL_VOLUME")) {
        tin_->rat = EQUAL_VOLUME;
    } else if (str == "EQUAL_AREA") {
        tin_->rat = EQUAL_AREA;
    }
    str = tm_get_option(INPUT_NP);
    shape_ = str;
    if (str.empty() || (str == "SPHEROID")) {
        tin_->np = SPHEROID;
    } else if (str == "CYLINDER") {
        tin_->np = CYLINDER;
    } else if (str == "CHEBYSHEV") {
        tin_->np = CHEBYSHEV;
    }
    str = tm_get_option(INPUT_NDISTR);
    distribution_ = str;
    if (str.empty() || (str == "MODIFIED_GAMMA")) {
        tin_->ndistr = MODIFIED_GAMMA;
    } else if (str == "LOGNORMAL") {
        tin_->ndistr = LOGNORMAL;
    } else if (str == "POWERLAW") {
        tin_->ndistr = POWERLAW;
    } else if (str == "GAMMA") {
        tin_->ndistr = GAMMA;
    } else if (str == "MODIFIED_POWERLAW") {
        tin_->ndistr = MODIFIED_POWERLAW;
    }
    tin_->lam = tm_get_option_doubles(INPUT_LAM, &num_bands_);
    tin_->eps_max = tm_get_option_double(INPUT_EPS_MAX);
    tin_->eps_min = tm_get_option_double(INPUT_EPS_MIN);
    tin_->eps_num = tm_get_option_int(INPUT_EPS_NUM);
    tin_->mrr_max = tm_get_option_double(INPUT_MRR_MAX);
    tin_->mrr_min = tm_get_option_double(INPUT_MRR_MIN);
    tin_->mrr_num = tm_get_option_int(INPUT_MRR_NUM);
    tin_->mri_max = tm_get_option_double(INPUT_MRI_MAX);
    tin_->mri_min = tm_get_option_double(INPUT_MRI_MIN);
    tin_->mri_num = tm_get_option_int(INPUT_MRI_NUM);
    tin_->b_max = tm_get_option_double(INPUT_B_MAX);
    tin_->b_min = tm_get_option_double(INPUT_B_MIN);
    tin_->b_num = tm_get_option_int(INPUT_B_NUM);
    tin_->ddelt = tm_get_option_double(INPUT_DDELT);
    tin_->axmax = tm_get_option_double(INPUT_AXMAX);
    tin_->gam = tm_get_option_double(INPUT_GAM);
    tin_->ndgs = tm_get_option_int(INPUT_NDGS);
    tin_->npnax = tm_get_option_int(INPUT_NPNAX);
    tin_->nkmax = tm_get_option_int(INPUT_NKMAX);
    tin_->npna = tm_get_option_int(INPUT_NPNA);
    tin_->ncoeff = tm_get_option_int(INPUT_NCOEFF);

    int nwl = num_bands_;
    int neps = tin_->eps_num;
    int nmrr = tin_->mrr_num;
    int nmri = tin_->mri_num;
    int nb =   tin_->b_num;
    int npax = tin_->npnax;
    int nang = tin_->npna;
    tout_ = new (tmatrix_out);
    tout_->wl.resize(boost::extents[nwl]);
    tout_->eps.resize(boost::extents[neps]);
    tout_->mrr.resize(boost::extents[nmrr]);
    tout_->mri.resize(boost::extents[nmri]);
    tout_->b.resize(boost::extents[nb]);
    tout_->a.resize(boost::extents[npax]);
    tout_->angles.resize(boost::extents[nang]);
    tout_->reff.resize(boost::extents[1][neps][nmrr][nmri][nb][npax]);
    tout_->veff.resize(boost::extents[1][neps][nmrr][nmri][nb][npax]);
    tout_->cext.resize(boost::extents[1][neps][nmrr][nmri][nb][npax]);
    tout_->cscat.resize(boost::extents[1][neps][nmrr][nmri][nb][npax]);
    tout_->walb.resize(boost::extents[1][neps][nmrr][nmri][nb][npax]);
    tout_->asymm.resize(boost::extents[1][neps][nmrr][nmri][nb][npax]);
    tout_->alpha.resize(boost::extents[neps][nmrr][nmri][nb][npax][NPL][ALPHA_COEFF]);
    tout_->beta.resize(boost::extents[neps][nmrr][nmri][nb][npax][NPL][BETA_COEFF]);
    tout_->f.resize(
            boost::extents[1][neps][nmrr][nmri][nb][npax][nang][STOKES6]);

    for (int iB = 0; iB < num_bands_; iB++) {
        tout_->wl[iB] = tin_->lam[iB] * MILLI;
    }
    if (tin_->eps_num < 2) {
        tout_->eps[0] = tin_->eps_max;
        tout_->eps[0] = (tout_->eps[0]==1.0) ? tout_->eps[0]+1.0e-6 : tout_->eps[0];
    } else {
        for (int i = 0; i < tin_->eps_num; i++) {
            tout_->eps[i] = tin_->eps_min +
                    i*(tin_->eps_max - tin_->eps_min)/(tin_->eps_num-1);
            tout_->eps[i] = (tout_->eps[i]==1.0) ? tout_->eps[i]+1.0e-6 : tout_->eps[i];
        }
    }
    if (tin_->mrr_num < 2) {
        tout_->mrr[0] = tin_->mrr_max;
    } else {
        for (int i = 0; i < tin_->mrr_num; i++) {
            tout_->mrr[i] = tin_->mrr_min +
                    i*(tin_->mrr_max - tin_->mrr_min)/(tin_->mrr_num-1);
        }
    }
    if (tin_->mri_num < 2) {
        tout_->mri[0] = tin_->mri_max;
    } else {
        for (int i = 0; i < tin_->mri_num; i++) {
            tout_->mri[i] = tin_->mri_min +
                    i*(tin_->mri_max - tin_->mri_min)/(tin_->mri_num-1);
        }
    }
    if (tin_->b_num < 2) {
        tout_->b[0] = tin_->b_max;
    } else {
        for (int i = 0; i < tin_->b_num; i++) {
            tout_->b[i] = tin_->b_min +
                    i*(tin_->b_max - tin_->b_min)/(tin_->b_num-1);
        }
    }
    if (tin_->npna < 2) {
        tout_->angles[0] = 0.0;
    } else {
        for (int i = 0; i < tin_->npna; i++) {
            tout_->angles[i] = i*180.0/(tin_->npna-1);
        }
    }

    status = create_output();

    return status;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Virtual function executes process algorithm.
 *
 *************************************************************************/

int TmProcess::compute_pace() {
    int status = TM_SUCCESS;

    for (int iB = 0; iB < NUM_PACE_BANDS; iB++) {

        double lam = dtdb_wavelengths[iB] * MILLI;

        calcrand_(&tin_->rat, &lam,
                &tin_->mrr_max, &tin_->mri_max,
                &tin_->eps_max, &tin_->np,
                &tin_->ddelt, &tin_->ndgs,
                &tin_->npnax, &tin_->axmax,
                &tin_->b_max, &tin_->gam,
                &tin_->ndistr, &tin_->nkmax,
                &tin_->npna, &tin_->ncoeff,
                &tout_->reff[iB][0][0][0][0][0],
                &tout_->veff[iB][0][0][0][0][0],
                &tout_->cext[iB][0][0][0][0][0],
                &tout_->cscat[iB][0][0][0][0][0],
                &tout_->walb[iB][0][0][0][0][0],
                &tout_->asymm[iB][0][0][0][0][0],
                &tout_->alpha[0][0][0][0][0][0][0],
                &tout_->beta[0][0][0][0][0][0][0],
                &tout_->f[iB][0][0][0][0][0][0][0],
                &tout_->a[0]);

        cout << "Wavelength = " << lam << endl;
    }

    return status;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Virtual function executes process algorithm.
 *
 *************************************************************************/

int TmProcess::compute_dtdb() {
    int status = TM_SUCCESS;

    for (int iB = 0; iB < num_bands_; iB++) {
        for (int iE = 0; iE < tin_->eps_num; iE++) {
            for (int iR = 0; iR < tin_->mrr_num; iR++) {
                for (int iI = 0; iI < tin_->mri_num; iI++) {
                    for (int iV = 0; iV < tin_->b_num; iV++) {

            double ddelt = tin_->ddelt;
            int status = calcrand_(&tin_->rat, &tout_->wl[iB],
                    &tout_->mrr[iR], &tout_->mri[iI],
                    &tout_->eps[iE], &tin_->np,
                    &ddelt, &tin_->ndgs,
                    &tin_->npnax, &tin_->axmax,
                    &tout_->b[iV], &tin_->gam,
                    &tin_->ndistr, &tin_->nkmax,
                    &tin_->npna, &tin_->ncoeff,
                    &tout_->reff[0][iE][iR][iI][iV][0],
                    &tout_->veff[0][iE][iR][iI][iV][0],
                    &tout_->cext[0][iE][iR][iI][iV][0],
                    &tout_->cscat[0][iE][iR][iI][iV][0],
                    &tout_->walb[0][iE][iR][iI][iV][0],
                    &tout_->asymm[0][iE][iR][iI][iV][0],
                    &tout_->alpha[iE][iR][iI][iV][0][0][0],
                    &tout_->beta[iE][iR][iI][iV][0][0][0],
                    &tout_->f[0][iE][iR][iI][iV][0][0][0],
                    &tout_->a[0]);

            if (status != TM_SUCCESS) {
                for( int iP=0; iP<tin_->npnax; iP++) {
                    tout_->reff[0][iE][iR][iI][iV][iP] = fillvalue;
                    tout_->veff[0][iE][iR][iI][iV][iP] = fillvalue;
                    tout_->cext[0][iE][iR][iI][iV][iP] = fillvalue;
                    tout_->cscat[0][iE][iR][iI][iV][iP] = fillvalue;
                    tout_->walb[0][iE][iR][iI][iV][iP] = fillvalue;
                    tout_->asymm[0][iE][iR][iI][iV][iP] = fillvalue;
                    for( int iA=0; iA<tin_->npna; iA++) {
                        for( int iS=0; iS<STOKES6; iS++) {
                            tout_->f[0][iE][iR][iI][iV][iP][iA][iS] = fillvalue;
                        }
                    }
                }
                cout << "\n Run failed for wl = " << tout_->wl[iB] <<
                        "  axial ratio = " << tout_->eps[iE] <<
                        "  mrr = " << tout_->mrr[iR] <<
                        "  mri = " << tout_->mri[iI] <<
                        "  var = " << tout_->b[iV]   <<
                        "  max_size = " << tin_->axmax << "\n" << endl;
            } else {
                cout << "\n Run complete for wl = " << tout_->wl[iB] <<
                        "  axial ratio = " << tout_->eps[iE] <<
                        "  mrr = " << tout_->mrr[iR] <<
                        "  mri = " << tout_->mri[iI] <<
                        "  var = " << tout_->b[iV]   <<
                        "  max_size = " << tin_->axmax << "\n" << endl;
            }
                    }
                }
            }
        }
        status = write_wavelength(iB);
    }

    return status;
}

/**************************************************************************
 * NAME: create_output()
 *
 * DESCRIPTION: Create NetCDF4 file containing 4x4 scattering matrices for
 * configured wavelength and scattering angle
 *
 *************************************************************************/

int TmProcess::create_output() {
    int status = TM_SUCCESS;

    output_filepath_ = tm_get_option(OUTPUT_NC4);
    if (output_filepath_.empty()) {
        cerr << "\nTmProcess:: Failure locating output file path.\n" << endl;
        return TM_FAIL;
    }
    size_t pos = output_filepath_.rfind("/");
    prod_name_ = output_filepath_.substr(pos + 1);
    NcFile* nc_output;
    try {
        nc_output = new NcFile(output_filepath_, NcFile::replace);
    } catch (NcException& e) {
        e.what();
        cerr
                << "\nTmProcess:: Failure creating product file: "
                        + output_filepath_ + "\n" << endl;
        return TM_FAIL;
    }

    write_global_attributes(nc_output);

    NcDim wave_dim = nc_output->addDim("dim_wavelength", num_bands_);
    NcDim eps_dim = nc_output->addDim("dim_axial_ratio", tin_->eps_num);
    NcDim mrr_dim = nc_output->addDim("dim_mrr", tin_->mrr_num);
    NcDim mri_dim = nc_output->addDim("dim_mri", tin_->mri_num);
    NcDim var_dim = nc_output->addDim("dim_variance", tin_->b_num);
    NcDim size_dim = nc_output->addDim("dim_radius", tin_->npnax);
    NcDim scat_dim = nc_output->addDim("dim_angles", tin_->npna);
    NcDim stokes_dim = nc_output->addDim("dim_stokes", STOKES6);
    vector<NcDim> scat_dims;
    scat_dims.push_back(wave_dim);
    scat_dims.push_back(eps_dim);
    scat_dims.push_back(mrr_dim);
    scat_dims.push_back(mri_dim);
    scat_dims.push_back(var_dim);
    scat_dims.push_back(size_dim);
    scat_dims.push_back(scat_dim);
    scat_dims.push_back(stokes_dim);
    vector<NcDim> size_dims;
    size_dims.push_back(wave_dim);
    size_dims.push_back(eps_dim);
    size_dims.push_back(mrr_dim);
    size_dims.push_back(mri_dim);
    size_dims.push_back(var_dim);
    size_dims.push_back(size_dim);

    NcVar var;

    var = nc_output->addVar("pts_wavelength", ncDouble, wave_dim);
    var.putVar(&tout_->wl[0]);
    var = nc_output->addVar("pts_axial_ratio", ncDouble, eps_dim);
    var.putVar(&tout_->eps[0]);
    var = nc_output->addVar("pts_mrr", ncDouble, mrr_dim);
    var.putVar(&tout_->mrr[0]);
    var = nc_output->addVar("pts_mri", ncDouble, mri_dim);
    var.putVar(&tout_->mri[0]);
    var = nc_output->addVar("pts_variance", ncDouble, var_dim);
    var.putVar(&tout_->b[0]);
    var = nc_output->addVar("pts_radius", ncDouble, size_dim);
    var = nc_output->addVar("pts_angles", ncDouble, scat_dim);
    var.putVar(&tout_->angles[0]);
    var = nc_output->addVar("cext", ncDouble, size_dims);
    var = nc_output->addVar("cscat", ncDouble, size_dims);
    var = nc_output->addVar("walb", ncDouble, size_dims);
    var = nc_output->addVar("asymmetry", ncDouble, size_dims);
    var = nc_output->addVar("reff", ncDouble, size_dims);
    var = nc_output->addVar("veff", ncDouble, size_dims);
    var = nc_output->addVar("F", ncDouble, scat_dims);
    var.putAtt("long_name", "Single scattering matrix parameters (F11,F22,F33,F44,F12,F34)");
    var.putAtt("dimensions", "In order, wavelength, axial_ratio, ior_real, ior_imag, size_spread, size_radius, scattering_angles, scattering matrix parameters");
    var.setFill(true, -999.9);

    delete nc_output;
    return status;
}

/**************************************************************************
 * NAME: write_wavelength()
 *
 * DESCRIPTION: Create NetCDF4 LUT containing 4x4 scattering matrices for
 * configured wavelength and scattering angle
 *
 *************************************************************************/

int TmProcess::write_wavelength( const int wl) {
    int status = TM_SUCCESS;

    output_filepath_ = tm_get_option(OUTPUT_NC4);
    if (output_filepath_.empty()) {
        cerr << "\nTmProcess:: Failure locating output file path.\n" << endl;
        return TM_FAIL;
    }
    size_t pos = output_filepath_.rfind("/");
    prod_name_ = output_filepath_.substr(pos + 1);
    NcFile* nc_output;
    try {
        nc_output = new NcFile(output_filepath_, NcFile::write);
    } catch (NcException& e) {
        e.what();
        cerr
                << "\nTmProcess:: Failure creating product file: "
                        + output_filepath_ + "\n" << endl;
        return TM_FAIL;
    }

    vector<size_t> start;
    start.push_back(wl);
    start.push_back(0);
    start.push_back(0);
    start.push_back(0);
    start.push_back(0);
    start.push_back(0);
    vector<size_t> count;
    count.push_back(1);
    count.push_back(tin_->eps_num);
    count.push_back(tin_->mrr_num);
    count.push_back(tin_->mri_num);
    count.push_back(tin_->b_num);
    count.push_back(tin_->npnax);

    NcVar var = nc_output->getVar("pts_radius");
    var.putVar(&tout_->a[0]);

    var = nc_output->getVar("cext");
    var.putVar(start, count, &tout_->cext[0][0][0][0][0][0]);

    var = nc_output->getVar("cscat");
    var.putVar(start, count, &tout_->cscat[0][0][0][0][0][0]);

    var = nc_output->getVar("walb");
    var.putVar(start, count, &tout_->walb[0][0][0][0][0][0]);

    var = nc_output->getVar("asymmetry");
    var.putVar(start, count, &tout_->asymm[0][0][0][0][0][0]);

    var = nc_output->getVar("reff");
    var.putVar(start, count, &tout_->reff[0][0][0][0][0][0]);

    var = nc_output->getVar("veff");
    var.putVar(start, count, &tout_->veff[0][0][0][0][0][0]);

    start.push_back(0);
    start.push_back(0);
    count.push_back(tin_->npna);
    count.push_back(STOKES6);
    var = nc_output->getVar("F");
    var.putVar(start, count, &tout_->f[0][0][0][0][0][0][0][0]);

    delete nc_output;
    return status;
}

/**************************************************************************
 * write_global_attributes()
 *
 * Write global attributes to specified netCDF file ID
 *
 **************************************************************************/

int TmProcess::write_global_attributes(NcFile* nc_output) {
// global attributes
    title_ = "Particle Scattering Characteristics";
    processing_version_ = "v0.0.0";
    Conventions_ = "CF-1.6";
    institution_ =
            "NASA Goddard Space Flight Center, Ocean Biology Group";
    license_ =
            "http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
    naming_authority_ = "gov.nasa.gsfc";

    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, 80, "%Y-%m-%dT%I:%M:%SZ", timeinfo);
    date_created_ = (string) buffer;
    nc_output->putAtt("processing_version", processing_version_);

    nc_output->putAtt("Conventions", Conventions_);
    nc_output->putAtt("institution", institution_);
    nc_output->putAtt("license", license_);
    nc_output->putAtt("naming_authority", naming_authority_);
    nc_output->putAtt("date_created", date_created_);
    nc_output->putAtt("history", history_);
    nc_output->putAtt("source", source_files_);
    string versionid_ = basename((char*) tm_get_option("VersionID").c_str());
    if (!versionid_.empty()) {
        nc_output->putAtt("VersionId", versionid_);
    }

    nc_output->putAtt("title", title_);
    nc_output->putAtt("product_name", prod_name_);
    nc_output->putAtt("particle_shape", shape_);
    nc_output->putAtt("size_distribution", distribution_);
    nc_output->putAtt("radius_type", rad_type_);

    return TM_SUCCESS;
}

