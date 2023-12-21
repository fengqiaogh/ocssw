
/**************************************************************************
*
* NAME: TmParamsReader.cpp
*
* DESCRIPTION: The TmParamsReader class sets up the items
*              that are in the program configuration file.
*              The class provides primitives to read the configurable
*              information and provide it to other classes as required.
*
*  Created on: April, 2018
*      Author: Sam Anderson, DT
*
**************************************************************************/

#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <libgen.h>

#include <TmParamsReader.h>
#include <TmConstants.h>

using namespace std;

//---------------------------------------------------------------------------
// Static constants
//---------------------------------------------------------------------------

const std::string INPUT_LAM         = "wavelength";
const std::string INPUT_EPS_MAX     = "axial_ratio_max";
const std::string INPUT_EPS_MIN     = "axial_ratio_min";
const std::string INPUT_EPS_NUM     = "axial_ratio_num";
const std::string INPUT_MRR_MAX     = "mrr_max";
const std::string INPUT_MRR_MIN     = "mrr_min";
const std::string INPUT_MRR_NUM     = "mrr_num";
const std::string INPUT_MRI_MAX     = "mri_max";
const std::string INPUT_MRI_MIN     = "mri_min";
const std::string INPUT_MRI_NUM     = "mri_num";
const std::string INPUT_B_MAX       = "var_max";
const std::string INPUT_B_MIN       = "var_min";
const std::string INPUT_B_NUM       = "var_num";
const std::string INPUT_NPNAX       = "rad_num";
const std::string INPUT_AXMAX       = "rad_max";
const std::string INPUT_NPNA        = "scattering_angles";
const std::string INPUT_DDELT       = "computational_accuracy";
const std::string INPUT_NDGS        = "division_points";
const std::string INPUT_NKMAX       = "quadrature_points";
const std::string INPUT_NCOEFF      = "expansion_coeffs";
const std::string INPUT_GAM         = "modified_gamma";
const std::string INPUT_RAT         = "radius_type";
const std::string INPUT_NP          = "shape_type";
const std::string INPUT_NDISTR      = "distribution_type";
const std::string OUTPUT_NC4        = "ofile";

// global variable to hold the command line parameters
static clo_optionList_t* tm_global_optionList = NULL;

void tm_set_optionList(clo_optionList_t* list) {
    tm_global_optionList = list;
}

clo_optionList_t* tm_get_optionList() {
    if(tm_global_optionList == NULL) {
        cerr << "tm_get_optionList: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    return tm_global_optionList;
}

std::string  tm_get_option(const std::string& name) {
    if(tm_global_optionList == NULL) {
        cerr << "tm_get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    clo_option_t* option = clo_findOption(tm_global_optionList, name.c_str());
    if(option == NULL || !clo_isOptionSet(option))
        return "";
    string result(clo_getOptionString(option));
    return result;
}

int  tm_get_option_int(const std::string& name) {
    if(tm_global_optionList == NULL) {
        cerr << "tm_get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    clo_option_t* option = clo_findOption(tm_global_optionList, name.c_str());
    if(option == NULL || !clo_isOptionSet(option))
        return -999;
    int result(clo_getOptionInt(option));
    return result;
}

double  tm_get_option_double(const std::string& name) {
    if(tm_global_optionList == NULL) {
        cerr << "tm_get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    clo_option_t* option = clo_findOption(tm_global_optionList, name.c_str());
    if(option == NULL || !clo_isOptionSet(option))
        return -999.9;
    float result(clo_getOptionDouble(option));
    return result;
}

double*  tm_get_option_doubles(const std::string& name, int *count) {
    if(tm_global_optionList == NULL) {
        cerr << "tm_get_option: optionList pointer needs to be set before accessing." << endl;
        exit(EXIT_FAILURE);
    }
    clo_option_t* option = clo_findOption(tm_global_optionList, name.c_str());
    if(option == NULL || !clo_isOptionSet(option))
        return nullptr;
    double* result(clo_getOptionDoubles(option, count));
    return result;
}


void tm_add_options(clo_optionList_t* list) {
    if(tm_global_optionList == NULL)
        tm_global_optionList = list;

    clo_addOption(list, OUTPUT_NC4.c_str(), CLO_TYPE_OFILE, NULL, "output t-matrix filename");
    clo_addOption(list, INPUT_RAT.c_str(), CLO_TYPE_STRING, NULL, "radius type");
    clo_addOption(list, INPUT_NP.c_str(), CLO_TYPE_STRING, NULL, "shape type");
    clo_addOption(list, INPUT_LAM.c_str(), CLO_TYPE_DOUBLE, NULL, "wavelengths");
    clo_addOption(list, INPUT_MRR_MAX.c_str(), CLO_TYPE_DOUBLE, NULL, "maximum real part index of refraction");
    clo_addOption(list, INPUT_MRR_MIN.c_str(), CLO_TYPE_DOUBLE, NULL, "minimum real part index of refraction");
    clo_addOption(list, INPUT_MRR_NUM.c_str(), CLO_TYPE_INT, NULL, "number of real part index of refraction");
    clo_addOption(list, INPUT_MRI_MAX.c_str(), CLO_TYPE_DOUBLE, NULL, "maximum imaginary part index of refraction");
    clo_addOption(list, INPUT_MRI_MIN.c_str(), CLO_TYPE_DOUBLE, NULL, "minimum imaginary part index of refraction");
    clo_addOption(list, INPUT_MRI_NUM.c_str(), CLO_TYPE_INT, NULL, "number of imaginary part index of refraction");
    clo_addOption(list, INPUT_EPS_MAX.c_str(), CLO_TYPE_DOUBLE, NULL, "maximum axial ratio");
    clo_addOption(list, INPUT_EPS_MIN.c_str(), CLO_TYPE_DOUBLE, NULL, "minimum axial ratio");
    clo_addOption(list, INPUT_EPS_NUM.c_str(), CLO_TYPE_INT, NULL, "number of axial ratio");
    clo_addOption(list, INPUT_B_MAX.c_str(), CLO_TYPE_DOUBLE, NULL, "maximum particle size variance");
    clo_addOption(list, INPUT_B_MIN.c_str(), CLO_TYPE_DOUBLE, NULL, "minimum particle size variance");
    clo_addOption(list, INPUT_B_NUM.c_str(), CLO_TYPE_INT, NULL, "number of particle size variance");
    clo_addOption(list, INPUT_AXMAX.c_str(), CLO_TYPE_DOUBLE, NULL, "maximum particle size radius");
    clo_addOption(list, INPUT_NPNAX.c_str(), CLO_TYPE_INT, NULL, "number of particle size radii");
    clo_addOption(list, INPUT_GAM.c_str(), CLO_TYPE_DOUBLE, NULL, "gamma distribution parameter");
    clo_addOption(list, INPUT_DDELT.c_str(), CLO_TYPE_DOUBLE, NULL, "computational accuracy");
    clo_addOption(list, INPUT_NDGS.c_str(), CLO_TYPE_INT, NULL, "division_points");
    clo_addOption(list, INPUT_NDISTR.c_str(), CLO_TYPE_STRING, NULL, "size distribution type");
    clo_addOption(list, INPUT_NKMAX.c_str(), CLO_TYPE_INT, NULL, "quadrature_points");
    clo_addOption(list, INPUT_NPNA.c_str(), CLO_TYPE_INT, NULL, "scattering_angles");
    clo_addOption(list, INPUT_NCOEFF.c_str(), CLO_TYPE_INT, NULL, "number of coefficients reported");

    clo_addOption(list, "pversion", CLO_TYPE_STRING, "Unspecified", "processing version string");

    string desc = "\n";
    desc += "This program can also accept a parameters file.  The par file option names take\n";
    desc += "precedence over the standard option parameters even if the standard option\n";
    desc += "is on the command line.  To over ride an option in a parameters file, use the \n";
    desc += "option name on the command line.\n\n";
    desc += "use \"-dump_options\" on the command line to see a list of ALL option names.\n";
    clo_addOption(list, "pcf_help", CLO_TYPE_HELP, NULL, desc.c_str());
}

void tm_copy_options() {
}

string tm_get_source() {
    vector<string> sourcesList;

    string source;
    
    for(vector<string>::iterator it = sourcesList.begin(); it < sourcesList.end(); it++) {
        string str = tm_get_option(*it);
        if(!str.empty()) {
            source.append(",");
            source.append(basename((char*)str.c_str()));
        }
    }
    return source;
}

string tm_get_history(int argc, char* argv[]) {
    string history = basename(argv[0]);
    for (int i=1; i<argc; i++) {
        history.append(" ");
        history.append(basename(argv[i]));
    }
    return history;
}

