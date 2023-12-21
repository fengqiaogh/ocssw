
/**************************************************************************
 *
 * NAME: TmParamsReader
 *
 * DESCRIPTION: Reads file path data from program
 * configuration file (pcf) and makes it available to applications via
 * various get commands.
 *
 *  Created on: April, 2018
 *      Author: Sam Anderson
 *
 **************************************************************************/

#ifndef _TmParamsReader_h_
#define _TmParamsReader_h_

#include <string>
#include <iostream>
#include <clo.h>

//
// Short name constants
//
extern const std::string INPUT_L1B;
extern const std::string OUTPUT_NC4;
extern const std::string NETCDF_LUT_PATH;

extern const std::string INPUT_RAT; // rat
extern const std::string INPUT_LAM; // lam
extern const std::string INPUT_MRR_MAX; // mrr maximum
extern const std::string INPUT_MRR_MIN; // mrr minimum
extern const std::string INPUT_MRR_NUM; // mrr number computed
extern const std::string INPUT_MRI_MAX; // mri maximum
extern const std::string INPUT_MRI_MIN; // mri minimum
extern const std::string INPUT_MRI_NUM; // mri number computed
extern const std::string INPUT_EPS_MAX; // eps maximum (axial ratio)
extern const std::string INPUT_EPS_MIN; // eps minimum (axial ratio)
extern const std::string INPUT_EPS_NUM; // eps number computed (axial ratio)
extern const std::string INPUT_B_MAX; // b maximum (size variance)
extern const std::string INPUT_B_MIN; // b minimum (size variance)
extern const std::string INPUT_B_NUM; // b number(size variance)
extern const std::string INPUT_NP; // np
extern const std::string INPUT_DDELT; // ddelt
extern const std::string INPUT_NDGS; // ndgs
extern const std::string INPUT_NPNAX; // npnax
extern const std::string INPUT_AXMAX; // axmax
extern const std::string INPUT_GAM; // gam
extern const std::string INPUT_NDISTR; // ndistr
extern const std::string INPUT_NKMAX; // nkmax
extern const std::string INPUT_NPNA; // npna (number of angles)
extern const std::string INPUT_NCOEFF; // ncoeff

void tm_set_optionList(clo_optionList_t* list);
clo_optionList_t* tm_get_optionList();
std::string tm_get_option(const std::string& name);
int  tm_get_option_int(const std::string& name);
double  tm_get_option_double(const std::string& name);
double*  tm_get_option_doubles(const std::string& name, int *count);

void tm_add_options(clo_optionList_t* list);
void tm_copy_options();
std::string tm_get_source();
std::string tm_get_history(int argc, char* argv[]);

#endif
