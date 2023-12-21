/* 
 * File:   tmatrix.h
 * Author: S. Anderson
 *
 * Created on April 12, 2018
 */

#ifndef tmatrix_H
#define tmatrix_H

#include <clo.h>
#include <TmParamsReader.h>

void tmatrix_init_options(clo_optionList_t* list, const char* softwareVersion);
void tmatrix_read_options(clo_optionList_t* list, int argc, char* argv[]);

#endif /* tmatrix_H */

