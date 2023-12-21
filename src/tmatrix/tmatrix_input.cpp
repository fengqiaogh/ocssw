#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <tmatrix.h>

static const char *help_optionKeys[] = {
        "-help",
        "-version",
        "-dump_options",
        "-dump_options_paramfile",
        "-dump_options_xmlfile",
        "par",
        "ofile",
        "verbose",
        "pversion",
        "pcf_help",
        NULL
};

// should the par file processing descend into other par files
static int enableFileDescending = 1;

/** CLO callback function for the "par" option.  Loads the parameter file
 into the list that is stored in option->cb_data */
void par_option_cb(struct clo_option_t *option) {
    if (enableFileDescending)
        clo_readFile((clo_optionList_t*) option->cb_data, option->valStr);
}

/** add all of the accepted command line options to list */
void tmatrix_init_options(clo_optionList_t* list, const char* softwareVersion) {
    char tmpStr[2048];
    clo_option_t* option;

    clo_setVersion2("tmatrix", softwareVersion);

    sprintf(tmpStr, "Usage: tmatrix argument-list\n\n");

    strcat(tmpStr, "  This program generates a T-Matrix product.\n\n");

    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the command line, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with command line over-riding.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    // add the par option and add the callback function
//    option = clo_addOption(list, "par", CLO_TYPE_STRING, NULL, "input parameter file");
//    option->cb_data = (void*) list;
//    option->cb = par_option_cb;

    tm_add_options(list);

    option = clo_addOption(list, "verbose", CLO_TYPE_BOOL, "false", "turn on verbose output");
    clo_addOptionAlias(option, "v");
        
    clo_setSelectOptionKeys((char**)help_optionKeys);
}


/* 
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - load the main program defaults file
    - load the command line (including specified par files)
    - re-load the command line disabling file descending so command
       line arguments will over ride

 */
void tmatrix_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        exit(EXIT_FAILURE);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);

    enableFileDescending = 1;

    // load program defaults
    sprintf(tmpStr, "%s/common/tmatrix_defaults.par", dataRoot);
    clo_readFile(list, tmpStr);

    // read all arguments including descending par files
    clo_readArgs(list, argc, argv);

    // if a file was descended, make sure the command args over-rides
    enableFileDescending = 0;
    clo_readArgs(list, argc, argv);
   
    // get sensor directory for this ifile
    //const char* sensorDir = getSensorDirectory(clo_getString(optionList, "ifile"));
    //const char* sensorDir = "viirsn";
        
    // load the sensor specific defaults file
    //sprintf(tmpStr, "%s/%s/instrument_defaults.par", dataRoot, sensorDir);
    //enableFileDescending = 1;
    //clo_readFile(list, tmpStr);
    
    // read all arguments including descending par files
    clo_readArgs(list, argc, argv);

    // enable the dump option
    clo_setEnableDumpOptions(1);

    // if a file was descended, make sure the command args over-rides
    enableFileDescending = 0;
    clo_readArgs(list, argc, argv);
    enableFileDescending = 1;

    // handle PCF file on command line
    if(clo_getPositionNumOptions(list) > 1) {
        fprintf(stderr, "-E- Too many command line parameters.  Only one par file name allowed.\n");
        exit(EXIT_FAILURE);
    }

    if(clo_getPositionNumOptions(list) == 1) {
        clo_readFile(list, clo_getPositionString(list, 0));
    }    
    
    tm_copy_options();
}
