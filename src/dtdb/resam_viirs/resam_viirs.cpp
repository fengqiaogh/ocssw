/**************************************************************************
 *
 * NAME: resam_viirs.cpp
 *
 * DESCRIPTION: Application to resample VIIRS GEO and resample and
 * fill bowtie regions of L1B files for easier L2 processing.
 *
 *
 * REFERENCES:
 *
 * REVISION HISTORY:
 *   DATE:       PR#         AUTHOR         Description
 *   --------    ------      --------       -----------------
 *   o2-17-2019				S. Anderson
 *
 * NOTES (MISCELLANEOUS) SECTION:
 * none
 *
 **************************************************************************/

#include <iostream>
#include <libgen.h>
#include <DDOptions.h>
#include "resam_viirs/RsViirs.h"

//-----------------------------------------------------------------------------
//
// Main Function
//
//-----------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    int status = RS_SUCCESS;

    clo_optionList_t* list = clo_createList();
    clo_setEnablePositionOptions(1);

    //char softwareVersion[200];
    //sprintf(softwareVersion, "%d.%d.%d-r%d", L3MAPGEN_VERSION_MAJOR,
    //L3MAPGEN_VERSION_MINOR, L3MAPGEN_VERSION_PATCH_LEVEL, SVN_REVISION);
    init_options(list, "3.1");
    if (argc == 1) {
        clo_printUsage(list);
        exit(1);
    }
    read_options(list, argc, argv);

    RsViirs* rs = new RsViirs();
    rs->history_ = get_history(argc, argv);
    rs->source_files_ = get_source();

    status = rs->process();
    if (status != RS_SUCCESS) {
        std::cerr << "Main:: Resampling processing failure for " << argv[1]
                << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "Main:: Resampling complete for "
            << get_option(INPUT_L1B) << std::endl;

    delete rs;

    return (status);
}

