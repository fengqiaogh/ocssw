/**************************************************************************
*
* NAME: generate_db_luts.cpp
*
* DESCRIPTION: Source file for generating netCDF4 LUTs from Deep blue algorithm
* LUTs.
*
* REFERENCES:
*
* REVISION HISTORY:
*   DATE:       PR#         AUTHOR         Description
*   --------    ------      --------       -----------------
*   05-02-2018				S. Anderson
*
* NOTES (MISCELLANEOUS) SECTION:
* none
*
**************************************************************************/

#include <genutils.h>
#include <DDOptions.h>
#include "deepblue/DbLutNetcdf.h"

//-----------------------------------------------------------------------------
//
// Main Function
//
//-----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    int status = 0;

    clo_optionList_t* list = clo_createList();
//    clo_setEnablePositionOptions(1);

    //char softwareVersion[200];
    //sprintf(softwareVersion, "%d.%d.%d-r%d", L3MAPGEN_VERSION_MAJOR, L3MAPGEN_VERSION_MINOR,
    //        L3MAPGEN_VERSION_PATCH_LEVEL, SVN_REVISION);
    init_options(list, "2.0");
    if(argc == 1) {
        clo_printUsage(list);
        exit(1);
    }
    read_options(list, argc, argv);

    if(clo_getBool(list, "verbose"))
        want_verbose = 1;
    else
        want_verbose = 0;

    DbLutNetcdf* generator = new DbLutNetcdf();
    generator->setHistory(get_history(argc, argv));
    
    status = generator->initialize();

   if ( status != DTDB_SUCCESS ) {
	   return(status);
   }

   status = generator->create_db_nc4_lut();

   if ( status != DTDB_SUCCESS ) {
	   std::cerr << "\nWARNING!!! Only "<< status <<" of 15 LUTs created" << std::endl;
	   status = DTDB_FAIL;
   }
   else {
	   std::cerr << "\nNetCDF4 Deep Blue LUTs created successfully" << std::endl;
   }

   delete generator;

   return(status);
}

