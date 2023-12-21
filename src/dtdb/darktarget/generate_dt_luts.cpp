/**************************************************************************
*
* NAME: generate_dt_luts.cpp
*
* DESCRIPTION: Source file for generating netCDF4 LUTs from Darktarget algorithm
* LUTs.
*
*
* REFERENCES:
*
* REVISION HISTORY:
*   DATE:       PR#         AUTHOR         Description
*   --------    ------      --------       -----------------
*   05-02-2017				S. Anderson
*
* NOTES (MISCELLANEOUS) SECTION:
* none
*
**************************************************************************/

#include <genutils.h>
#include <DDOptions.h>
#include "darktarget/DtLutNetcdf.h"

//-----------------------------------------------------------------------------
//
// Main Function
//
//-----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    int status = 0;

    clo_optionList_t* list = clo_createList();
    clo_setEnablePositionOptions(1);

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

    DtLutNetcdf* generator = new DtLutNetcdf();
    generator->setHistory(get_history(argc, argv));
    
    status = generator->initialize();

   if ( status != DD_SUCCESS ) {
	   return(status);
   }

   status = generator->create_dt_nc4_lut();

   if ( status != DD_SUCCESS ) {
	   std::cerr << "WARNING!!! NetCDF4 Darktarget Aerosol LUT NOT created" << std::endl;
   }
   else {
	   std::cerr << "NetCDF4 Darktarget Aerosol LUT created" << std::endl;
   }

   delete generator;

   std::cerr << "\nAll NetCDF4 LUTs created. \n" << std::endl;

   return(status);
}

