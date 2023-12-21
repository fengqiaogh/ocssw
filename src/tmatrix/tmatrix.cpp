/**************************************************************************
*
* NAME: tmatrix.cpp
*
* DESCRIPTION: Source file for T-Matrix calculation
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

#include <tmatrix.h>

#include <genutils.h>
#include <TmParamsReader.h>
#include <TmProcess.h>

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
    tmatrix_init_options(list, "2.0");
    if(argc == 1) {
        clo_printUsage(list);
        exit(1);
    }
    tmatrix_read_options(list, argc, argv);

    if(clo_getBool(list, "verbose"))
        want_verbose = 1;
    else
        want_verbose = 0;

    TmProcess* generator = new TmProcess();
    generator->setHistory(tm_get_history(argc, argv));
    
    status = generator->initialize();

    if ( status != TM_SUCCESS ) {
        return(status);
    }

    status = generator->compute_dtdb();

    if ( status != TM_SUCCESS ) {
        return(status);
    }

//    status = generator->write_scatter_lut();

   if ( status != TM_SUCCESS ) {
	   std::cerr << "WARNING!!! NetCDF4 Scattering Matrix LUT NOT created" << std::endl;
   }
   else {
	   std::cerr << "NetCDF4 Scattering Matrix LUT created" << std::endl;
   }

   delete generator;

   std::cerr << "\nAll NetCDF4 LUTs created. \n" << std::endl;

   return(status);
}

