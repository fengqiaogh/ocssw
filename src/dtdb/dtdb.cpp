/**************************************************************************
*
* NAME: dtdb.cpp
*
* DESCRIPTION: Source file for deep blue and dark target executable.
*
* REFERENCES:
*
* REVISION HISTORY:
*   DATE:       PR#         AUTHOR         Description
*   --------    ------      --------       -----------------
*   01-12-2020				S. Anderson
*
* NOTES (MISCELLANEOUS) SECTION:
* none
*
**************************************************************************/

#include <genutils.h>
#include <DDOptions.h>
#include <DDProcess.h>

//-----------------------------------------------------------------------------
//
// Main Function
//
//-----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    int status = 0;

    clo_optionList_t* list = clo_createList();
    //clo_setEnablePositionOptions(1);

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
        
    DDProcess* gran = new DDProcess();
    gran->history_ = get_history(argc, argv);
    gran->source_files_ = get_source();

	status = gran->initialize();
	if (status != DTDB_SUCCESS) {
	   delete gran;
	   std::cerr << "Main:: Initialization failure" << std::endl;
	   exit(EXIT_FAILURE);
	}

	status = gran->process();
	if (status != DTDB_SUCCESS) {
	  delete gran;
	  std::cerr << "Main:: Processing failure" << std::endl;
	  exit(EXIT_FAILURE);
	}

	std::cerr << "Main:: Processing complete for  "
			<< get_option("ofile") << "\n" << std::endl;

    delete gran;

    return(status);
}

