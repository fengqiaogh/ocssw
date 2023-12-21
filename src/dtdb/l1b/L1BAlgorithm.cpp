/*******************************************************************************
 *
 * NAME: ViirsAlg.cpp
 *
 * DESCRIPTION: Object class that provides data structures and processes that
 * compute cloud masks for a given DtDDProcess object class.
 *
 *  Created on: November 3, 2016
 *      Author: Sam Anderson, DT
 *
 *  Modified:
 *
 *******************************************************************************/

#include <L1BAlgorithm.h>

#include <new>    		// nothrow
#include <algorithm>    // std::sort
#include <iostream>     // std::cout
#include <functional>   // std::bind
#include <vector>

using namespace std;


/**************************************************************************
 * NAME: OciAlg()
 *
 * DESCRIPTION: Class Constructor/Destructor
 *
 *************************************************************************/

OciAlg::OciAlg() {
}

OciAlg::~OciAlg() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * process operations.
 *
 *************************************************************************/

int OciAlg::initialize( map<string, ddata*> imap ) {
    int status = DT_SUCCESS;
	return status;
}

/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Virtual function executes process algorithm.
 *
 *************************************************************************/

map<string, ddata*> OciAlg::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap) {
	map<string, ddata*> omap;
	return omap;
}

/**************************************************************************
 * NAME: ViirsAlg()
 *
 * DESCRIPTION: Class Constructor/Destructor
 *
 *************************************************************************/

ViirsAlg::ViirsAlg() {
}

ViirsAlg::~ViirsAlg() {
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Virtual function initializes data and object classes for
 * process operations.
 *
 *************************************************************************/

int ViirsAlg::initialize( map<string, ddata*> imap ) {
    int status = DT_SUCCESS;
	return status;
}


/**************************************************************************
 * NAME: compute()
 *
 * DESCRIPTION: Virtual function executes process algorithm.
 *
 *************************************************************************/

map<string, ddata*> ViirsAlg::process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap) {
	map<string, ddata*> omap;
	return omap;
}


