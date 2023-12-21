/*
 * L1BAlgorithm.h
 *
 *  Created on: April 3, 2021
 *      Author: ssander2
 */

#ifndef INCLUDE_L1BALGORITHM_H_
#define INCLUDE_L1BALGORITHM_H_

#include <DDAlgorithm.h>


class OciAlg : public DDAlgorithm
{
public:
/**
 *  Class constructor/destructor
 */
	OciAlg();
	~OciAlg ();

/**
 *  Initialize Input data
 */
	int initialize( map<string, ddata*> imap );

/**
 *  Compute viirs to oci interpolated reflectances
 */

	map<string, ddata*> process(vector<size_t> start, vector<size_t> count,
			map<string, ddata*> imap);

};


class ViirsAlg : public DDAlgorithm
{
public:

/**
 *  Class constructor/destructor
 */
	ViirsAlg();
	~ViirsAlg ();

/**
 *  Initialize Input data
 */
	int initialize( map<string, ddata*> imap );

/**
 *  Return VIIRS reflectances
 */
map<string, ddata*> process(vector<size_t> start, vector<size_t> count,
		map<string, ddata*> imap);

};

#endif
