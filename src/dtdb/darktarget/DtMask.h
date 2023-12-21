/*
 * DtMask.h
 *
 *  Created on: Nov 3, 2016
 *      Author: ssander2
 */

#ifndef INCLUDE_DTMASK_H_
#define INCLUDE_DTMASK_H_

#include "darktarget/DtAlgorithm.h"

class DtAlgorithm;
class DtAlgLand;
class DtAlgOcean;

class DtCloudMaskOcean
{
public:
    // algorithm thresholds
    static constexpr float M01_STDV_THOLD  = 0.0025;
    static constexpr float M08_STDV_THOLD  = 0.0025;
    static constexpr float M08_HILAT_STDV_THOLD  = 0.001;
    static constexpr float M09_THOLD   = 0.004;
    static constexpr float M03_THOLD   = 0.11;

    static constexpr float THRHLD550_STD = 0.0025;
    static constexpr float THRHLD470 = 0.4;
    static constexpr float THRHLD138 = 0.03;
    static constexpr float THRHLD138_LOW = 0.005;
    static constexpr float THRHLD_DUST = 0.75;
    static constexpr float THRHLD_CIRRUS = 0.3;

/**
 *  Class constructor/destructor
 */
	DtCloudMaskOcean();
	DtCloudMaskOcean(DtAlgorithm* process);
	~DtCloudMaskOcean();

/**
 *  Initialize Input data
 */
	int initialize();

/**
 *  Compute the mask for one pixel
 */
    int compute( unsigned char& mask );

protected:

    DtAlgorithm* p_;

};


class DtCloudMaskLand
{
public:

	// algorithm thresholds
    static constexpr float THRHLD1380_STD = 0.003;
    static constexpr float THRHLD1380 = 0.025;
    static constexpr float THRHLD1380_LOW = 0.01;
    static constexpr float THRHLD470_STD = 0.0025;
    static constexpr float THRHLD470 = 0.4;

    unsigned char snowmask_;

/**
 *  Class constructor/destructor
 */
	DtCloudMaskLand();
	DtCloudMaskLand(DtAlgorithm* process);
	~DtCloudMaskLand ();

/**
 *  Initialize Input data
 */
	int initialize();

/**
 *  Compute the mask for one pixel
 */
    int compute( unsigned char& mask );

/**
 *  Compute snow mask
 */
    int compute_snowmask();

    unsigned char mask_cirrus_;

protected:

    DtAlgorithm* p_;
};


class DtSedimentMask
{
public:

/**
 *  Class constructor
 */
    DtSedimentMask();

/**
 *  Class constructor
 */
    DtSedimentMask(DtAlgorithm* process);

/**
 *  Class destructor
 */
    ~DtSedimentMask ();

/**
 *  Initialize Input data
 */
    int initialize();

/**
 *  Compute the mask for a single pixel
 */
    int compute(short &mask);

protected:

    DtAlgorithm* p_;
};


#endif /* INCLUDE_DTMASK_H_ */
