/*
 * CldMask.h
 *
 *  Created on: April 3, 2021
 *      Author: ssander2
 */

#ifndef INCLUDE_CLDMASK_H_
#define INCLUDE_CLDMASK_H_

class DDAlgorithm;

class CldMaskOcean
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
	CldMaskOcean();
	CldMaskOcean(DDAlgorithm* process);
	~CldMaskOcean();

/**
 *  Initialize Input data
 */
	int initialize();

/**
 *  Compute the mask for one pixel
 */
    int compute( unsigned char& mask );

protected:

    DDAlgorithm* p_;

};


class CldMaskLand
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
	CldMaskLand();
	CldMaskLand(DDAlgorithm* process);
	~CldMaskLand ();

/**
 *  Initialize Input data
 */
	int initialize();

/**
 *  Compute the mask for one pixel
 */
    int compute( unsigned char& mask );

    unsigned char mask_cirrus_;

protected:

    DDAlgorithm* p_;
};


#endif /* INCLUDE_CLDMASK_H_ */
