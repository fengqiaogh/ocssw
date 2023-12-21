/*
 * DbMaskOCI.h
 *
 *  Created on: Jul 15, 2018
 *      Author: ssander2
 */

#ifndef INCLUDE_DBMASKOCI_H_
#define INCLUDE_DBMASKOCI_H_

#include "deepblue/DbProcess.h"

class DbProcessLand;
class DbProcessOcean;

class DbCloudMaskLandOCI
{
public:

/**
 *  Class constructor / destructor
 */
    DbCloudMaskLandOCI();
    DbCloudMaskLandOCI(DbProcessLand* proc);
    ~DbCloudMaskLandOCI();

/**
 *  Compute the 2 stages of cloud mask calculations
 */
    int compute_1();
    int compute_2();

    int compute_1( const int iy, const int ix, short& mask,
                   short& snow1, short& snow2);
    int compute_2( const int iy, const int ix, short& mask,
                   const short snow2 );

protected:
    DbProcessLand* p_;

};

class DbSmokeMaskOCI
{
public:

/**
 *  Class constructor / destructor
 */
    DbSmokeMaskOCI();

    DbSmokeMaskOCI(DbProcessLand* proc);

    ~DbSmokeMaskOCI ();

/**
 *  Compute the cloud mask
 */
    int compute();

    int compute( const int iy, const int ix, short& mask );

protected:

    DbProcessLand* p_;

};


class DbPyrocbMaskOCI : public Mask
{
public:

/**
 *  Class constructor / destructor
 */
    DbPyrocbMaskOCI();

    DbPyrocbMaskOCI(DbProcessLand* proc);

    ~DbPyrocbMaskOCI ();

/**
 *  Compute the cloud mask
 */
    int compute( const int iy, const int ix, short& mask );

protected:
    DbProcessLand* p_;

};


class DbHighAltSmokeMaskOCI : public Mask
{
public:

/**
 *  Class constructor / destructor
 */
    DbHighAltSmokeMaskOCI();
    DbHighAltSmokeMaskOCI(DbProcessLand* proc);
    ~DbHighAltSmokeMaskOCI ();

/**
 *  Compute the cloud mask
 */
    int compute( const int iy, const int ix, short& mask );

protected:

    DbProcessLand* p_;
};

class DbCloudMaskOceanOCI
{
public:
    static constexpr float M01_STDV_THOLD  = 0.0025;
    static constexpr float M08_STDV_THOLD  = 0.0025;
    static constexpr float M08_HILAT_STDV_THOLD  = 0.001;
    static constexpr float M09_THOLD   = 0.004;
    static constexpr float M03_THOLD   = 0.11;

/**
 *  Class constructor / destructor
 */
    DbCloudMaskOceanOCI();
    DbCloudMaskOceanOCI(DbProcessOcean* proc);
    ~DbCloudMaskOceanOCI ();

/**
 *  Compute the cloud mask
 */

    int compute( const int iy, const int ix, short& mask );

protected:

    DbProcessOcean* p_;

};

#endif /* INCLUDE_DBMASK_H_ */
