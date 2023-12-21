/*
 * DbMask.h
 *
 *  Created on: Jul 15, 2018
 *      Author: ssander2
 */

#ifndef INCLUDE_DBMASK_H_
#define INCLUDE_DBMASK_H_

#include "deepblue/DbAlgorithm.h"

class DDProcess;
class DbAlgLand;
class DbAlgOcean;

class DbCloudMaskLand
{
public:

/**
 *  Class constructor / destructor
 */
    DbCloudMaskLand();

    DbCloudMaskLand(DbAlgLand* proc);
    DbCloudMaskLand(DDProcess* granule, DbAlgLand* proc);

    ~DbCloudMaskLand();

/**
 *  Compute the 2 stages of cloud mask calculations
 */
    int compute_1( const size_t iy, const size_t ix, short& mask,
                   short& snow1, short& snow2);
    int compute_2( const size_t iy, const size_t ix, short& mask,
                   const short snow2 );

protected:

    DbAlgLand* p_;
};

class DbSmokeMask
{
public:

/**
 *  Class constructor / destructor
 */
    DbSmokeMask();
    DbSmokeMask(DbAlgLand* proc);
    ~DbSmokeMask ();

/**
 *  Compute the cloud mask
 */
    int compute( const size_t iy, const size_t ix, short& mask );

protected:

    DbAlgLand* p_;
};


class DbPyrocbMask
{
public:

/**
 *  Class constructor / destructor
 */
    DbPyrocbMask();
    DbPyrocbMask(DbAlgLand* proc);
    ~DbPyrocbMask ();

/**
 *  Compute the cloud mask
 */
    int compute( const size_t iy, const size_t ix, short& mask );

protected:

    DbAlgLand* p_;

};


class DbHighAltSmokeMask
{
public:

/**
 *  Class constructor / destructor
 */
    DbHighAltSmokeMask();
    DbHighAltSmokeMask(DbAlgLand* proc);
    ~DbHighAltSmokeMask ();

/**
 *  Compute the cloud mask
 */

    int compute( const size_t iy, const size_t ix, short& mask );

protected:

    DbAlgLand* p_;
};

class DbCloudMaskOcean
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
    DbCloudMaskOcean();
    DbCloudMaskOcean(DbAlgOcean* process);
    ~DbCloudMaskOcean ();

/**
 *  Compute the cloud mask
 */

    int compute( short& mask );

protected:

    DbAlgOcean* p_;

};

#endif /* INCLUDE_DBMASK_H_ */
