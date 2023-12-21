/*
 * DbAlgorithm.h
 *
 *  Created on: May 2, 2018
 *      Author: Sam Anderson
 */

#ifndef INCLUDE_DBALGORITHM_H_
#define INCLUDE_DBALGORITHM_H_

#include <DDAlgorithm.h>
#include <DDProcess.h>
#include <DDOptions.h>
#include "deepblue/DbLutNetcdf.h"
#include "deepblue/DbMask.h"

using namespace std;

enum OBANDS_ENUM {
    O488,           // m03
    O550,           // m04
    O670,           // m05
    O865,           // m07
    O1240,          // m08
    O1610,          // m10
    O2250,          // m11
    OMAXBAND
};

enum NCBANDS_ENUM {
    NC412,           // m01
    NC488,           // m03
    NC670,           // m05
    NCMAXBAND
};

enum LBANDS_ENUM {
    OL412,           // m01
    OL488,           // m03
    OL670,           // m05
    LMAXBAND
};

enum OMODEL_ENUM {
    DUST,
    FINE,
    MARI,
    MIX,
    BEST,
    NMODEL
};


class DbCloudMaskOcean;
class DbCloudMaskLand;
class DbSmokeMask;
class DbHighAltSmokeMask;
class DbPyrocbMask;

class DbAlgorithm : public DDAlgorithm
{
	friend class DbCloudMaskLand;
	friend class DbCloudMaskOcean;
	friend class DbSmokeMask;
	friend class DbHighAltSmokeMask;
	friend class DbPyrocbMask;
public:

// algorithm constants

    static constexpr int No_byte = 6;
    static constexpr int Fmax = 150;
    static constexpr int Lmax = 210;
    static constexpr int Fbmax = 1500;
    static constexpr int Lbmax = 2100;
    static constexpr int No_byte_O = 5;
    static constexpr float delta = -0.000001;
    static constexpr float SolarZenithAngleZEPS = 84.000001;
    static constexpr float NDVI1_CUTOFF = 0.18;
    static constexpr float NDVI2_CUTOFF = 0.35;
    static const float xzlog[10];
    static const float xlog[8];
    static const float htab[8];
    static const float ttab[8];
    static const float ptab[8];
    static const float gtab[8];

/**
 *  Class constructor
 */
	DbAlgorithm();

/**
 *  Class destructor
 */
	virtual ~DbAlgorithm ();

/**
 *  Initialize Input data
 */
	virtual int initialize( map<string, ddata*> imap );

/**
 *  Return a list of product names generated by the algorithm
 */
	virtual vector<string>  get_products() { return {};};

/**
 *  Compute aerosol deep blue algorithm
 */
	virtual map<string, ddata*> process(vector<size_t> start, vector<size_t> count,
			map<string, ddata*> imap);


/**
 *  Compute gas corrections for one pixel
 */
    int compute_gas_correction();

    int   ler_start_[2];
    int   ler_edge_[2];
    int   dateline_;
    float cphi_;
    float cphir_;
    float phs_;
    float phsr_;

    float NC_[DB_NC_BANDS];

    float   scatter_angle_;
    float   glint_angle_;
    float   glint_refl_;
    float   ndvi_;

    dbTablesLUT*            mt_lut_;
    dbSurfacePressureLUT*   sp_lut_;
    dbSurfCoeffLimited*     scl_lut_;
    dbGeozoneLUT*           gz_lut_;
    dbViirsSurfReflLimited* vsr_lut_;
    dbModisSurfReflLimited* msr_lut_;
    dbOceanAerosolLUMA*     oaLut_[NDBMDL];
    dbBathymetryLUT*        bath_lut_;
    dbChlLUT*               chl_lut_;

//    dbLUT*                  db_lut_;
//    dbNvalxLUT*             nv_lut_;
//    dbViirsSwirVsVisLUT*    vsv_lut_;
//    dbRayleighLUT*          ra_lut_;
//    dbLandcoverLUT*         lc_lut_;
//    dbModisCorrectionsLUT*  mc_lut_;
//    dbModisSwirVsVisLUT*    msv_lut_;


protected:

/**
 *  Read LUT files into their respective data structures.
 */
    int initialize_LUT_data( map<string, ddata*> imap );

/**
 *  Perform binary search of array y, such that x lies between
 *  y[i] and y[i+1]
 */
    int locate(int size, float y[], float x, int& status);

/**
 *  Compute pressure from elevation data
 */
    int compute_pressure (float height, float& sigma, float& ps, float& theta);

/**
 *  Compute glint reflectance
 */
    int compute_glint_refl(float& glint_angle);

/**
 *  Compute scatter angle based on geometry
 */
    int compute_scatter_angle(float& scat_angle);

};


class DbAlgOcean : public DbAlgorithm
{
public:

/**
 *  Class constructor
 */
    DbAlgOcean();

/**
 *  Class destructor
 */
    ~DbAlgOcean ();

/**
 *  Initialize Input data
 */
	int initialize( map<string, ddata*> imap );

/**
 *  Return a list of product names generated by the algorithm
 */
	vector<string>  get_products() { return {"cloud_mask", "quality", "aerosol_type",
		"l2_flags", "scattang", "chlorophyll", "fmf_550", "angstrom", "aot_380",
		"aot_490", "aot_550", "aot_670","aot_865", "aot_1240", "aot_1610", "aot_2250"};};

/**
 *  Compute ocean aerosol Deep Blue algorithm
 */
	map<string, ddata*> process(vector<size_t> start, vector<size_t> count,
			map<string, ddata*> imap);

	struct osOut {
        float aot550;
        float chl;
        float fmf;
        float ae;
        float ss;
        float aot[NOWL];
        short alg_flag;
        short model_flag;
    };

    osOut oOut_[NDBMDL+1];

    float chl_;
    float bathy_;
    short mask_cm_;
    short mask_cm_osi_;

    DbCloudMaskOcean* cm_;

protected:

/**
 *  Read LUT files into their respective data structures.
 */
    int initialize_LUT_data();

/**
 *  Calculate glint reflectance
 */
    float calc_glint_refl(size_t iy, size_t ix, int& status);

/**
 *  Run inversion algorithm for the selected pixel
 */
    int run_inversion(size_t iy, size_t ix, dbOceanAerosolLUMA* lut, osOut* iout);

 /**
 *  Compute linear fit to data in x and y arrays
 */
    int linfit(int size, float x[], float y[], float r[]);

/**
 *  Compute turbidity residual
 */
    float calc_turbid_residual(float sza, float r488,
            float r1240, float r1600, float r2250, float r550, int& status );

/**
 *  Set output to fill values
 */
    int set_fill_out();
};


class DbAlgLand : public DbAlgorithm
{
public:

/**
 *  Class constructor
 */
    DbAlgLand();

/**
 *  Class destructor
 */
    ~DbAlgLand ();

/**
 *  Initialize Input data
 */
	int initialize( map<string, ddata*> imap );

/**
 *  Return a list of product names generated by the algorithm
 */
	vector<string>  get_products() { return {"cloud_mask", "quality", "aerosol_type",
		"l2_flags", "scattang", "fmf_550", "angstrom", "aot_380",
		"aot_410", "aot_550", "aot_490", "aot_670"};};

/**
 */
	map<string, ddata*> process(vector<size_t> start, vector<size_t> count,
			map<string, ddata*> imap);

/**
 *  Public variables
 */
    int algflg_;
    float cofs_[16];

    float   sca_;
    float   gla_;
    float   amf_;
    float   ler412_;
    float   ler488_;
    float   ler670_;
    float   qdf412_;
    float   qdf488_;
    float   qdf670_;
    float   sr412_;
    float   sr488_;
    float   sr670_;
    float   btd8_;
    float   btd11_;
    float   dstar_;

    struct lsOut {
        float aot550;
        float ae;
        float ndvi;
        float aot[NLWL];
        float ssa[NLWL];
        float sr[NLWL];
        short sfc_type;
        short alg_flag;
        short aerosol_type;
    };

    lsOut lOut_;

    DbCloudMaskLand*        cm_;
    DbSmokeMask*            smoke_;
    DbHighAltSmokeMask*     ha_smoke_;
    DbPyrocbMask*           pyrocb_;

    short mask_cm_;
    short mask_smoke_;
    short mask_ha_smoke_;
    short mask_pyrocb_;

protected:

    float densol_[7][4];
    float denscn_[5][4];
    int   isnow_ = 0;
    float pcloud_ = 0.7;
    float pteran_ = 0;

/**
 *  Read LUT files into their respective data structures.
 */
    int initialize_LUT_data();
    int initialize_LUT_data( map<string, ddata*> imap );

/**
 *  Compute stdv
 */
    int compute_stdv();

/**
 *  Compute D* and sr670
 */
    int compute_dstar(const size_t iy, const size_t ix);

/**
 *  Compute Surface Reflectance
 */
    int compute_sr(const size_t iy, const size_t ix);

/**
 *  Compute Lambertian Equivalent Reflectance
 */
    int compute_ler(const size_t iy, const size_t ix);

/**
 *  interx:: Calculates intensity terms used in computation of
 *  reflecitivity
 */
    int interx(int i1, int i2, rhot_band w, float& ezero,
                float& tr, float& sb);

/**
 *  Set output to fill values
 */
    int set_fill_out();

};

#endif /* INCLUDE_DBALGORITHM_H_ */