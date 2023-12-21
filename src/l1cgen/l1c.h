/*
 * File:   l1c.h
 * Author: mmontes
 *
 * Created on November 4, 2020, 8:45 AM
 * //  last version 8/15/2022
 */

#ifndef L1C_H
#define L1C_H

#include <stdio.h>
#include <string>
#include <vector>
#include <filetype.h>
#include "l1c_filehandle.h"
#include "l1c_str.h"
#include "l2_str.h"
#include "l1c_input.h"
#include "hawkeye_methods.h"
#include <boost/assign/list_of.hpp> // for 'list_of()'
#include <boost/assert.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/foreach.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include<netcdf>

#define READ  0
#define WRITE 1

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> Point_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::box<Point_t> Box_t;


namespace l1c {

class L1C {
protected:

public:
    //methods--------------------------------------
    L1C();
    virtual ~L1C();

    virtual int32_t load_l1c_filehandle4(l1c_filehandle *l1cfile,L1C_input *l1cinput);
    virtual int32_t ect_sf2(const char *filename,L1C_input *l1cinput,l1c_filehandle *l1cfile);
    virtual int32_t ect_swt(int swt,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double *latswt,double *lonswt,float*tcross,float*loncross);
    virtual int32_t mov_SOCEA(l1c_filehandle* l1cfile, L1C_input* l1cinput, double* tcross, int16_t* file_id, int16_t* swtd_id, int16_t* nfiles_swath, double* ect_swtd, int16_t* tod, int16_t* orbdir, float* mgv_swath,size_t *num_scans_swt,double *time_swt,float **pos_swt,float **vel_swt);
    virtual int32_t time_swt2(int swtd,l1c_filehandle *l1cfile,L1C_input *l1cinput,double *ect_d,int16_t *swtdid,int16_t *fileid,int16_t *nfiles_swt,float *mgv_swt,double *time_mgv);   
    virtual int32_t swtime_swt2(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv);
    virtual int32_t swtime_swt2_segment(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv);
    virtual int32_t swtime_swt3(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv);
    virtual int32_t swath_latlon(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,int16_t* swtd_id, int16_t* file_id, int16_t* nfiles_swt,float** lat_tot, float** lon_tot);
    virtual int32_t calc_biny(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv);
    virtual int32_t calc_biny_dist(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,float* lati2, float* loni2);
    virtual int32_t azmean_swt3(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,float **lat_tot,float **lon_tot);
    virtual int32_t across_gridlines_l1c2(int swtd,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t *swtd_id,int16_t *file_id,int16_t *nfiles_swt,float *lati3,float *loni3, float **lat_cgd,float **lon_cgd,float *az_east);
    virtual int32_t write_L1C_granule2(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,double *tmgv,float** lat_gd, float** lon_gd,float** alt_gd);
    virtual int32_t open_l1atol1c3(L1C_input *l1cinput,l1c_filehandle *l1cfile);
    virtual bool sbs2_l1c(L1C_input *l1cinput,int32_t ydim,int32_t xdim,float **alat, short **alat_index,float latpix,float lonpix,float **lon_gd,short *erow, short *ecol);
    virtual bool sbs2_l1c2(L1C_input *l1cinput,int32_t ydim,int32_t xdim,float **alat, short **alat_index,float latpix,float lonpix,float **lon_gd,short *erow, short *ecol);
    virtual bool sbs2_l1c4(L1C_input *l1cinput,int32_t ydim,int32_t xdim,float **alat, short **alat_index,float latpix,float lonpix,float **lat_gd,float **lon_gd,short *erow, short *ecol);
    virtual bool search_rc_l1c3(L1C_input* l1cinput, l1c_filehandle* l1cfile,int32_t ydim,int32_t xdim,float latpix,float lonpix,float **lat_gd,float **lon_gd,short *erow,short *ecol);
    virtual int search_rc_l1c5(L1C_input* l1cinput, l1c_filehandle* l1cfile, l1c_str *l1cstr,short **gdindex);//line by line 
    virtual int32_t create_SOCEA2(int swtd,L1C_input* l1cinput, l1c_filehandle* l1cfile,float** lat_gd, float **lon_gd,float **altitude,double *tswt);
    virtual int32_t binL1C_wgranule_aw3(int swtd,l1c_filehandle *l1cfile, L1C_input *l1cinput,l1c_str *l1cstr,float **Ltfracsum,float **areabinsum,float **nobs_perbin,size_t sline);         
    virtual int32_t openL1Cgrid(l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput);
    virtual int32_t openL1Cgrid3(int swtd,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t* swtd_id,int16_t* file_id, int16_t* nfiles_swt);
    virtual int32_t binL1C_sbs_line3(int swtd,L1C *l1c,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t* swtd_id,int16_t* file_id, int16_t* nfiles_swt,float ****binLt,int ****bincount,float ****binLt_pol,int ****bincount_pol,size_t recnums,int granid);  
    virtual int32_t binL1C_sbs_line_l2(L1C *l1c,l2_str *l2str,l1c_filehandle *l1cfile,L1C_input *l1cinput,float ****binmean_prod,int ****bincount,size_t sline,int granid);
    virtual int32_t xy_pixsize_sf4(const char*ptstr,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,float **pix_size_u,float **pix_size_v,float **Ltfracsum,float **areabinsum,float **nobs_perbin,float ****binLt,int ****bincount,size_t recnums);
    virtual int32_t pix_corners4_l1c(l1c_filehandle *l1cfile,L1C_input *l1cinput,float dist_u,float dist_v,float azpix, int32_t scanline,int32_t pix,float pixlat,float pixlon,float pixLt,float **lat_asort,short **index_xy,float **lat_gd,float **lon_gd,double areaFracBox[3][3],float **Ltfracsum,float **areabinsum,float **nobs_perbin);
    virtual  int32_t gwindowTopix_l1c2(l1c_filehandle* l1cfile, L1C_input* l1cinput, short row, short col, double** latcornBox, double** loncornBox);
    virtual bool binIntersectsPix4corn4_l1c2(l1c_filehandle *l1cfile,L1C_input *l1cinput,short row, short col, float **lat_gd,float **lon_gd,Polygon_t &pixelPoly, double areaFracBox[3][3],double areabinBox[3][3]);
    virtual int32_t  l1b_cloud_correct(L1C_input *l1cinput,l1c_filehandle *l1cfile, netCDF::NcFile *nc_l1cgrid);
    virtual int32_t  l1c_cloud_correct(L1C_input *l1cinput,l1c_filehandle *l1cfile);

    //---------------------------------------------
    //global attributes----
    std::string l1b_name;//just the input l1b filenamed
    size_t sd_id;


    //more L1C input vars
    size_t l1c_pflag;//L1C processing 
    std::vector<std::string> cust_l1cprod;//list of L1c products to be included, flexible approach


    file_type format;//file type, netcdf4 or hdf4/5
    //size_t mode;//l1c processing mode ---
    size_t length;//data block
    size_t sensorID;//original int32 changed to size_t type
    size_t subsensorID;
    float res_spat;//spatial resolution in km
    float res_spec;//spectral resolution in nm

    //dimensions--
    size_t ndets;
    size_t nscan;
    size_t n_views;//sensor views
    size_t npols; //polarization states
    size_t nbands;//number of total bands
    size_t nband_blue;//this includes uv + visible bands
    size_t nband_red;//this inc
    size_t nband_swir;
    size_t npix;//num

 //sensor characteritics
    float *views;//indexes are not defined
    size_t pols[3];//[1 0 0] non-pol, [1 1 0]: non+cross, [1 1 1]: non+cross+parall
    float *bbands;//array with bands wavelengths
    float *rbands;
    float *swirbands;


    //navigation attributes
    //fred attr here
    size_t orbit_number;
    size_t orb_dir;//asc or desc orbit
    float orbit_node_lon;//long at which is crossing the equator asc or desc

    //geolocation attributes
    size_t terrain_corrected;
    size_t cloud_corrected;
    float *cloud_height; //should be in geolocation group nc file  

    //calibration attributes
    float *Fobar;

    //projection attributes-
    //std::string proj_type;
    size_t projection;
    float grid_resolution;//grid resolution in km
    //these params go to proj class

    //multi  attributes (view, pol, bands)
    float *view_agg; //views to be aggregated for later products such as vsfs etc
    float *pol_agg;//polarization states to be aggregated for post-processing products, linear depolarization ratio etc
    float *band_agg;//specific bands for future merged products
    bool overlap_vflag;//tells if we want merged views
    bool overlap_pflag;//tells if we want merged polarizations
    bool overlap_bflag;//tells if we want merged spectral bands 
    //uncertainty params for merged l1c products
    size_t unc_meth;//0: no error calculation, 1: propagation, 2: Monte Carlo
    float unc_thres_v; //uncertainity threshold of angular merged products as %
    float unc_thres_p;//same but for polarization
    float unc_thres_b;//same but for spectral bands
    //ancill info, requested products?
    //float *cloud_h;//this is cloud_height   
};


//protos---
    int32_t load_l1c_filehandle4(l1c_filehandle *l1cfile,L1C_input *l1cinput);
    int32_t ect_sf2(const char *filename,L1C_input *l1cinput,l1c_filehandle *l1cfile);
    int32_t ect_swt(int swt,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double *latswt,double *lonswt,float*tcross,float*loncross);
    int32_t mov_SOCEA(l1c_filehandle* l1cfile, L1C_input* l1cinput, double* tcross, int16_t* file_id, int16_t* swtd_id, int16_t* nfiles_swath, double* ect_swtd, int16_t* tod, int16_t* orbdir, float* mgv_swath,size_t *num_scans_swt,double *time_swt,float **pos_swt,float **vel_swt);
    int32_t time_swt2(int swtd,l1c_filehandle *l1cfile,L1C_input *l1cinput,double *ect_d,int16_t *swtdid,int16_t *fileid,int16_t *nfiles_swt,float *mgv_swt,double *time_mgv);
    int32_t swtime_swt2(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv);
    int32_t swtime_swt2_segment(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv);
    int32_t swtime_swt3(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv);
    int32_t swath_latlon(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,int16_t* swtd_id, int16_t* file_id, int16_t* nfiles_swt,float** lat_tot, float** lon_tot);
    int32_t calc_biny(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv);
    int32_t calc_biny_dist(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,float* lati2, float* loni2);
    int32_t azmean_swt3(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,float **lat_tot,float **lon_tot);
    int32_t across_gridlines_l1c2(int swtd,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t *swtd_id,int16_t *file_id,int16_t *nfiles_swt,float *lati3,float *loni3, float **lat_cgd,float **lon_cgd,float *az_east);
    int32_t write_L1C_granule2(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,double *tmgv,float** lat_gd, float** lon_gd,float** alt_gd);
    int32_t open_l1atol1c3(L1C_input *l1cinput,l1c_filehandle *l1cfile);
    bool sbs2_l1c(L1C_input *l1cinput,int32_t ydim,int32_t xdim,float **alat, short **alat_index,float latpix,float lonpix,float **lon_gd,short *erow, short *ecol);
    bool sbs2_l1c2(L1C_input *l1cinput,int32_t ydim,int32_t xdim,float **alat, short **alat_index,float latpix,float lonpix,float **lon_gd,short *erow, short *ecol);
    bool sbs2_l1c4(L1C_input *l1cinput,int32_t ydim,int32_t xdim,float **alat, short **alat_index,float latpix,float lonpix,float **lat_gd,float **lon_gd,short *erow, short *ecol);
    bool search_rc_l1c3(L1C_input* l1cinput, l1c_filehandle* l1cfile,int32_t ydim,int32_t xdim,float latpix,float lonpix,float **lat_gd,float **lon_gd,short *erow,short *ecol);
    int search_rc_l1c5(L1C_input* l1cinput, l1c_filehandle* l1cfile, l1c_str *l1cstr,short **gdindex);//line by line 
    int32_t create_SOCEA2(int swtd,L1C_input* l1cinput, l1c_filehandle* l1cfile,float** lat_gd, float **lon_gd,float **altitude,double *tswt);
    int32_t binL1C_wgranule_aw3(int swtd,l1c_filehandle *l1cfile, L1C_input *l1cinput,l1c_str *l1cstr,float **Ltfracsum,float **areabinsum,float **nobs_perbin,size_t sline);
    int32_t openL1Cgrid(l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput);
    int32_t openL1Cgrid3(int swtd,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t* swtd_id,int16_t* file_id, int16_t* nfiles_swt);
    int32_t binL1C_sbs_line3(int swtd,L1C *l1c,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t* swtd_id,int16_t* file_id, int16_t* nfiles_swt,float ****binLt,int ****bincount,float ****binLt_pol,int ****bincount_pol,size_t recnums,int granid);
    int32_t binL1C_sbs_line_l2(L1C *l1c,l2_str *l2str,l1c_filehandle *l1cfile,L1C_input *l1cinput,float ****binmean_prod,int ****bincount,size_t sline,int granid);
    int32_t xy_pixsize_sf4(const char*ptstr,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,float **pix_size_u,float **pix_size_v,float **Ltfracsum,float **areabinsum,float **nobs_perbin,float ****binLt,int ****bincount,size_t recnums);
    int32_t pix_corners4_l1c(l1c_filehandle *l1cfile,L1C_input *l1cinput,float dist_u,float dist_v,float azpix, int32_t scanline,int32_t pix,float pixlat,float pixlon,float pixLt,float **lat_asort,short **index_xy,float **lat_gd,float **lon_gd,double areaFracBox[3][3],float **Ltfracsum,float **areabinsum,float **nobs_perbin);
    int32_t gwindowTopix_l1c2(l1c_filehandle* l1cfile, L1C_input* l1cinput, short row, short col, double** latcornBox, double** loncornBox);
    bool binIntersectsPix4corn4_l1c2(l1c_filehandle *l1cfile,L1C_input *l1cinput,short row, short col, float **lat_gd,float **lon_gd,Polygon_t &pixelPoly, double areaFracBox[3][3],double areabinBox[3][3]);
    int32_t  l1b_cloud_correct(L1C_input *l1cinput,l1c_filehandle *l1cfile,netCDF::NcFile *nc_l1cgrid);
    int32_t  l1c_cloud_correct(L1C_input *l1cinput,l1c_filehandle *l1cfile);


}  //end namespace
#endif
