#ifndef L1C_LATLONGRID_H
#define L1C_LATLONGRID_H

#include "filehandle.h"
#include "l1.h"
#include <netcdf.h>
#include <unistd.h>
#include <nc4utils.h>
#include <netcdf>


#ifdef __cplusplus
extern "C" {
#endif

class bin_str{
   protected:

   public:
   bin_str();
   virtual ~bin_str();    
   virtual  int alloc_bin(bin_str* binstr);    
   virtual  int close_bin(bin_str* binstr);

 //global attributes
   std::string version;
   std::string history;

  //dims
   int16_t num_gridlines;
   int16_t nbinx;
   int16_t nviews;
   int16_t nbands;
   int32_t num_pixels;
   int32_t nscans;

//geo
   float **lat_gd;
   float **lon_gd;  
   double ***time_l1b;
   double *time_gd; 
   float **alt;

//binned vars
   short **nrec_2D;//row/col
   short ***nrec_3D;//row/col/view
   short ***nrec_3D_view;
   float ****nrec_4D_band;//row/col/view/bands
   float **alt_mean;
   float **alt_rmse;
   short **alt_2D;
   short **alt_diff2;//row/col
   float ***suna_3D;
   float ***sunz_3D;
   float ***sena_3D;
   float ***senz_3D;
   float ***sca_3D;//row/col/view
   float ****QC_bitwise_4D;//row/col/view/bands
   float ***QC_4D;//row/col/view
   float ****I_4D;//row/col/view/bands
   float ****I_noise_4D;//#pixels

    //OCIS line by line
   short *obs_per_view;
   float *QC_bitwise;
   float *QC;
   float *I;
   float *I_noise;

   file_format format;

   filehandle *l1file;
   std::string full_l1cgrid;

   short fillval1=-32768;
   float fillval2=-999.;

   int32_t inpix;
   int32_t outpix;


   double tini_l1c;
   double tend_l1c;
   double tini_l1b;
   double tend_l1b;


};    

int close_bin(bin_str* binstr);
int alloc_bin(bin_str* binstr);

int meta_l1c_grid(char* gridname,int16_t num_gridlines,netCDF::NcFile* nc_output);
int meta_l1c_global(char* gridname,int16_t num_gridlines,netCDF::NcFile* nc_output);
int meta_l1c_full(filehandle *l1file,bin_str *binstr,const char *l1c_grid,netCDF::NcFile *nc_output);
int meta_l1c_bin(filehandle *l1file,bin_str *binstr,netCDF::NcFile *nc_output);
int meta_l1c_altvar(bin_str *binstr,netCDF::NcFile *nc_output);//computes rmse for height

int open_l1c(const char *l1c_grid,size_t *ybins, size_t *xbins,float **lat_gd,float **lon_gd);    
int search_l1c(filehandle* l1file, l1str *l1rec,bin_str *binstr,short **gdindex);
int search2_l1c(size_t ybins,size_t xbins,float lat,float lon,float **lat_gd,float **lon_gd,short *row,short *col);
int bin_l1c(filehandle* l1file,l1str *l1rec,bin_str *binstr,short **gdindex,netCDF::NcFile *nc_output,double scantime);
int bintime_l1c(filehandle* l1file,l1str *l1rec,bin_str *binstr,short **gdindex,double scantime, netCDF::NcFile *nc_output);

int check_swath_time(filehandle* l1file,const char *l1c_grid,bin_str *binstr);//check if chosen l1c grid is correct for the chosen l1b granule used for binning and FULL l1c 
int rmse_l1c_alt(filehandle *l1file,bin_str *binstr, l1str *l1rec,short **gdindex,double scantime);//computes diff squared for rmse for height
#ifdef __cplusplus
} // extern "C"
#endif

#endif
