//**********************************************************
//library for searching row/col from lat/lon of L1C file
//  Originally created by Martin Montes on 11/4/2022
//  last version 3/13/2023
//**********************************************************
#include <allocate2d.h>
#include <allocate3d.h>
#include "allocate4d.h"
#include <iostream>
#include <chrono>
#include <sys/stat.h>
#include <genutils.h>
#include <iostream>
#include "l1c_latlongrid.h"
#include<netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


bin_str::bin_str(){

 //global attrs
   version="";
   history="";

 //dimensions
   nviews=-1;   
   nbands=-1;
   num_gridlines=-1;
   nbinx=-1;
   num_pixels=-1;
   nscans=-1;
//geo
   lat_gd=nullptr;
   lon_gd=nullptr;
   time_l1b=nullptr;
   time_gd=nullptr;
   alt=nullptr;

 //binned vars
   nrec_2D=nullptr;//row/col
   nrec_3D=nullptr;
   nrec_3D_view=nullptr;
   nrec_4D_band=nullptr;//row/col/view
   alt_mean=nullptr;
   alt_rmse=nullptr;
   alt_2D=nullptr;
   alt_diff2=nullptr;//row/col
   sca_3D=nullptr;//row/col/view
   QC_bitwise_4D=nullptr;//row/col/view/bands
   QC_4D=nullptr;//row/col/view
   I_4D=nullptr;//row/col/view/bands
   I_noise_4D=nullptr;//#pixels

    //OCIS line by line
   obs_per_view=nullptr;
   QC_bitwise=nullptr;
   QC=nullptr;
   I=nullptr;
   I_noise=nullptr;

   l1file=nullptr;
   full_l1cgrid="";

   fillval1=-32768;
   fillval2=-999.;

   inpix=0;
   outpix=0;

   tini_l1c=-1;
   tend_l1c=-1;
   tini_l1b=-1;
   tend_l1b=-1;
}

bin_str::~bin_str(){

}    




int bin_str::close_bin(bin_str* binl1c){

   if(binl1c->time_l1b!=nullptr)
    delete [] (binl1c->time_l1b);
   binl1c->time_l1b=nullptr;
   if(binl1c->time_gd!=nullptr)
    delete [] (binl1c->time_gd);
   binl1c->time_gd=nullptr;
   if(binl1c->alt!=nullptr)
    delete [] (binl1c->alt);
   binl1c->alt=nullptr;
   if(binl1c->lat_gd!=nullptr)
    delete [] (binl1c->lat_gd);
   binl1c->lat_gd=nullptr;
   if(binl1c->lon_gd!=nullptr)
    delete [] (binl1c->lon_gd);
   binl1c->lon_gd=nullptr;
   if(binl1c->alt_rmse!=nullptr)
    delete [] (binl1c->alt_rmse);
   binl1c->alt_rmse=nullptr;
   if(binl1c->alt_mean!=nullptr)
    delete [] (binl1c->alt_mean);
   binl1c->alt_mean=nullptr;
   if(binl1c->alt_2D!=nullptr)
    delete [] (binl1c->alt_2D);
   binl1c->alt_2D=nullptr;
   if(binl1c->alt_diff2!=nullptr)
    delete [] (binl1c->alt_diff2);
   binl1c->alt_diff2=nullptr;
   if(binl1c->suna_3D!=nullptr)
    delete [] (binl1c->suna_3D);
   binl1c->suna_3D=nullptr;
   if(binl1c->sunz_3D!=nullptr)
    delete [] (binl1c->sunz_3D);
   binl1c->sunz_3D=nullptr;
   if(binl1c->sena_3D!=nullptr)
    delete [] (binl1c->sena_3D);
   binl1c->sena_3D=nullptr;
   if(binl1c->senz_3D!=nullptr)
    delete [] (binl1c->senz_3D);
   binl1c->senz_3D=nullptr;
   if(binl1c->sca_3D!=nullptr)
    delete [] (binl1c->sca_3D);
   binl1c->sca_3D=nullptr;
   if(binl1c->nrec_2D!=nullptr)
    delete [] (binl1c->nrec_2D);
   binl1c->nrec_2D=nullptr;
   if(binl1c->nrec_3D!=nullptr)
    delete [] (binl1c->nrec_3D);
   binl1c->nrec_3D=nullptr;
   if(binl1c->nrec_3D_view!=nullptr)
    delete [] (binl1c->nrec_3D_view);
   binl1c->nrec_3D_view=nullptr;
   if(binl1c->nrec_4D_band!=nullptr)
    delete [] (binl1c->nrec_4D_band);
   binl1c->nrec_4D_band=nullptr;
   if(binl1c->I_4D!=nullptr)
    delete [] (binl1c->I_4D);
   binl1c->I_4D=nullptr;
   if(binl1c->I_noise_4D!=nullptr)
    delete [] (binl1c->I_noise_4D);
   binl1c->I_noise_4D=nullptr;

return 0;

}    


int bin_str::alloc_bin(bin_str* binl1c){

   cout<<"allocating/init high order dimensional arrays for L1C binning----------------------------"<<endl;
   cout<<"nviews.."<<binl1c->nviews<<"nbands.."<<binl1c->nbands<<"#gridlines.."<<binl1c->num_gridlines<<"nbinx.."<<binl1c->nbinx<<endl;

   binl1c->nrec_2D=allocate2d_short(binl1c->num_gridlines,binl1c->nbinx);
   binl1c->alt = allocate2d_float(binl1c->num_gridlines,binl1c->nbinx);
   binl1c->alt_rmse = allocate2d_float(binl1c->num_gridlines,binl1c->nbinx);
   binl1c->alt_mean = allocate2d_float(binl1c->num_gridlines,binl1c->nbinx);
   binl1c->alt_2D = allocate2d_short(binl1c->num_gridlines,binl1c->nbinx);
   binl1c->alt_diff2 = allocate2d_short(binl1c->num_gridlines,binl1c->nbinx);
   binl1c->lat_gd = allocate2d_float(binl1c->num_gridlines,binl1c->nbinx);
   binl1c->lon_gd = allocate2d_float(binl1c->num_gridlines,binl1c->nbinx);

   binl1c->time_gd=(double*)calloc(binl1c->num_gridlines,sizeof(double));
   binl1c->time_l1b = allocate3d_double(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   binl1c->suna_3D = allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   binl1c->sunz_3D= allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   binl1c->sena_3D = allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   binl1c->senz_3D = allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   binl1c->sca_3D = allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   binl1c->nrec_3D=allocate3d_short(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   binl1c->nrec_3D_view=allocate3d_short(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);

   binl1c->nrec_4D_band=allocate4d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews,binl1c->nbands);
   binl1c->I_4D=allocate4d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews,binl1c->nbands);
   binl1c->I_noise_4D=allocate4d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews,binl1c->nbands);


//init
   for(int i=0;i<binl1c->num_gridlines;i++){
      binl1c->time_gd[i]=0;
       for(int j=0;j<binl1c->nbinx;j++){
           binl1c->nrec_2D[i][j]=0;
           binl1c->alt[i][j]=0.;
           binl1c->alt_rmse[i][j]=0.;
           binl1c->alt_mean[i][j]=0.;
           binl1c->alt_2D[i][j]=0;
           binl1c->alt_diff2[i][j]=0;
           binl1c->lat_gd[i][j]=0.;
           binl1c->lon_gd[i][j]=0.;

          for(int v=0;v<binl1c->nviews;v++){
               binl1c->time_l1b[i][j][v]=0.;
               binl1c->suna_3D[i][j][v]=0.;
               binl1c->sunz_3D[i][j][v]=0.;
               binl1c->sena_3D[i][j][v]=0.;
               binl1c->senz_3D[i][j][v]=0.;
               binl1c->sca_3D[i][j][v]=0.;
               binl1c->nrec_3D[i][j][v]=0;
               binl1c->nrec_3D_view[i][j][v]=0;
            for(int sb=0;sb<binl1c->nbands;sb++){
               binl1c->nrec_4D_band[i][j][v][sb]=0.;
               binl1c->I_4D[i][j][v][sb]=0.;
               binl1c->I_noise_4D[i][j][v][sb]=0.;
              }}}}
   
return 0;
}    


int rmse_l1c_alt(filehandle *l1file,bin_str *binl1c, l1str *l1rec,short **gdindex,double scantime){
   short gd_row=-1,gd_col=-1;
              int npix=l1file->npix;
  
   for(int pix=0;pix<npix;pix++){
        gd_row=gdindex[pix][0]-1;
        gd_col=gdindex[pix][1]-1;


          if(scantime>=binl1c->tini_l1c && scantime<=binl1c->tend_l1c && gd_row>=0 && gd_col>=0 && l1rec->lat[pix]!=binl1c->fillval1 && l1rec->lon[pix]!=binl1c->fillval2)



        //      binl1c->alt_diff2[gd_row][gd_col]=binl1c->alt_diff2[gd_row][gd_col]+(l1rec->height[pix]-binl1c->alt_mean[gd_row][gd_col])*(l1rec->height[pix]-binl1c->alt_mean[gd_row][gd_col]);
     //         binl1c->alt_diff2[gd_row][gd_col]=binl1c->alt_diff2[gd_row][gd_col]+sqrt((rand()%10000-binl1c->alt_mean[gd_row][gd_col])*(rand()%10000-binl1c->alt_mean[gd_row][gd_col]));   
              binl1c->alt_diff2[gd_row][gd_col]=binl1c->alt_diff2[gd_row][gd_col]+(10000-binl1c->alt_mean[gd_row][gd_col])*(10000-binl1c->alt_mean[gd_row][gd_col]);
          
   }
                   
return 0;
}


int meta_l1c_altvar(bin_str *binl1c,NcFile* nc_output){
   float term;

   cout<<"adding alt_rmse to L1C file..."<<endl;

   NcGroup geo_grp=nc_output->getGroup("geolocation_data");
   NcVar v1=geo_grp.getVar("altitude_variability");

   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
           term=binl1c->alt_diff2[i][j];
          if(term>=0 && binl1c->nrec_2D[i][j]>0){
           binl1c->alt_rmse[i][j]=sqrt(binl1c->alt_diff2[i][j]/binl1c->nrec_2D[i][j]);       
          }
          else binl1c->alt_rmse[i][j]=binl1c->fillval1; 

         }}
   v1.putVar(&binl1c->alt_rmse[0][0]);

  

    return 0;
}



int meta_l1c_bin(filehandle *l1file,bin_str *binl1c,NcFile *nc_output){

   cout<<"adding final binned vars to ...."<<binl1c->full_l1cgrid<<endl;

   NcDim ydim=nc_output->getDim("bins_along_track"); 
   NcDim xdim=nc_output->getDim("bins_across_track");
   NcDim vdim=nc_output->getDim("number_of_views");

   std::vector<NcDim> dimvec3;
   dimvec3.push_back(ydim);
   dimvec3.push_back(xdim);
   dimvec3.push_back(vdim);

/*  
  //mean height---
  temp=allocate2d_float(binl1c->num_gridlines,binl1c->nbinx);
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
                if(binl1c->nrec_2D[i][j]>0){
                   binl1c->alt_mean[i][j]=binl1c->alt_2D[i][j]/binl1c->nrec_2D[i][j];
                }

          }}

   v1=geo_grp.getVar("altitude");
   v1.putVar(&temp[0][0]);
   delete [] (temp);
*/

//time_offsets

   double ***time_offsets=nullptr, toff,toff_fill,toff_min,toff_max;

//   time_mean=allocate3d_double(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   time_offsets=allocate3d_double(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);

      for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
            time_offsets[i][j][v]=0;
          }}}

   NcGroup geo_grp=nc_output->getGroup("geolocation_data");
   NcGroup ba_grp=nc_output->getGroup("bin_attributes");
   NcVar v1=ba_grp.getVar("nadir_view_time");
   v1.getVar(&binl1c->time_gd[0]);


   v1=ba_grp.getVar("view_time_offsets");
   NcVarAtt a1=v1.getAtt("_FillValue");//root group
   a1.getValues(&toff_fill);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&toff_min);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&toff_max);

//sensor azimuth---   
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
                if(binl1c->nrec_3D[i][j][v]>0){

                  toff=binl1c->time_gd[i]-(binl1c->time_l1b[i][j][v]/binl1c->nrec_3D[i][j][v]);   
                  if(toff>=toff_min && toff<=toff_max)
                  { 
                    time_offsets[i][j][v]=toff;
           //         if(v==1) cout<<"toff_fill.."<<toff_fill<<"toff.."<<toff<<endl;
                  }
                  else time_offsets[i][j][v]=toff_fill;
                }
               else{
                    time_offsets[i][j][v]=toff_fill;
         //           if(v==1) cout<<"toff_fill.."<<toff_fill<<endl; 
                     }             
          }}}

   v1=ba_grp.getVar("view_time_offsets");
   v1.putVar(&time_offsets[0][0][0]);



   delete [] (time_offsets);

 
   v1=geo_grp.getVar("sensor_azimuth");

   float ***temp=nullptr;

   temp=allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);

//sensor azimuth---   
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
                if(binl1c->nrec_3D[i][j][v]>0){
                 temp[i][j][v]=binl1c->sena_3D[i][j][v]/binl1c->nrec_3D[i][j][v];                
                }
                else{
                  temp[i][j][v]=binl1c->fillval1;
                }
              
          }}}

   v1.putVar(&temp[0][0][0]);
   delete [] (temp); 

//sensor zenith---   
   temp=allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
                if(binl1c->nrec_3D[i][j][v]>0){
                 temp[i][j][v]=binl1c->senz_3D[i][j][v]/binl1c->nrec_3D[i][j][v];
                }
                else{
                  temp[i][j][v]=binl1c->fillval1;
                }
             
          }}}

   v1=geo_grp.getVar("sensor_zenith");
   v1.putVar(&temp[0][0][0]);
   delete [] (temp);



   //solar azimuth---
   temp=allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
                if(binl1c->nrec_3D[i][j][v]>0){
                 temp[i][j][v]=binl1c->suna_3D[i][j][v]/binl1c->nrec_3D[i][j][v];
                }
                else{
                  temp[i][j][v]=binl1c->fillval1;
                }
            
          }}}

   v1=geo_grp.getVar("solar_azimuth");
   v1.putVar(&temp[0][0][0]);
   delete [] (temp);

   //sensor zenith---
   temp=allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
                if(binl1c->nrec_3D[i][j][v]>0){
                 temp[i][j][v]=binl1c->sunz_3D[i][j][v]/binl1c->nrec_3D[i][j][v];
                }
                else{
                  temp[i][j][v]=binl1c->fillval1;
                }
           
          }}}

   v1=geo_grp.getVar("solar_zenith");
   v1.putVar(&temp[0][0][0]);
   delete [] (temp);

 //scattering angle---
   temp=allocate3d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews);
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
                if(binl1c->nrec_3D[i][j][v]>0){
                 temp[i][j][v]=binl1c->sca_3D[i][j][v]/binl1c->nrec_3D[i][j][v];             
                }
                else{
                  temp[i][j][v]=binl1c->fillval1;
                }
          
          }}}

   v1=geo_grp.getVar("scattering_angle");
   v1.putVar(&temp[0][0][0]);
   delete [] (temp);



   NcGroup od_grp=nc_output->getGroup("observation_data");
   v1=od_grp.getVar("obs_per_view"); //ONLY 1 BAND!!!!!! NO QUALITY CONTROL,constrained no fillvalues, and withing min/max Lt

   v1.putVar(&binl1c->nrec_3D[0][0][0]);



//I/reflectance
   float ****temp2=allocate4d_float(binl1c->num_gridlines,binl1c->nbinx,binl1c->nviews,binl1c->nbands);
   for(int i=0;i<binl1c->num_gridlines;i++){
      for(int j=0;j<binl1c->nbinx;j++){
          for(int v=0;v<binl1c->nviews;v++){
               for(int sb=0;sb<binl1c->nbands;sb++){
                if(binl1c->nrec_4D_band[i][j][v][sb]>0)
                {
                  temp2[i][j][v][sb]=binl1c->I_4D[i][j][v][sb]/binl1c->nrec_4D_band[i][j][v][sb];
/*                    if(v==1)
                    {
                        cout<<"row."<<i+1<<"col.."<<j+1<<"MEAN reflectance..band #1 OCIS = "<<temp2[i][j][1][0]<<endl;
                    }       
*/                    
                }
                else{
                  temp2[i][j][v][sb]=binl1c->fillval2;
                }
         
          }}}}

   v1=od_grp.getVar("I");
   v1.putVar(&temp2[0][0][0][0]);
   delete [] (temp2);

return 0;    
}


int meta_l1c_full(filehandle *l1file,bin_str *binl1c,const char *l1c_grid,NcFile *nc_output){
   string senstr,GATT_VAL1,prodstr,binstr,ifile_str;
   string date_created,y_create,m_create,d_create,t_create;
   int NVIEWS,NBANDS;//NBANDS_POL;
   int32_t xbins;//,ybins;

   char *ifile_char=strdup(l1c_grid);
   file_format format = getFormat(ifile_char);
   if(format.type==FT_HKT || format.type==FT_L1C) format.type=FT_OCIS; //if HKT then FT_OCIS----


   if(format.type==FT_SPEXONE){
                         senstr = "SPEXONE";
                         GATT_VAL1="PACE SPEXone Level-1C Data";
                         binstr="12";
                         xbins=25;
                         NVIEWS=5;
                         NBANDS=400;
                     //    NBANDS_POL=50;
                     }
   else if(format.type==FT_HARP2){
                         senstr = "HARP2";
                         GATT_VAL1="PACE HARP2 Level-1C Data";
                         binstr="228";
                         xbins=457;
                         NVIEWS=90;
                         NBANDS=1;
                     //    NBANDS_POL=1;
                     }
   else if(format.type==FT_OCIS || format.type==FT_OCIL1B){
                         senstr = "OCI";
                         GATT_VAL1="PACE OCI Level-1C Data";
                         binstr="259";
                         xbins=519;
                         NVIEWS=2;
                         NBANDS=239;//249 originally no polarization bands
                     }    
   else//OCIS
                     {
                         senstr = "OCI";
                         GATT_VAL1="PACE OCI Level-1C Data";
                         binstr="259";
                         xbins=519;
                         NVIEWS=2;
                         NBANDS=239;//249 originally no polarization bands
                     }

    //l1c grid
   string l1c_str=l1c_grid;

   NcFile* nc_l1cgrid;
   try {
        nc_l1cgrid = new NcFile(l1c_grid, NcFile::read);
        }
   catch (NcException& e) {
          e.what();
          cerr << "l1cgen l1c_pflag= 8:: Failure reading L1C grid: "
          + l1c_str << endl;
          exit(1);
   }

   string file_str=l1c_str.substr(5,15);
   string l1c_full_str="PACE_"+senstr+"."+file_str+"Z.L1C.5.2km.nc";
   binl1c->full_l1cgrid=l1c_full_str;
 
   //l1c_grid--------
   NcDim yd=nc_l1cgrid->getDim("bins_along_track");
   int32_t num_gridlines=yd.getSize();
   binl1c->nviews=NVIEWS;
   binl1c->nbands=NBANDS;
   binl1c->num_gridlines=num_gridlines;
   binl1c->nbinx=xbins;

//alloc binned vars in binl1c class
   binl1c->alloc_bin(binl1c);

   NcDim vdim=nc_l1cgrid->getDim("number_of_views"); 
   NcDim idim=nc_l1cgrid->getDim("intensity_bands_per_view"); 

   meta_l1c_global(ifile_char,num_gridlines,nc_output);
 
//adding fill values of binned variables----
//sena,senz,sola,solz,sca, I
   NcGroup geo_grp=nc_output->getGroup("geolocation_data");
   NcVar v1=geo_grp.getVar("sensor_azimuth");
   short value;
   NcVarAtt a1=v1.getAtt("_FillValue");//root group
   a1.getValues(&value);

   binl1c->fillval1=value;//-32768s 

   NcGroup od_grp=nc_output->getGroup("observation_data");
   v1=od_grp.getVar("I");
   float value2;
   a1=v1.getAtt("_FillValue");//root group
   a1.getValues(&value2);

   binl1c->fillval2=value2;//-999.f ;

//global attributes
   nc_output->putAtt("processing_version",binl1c->version);
   nc_output->putAtt("history",binl1c->history);   
   string name;
   nc_output->putAtt("product_name",binl1c->full_l1cgrid);
   NcGroupAtt i1=nc_l1cgrid->getAtt("startDirection");//NcGroupAtt is a global attr!!
   i1.getValues(name);
   nc_output->putAtt("startDirection",name);
   i1=nc_l1cgrid->getAtt("endDirection");
   i1.getValues(name);
   nc_output->putAtt("endDirection",name);
   i1=nc_l1cgrid->getAtt("time_coverage_start");
   i1.getValues(name);
   nc_output->putAtt("time_coverage_start",name);
   i1=nc_l1cgrid->getAtt("time_coverage_end");
   i1.getValues(name);
   nc_output->putAtt("time_coverage_end",name);
   
   double *time_nad = (double*)calloc(num_gridlines,sizeof(double));
   NcGroup ba_grp=nc_l1cgrid->getGroup("bin_attributes");
   v1=ba_grp.getVar("nadir_view_time");
   v1.getVar(time_nad);
   ba_grp=nc_output->getGroup("bin_attributes");
   v1=ba_grp.getVar("nadir_view_time");
   v1.putVar(&time_nad[0]);
   if(time_nad!=nullptr)
     delete [] (time_nad);
   time_nad=nullptr;
 
   geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   v1=geo_grp.getVar("latitude");
   v1.getVar(&binl1c->lat_gd[0][0]);
   geo_grp=nc_output->getGroup("geolocation_data");
   v1=geo_grp.getVar("latitude");
   v1.putVar(&binl1c->lat_gd[0][0]);
   
   geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   v1=geo_grp.getVar("longitude");
   v1.getVar(&binl1c->lon_gd[0][0]);
   geo_grp=nc_output->getGroup("geolocation_data");
   v1=geo_grp.getVar("longitude");
   v1.putVar(&binl1c->lon_gd[0][0]);

   geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   v1=geo_grp.getVar("altitude");
   v1.getVar(&binl1c->alt[0][0]);
   geo_grp=nc_output->getGroup("geolocation_data");
   v1=geo_grp.getVar("altitude");
   v1.putVar(&binl1c->alt[0][0]);

 
     //GRING----
   geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   NcGroupAtt ga1=geo_grp.getAtt("GRingPointLatitude");
   size_t dp=ga1.getAttLength(); 
   float *latarr=(float*)calloc(dp,sizeof(float));
   ga1.getValues(&latarr[0]);
   geo_grp=nc_output->getGroup("geolocation_data");
   geo_grp.putAtt("GRingPointLatitude",ncFloat,dp,latarr);

   geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   float *lonarr=(float*)calloc(dp,sizeof(float));
   ga1=geo_grp.getAtt("GRingPointLongitude");
   ga1.getValues(&lonarr[0]);
   geo_grp=nc_output->getGroup("geolocation_data");
   geo_grp.putAtt("GRingPointLongitude",ncFloat,dp,lonarr);

   geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   int *narr=(int*)calloc(dp,sizeof(int));
   ga1=geo_grp.getAtt("GRingPointSequenceNo");
   ga1.getValues(&narr[0]);
   geo_grp=nc_output->getGroup("geolocation_data");
   geo_grp.putAtt("GRingPointSequenceNo",ncInt,dp,narr);

   if (latarr!=nullptr)
      delete [] (latarr);
   if (lonarr!=nullptr)
      delete [] (lonarr);
   if (narr!=nullptr)
      delete [] (narr);

   std::vector<NcDim> dimvec2_rad;
   dimvec2_rad.push_back(vdim);
   dimvec2_rad.push_back(idim);

   int32_t *iwave=l1file->iwave;
   float *Fobar=l1file->Fobar;
   float *fwhm=l1file->fwave;

   float fwave_view[NVIEWS][NBANDS],fwhm_view[NVIEWS][NBANDS],Fobar_view[NVIEWS][NBANDS];

   for(int i=0;i<NVIEWS;i++){
   for(int j=0;j<NBANDS;j++){
       fwave_view[i][j]=iwave[j];
       fwhm_view[i][j]=fwhm[j];
       Fobar_view[i][j]=Fobar[j];
   }}

   //intensity_wavelengths, bandpasses and F0---
   NcGroup svb_grp=nc_output->getGroup("sensor_views_bands");
   v1=svb_grp.addVar("intensity_wavelengths",ncFloat,dimvec2_rad);//this is infull meta
   string longName="Intensity field center wavelengths at each view";
   v1.putAtt( "long_name",longName);
   string units="nm";
   v1.putAtt( "units",units);
   float FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   float valid_min=320;
   v1.putAtt("valid_min",ncFloat,valid_min);
   float valid_max=2260;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(fwave_view);

   //intensity_wavelengths, bandpasses and F0---
   v1=svb_grp.addVar("intensity_bandpasses",ncFloat,dimvec2_rad);
   longName="Intensity field bandpasses at each view";
   v1.putAtt( "long_name",longName);
   units="nm";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=2.5;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=100;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(fwhm_view);

   v1=svb_grp.addVar("intensity_F0",ncFloat,dimvec2_rad);
   longName="Intensity band solar irradiance";
   v1.putAtt( "long_name",longName);
   units="W m^-2 µm^-1";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=0;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=900;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(Fobar_view);

   nc_l1cgrid->close();


return 0;
}





int meta_l1c_grid(char* gridname,int16_t num_gridlines,NcFile* nc_output){
   string senstr,GATT_VAL1,prodstr,binstr,ifile_str;
   string date_created,y_create,m_create,d_create,t_create;
   int NVIEWS,NBANDS,NBANDS_POL;
   int32_t xbins,ybins;

   ybins=num_gridlines;
   prodstr=string(gridname);


   file_format format=getFormat(gridname);
   if(format.type==FT_HKT || format.type==FT_L1C) format.type=FT_OCIS; 

   if(format.type==FT_SPEXONE){
                         senstr = "SPEXONE";
                         GATT_VAL1="PACE SPEXone Level-1C Data";
                         binstr="12";
                         xbins=25;
                         NVIEWS=5;
                         NBANDS=400;
                         NBANDS_POL=50;
                     }
   else if(format.type==FT_HARP2){
                         senstr = "HARP2";
                         GATT_VAL1="PACE HARP2 Level-1C Data";
                         binstr="228";
                         xbins=457;
                         NVIEWS=90;
                         NBANDS=1;
                         NBANDS_POL=1;
                     }
   else if(format.type==FT_OCIS || format.type==FT_OCIL1B){
                         senstr = "OCI";
                         GATT_VAL1="PACE OCI Level-1C Data";
                         binstr="259";
                         xbins=519;
                         NVIEWS=2;
                         NBANDS=239;//249 originally no polarization bands
                     }
   else//OCIS
                     {
                         senstr = "OCI";
                         GATT_VAL1="PACE OCI Level-1C Data";
                         binstr="259";
                         xbins=519;
                         NVIEWS=2;
                         NBANDS=239;//249 ORIGINALLY no polarization bands
                     }
   
 // creation date---
// Current date/time based on current system
   time_t now = time(0);
  // Convert now to tm struct for UTC
   tm* gmtm = gmtime(&now);
   if (gmtm != NULL) {

    }
   else {
    cerr << "Failed to get the UTC date and time" << endl;
    return EXIT_FAILURE;
   }

   y_create=std::to_string(1900 + gmtm->tm_year);
   m_create=std::to_string(1 + gmtm->tm_mon);
   d_create=std::to_string(gmtm->tm_mday);
   t_create=std::to_string(gmtm->tm_hour) +":"+std::to_string(gmtm->tm_min)+":"+std::to_string(gmtm->tm_sec);
   date_created=y_create+"-"+m_create+"-"+d_create+"T"+t_create+"Z";

//global attributes---

   nc_output->putAtt("title", GATT_VAL1);
   nc_output->putAtt("instrument", senstr);
   nc_output->putAtt("processing_version","");
   nc_output->putAtt("Conventions","CF-1.8 ACDD-1.3");
   nc_output->putAtt("institution","NASA Goddard Space Flight Center, Ocean Biology Processing Group");
   nc_output->putAtt("license","http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/");
   nc_output->putAtt("naming_authority","gov.nasa.gsfc.sci.oceancolor");
   nc_output->putAtt("keywords_vocabulary","NASA Global Change Master Directory (GCMD) Science Keywords");
   nc_output->putAtt("stdname_vocabulary","NetCDF Climate and Forecast (CF) Metadata Convention");
   nc_output->putAtt("creator_name","NASA/GSFC");
   nc_output->putAtt("creator_email","data@oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("creator_url","http://oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("project","PACE Project");
   nc_output->putAtt("publisher_name","NASA/GSFC");
   nc_output->putAtt("publisher_email","data@oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("publisher_url","http://oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("processing_level","L1C");
   nc_output->putAtt("cdm_data_type","swath");
   nc_output->putAtt("orbit_number","xxx");
   nc_output->putAtt("history","");
   nc_output->putAtt("CDL_version_date","2021-09-10");
   nc_output->putAtt("product_name",prodstr);
   nc_output->putAtt("startDirection","");
   nc_output->putAtt("endDirection","");
   nc_output->putAtt("time_coverage_start","");
   nc_output->putAtt("time_coverage_end","");
   nc_output->putAtt("date_created",date_created);
   nc_output->putAtt("sun_earth_distance","0.990849042172323");
   nc_output->putAtt("terrain_data_source","");
   nc_output->putAtt("spectral_response_function","");
   nc_output->putAtt("systematic_uncertainty_model","");
   nc_output->putAtt("nadir_bin",binstr);
   nc_output->putAtt("bin_size_at_nadir","5.2km2");

   NcDim ydim=nc_output->addDim("bins_along_track",ybins);
   NcDim xdim=nc_output->addDim("bins_across_track",xbins);
   NcDim vdim=nc_output->addDim("number_of_views", NVIEWS);
   NcDim idim=nc_output->addDim("intensity_bands_per_view",NBANDS);

//dims
   std::vector<NcDim> dimvec2_geo;
   dimvec2_geo.push_back(ydim);
   dimvec2_geo.push_back(xdim);
   std::vector<NcDim> dimvec2_rad;
   dimvec2_rad.push_back(vdim);
   dimvec2_rad.push_back(idim);
   std::vector<NcDim> dimvec3;
   dimvec3.push_back(ydim);
   dimvec3.push_back(xdim);
   dimvec3.push_back(vdim);
   std::vector<NcDim> dimvec4;
   dimvec4.push_back(ydim);
   dimvec4.push_back(xdim);
   dimvec4.push_back(vdim);
   dimvec4.push_back(idim);

//groups
   NcGroup ba_grp = nc_output->addGroup("bin_attributes");
   NcGroup geo_grp = nc_output->addGroup("geolocation_data");

    //vars
 /*  NcVar v2=svb_grp.addVar("intensity_wavelengths",ncFloat,dimvec2_rad);//this is infull meta
   longName="Intensity field center wavelengths at each view";
   v1.putAtt( "long_name",longName);
   units="nm";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=320;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=2260;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(fwave_view);
*/
   if(format.type==FT_HARP2 || format.type==FT_SPEXONE){

       NcDim pdim=nc_output->addDim("polarization_bands_per_view",NBANDS_POL);
       std::vector<NcDim> dimvec2b_rad;
       dimvec2b_rad.push_back(vdim);
       dimvec2b_rad.push_back(pdim);
       std::vector<NcDim> dimvec4b;
       dimvec4b.push_back(ydim);
       dimvec4b.push_back(xdim);
       dimvec4b.push_back(vdim);
       dimvec4b.push_back(pdim);

        //grp svb
 /*      float fwave_view_pol[NVIEWS][NBANDS_POL],fwhm_view_pol[NVIEWS][NBANDS_POL],Fobar_view_pol[NVIEWS][NBANDS_POL];

       for(int i=0;i<NVIEWS;i++){
         for(int j=0;j<NBANDS_POL;j++){
             fwave_view_pol[i][j]=fwave[j];
             fwhm_view_pol[i][j]=fwhm[j];
             Fobar_view_pol[i][j]=Fobar[j];
       }}
       v1=svb_grp.addVar("polarization_wavelengths",ncFloat,dimvec2b_rad);
       string longName2="Polarization field wavelengths at each view";
       v1.putAtt( "long_name",longName2);
       string units2="nm";
       v1.putAtt( "units",units2);
       float FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       float valid_min_p=320;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       float valid_max_p=2260;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
       v1.putVar(fwave_view_pol);

       v1=svb_grp.addVar("polarization_bandpasses",ncFloat,dimvec2b_rad);
       longName2="Polarization field bandpasses at each view";
       v1.putAtt( "long_name",longName2);
       units2="nm";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=2.5;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=100;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
       v1.putVar(fwhm_view_pol);

       v1=svb_grp.addVar("polarization_F0",ncFloat,dimvec2b_rad);
       longName2="Polarization band solar irradiance";

       v1.putAtt( "long_name",longName2);
       units2="W m^-2 µm^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=0;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=900;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
       v1.putVar(Fobar_view_pol);
*/

   if(format.type==FT_HARP2){
       //grp od_grp
 /*      v1 = od_grp.addVar("Q", ncFloat,dimvec4b);
       longName2="Q Stokes vector component";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=-800;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);

       v1 = od_grp.addVar("Q_noise", ncFloat,dimvec4b);
       longName2="Random noise of Q in bin";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=0;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);

       v1 = od_grp.addVar("U", ncFloat,dimvec4b);
       longName2="U Stokes vector component";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=-800;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);

       v1 = od_grp.addVar("U_noise", ncFloat,dimvec4b);
       longName2="Random noise of U in bin";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=0;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
*/ 
     }

/*   v1 = od_grp.addVar("DOLP", ncFloat,dimvec4b);
   longName2="Degree of linear polarization";
   v1.putAtt( "long_name",longName2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=1;
   v1.putAtt("valid_max",ncFloat,valid_max_p);

   v1 = od_grp.addVar("DOLP_noise", ncFloat,dimvec4b);
   longName2="Random noise of DOLP in bin";

   v1.putAtt( "long_name",longName2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=1;
   v1.putAtt("valid_max",ncFloat,valid_max_p);

   v1 = od_grp.addVar("AOLP", ncFloat,dimvec4b);
   longName2="Angle of linear polarization";
   v1.putAtt( "long_name",longName2);
   units2="degrees";
   v1.putAtt("units",units2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=360;
   v1.putAtt("valid_max",ncFloat,valid_max_p);

   v1 = od_grp.addVar("AOLP_noise", ncFloat,dimvec4b);
   longName2="Random noise of AOLP in bin";
   v1.putAtt( "long_name",longName2);
   units2="degrees";
   v1.putAtt("units",units2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=360;
   v1.putAtt("valid_max",ncFloat,valid_max_p);
*/
    }
/*   v1=svb_grp.addVar("intensity_bandpasses",ncFloat,dimvec2_rad);
   longName="Intensity field bandpasses at each view";
   v1.putAtt( "long_name",longName);
   units="nm";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=2.5;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=100;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(fwhm_view);

   v1=svb_grp.addVar("intensity_F0",ncFloat,dimvec2_rad);
   longName="Intensity band solar irradiance";
   v1.putAtt( "long_name",longName);
   units="W m^-2 µm^-1";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=0;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=900;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(Fobar_view);
*/
   NcVar v1=ba_grp.addVar("nadir_view_time",ncDouble,ydim);
   double FillValue2;
   string longName="Time bin was viewed at nadir view";
   v1.putAtt( "long_name",longName);
   string units="seconds";
   v1.putAtt( "units",units);
   FillValue2=-999.;
   v1.putAtt( "_FillValue",ncDouble,FillValue2);
   double valid_min_d=0;
   v1.putAtt("valid_min",ncDouble,valid_min_d);
   double valid_max_d=86400.999;
   v1.putAtt("valid_max",ncDouble,valid_max_d);

   v1 = geo_grp.addVar("latitude", ncFloat,dimvec2_geo);
   longName="Latitudes of bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees_north";
   v1.putAtt( "units",units);
   float FillValue=-999.;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   float valid_min=-90;
   v1.putAtt("valid_min",ncFloat,valid_min);
   float valid_max=90;
   v1.putAtt("valid_max",ncFloat,valid_max);

   v1 = geo_grp.addVar("longitude", ncFloat,dimvec2_geo);
   longName="Longitudes of bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees_east";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=-180;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=180;
   v1.putAtt("valid_max",ncFloat,valid_max);

   v1 = geo_grp.addVar("altitude", ncFloat,dimvec2_geo);
   longName="Altitude at bin locations";
   v1.putAtt( "long_name",longName);
   units="meters";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   short valid_min_s=-1000;
   v1.putAtt("valid_min",ncFloat,valid_min_s);
   short valid_max_s=10000;
   v1.putAtt("valid_max",ncFloat,valid_max_s);
   float scale_factor=1;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   float add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);
/*
   v1 = geo_grp.addVar("altitude_variability", ncShort,dimvec2_geo);
   longName="RMS variability of altitude at bin locations";
   v1.putAtt( "long_name",longName);
   units="meters";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncShort,FillValue);
   valid_min_s=0;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=1000;
   v1.putAtt("valid_max",ncShort,valid_max_s);
   scale_factor=1;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);
*/



return 0;
}







int meta_l1c_global(char* gridname,int16_t num_gridlines,NcFile* nc_output){
   string senstr,GATT_VAL1,prodstr,binstr,ifile_str;
   string date_created,y_create,m_create,d_create,t_create;
   int NVIEWS,NBANDS,NBANDS_POL;
   int32_t xbins,ybins;

   ybins=num_gridlines;
   prodstr=string(gridname);

   file_format format=getFormat(gridname);
   if(format.type==FT_HKT || format.type==FT_L1C) format.type=FT_OCIS; 

   if(format.type==FT_SPEXONE){
                         senstr = "SPEXONE";
                         GATT_VAL1="PACE SPEXone Level-1C Data";
                         binstr="12";
                         xbins=25;
                         NVIEWS=5;
                         NBANDS=400;
                         NBANDS_POL=50;
                     }
   else if(format.type==FT_HARP2){
                         senstr = "HARP2";
                         GATT_VAL1="PACE HARP2 Level-1C Data";
                         binstr="228";
                         xbins=457;
                         NVIEWS=90;
                         NBANDS=1;
                         NBANDS_POL=1;
                     }
   else if(format.type==FT_OCIS || format.type==FT_OCIL1B){
                         senstr = "OCI";
                         GATT_VAL1="PACE OCI Level-1C Data";
                         binstr="259";
                         xbins=519;
                         NVIEWS=2;
                         NBANDS=239;//249 originally no polarization bands
                     }
   else//OCIS
                     {
                         senstr = "OCI";
                         GATT_VAL1="PACE OCI Level-1C Data";
                         binstr="259";
                         xbins=519;
                         NVIEWS=2;
                         NBANDS=239;//249 ORIGINALLY no polarization bands
                     }


   //views
   float views[NVIEWS];
   if(format.type==FT_SPEXONE){
      views[0]=-58;
      views[1]=-20;
      views[2]=0;
      views[3]=20;
      views[4]=58;
   }
   else if (format.type==FT_OCIS){
      views[0]=-20;//OCI
      views[1]=20;
   }
   else if (format.type==FT_HARP2){
    cout<<"# views TBD.....ERROR.."<<endl;
    exit(1);
   }
   else{//OCI default in HKT---
    views[0]=-20;//OCI
    views[1]=20;
   }    
 // creation date---
// Current date/time based on current system
   time_t now = time(0);
  // Convert now to tm struct for UTC
   tm* gmtm = gmtime(&now);
   if (gmtm != NULL) {
  //   cout << "The UTC date and time is: " << asctime(gmtm) << endl;
   }
   else {
    cerr << "Failed to get the UTC date and time" << endl;
    return EXIT_FAILURE;
   }

   y_create=std::to_string(1900 + gmtm->tm_year);
   m_create=std::to_string(1 + gmtm->tm_mon);
   d_create=std::to_string(gmtm->tm_mday);
   t_create=std::to_string(gmtm->tm_hour) +":"+std::to_string(gmtm->tm_min)+":"+std::to_string(gmtm->tm_sec);
   date_created=y_create+"-"+m_create+"-"+d_create+"T"+t_create+"Z";

//global attributes---

   nc_output->putAtt("title", GATT_VAL1);
   nc_output->putAtt("instrument", senstr);
   nc_output->putAtt("processing_version","");
   nc_output->putAtt("Conventions","CF-1.8 ACDD-1.3");
   nc_output->putAtt("institution","NASA Goddard Space Flight Center, Ocean Biology Processing Group");
   nc_output->putAtt("license","http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/");
   nc_output->putAtt("naming_authority","gov.nasa.gsfc.sci.oceancolor");
   nc_output->putAtt("keywords_vocabulary","NASA Global Change Master Directory (GCMD) Science Keywords");
   nc_output->putAtt("stdname_vocabulary","NetCDF Climate and Forecast (CF) Metadata Convention");
   nc_output->putAtt("creator_name","NASA/GSFC");
   nc_output->putAtt("creator_email","data@oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("creator_url","http://oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("project","PACE Project");
   nc_output->putAtt("publisher_name","NASA/GSFC");
   nc_output->putAtt("publisher_email","data@oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("publisher_url","http://oceancolor.gsfc.nasa.gov");
   nc_output->putAtt("processing_level","L1C");
   nc_output->putAtt("cdm_data_type","swath");
   nc_output->putAtt("orbit_number","xxx");
   nc_output->putAtt("history","");
   nc_output->putAtt("CDL_version_date","2021-09-10");
   nc_output->putAtt("product_name",prodstr);
   nc_output->putAtt("startDirection","");
   nc_output->putAtt("endDirection","");
   nc_output->putAtt("time_coverage_start","");
   nc_output->putAtt("time_coverage_end","");
   nc_output->putAtt("date_created",date_created);
   nc_output->putAtt("sun_earth_distance","0.990849042172323");
   nc_output->putAtt("terrain_data_source","");
   nc_output->putAtt("spectral_response_function","");
   nc_output->putAtt("systematic_uncertainty_model","");
   nc_output->putAtt("nadir_bin",binstr);
   nc_output->putAtt("bin_size_at_nadir","5.2km2");

   NcDim ydim=nc_output->addDim("bins_along_track",ybins);
   NcDim xdim=nc_output->addDim("bins_across_track",xbins);
   NcDim vdim=nc_output->addDim("number_of_views", NVIEWS);
   NcDim idim=nc_output->addDim("intensity_bands_per_view",NBANDS);

//dims
   std::vector<NcDim> dimvec2_geo;
   dimvec2_geo.push_back(ydim);
   dimvec2_geo.push_back(xdim);
   std::vector<NcDim> dimvec2_rad;
   dimvec2_rad.push_back(vdim);
   dimvec2_rad.push_back(idim);
   std::vector<NcDim> dimvec3;
   dimvec3.push_back(ydim);
   dimvec3.push_back(xdim);
   dimvec3.push_back(vdim);
   std::vector<NcDim> dimvec4;
   dimvec4.push_back(ydim);
   dimvec4.push_back(xdim);
   dimvec4.push_back(vdim);
   dimvec4.push_back(idim);

//groups
   NcGroup svb_grp = nc_output->addGroup("sensor_views_bands");
   NcGroup ba_grp = nc_output->addGroup("bin_attributes");
   NcGroup geo_grp = nc_output->addGroup("geolocation_data");
   NcGroup od_grp = nc_output->addGroup("observation_data");
    //vars
   NcVar v1=svb_grp.addVar("view_angles",ncFloat,vdim);
   string longName="Along-track view angles for sensor";
   v1.putAtt( "long_name",longName);
   string units="degrees";
   v1.putAtt( "units",units);
   float FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   float valid_min=-89;
   v1.putAtt("valid_min",ncFloat,valid_min);
   float valid_max=89;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(views);    

 /*  NcVar v2=svb_grp.addVar("intensity_wavelengths",ncFloat,dimvec2_rad);//this is infull meta
   longName="Intensity field center wavelengths at each view";
   v1.putAtt( "long_name",longName);
   units="nm";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=320;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=2260;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(fwave_view);
*/
   if(format.type==FT_HARP2 || format.type==FT_SPEXONE){

       NcDim pdim=nc_output->addDim("polarization_bands_per_view",NBANDS_POL);
       std::vector<NcDim> dimvec2b_rad;
       dimvec2b_rad.push_back(vdim);
       dimvec2b_rad.push_back(pdim);
       std::vector<NcDim> dimvec4b;
       dimvec4b.push_back(ydim);
       dimvec4b.push_back(xdim);
       dimvec4b.push_back(vdim);
       dimvec4b.push_back(pdim);

        //grp svb
 /*      float fwave_view_pol[NVIEWS][NBANDS_POL],fwhm_view_pol[NVIEWS][NBANDS_POL],Fobar_view_pol[NVIEWS][NBANDS_POL];

       for(int i=0;i<NVIEWS;i++){
         for(int j=0;j<NBANDS_POL;j++){
             fwave_view_pol[i][j]=fwave[j];
             fwhm_view_pol[i][j]=fwhm[j];
             Fobar_view_pol[i][j]=Fobar[j];
       }}
       v1=svb_grp.addVar("polarization_wavelengths",ncFloat,dimvec2b_rad);
       string longName2="Polarization field wavelengths at each view";
       v1.putAtt( "long_name",longName2);
       string units2="nm";
       v1.putAtt( "units",units2);
       float FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       float valid_min_p=320;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       float valid_max_p=2260;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
       v1.putVar(fwave_view_pol);

       v1=svb_grp.addVar("polarization_bandpasses",ncFloat,dimvec2b_rad);
       longName2="Polarization field bandpasses at each view";
       v1.putAtt( "long_name",longName2);
       units2="nm";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=2.5;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=100;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
       v1.putVar(fwhm_view_pol);

       v1=svb_grp.addVar("polarization_F0",ncFloat,dimvec2b_rad);
       longName2="Polarization band solar irradiance";

       v1.putAtt( "long_name",longName2);
       units2="W m^-2 µm^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=0;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=900;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
       v1.putVar(Fobar_view_pol);
*/

   if(format.type==FT_HARP2){
       //grp od_grp
 /*      v1 = od_grp.addVar("Q", ncFloat,dimvec4b);
       longName2="Q Stokes vector component";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=-800;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);

       v1 = od_grp.addVar("Q_noise", ncFloat,dimvec4b);
       longName2="Random noise of Q in bin";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=0;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);

       v1 = od_grp.addVar("U", ncFloat,dimvec4b);
       longName2="U Stokes vector component";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=-800;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);

       v1 = od_grp.addVar("U_noise", ncFloat,dimvec4b);
       longName2="Random noise of U in bin";
       v1.putAtt( "long_name",longName2);
       units2="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units2);
       FillValue_p=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue_p);
       valid_min_p=0;
       v1.putAtt("valid_min",ncFloat,valid_min_p);
       valid_max_p=800;
       v1.putAtt("valid_max",ncFloat,valid_max_p);
*/ 
     }

/*   v1 = od_grp.addVar("DOLP", ncFloat,dimvec4b);
   longName2="Degree of linear polarization";
   v1.putAtt( "long_name",longName2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=1;
   v1.putAtt("valid_max",ncFloat,valid_max_p);

   v1 = od_grp.addVar("DOLP_noise", ncFloat,dimvec4b);
   longName2="Random noise of DOLP in bin";

   v1.putAtt( "long_name",longName2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=1;
   v1.putAtt("valid_max",ncFloat,valid_max_p);

   v1 = od_grp.addVar("AOLP", ncFloat,dimvec4b);
   longName2="Angle of linear polarization";
   v1.putAtt( "long_name",longName2);
   units2="degrees";
   v1.putAtt("units",units2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=360;
   v1.putAtt("valid_max",ncFloat,valid_max_p);

   v1 = od_grp.addVar("AOLP_noise", ncFloat,dimvec4b);
   longName2="Random noise of AOLP in bin";
   v1.putAtt( "long_name",longName2);
   units2="degrees";
   v1.putAtt("units",units2);
   FillValue_p=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue_p);
   valid_min_p=0;
   v1.putAtt("valid_min",ncFloat,valid_min_p);
   valid_max_p=360;
   v1.putAtt("valid_max",ncFloat,valid_max_p);
*/
    }
/*   v1=svb_grp.addVar("intensity_bandpasses",ncFloat,dimvec2_rad);
   longName="Intensity field bandpasses at each view";
   v1.putAtt( "long_name",longName);
   units="nm";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=2.5;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=100;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(fwhm_view);

   v1=svb_grp.addVar("intensity_F0",ncFloat,dimvec2_rad);
   longName="Intensity band solar irradiance";
   v1.putAtt( "long_name",longName);
   units="W m^-2 µm^-1";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=0;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=900;
   v1.putAtt("valid_max",ncFloat,valid_max);
   v1.putVar(Fobar_view);
*/
   v1=ba_grp.addVar("nadir_view_time",ncDouble,ydim);
   float FillValue2;
   longName="Time bin was viewed at nadir view";
   v1.putAtt( "long_name",longName);
   units="seconds";
   v1.putAtt( "units",units);
   FillValue2=-999.;
   v1.putAtt( "_FillValue",ncDouble,FillValue2);
   double valid_min_d=0;
   v1.putAtt("valid_min",ncDouble,valid_min_d);
   double valid_max_d=86400.999;
   v1.putAtt("valid_max",ncDouble,valid_max_d);

   v1=ba_grp.addVar("view_time_offsets",ncDouble,dimvec3);
   longName="Time offsets of views from nadir view";
   v1.putAtt( "long_name",longName);
   units="seconds";
   v1.putAtt( "units",units);
   FillValue2=-999;
   v1.putAtt( "_FillValue",ncDouble,FillValue2);
   valid_min_d=-200;
   v1.putAtt("valid_min",ncDouble,valid_min_d);
   valid_max_d=200;
   v1.putAtt("valid_max",ncDouble,valid_max_d);

   v1 = geo_grp.addVar("latitude", ncFloat,dimvec2_geo);
   longName="Latitudes of bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees_north";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=-90;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=90;
   v1.putAtt("valid_max",ncFloat,valid_max);

   v1 = geo_grp.addVar("longitude", ncFloat,dimvec2_geo);
   longName="Longitudes of bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees_east";
   v1.putAtt( "units",units);
   FillValue=-999;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   valid_min=-180;
   v1.putAtt("valid_min",ncFloat,valid_min);
   valid_max=180;
   v1.putAtt("valid_max",ncFloat,valid_max);

   v1 = geo_grp.addVar("altitude", ncFloat,dimvec2_geo);
   longName="Altitude at bin locations";
   v1.putAtt( "long_name",longName);
   units="meters";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncFloat,FillValue);
   short valid_min_s=-1000;
   v1.putAtt("valid_min",ncFloat,valid_min_s);
   short valid_max_s=10000;
   v1.putAtt("valid_max",ncFloat,valid_max_s);
   float scale_factor=1;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   float add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);

   v1 = geo_grp.addVar("altitude_variability", ncShort,dimvec2_geo);
   longName="RMS variability of altitude at bin locations";
   v1.putAtt( "long_name",longName);
   units="meters";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncShort,FillValue);
   valid_min_s=0;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=1000;
   v1.putAtt("valid_max",ncShort,valid_max_s);
   scale_factor=1;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);

   v1 = geo_grp.addVar("sensor_azimuth", ncShort,dimvec3);
   longName="Sensor azimuth angles at bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees from north";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncShort,FillValue);
   valid_min_s=-18000;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=18000;
   v1.putAtt("valid_max",ncShort,valid_max_s);
   scale_factor=0.01;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);

   v1 = geo_grp.addVar("sensor_zenith", ncShort,dimvec3);
   longName="Sensor zenith angles at bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncShort,FillValue);
   valid_min_s=0;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=18000;
   v1.putAtt("valid_max",ncShort,valid_max_s);
   scale_factor=0.01;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);

   v1 = geo_grp.addVar("solar_azimuth", ncShort,dimvec3);
   longName="Solar azimuth angle at bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees from north";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncShort,FillValue);
   valid_min_s=-18000;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=18000;
   v1.putAtt("valid_max",ncShort,valid_max_s);
   scale_factor=0.01;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);

   v1 = geo_grp.addVar("solar_zenith", ncShort,dimvec3);
   longName="Solar zenith angle at bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncShort,FillValue);
   valid_min_s=0;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=18000;
   v1.putAtt("valid_max",ncShort,valid_max_s);
   scale_factor=0.01;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);

   v1 = geo_grp.addVar("scattering_angle", ncShort,dimvec3);
   longName="Scattering angle at bin locations";
   v1.putAtt( "long_name",longName);
   units="degrees";
   v1.putAtt( "units",units);
   FillValue=-32768;
   v1.putAtt( "_FillValue",ncShort,FillValue);
   valid_min_s=0;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=18000;
   v1.putAtt("valid_max",ncShort,valid_max_s);
   scale_factor=0.01;
   v1.putAtt("scale_factor",ncFloat,scale_factor);
   add_offset=0;
   v1.putAtt("add_offset",ncFloat,add_offset);

   v1 = od_grp.addVar("obs_per_view", ncShort,dimvec3);
   longName="Observations contributing to bin from each view";
   v1.putAtt( "long_name",longName);
   valid_min_s=0;
   v1.putAtt("valid_min",ncShort,valid_min_s);
   valid_max_s=999;
   v1.putAtt("valid_max",ncShort,valid_max_s);

   if(format.type==FT_OCIL1B || format.type==FT_OCIS || format.type==FT_HARP2){   
       uint8_t FillValue3,valid_min2,valid_max2;  
       v1 = od_grp.addVar("QC_bitwise", ncUbyte,dimvec4);
       longName="bitwise quality indicator";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue3=255;
       v1.putAtt( "_FillValue",ncUbyte,FillValue3);
       valid_min2=0;
       v1.putAtt("valid_min",ncUbyte,valid_min2);
       valid_max2=255;
       v1.putAtt("valid_max",ncUbyte,valid_max2);

       v1 = od_grp.addVar("QC", ncUbyte,dimvec4);
       longName="quality indicator";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue3=255;
       v1.putAtt( "_FillValue",ncUbyte,FillValue3);
       valid_min2=0;
       v1.putAtt("valid_min",ncUbyte,valid_min2);
       valid_max2=10;
       v1.putAtt("valid_max",ncUbyte,valid_max2);

       v1 = od_grp.addVar("I", ncFloat,dimvec4);
  //     longName="I Stokes vector component";
       longName="Reflectance";
       v1.putAtt( "long_name",longName);
   //    units="W m^-2 sr^-1 um^-1";
       units="dimensionless";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=0;
       v1.putAtt("valid_min",ncFloat,valid_min);
//       valid_max=999;
       valid_max=1.3;       
       v1.putAtt("valid_max",ncFloat,valid_max);

       v1 = od_grp.addVar("I_noise", ncFloat,dimvec4);
       longName="Random noise of I in bin";
       v1.putAtt( "long_name",longName);
       units="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=0;
       v1.putAtt("valid_min",ncFloat,valid_min);
       valid_max=800;
       v1.putAtt("valid_max",ncFloat,valid_max);
       
     }
   else if(format.type==FT_SPEXONE)//spexone
    {
  /*     v1 = od_grp.addVar("QC_polsample_bitwise", ncUbyte,dimvec4);
       longName="bitwise quality indicator";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue3=255;
       v1.putAtt( "_FillValue",ncUbyte,FillValue3);
       valid_min2=0;
       v1.putAtt("valid_min",ncUbyte,valid_min2);
       valid_max2=255;
       v1.putAtt("valid_max",ncUbyte,valid_max2);

       v1 = od_grp.addVar("QC_polsample", ncUbyte,dimvec4);
       longName="quality indicator";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue3=255;
       v1.putAtt( "_FillValue",ncUbyte,FillValue3);
       valid_min2=0;
       v1.putAtt("valid_min",ncUbyte,valid_min2);
       valid_max2=10;
       v1.putAtt("valid_max",ncUbyte,valid_max2);

       v1 = od_grp.addVar("I_polsample", ncFloat,dimvec4);
       longName="I Stokes vector component at polarimeter spectral sampling";
       v1.putAtt( "long_name",longName);
       units="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=0;
       v1.putAtt("valid_min",ncFloat,valid_min);
       valid_max=999;
       v1.putAtt("valid_max",ncFloat,valid_max);

       v1 = od_grp.addVar("I_polsample_noise", ncFloat,dimvec4);
       longName="Random noise of I_polsample in bin";
       v1.putAtt( "long_name",longName);
       units="W m^-2 sr^-1 um^-1";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=0;
       v1.putAtt("valid_min",ncFloat,valid_min);
       valid_max=800;
       v1.putAtt("valid_max",ncFloat,valid_max);

       v1 = od_grp.addVar("Q_over_I", ncFloat,dimvec4);
       longName="Q over I (little q) Stokes vector component";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=-800;
       v1.putAtt("valid_min",ncFloat,valid_min);
       valid_max=800;
       v1.putAtt("valid_max",ncFloat,valid_max);

       v1 = od_grp.addVar("Q_over_I_noise", ncFloat,dimvec4);
       longName="Random noise of Q over I in bin";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=0;
       v1.putAtt("valid_min",ncFloat,valid_min);
       valid_max=800;
       v1.putAtt("valid_max",ncFloat,valid_max);

       v1 = od_grp.addVar("U_over_I", ncFloat,dimvec4);
       longName="U over I (little q) Stokes vector component";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=-800;
       v1.putAtt("valid_min",ncFloat,valid_min);
       valid_max=800;
       v1.putAtt("valid_max",ncFloat,valid_max);

       v1 = od_grp.addVar("U_over_I_noise", ncFloat,dimvec4);
       longName="Random noise of U over I in bin";
       v1.putAtt( "long_name",longName);
       units="unitless";
       v1.putAtt( "units",units);
       FillValue=-999;
       v1.putAtt( "_FillValue",ncFloat,FillValue);
       valid_min=0;
       v1.putAtt("valid_min",ncFloat,valid_min);
       valid_max=800;
       v1.putAtt("valid_max",ncFloat,valid_max);
*/
    }



return 0;
}


int  bintime_l1c(filehandle* l1file,l1str *l1rec,bin_str *binstr,short **gdindex,double scantime, NcFile *nc_output)
{ 
   short gd_row=-1,gd_col=-1;
   int npix=l1file->npix;           
   int view;
   int nadpix=(npix-1)/2;
   double timefill,mintime,maxtime;
 
   NcGroup ba_grp=nc_output->getGroup("bin_attributes");
   NcVar v1=ba_grp.getVar("nadir_view_time");    
   NcVarAtt a1=v1.getAtt("_FillValue");//root group
   a1.getValues(&timefill);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&mintime);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxtime);

   //ASSUMING same time for each line!!
   for(int pix=0;pix<npix;pix++)
                  {
                      gd_row=gdindex[pix][0]-1;
                      gd_col=gdindex[pix][1]-1;

            //1 view at the time for OCI!!!!
                    if(l1rec->senz[pix]!=binstr->fillval1 && l1rec->senz[nadpix]<0){
                       view=0;
                    }
                    else if(l1rec->senz[pix]!=binstr->fillval1 && l1rec->senz[nadpix]>0){
                        view=1;
                    }
                 //Cumulative values---

                   if(scantime>=binstr->tini_l1c && scantime<=binstr->tend_l1c && gd_row>=0 && gd_col>=0 && gd_row<=binstr->num_gridlines-1 && gd_col<=binstr->nbinx-1 && scantime!=timefill && scantime>=mintime &&scantime<=maxtime){ //and row<=ybin and col<xbin, I need time counter 
                      binstr->time_l1b[gd_row][gd_col][view] = binstr->time_l1b[gd_row][gd_col][view] + scantime;                                        
                   }              

                  }   


return 0;
}    


int bin_l1c(filehandle *l1file,l1str *l1rec,bin_str *binl1c,short **gdindex,NcFile *nc_output, double scantime){
   int32_t ibp=0,sb=0;
   short gd_row=-1,gd_col=-1;
   int npix=l1file->npix;
   int nbands=l1file->nbands;
   int view;
   int nadpix=(npix-1)/2;
   float term1,term2,term3,sca_pix,cosu,cose;
   float scale_all[6],offset_all[6];
   short minval1_all[6],maxval1_all[6];
   float minval2_all[6],maxval2_all[6];

   NcGroup geo_grp=nc_output->getGroup("geolocation_data");
   NcGroup od_grp=nc_output->getGroup("observation_data");
   //altitude
   NcVar v1=geo_grp.getVar("altitude");
   NcVarAtt a1=v1.getAtt("scale_factor");//root group
   a1.getValues(&scale_all[0]);
   a1=v1.getAtt("add_offset");//root group
   a1.getValues(&offset_all[0]);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&minval1_all[0]);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxval1_all[0]);

  //sensor azimuth
   v1=geo_grp.getVar("sensor_azimuth");
   a1=v1.getAtt("scale_factor");//root group
   a1.getValues(&scale_all[1]);
   a1=v1.getAtt("add_offset");//root group
   a1.getValues(&offset_all[1]);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&minval1_all[1]);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxval1_all[1]);

   //sensor zenith
   v1=geo_grp.getVar("sensor_zenith");
   a1=v1.getAtt("scale_factor");//root group
   a1.getValues(&scale_all[2]);
   a1=v1.getAtt("add_offset");//root group
   a1.getValues(&offset_all[2]);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&minval1_all[2]);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxval1_all[2]);

   minval1_all[2]=(minval1_all[2]-offset_all[2])/scale_all[2];
   maxval1_all[2]=(maxval1_all[2]-offset_all[2])/scale_all[2];

   //solar zenith
   v1=geo_grp.getVar("solar_azimuth");
   a1=v1.getAtt("scale_factor");//root group
   a1.getValues(&scale_all[3]);
   a1=v1.getAtt("add_offset");//root group
   a1.getValues(&offset_all[3]);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&minval1_all[3]);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxval1_all[3]);

   //solar zenith
   v1=geo_grp.getVar("solar_zenith");
   a1=v1.getAtt("scale_factor");//root group
   a1.getValues(&scale_all[4]);
   a1=v1.getAtt("add_offset");//root group
   a1.getValues(&offset_all[4]);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&minval1_all[4]);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxval1_all[4]);

   //solar zenith
   v1=geo_grp.getVar("scattering_angle");
   a1=v1.getAtt("scale_factor");//root group
   a1.getValues(&scale_all[5]);
   a1=v1.getAtt("add_offset");//root group
   a1.getValues(&offset_all[5]);
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&minval1_all[5]);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxval1_all[5]);

   //I
   v1=od_grp.getVar("I");                             
   a1=v1.getAtt("valid_min");//root group
   a1.getValues(&minval2_all[0]);
   a1=v1.getAtt("valid_max");//root group
   a1.getValues(&maxval2_all[0]);
    
   for(int pix=0;pix<npix;pix++){                  
                      gd_row=gdindex[pix][0]-1;
                      gd_col=gdindex[pix][1]-1;
                   

             //2-D arrays---needing to add fillvalue != for l1rec->alt[pix] 
                    if(scantime>=binl1c->tini_l1c && scantime<=binl1c->tend_l1c && gd_row>=0 && gd_col>=0 && gd_row<=binl1c->num_gridlines-1 && gd_col<=binl1c->nbinx-1 && l1rec->lat[pix]!=binl1c->fillval2 && l1rec->lon[pix]!=binl1c->fillval2)
                    {
                           if((rand()%10000-1000-offset_all[0])/scale_all[0]!=binl1c->fillval1)//ADD MIN/MAX CONSTRAINT
                           {
                             binl1c->nrec_2D[gd_row][gd_col] += 1;                              
                             binl1c->alt_2D[gd_row][gd_col] = binl1c->alt_2D[gd_row][gd_col]+(rand()%10000-1000-offset_all[0])/scale_all[0];//l1rec->height[pix];
               //              binl1c->alt_diff2[gd_row][gd_col] = binl1c->alt_diff2[gd_row][gd_col]+fudge_factor;//altitude reader in l1???? + l1rec->sca[pix];
                           }
                                  
                      binl1c->inpix++;              
                    }                    
                    else 
                    {                                                   
                      binl1c->outpix++;     
                    }

                  //1 view at the time for OCI!!!!
                    if(l1rec->senz[pix]!=binl1c->fillval1 && l1rec->senz[nadpix]<0){
                        view=0;
                    }
                    else if(l1rec->senz[pix]!=binl1c->fillval1 && l1rec->senz[nadpix]>0){
                        view=1;
                    }
                                       
                 //Cumulative values---

                   if(scantime>=binl1c->tini_l1c && scantime<=binl1c->tend_l1c && gd_row>=0 && gd_col>=0 && gd_row<=binl1c->num_gridlines-1 && gd_col<=binl1c->nbinx-1 && l1rec->senz[pix]!=binl1c->fillval1 && l1rec->senz[pix]>=minval1_all[2] && l1rec->senz[pix]<=maxval1_all[2]){
                      binl1c->nrec_3D[gd_row][gd_col][view] += 1; 
                      binl1c->sena_3D[gd_row][gd_col][view] = binl1c->sena_3D[gd_row][gd_col][view] + (l1rec->sena[pix]-offset_all[1])/scale_all[1];
                      binl1c->senz_3D[gd_row][gd_col][view] = binl1c->senz_3D[gd_row][gd_col][view] + (l1rec->senz[pix]-offset_all[2])/scale_all[2];  //sensor zenith goes from 0 to 180
                      binl1c->suna_3D[gd_row][gd_col][view] = binl1c->suna_3D[gd_row][gd_col][view] + (l1rec->sola[pix]-offset_all[3])/scale_all[3];
                      binl1c->sunz_3D[gd_row][gd_col][view] = binl1c->sunz_3D[gd_row][gd_col][view] + (l1rec->solz[pix]-offset_all[4])/scale_all[4];

                      cose=cos((l1rec->senz[pix]+180)*M_PI/180);
                      cosu=cos((l1rec->solz[pix])*M_PI/180);   
                      term1=cose*cosu;
                      term2=sqrt((1-cose*cose))*sqrt((1-cosu*cosu));
                      term3=cos((l1rec->sena[pix]+180-l1rec->sola[pix])*M_PI/180);
                      sca_pix=acos(term1+term2*term3)*180/M_PI;                    
                      binl1c->sca_3D[gd_row][gd_col][view] = binl1c->sca_3D[gd_row][gd_col][view]+(sca_pix-offset_all[5])/scale_all[5];  
                     }
                  
                     ibp=pix*nbands;
  
                     for(sb=0;sb<nbands;sb++){

                           if (scantime>=binl1c->tini_l1c && scantime<=binl1c->tend_l1c && gd_row>=0 && gd_col>=0 && gd_row<=binl1c->num_gridlines-1 && gd_col<=binl1c->nbinx-1 && l1rec->Lt[ibp]!=binl1c->fillval2 && l1rec->Lt[ibp]>=minval2_all[0] && l1rec->Lt[ibp]<=maxval2_all[0]) 
                           {
                               binl1c->I_4D[gd_row][gd_col][view][sb] = binl1c->I_4D[gd_row][gd_col][view][sb] + l1rec->Lt[ibp];
                               binl1c->nrec_4D_band[gd_row][gd_col][view][sb] += 1;
                               binl1c->nrec_3D_view[gd_row][gd_col][view] += 1;

                               if(gd_row>binl1c->num_gridlines-1 || gd_col>binl1c->nbinx-1)
                               {
                                cout<<"ERROR IN BINNING Lt/rhot --gd_row/gd-col outside the grid limits"<<endl;
                                cout<<"row.."<<gd_row<<"gd_col.."<<gd_col<<endl; 
                                
                               }
                              //missing QC an QC_bitwise variables
                           }     
                           ibp++;                       
                          }
                   }
                  
return 0;                       
}



int check_swath_time(filehandle* l1file,const char *l1c_grid,bin_str *binl1c){
   string l1c_str=l1c_grid;
   string l1b_str=l1file->name;

   NcFile* nc_l1cgrid;
   try {
        nc_l1cgrid = new NcFile(l1c_grid, NcFile::read);
        }
   catch (NcException& e) {
          e.what();
          cerr << "l1cgen l1c_pflag= 8:: Failure reading L1C grid: "
          + l1c_str << endl;
          exit(1);
   }


   NcDim yd=nc_l1cgrid->getDim("bins_along_track");
   int32_t num_gridlines=yd.getSize();
   double *time_nad_l1c = (double*)calloc(num_gridlines,sizeof(double));
   NcGroup ba_grp=nc_l1cgrid->getGroup("bin_attributes");
   NcVar v1=ba_grp.getVar("nadir_view_time");
   v1.getVar(time_nad_l1c);

   NcFile* nc_l1b;
   try {
        nc_l1b = new NcFile(l1file->name, NcFile::read);
        }
   catch (NcException& e) {
          e.what();
          cerr << "l1cgen l1c_pflag= 8:: Failure reading L1B file: "
          + l1b_str << endl;
          exit(1);
   }

    //grab nadir time for the grid
   yd=nc_l1b->getDim("number_of_scans");
   int32_t nscans=yd.getSize();

   binl1c->nscans=nscans;

   double *time_l1b = (double*)calloc(nscans,sizeof(double));
   NcGroup sla_grp=nc_l1b->getGroup("scan_line_attributes");
   v1=sla_grp.getVar("time");
   v1.getVar(time_l1b);


//check time limits
   double  tini_l1c=time_nad_l1c[0];
   double tend_l1c=time_nad_l1c[num_gridlines-1];
   double tini_l1b=time_l1b[0];
   double tend_l1b=time_l1b[nscans-1];  


   binl1c->tini_l1c=tini_l1c;
   binl1c->tend_l1c=tend_l1c;
   binl1c->tini_l1b=tini_l1b;
   binl1c->tend_l1b=tend_l1b;

   int time_flag=-1;
   if(tini_l1b>=tini_l1c && tini_l1b<=tend_l1c || tend_l1b>=tini_l1c && tend_l1b<=tend_l1c)
    {
      time_flag=0;

      cout<<"OK--L1B granule is INSIDE of the L1C grid--"<<endl;                                     
    }    
   else
    {
   time_flag=1;
   cout<<"ERROR--L1B granule is outside of the L1C grid--"<<endl;
   exit(1);
    }


   if(time_nad_l1c!=nullptr)
     delete [] (time_nad_l1c);
   time_nad_l1c=nullptr;

   if(time_l1b!=nullptr)
     delete [] (time_l1b);
   time_l1b=nullptr;

   nc_l1cgrid->close();
   nc_l1b->close();

return time_flag;
}



int open_l1c(const char *ifile_l1c,size_t *ybins,size_t *xbins,float **lat_gd,float **lon_gd){   
   std::string str;
   std::string ifile_str;
   string gridname, azeast_name;
   std::string fname_out, pathstr, senstr, monstr, daystr, yearstr, prodstr, gdstr, swtstr, swtnum, extstr, granstr,timestr,azstr,missionstr,ofilestr;


   cout<<"Opening L1C grid........................................................................."<<endl;
   string l1c_str(ifile_l1c);

   NcFile* nc_l1cgrid;

           try {
                 nc_l1cgrid = new NcFile(ifile_l1c, NcFile::read);
               }
           catch (NcException& e) {
                   e.what();
                 cerr << "l1cgen :: Failure reading L1C grid: "
                  + l1c_str << endl;
                  exit(1);
                 }
   NcGroup geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   NcVar v1=geo_grp.getVar("latitude");
   v1.getVar(&lat_gd[0][0]);
   v1=geo_grp.getVar("longitude");
   v1.getVar(&lon_gd[0][0]);



   nc_l1cgrid->close();

    return 0;
}





//fred/me implementation
int search_l1c(filehandle* l1file, l1str *l1rec,bin_str *binl1c,short **gdindex){
   int32_t num_gridlines, nbinx;
   float gnvm;
   float gnvec[3],gvec[3],bvec[3];
   int flag_out=-1;
   size_t pix;
   int32_t i;
   size_t inpix=0;
   size_t num_pixels;
   int irow=-1,col=-1;
   float bmcm,bm=100;
   double db;
   float c1,c2,c3;
   double fudge=0.00001,dotprod,dot_firstline,dot_lastline;   
   flag_out=-1;
        
   num_gridlines=binl1c->num_gridlines;
   nbinx=binl1c->nbinx;
   num_pixels = l1file->npix;

   db = (5.2)/6371/2; //Half of bin size in radians, resolution in km

    //big loop
      //dot product
   for(pix=0;pix<num_pixels;pix++)
        {
             bvec[0]=cos(l1rec->lon[pix]*M_PI/180)*cos(l1rec->lat[pix]*M_PI/180);
             bvec[1]=sin(l1rec->lon[pix]*M_PI/180)*cos(l1rec->lat[pix]*M_PI/180);
             bvec[2]=sin(l1rec->lat[pix]*M_PI/180);

          for(i=0;i<num_gridlines;i++)
          {
           if(l1rec->lat[pix]>90 || l1rec->lat[pix]<-90 || l1rec->lon[pix]<-180 || l1rec->lon[pix]>180){cout<<"latitude longitude pixel out of the boundaries.."<<endl; exit(1);}
              //normal vectors for L1C rows
           if(binl1c->lat_gd[i][nbinx-1]>90 || binl1c->lat_gd[i][nbinx-1]<-90 || binl1c->lon_gd[i][nbinx-1]<-180 || binl1c->lon_gd[i][nbinx-1]>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}
           if(binl1c->lat_gd[i][0]>90 || binl1c->lat_gd[i][0]<-90 || binl1c->lon_gd[i][0]<-180 || binl1c->lon_gd[i][0]>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}

           gnvec[0] = sin(binl1c->lon_gd[i][nbinx-1]*M_PI/180)*cos(binl1c->lat_gd[i][nbinx-1]*M_PI/180)*sin(binl1c->lat_gd[i][0]*M_PI/180) - sin(binl1c->lat_gd[i][nbinx-1]*M_PI/180)*sin(binl1c->lon_gd[i][0]*M_PI/180)*cos(binl1c->lat_gd[i][0]*M_PI/180);
           gnvec[1] = sin(binl1c->lat_gd[i][nbinx-1]*M_PI/180)*cos(binl1c->lon_gd[i][0]*M_PI/180)*cos(binl1c->lat_gd[i][0]*M_PI/180) - cos(binl1c->lon_gd[i][nbinx-1]*M_PI/180)*cos(binl1c->lat_gd[i][nbinx-1]*M_PI/180)*sin(binl1c->lat_gd[i][0]*M_PI/180);
           gnvec[2] = cos(binl1c->lon_gd[i][nbinx-1]*M_PI/180)*cos(binl1c->lat_gd[i][nbinx-1]*M_PI/180)*sin(binl1c->lon_gd[i][0]*M_PI/180)*cos(binl1c->lat_gd[i][0]*M_PI/180) - sin(binl1c->lon_gd[i][nbinx-1]*M_PI/180)*cos(binl1c->lat_gd[i][nbinx-1]*M_PI/180)*cos(binl1c->lon_gd[i][0]*M_PI/180)*cos(binl1c->lat_gd[i][0]*M_PI/180);


         //vector norm
           gnvm=sqrt(gnvec[0]*gnvec[0]+gnvec[1]*gnvec[1]+gnvec[2]*gnvec[2]);
           if(isnan(gnvm)==1){cout<<"NAN value for gnvm.."<<endl; exit(1);}
           if(gnvm==0){ cout<<"ERROR gnvm == 0--- WE CANT NORMALIZE..."<<endl;exit(1);}
              //normalization
           gnvec[0] = gnvec[0]/gnvm;
           gnvec[1] = gnvec[1]/gnvm;
           gnvec[2] = gnvec[2]/gnvm;

        //for each pixels      
          //dot prod, orbital normaliz and transposed by pixel vector   
           dotprod=gnvec[0]*bvec[0]+gnvec[1]*bvec[1]+gnvec[2]*bvec[2];                  

           if(i==0){
               dot_firstline=dotprod;
                }
           if(i==num_gridlines-1){
               dot_lastline=dotprod;
                }   


           if(dotprod-fudge<=db && dotprod+fudge>-db)
              {
            gdindex[pix][0]=i+1;//first found
              }

          } //end lines
       
         //for each pixels
          if(dot_firstline<=db && dot_lastline>-db)
              {
               inpix++;

            //find column
                     irow=gdindex[pix][0]-1;
                     if(irow<0){cout<<"ERROR icol in search_l1c..."<<"at pix#.."<<pix+1<<"and irow#.."<<irow<<endl;exit(1);}

               for(int j=0;j<nbinx;j++)
                 {
                 gvec[0]=cos(binl1c->lon_gd[irow][j]*M_PI/180)*cos(binl1c->lat_gd[irow][j]*M_PI/180);
                 gvec[1]=sin(binl1c->lon_gd[irow][j]*M_PI/180)*cos(binl1c->lat_gd[irow][j]*M_PI/180);
                 gvec[2]=sin(binl1c->lat_gd[irow][j]*M_PI/180);

                 c1=bvec[0]-gvec[0];
                 c2=bvec[1]-gvec[1];
                 c3=bvec[2]-gvec[2];

                 bmcm = sqrt(c1*c1+c2*c2+c3*c3);////bmcm only checked one pixel!!I
                 if(bmcm<bm){
                    bm=bmcm;
                    col=j+1;
                       }
                 } 
               if(col<1){cout<<"ERROR col in search_l1c..."<<"at pix#.."<<pix+1<<"and row#.."<<irow+1<<endl;exit(1);}

               gdindex[pix][1]=col;
               if(irow==num_gridlines-1 && l1rec->lat[pix] < 85.0) {
                  gdindex[pix][0] = -1;
                  gdindex[pix][1] = -1;
               }
               bm=100;
               col=-1;
           }
           else
           {
            gdindex[pix][0]=-1;//actually these are row/col not irow/icol like in my searching algos!!!
            gdindex[pix][1]=-1;
           }

        }//end pixels

       flag_out=-1;
       for(pix=0;pix<num_pixels;pix++){
               if(gdindex[pix][0]>0 && gdindex[pix][1]>0){
                         flag_out=0;
               }

               else{
               //    cout<<"THIS LINE WILL BE SKIPPED -- NO PIXELS BINNED.............."<<endl;
               } 
       }
                    
    if(flag_out==0)
       return 0;
    else
       return 1;
}

//fred/richard implemtation
//one pixel at the time----
int search2_l1c(size_t ybins,size_t xbins,float lat,float lon,float **lat_gd,float **lon_gd,short *row,short *col){
   int32_t num_gridlines, nbinx,ic;
   short gdrow=-1,gdcol=-1;
   float pvec[3],gnvm,db,pdotgn_first,pdotgn_last;
   float **gnvec=nullptr,**cnvec=nullptr,**dc=nullptr;;
   int flag_out=-1;
   int32_t i;
   size_t inpix=0,outpix=0;
   float dotprod,dotprod2,dcm;
   flag_out=-1;

   num_gridlines=ybins;
   nbinx=xbins;

   db = (5.2)/6371/2; //Half of bin size in radians, resolution in km

   //fred/richard implemtation
   gnvec=allocate2d_float(3,num_gridlines);
   cnvec=allocate2d_float(3,num_gridlines);
   dc=allocate2d_float(3,num_gridlines);

   ic = nbinx/2;


   for(i=0;i<num_gridlines;i++){
          //normal vectors for L1C rows
             if(lat_gd[i][nbinx-1]>90 || lat_gd[i][nbinx-1]<-90 || lon_gd[i][nbinx-1]<-180 || lon_gd[i][nbinx-1]>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}
             if(lat_gd[i][0]>90 || lat_gd[i][0]<-90 || lon_gd[i][0]<-180 || lon_gd[i][0]>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}
         //compute normal vector for each rowgrid----
              gnvec[0][i] = sin(lon_gd[i][nbinx-1]*M_PI/180)*cos(lat_gd[i][nbinx-1]*M_PI/180)*sin(lat_gd[i][0]*M_PI/180) - sin(lat_gd[i][nbinx-1]*M_PI/180)*sin(lon_gd[i][0]*M_PI/180)*cos(lat_gd[i][0]*M_PI/180);
              gnvec[1][i] = sin(lat_gd[i][nbinx-1]*M_PI/180)*cos(lon_gd[i][0]*M_PI/180)*cos(lat_gd[i][0]*M_PI/180) - cos(lon_gd[i][nbinx-1]*M_PI/180)*cos(lat_gd[i][nbinx-1]*M_PI/180)*sin(lat_gd[i][0]*M_PI/180);
              gnvec[2][i] = cos(lon_gd[i][nbinx-1]*M_PI/180)*cos(lat_gd[i][nbinx-1]*M_PI/180)*sin(lon_gd[i][0]*M_PI/180)*cos(lat_gd[i][0]*M_PI/180) - sin(lon_gd[i][nbinx-1]*M_PI/180)*cos(lat_gd[i][nbinx-1]*M_PI/180)*cos(lon_gd[i][0]*M_PI/180)*cos(lat_gd[i][0]*M_PI/180);

              gnvm=sqrt(gnvec[0][i]*gnvec[0][i]+gnvec[1][i]*gnvec[1][i]+gnvec[2][i]*gnvec[2][i]);
              if(isnan(gnvm)==1){cout<<"NAN value for gnvm.."<<endl; exit(1);}
              if(gnvm==0){ cout<<"ERROR gnvm == 0--- WE CANT NORMALIZE..."<<endl;exit(1);}
              //normalization
              gnvec[0][i] = gnvec[0][i]/gnvm;
              gnvec[1][i] = gnvec[1][i]/gnvm;
              gnvec[2][i] = gnvec[2][i]/gnvm;

         // Compute normals to center columns
              cnvec[0][i] = sin(lon_gd[i][ic]*M_PI/180)*cos(lat_gd[i][ic]*M_PI/180)*gnvec[2][i] - sin(lat_gd[i][ic]*M_PI/180)*gnvec[1][i];
              cnvec[1][i] = sin(lat_gd[i][ic]*M_PI/180)*gnvec[0][i] - cos(lon_gd[i][ic]*M_PI/180)*cos(lat_gd[i][ic]*M_PI/180)*gnvec[2][i];
              cnvec[2][i] = cos(lon_gd[i][ic]*M_PI/180)*cos(lat_gd[i][ic]*M_PI/180)*gnvec[1][i] - sin(lon_gd[i][ic]*M_PI/180)*cos(lat_gd[i][ic]*M_PI/180)*gnvec[0][i];

            //; Compute grid row nadir resolution
              dc[0][i] = cos(lon_gd[i][ic+1]*M_PI/180)*cos(lat_gd[i][ic+1]*M_PI/180)- cos(lon_gd[i][ic]*M_PI/180)*cos(lat_gd[i][ic]*M_PI/180);
              dc[1][i] = sin(lon_gd[i][ic+1]*M_PI/180)*cos(lat_gd[i][ic+1]*M_PI/180)- sin(lon_gd[i][ic]*M_PI/180)*cos(lat_gd[i][ic]*M_PI/180);
              dc[2][i] = sin(lat_gd[i][ic+1]*M_PI/180)- sin(lat_gd[i][ic]*M_PI/180);
           }

      //dot product
//        for(pix=0;pix<num_pixels;pix++){
   for(i=0;i<num_gridlines;i++){
             if(lat>90 || lat<-90 || lon<-180 || lon>180){cout<<"latitude longitude pixel out of the boundaries.."<<endl; exit(1);}

            pvec[0]=cos(lon*M_PI/180)*cos(lat*M_PI/180);
            pvec[1]=sin(lon*M_PI/180)*cos(lat*M_PI/180);
            pvec[2]=sin(lat*M_PI/180);

            dotprod=pvec[0]*gnvec[0][i]+pvec[1]*gnvec[1][i]+pvec[2]*gnvec[2][i];


            if(i==0){
                 pdotgn_first=dotprod;//first gridline
                 }
            if(i==num_gridlines-1){
                  pdotgn_last=dotprod;//last line
                 }

           if(pdotgn_first<=db && pdotgn_last>-db && gdrow<0 && gdcol<0){
 //           cout<<"this pix is inside the L1C grid.."<<endl;
              }

           //check row for pixel j
           if(dotprod<=db && dotprod>-db && gdrow<0){
                 gdrow=i;
               }

           if(gdrow>=0){
              dotprod2=pvec[0]*cnvec[0][gdrow]+pvec[1]*cnvec[1][gdrow]+pvec[2]*cnvec[2][gdrow];
              dcm=sqrt(dc[0][gdrow]*dc[0][gdrow]+dc[1][gdrow]*dc[1][gdrow]+dc[2][gdrow]*dc[2][gdrow]);
              gdcol = ic + dotprod2/dcm + .5;
              if(gdrow>=0 && gdcol>=0) {
                 flag_out=0;
                 inpix++;}
              else{
                 flag_out=1;
                 outpix++; }
                    }
               }
           *row=gdrow;
           *col=gdcol;

           gdrow=-1;
           gdcol=-1;


   //bmcm only checked one pixel!!I
   if(dc!=nullptr)
        delete [] (dc);
   if(gnvec!=nullptr)
        delete [] (gnvec);
   if(cnvec!=nullptr)
        delete [] (cnvec);

   if(flag_out==0)
       return 0;
   else
       return 1;
}


