
//**********************************************************
//l1cgen firstly defined as mainL1C.cpp
//  Originally created by Martin Montes on 3/27/2021 (first version)
//  last version 1/6/2023
//**********************************************************
#include <allocate2d.h>
#include <allocate3d.h>
#include "allocate4d.h"
#include <iostream>
#include "l1c_input.h"
#include "l1c_filehandle.h"
#include "l1c_str.h"
#include "l1c.h"
#include "hawkeye_methods.h"
#include <chrono>
#include <sys/stat.h>
#include "l2_str.h"

#include <filehandle.h>
#include <l1.h>
#include <genutils.h>

#include <l1c_latlongrid.h>
#include <netcdf>


//#include <l1b_misr.h>

using namespace std;
using namespace l1c; 
using namespace std::chrono;
using namespace netCDF;
using namespace netCDF::exceptions;

#define PROGRAM_NAME "l1cgen"
#define VERSION "5.16 3/8/2023"

//function to calculate cross product of two vectors, 3 components each vector
   void cross_product_double(double vector_a[], double vector_b[], double temp[]) {
        temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
        temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
        temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
        }

    double cross_product_norm_double(double vector_a[], double vector_b[]) {
        double temp[3],nvec;
        temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
        temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
        temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

        nvec=sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
        return nvec;
        }

int main(int argc, char** argv)
{
   int status;  
   std::string str;
   const char* ptstr;
   double etime;
   string ifile_str;
   char* ifile_char;
  //------------------------------------------------------------------
   auto start = high_resolution_clock::now();
   int swtd = 1;
   size_t num_scans, num_pixels,num_views;
   double* time_mgv = nullptr,*time_mgv_off=nullptr;
   float ** lat_gd = nullptr, ** lon_gd = nullptr,**alt=nullptr;
   int32_t num_gridlines;
   int nfiles;
   double* tcross = nullptr, * ect_d = nullptr;
   float* mgvswt = nullptr;
   int16_t* tiod = nullptr, * odir = nullptr, * fileid = nullptr, * swtdid = nullptr, * nfiles_swt = nullptr;
   float** pix_size_u, ** pix_size_v;
   float** Ltfracsum = nullptr, ** areabinsum = nullptr, ** nobs_perbin = nullptr;
   vector<pair<float, int>> vp;
//    char* diroutname="out/";
//    int check;
   size_t sline;
   static float Re = 6378.137;

   L1C_input l1cinput;
   l1c_filehandle l1cfile;
   l1c_str l1cstr;
   l2_str l2str;
   L1C* ptl1c = new L1C(); 
   filehandle l1file;
   l1str l1rec;
//-------------------------------------------------------------------
   if (argc == 1)
    {
        l1cinput.l1c_usage(PROGRAM_NAME, VERSION);
        return 1;
    }

   for (int i = 0; i < argc; i++)
    {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0))
        {
            l1cinput.l1c_usage(PROGRAM_NAME, VERSION);
        }
    }

   cout<<"PROGRAM_NAME:"<<PROGRAM_NAME<<"\n"<<"VERSION:"<<VERSION<<endl;
   l1cfile.version=VERSION;
    //grab CLI info
   l1cinput.l1c_inputmain(argc, argv, &l1cinput, &l1cfile, PROGRAM_NAME, VERSION);

   ifile_str = l1cinput.files[0];
   ifile_char = &ifile_str[0];
   file_format format = getFormat(ifile_char);
   l1cfile.format=format.type;
    //create output file dir
/*    check=mkdir(diroutname,0777);//add error here in case the dir exists
        if (check == -1)
        // cerr << "Error :  " << strerror(errno) << endl;
        cout<<"the out directory was already created-------------------------------"<<endl;
        else
         cout << "Output Directory created..."<<diroutname<<endl;
*/
   if (format.type != FT_INVALID)
    {


            if (ptl1c->load_l1c_filehandle4(&l1cfile, &l1cinput) != 0)
            {
                printf("-E- %s: Error loading %sl1c filehandle.\n", argv[0], l1cfile.l1b_name.c_str());
                exit(1);
            }
                 
            //******************************************************************************************************************************
                 //LEGACY OPTION 1 (L1C grid generation from l1b files)
              //**********************************************************************************************************************
            if (l1cfile.l1c_pflag == 1)
            { 
                size_t sfiles=0;
               double tcross1=-1.; 
               size_t nscans_swt;
               double *time_swt=nullptr;
               float **pos_swt=nullptr,**vel_swt=nullptr;    
               double rl2,pos_norm,clatg,clatg2,fe=1/298.257;
               double v1[3],v2[3],vecout[3],orbnorm[3],nvec,vi,toff;
               double oangle,G[3],glat,glon,gnorm,rem=6371,omf2,omf2p,pxy,temp;
               string senstr;
               double ***time_off=nullptr;
               double rl;

               

                cout << "producing L1C grid..........from l1b files............................................" << endl;

          

                nfiles = l1cfile.ifiles.size();
                if(l1cinput.sensor==34){
                 senstr="SPEXone";
                 l1cfile.nbinx=25;
                 l1cfile.n_views=5;
                       }
                else if (l1cinput.sensor==30){
                 senstr="OCI";
                 l1cfile.nbinx=519;
                 l1cfile.n_views=2;
                    }
                else if (l1cinput.sensor==35){
                 senstr="HARP2";
                 l1cfile.nbinx=457;
                 l1cfile.n_views=90;
                       }
                else if (l1cinput.sensor==28){
                    senstr="MISR";
                    l1cfile.nbinx=81;
                       }
                else{cout<<"sensor by default is OCI option 2....."<<endl;
                 senstr="OCI";
                 l1cfile.nbinx=519;
                 l1cfile.n_views=2;
                    }

                num_views=l1cfile.n_views;

                 try{
                  int result;
                  result=l1cfile.selgran[0];
                  if(result<0){
                cout << "number of files to be processed..." << l1cfile.ifiles.size() << endl;
                        while(l1cfile.selgran[sfiles]<0 && sfiles<l1cfile.ifiles.size()){
                             l1cfile.selgran[sfiles]=sfiles+1;
                             cout<<"selected granule #..........................."<<l1cfile.selgran[sfiles]<<endl;
                             sfiles++;                 
                        }}
                  else{
                  throw(result);
                  }}

                catch(int e){
                   if(e>=0)//process all granules
                    {
                     cout<<"ERROR selgran>=0 ..all granules for the swath are used to produce the L1C grid, not specific ones.."<<"l1cfile.selgran[0]<=0...."<<e<<endl;
                     cout<<"DO NOT USE selgran option with l1c_pflag==1"<<endl;
                     exit(1);
                     }
                  }
               

//MISR sensor
//---------------------------------------------------------------------
               if(l1cfile.format == FT_MISR){
                      filehandle_init(&l1file);
                      if(status = l1cstr.openl1b_misr_l1c(&l1cstr, &l1cfile, &l1file)>1){
                         printf("-E- %d: status=2  error opening..L1B file for sensor...............\n", format.type);
                         exit(1);
                      }
                   }

//---------------------------------------------------------------------

                int16_t c = 0;
                tcross = (double*)calloc(nfiles,sizeof(double));

                                                                   
                for (int i = 0; i < nfiles; i++)
                {
                    str = l1cfile.ifiles[i];
                    ptstr = str.c_str();
         
               
                    if (status = ptl1c->ect_sf2(ptstr, &l1cinput,&l1cfile) > 1)
                    {
                        printf("-E- %s: status=2  error computing equatorial crossing time%f for file%s.................\n", argv[0], tcross[i], ptstr);
                        exit(1);
                    }
                    tcross[i] = l1cfile.eqt;
                    if (l1cfile.eqt > 0)
                    {
                        c++;
                        tcross1=l1cfile.eqt;
                    }
                        
                } 
               

                cout << "number of eq crossings per day..." << c << endl;
                l1cfile.nswath = c;

                fileid = (int16_t*)calloc(nfiles, sizeof(int16_t));
                swtdid = (int16_t*)calloc(nfiles,sizeof(int16_t));
                ect_d = (double*)calloc(nfiles,sizeof(double));
                nfiles_swt = (int16_t*)calloc(l1cfile.nswath,sizeof(int16_t));
                tiod = (int16_t*)calloc(l1cfile.nswath,sizeof(int16_t));
                odir = (int16_t*)calloc(l1cfile.nswath,sizeof(int16_t));
                mgvswt = (float*)calloc(l1cfile.nswath,sizeof(float)); //c
                
                num_scans = l1cfile.nscan;


                time_swt = (double*)calloc(nfiles*num_scans,sizeof(double));
                pos_swt = allocate2d_float(nfiles*num_scans, 3);
                vel_swt = allocate2d_float(nfiles*num_scans, 3);
           
                ptl1c->mov_SOCEA(&l1cfile, &l1cinput, tcross, fileid, swtdid, nfiles_swt, ect_d, tiod, odir, mgvswt,&nscans_swt,time_swt,pos_swt,vel_swt);

                if (tcross != nullptr)
                delete [](tcross);
                if (tiod != nullptr)
                delete [](tiod);
                if (odir != nullptr)
                delete [](odir);

                num_pixels = l1cfile.npix;                
                ptl1c->calc_biny(swtd,&l1cinput,&l1cfile,nscans_swt,time_swt,tcross1,mgvswt[0]);
              
                num_gridlines = l1cfile.num_gridlines;

                time_mgv = (double*)calloc(num_gridlines,sizeof(double)); //estimate of number of grids, assuming 4000 gridlines MAX INITIALLY!!...
                time_mgv_off = (double*)calloc(num_gridlines,sizeof(double)); 


                cout<<"tcross1...."<<tcross1<<"mean ground velocity #1..."<<mgvswt[0]<<"num_gridlines..."<<num_gridlines<<endl;
  

                ptl1c->swtime_swt3(1,&l1cinput,&l1cfile,nscans_swt,time_swt,tcross1,mgvswt[0],time_mgv);
                                                                 
                float mgv=mgvswt[0];

                if (mgvswt != nullptr)
                delete [](mgvswt);
                if (ect_d != nullptr)
                delete [](ect_d);
                         
                cout<<"number of across bins L1C grid...#.."<<l1cfile.nbinx<<"num of views.."<<num_views<<endl;

                num_gridlines=l1cfile.num_gridlines;
                lat_gd = allocate2d_float(num_gridlines, l1cfile.nbinx);
                lon_gd = allocate2d_float(num_gridlines, l1cfile.nbinx);
                alt=allocate2d_float(num_gridlines, l1cfile.nbinx);
                time_off = allocate3d_double(num_gridlines, l1cfile.nbinx,num_views);


                orb_array* posgrid = new orb_array[num_gridlines]();
                orb_array* velgrid = new orb_array[num_gridlines]();

                orb_array* pos_tot = new orb_array[num_scans*nfiles]();
                orb_array* vel_tot = new orb_array[num_scans*nfiles]();


                for(size_t i=0;i<num_scans*nfiles;i++){
                     for(size_t j=0;j<3;j++){
                       pos_tot[i][j]=pos_swt[i][j];
                       vel_tot[i][j]=vel_swt[i][j];
                     }                                 
                }


                if (pos_swt != nullptr)
                   delete [](pos_swt);
                if (vel_swt != nullptr)
                  delete [](vel_swt);


                orb_interp(nscans_swt,num_gridlines,time_swt,pos_tot,vel_tot,time_mgv,posgrid,velgrid);

                for(int i=0;i<num_gridlines;i++){
                          pos_norm=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1]+posgrid[i][2]*posgrid[i][2]);
                          clatg2=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1])/pos_norm;
                          rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                          v2[0]=velgrid[i][0]*rl2/pos_norm;
                          v2[1]=velgrid[i][1]*rl2/pos_norm;
                          v2[2]=velgrid[i][2]*rl2/pos_norm;

                          vi=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])*1000;//ground velocity in m/s
                          toff=vi/mgv;

                          time_mgv_off[i]=time_mgv[i]+toff;
                            }

                 cout<<"# gridlines after time interpolation correction...."<<num_gridlines<<endl;


                if (posgrid != nullptr)
                     delete[](posgrid);

                if (velgrid != nullptr)
                     delete[](velgrid); 


                orb_array* posgrid2 = new orb_array[num_gridlines]();
                orb_array* velgrid2 = new orb_array[num_gridlines]();
            
               
                orb_interp(nscans_swt, num_gridlines,time_swt,pos_tot,vel_tot,time_mgv_off, posgrid2,velgrid2);


                if (time_mgv != nullptr)
                    delete [](time_mgv);
                if (time_mgv_off != nullptr)
                    delete [](time_mgv_off);
                if (time_swt != nullptr)
                 delete [](time_swt);


                if (pos_tot != nullptr)
                 delete [](pos_tot);
                if (vel_tot != nullptr)
                 delete [](vel_tot);

                 //angle subsat track----
                omf2=(1-fe)*(1-fe);

                  for(int i=0;i<num_gridlines;i++){
                        pos_norm=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1]+posgrid2[i][2]*posgrid2[i][2]);
                        clatg2=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1])/pos_norm;
                        rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));
               //ground pos
                        v1[0]=(posgrid2[i][0])*rl2/pos_norm;
                        v1[1]=(posgrid2[i][1])*rl2/pos_norm;
                        v1[2]=(posgrid2[i][2])*rl2/pos_norm;
               //ground vel
                        v2[0]=(velgrid2[i][0])*rl2/pos_norm;
                        v2[1]=(velgrid2[i][1])*rl2/pos_norm;
                        v2[2]=(velgrid2[i][2])*rl2/pos_norm;

                        cross_product_double(v1,v2,vecout);
                        nvec=cross_product_norm_double(v1,v2);

                        orbnorm[0]=vecout[0]/nvec;
                        orbnorm[1]=vecout[1]/nvec;
                        orbnorm[2]=vecout[2]/nvec;

                   for(int j=0;j<l1cfile.nbinx;j++){

                      pos_norm=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);

                      oangle=asin((j-(l1cfile.nbinx-1)/2)*5.2/pos_norm);

                  //Geocentric vector
                      G[0]=v1[0]*cos(oangle)-orbnorm[0]*pos_norm*sin(oangle);
                      G[1]=v1[1]*cos(oangle)-orbnorm[1]*pos_norm*sin(oangle);
                      G[2]=v1[2]*cos(oangle)-orbnorm[2]*pos_norm*sin(oangle);

                      glon=atan2(G[1],G[0])*180/M_PI;

                      gnorm = sqrt(G[0]*G[0]+G[1]*G[1]+G[2]*G[2]);
                      omf2p = (omf2*rem + gnorm - rem)/gnorm;
                      pxy = G[0]*G[0]+G[1]*G[1];
                      temp = sqrt(G[2]*G[2] + omf2p*omf2p*pxy);
                      glat=asin(G[2]/temp)*180/M_PI;

                      lat_gd[i][j]=glat;
                      lon_gd[i][j]=glon;

                      for(size_t v=0;v<num_views;v++){
                         time_off[i][j][v]=time_mgv_off[i]-tcross1;
                      }

                      //altitude
                      clatg = cos(atan(omf2*tan(glat*M_PI/180.)));
                      rl = rem*(1.-fe)/sqrt(1.-(2.-fe)*fe*clatg*clatg);
                      alt[i][j] = (gnorm - rl)*1000;//in meters
                   }

                  }

                ptl1c->create_SOCEA2(1,&l1cinput,&l1cfile,lat_gd,lon_gd,alt,time_mgv_off);//THIS IS SWATH PROCESSING 


                if (posgrid2 != nullptr)
                     delete[](posgrid2);

                if (velgrid2 != nullptr)
                     delete[](velgrid2);
                
                if (lat_gd != nullptr)
                    delete [](lat_gd);
                if (lon_gd != nullptr)
                    delete [](lon_gd);
                if (fileid != nullptr)
                delete [](fileid);
                if (swtdid != nullptr)
                delete [](swtdid);
                if (nfiles_swt != nullptr)
                delete [](nfiles_swt);
                if (ptl1c != nullptr)
                delete ptl1c;
                if (time_off != nullptr)
                delete [] (time_off);
                if (alt != nullptr)
                delete [] (alt);
 
                auto stop1 = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop1 - start);
                etime = duration.count();
                etime /= (1000000 * 60);
                cout<<"done processing options #1..(creating L1C grid) in "<<etime<<"minutes..."<<endl;          
            }
 
          
 //OPTION 2 (binning SBS)
  //**********************************************************************************************************************************************
            else if (l1cfile.l1c_pflag == 2)
            {                 
            float ****binLt=nullptr,****binLt2;
            int ****bincount=nullptr,****bincount2;
            int sfiles=0,fi=0;

                cout<<"------------------------------SBS/SOCEA discrete binning methods--------------------"<<endl;
                cout<<"processing legacy sensor..."<<"l1cgen option #...."<<l1cfile.l1c_pflag<<endl;
                cout << "loading L1C grid and processing L1C granule............................................................................." << endl;
  
                nfiles=l1cfile.ifiles.size();
                
                try{
                  int result;
                  result=l1cfile.selgran[0];
                  if(result>0){
                    while(l1cfile.selgran[sfiles]>0){
                       cout<<"selected granule #..........................."<<l1cfile.selgran[sfiles]<<endl;
                       sfiles++;
                     }}
                else{ 
                  throw(result);
                  }}
                catch(int e){
                   if(e<0)
                    {
                     cout<<"all L1B granules will be processed...."<<nfiles<<endl;
                     cout<<"the number of files must be <=10!!..."<<endl;
                     sfiles=nfiles;
                     for(int j=0;j<sfiles;j++) l1cfile.selgran[j]=j+1; //the number of files must be <=10
                     }
                    if(e==0)
                    {
                     cout<<"wrong granule #..."<<e<<"must be >0..."<<endl;
                     exit(1);
                     }
                  }

                cout<<"number of granules to be processed...."<<sfiles<<endl;
                
                l1cfile.nswath = 1;

                swtdid = (int16_t*)calloc(nfiles , sizeof(int16_t));
                odir = (int16_t*)calloc(l1cfile.nswath , sizeof(int16_t));
                fileid = (int16_t*)calloc(nfiles, sizeof(int16_t));
                nfiles_swt = (int16_t*)calloc(l1cfile.nswath , sizeof(int16_t));

                nfiles_swt[0]=nfiles;
                odir[0]=1;
                int swtd=1;
                swtdid[0]=1; 

                l1cfile.nbinx=519;

                ptl1c->openL1Cgrid3(swtd,&l1cstr,&l1cfile,&l1cinput,swtdid,fileid,nfiles_swt); 

                cout<<"OK after sorting latitudes..."<<endl; 

                for(int j=0;j<sfiles;j++){
                 fi=l1cfile.selgran[j];
                 if(l1cfile.ifiles.size()>1){
                       l1cfile.l1b_name=l1cfile.ifiles[fi-1];
                   }               
                 cout<<"binning granule....................................................................................."<<fi<<"filename......................................................#"<<l1cfile.l1b_name<<endl;        

                 if(l1cfile.format == FT_SPEXONE){
                     if(status = l1cstr.openl1b_spex_l1c(&l1cstr, &l1cfile)>1){
                        printf("-E- %d: status=2  error opening..L1B file for sensor...............\n", format.type);
                        exit(1);
                    }
                  }
                 else if(l1cfile.format == FT_HARP2){
                      if(status = l1cstr.openl1b_harp2_l1c(&l1cstr, &l1cfile)>1){
                         printf("-E- %d: status=2  error opening..L1B file for sensor...............\n", format.type);
                         exit(1);
                     }
                 }    
                 else{//OCI
                     if(status = l1cstr.openl1b_ocis_l1c(&l1cstr, &l1cfile)>1){
                         printf("-E- %d: status=2  error opening..L1B file for sensor...............\n", format.type);
                         exit(1);
                     }
                   }
                      num_scans = l1cfile.nscan;
                      cout<<"#gridlines.."<<l1cfile.num_gridlines<<"nbinx..."<<l1cfile.nbinx<<endl;



                   if(l1cfile.format == FT_SPEXONE || l1cfile.format ==FT_HARP2){
               //intensity bands        
                      binLt = allocate4d_float(l1cfile.n_views,l1cfile.nband_view,l1cfile.num_gridlines,l1cfile.nbinx);
                      bincount = allocate4d_int(l1cfile.n_views,l1cfile.nband_view,l1cfile.num_gridlines,l1cfile.nbinx);
               //polarized bands
                      binLt2 = allocate4d_float(l1cfile.n_views,l1cfile.npol_band_view,l1cfile.num_gridlines,l1cfile.nbinx);
                      bincount2 = allocate4d_int(l1cfile.n_views,l1cfile.npol_band_view,l1cfile.num_gridlines,l1cfile.nbinx);
                      }
                   else{
                      binLt = allocate4d_float(l1cfile.n_views,l1cfile.nbands,l1cfile.num_gridlines,l1cfile.nbinx);
                      bincount = allocate4d_int(l1cfile.n_views,l1cfile.nbands,l1cfile.num_gridlines,l1cfile.nbinx);
                      }



                      for(size_t v=0;v<l1cfile.n_views;v++){
                           for(size_t b=0;b<l1cfile.nbands;b++){
                                for(int gd=0;gd<l1cfile.num_gridlines;gd++){
                                     for(int xbin=0;xbin<l1cfile.nbinx;xbin++){
                                      binLt[v][b][gd][xbin]=0;
                                      bincount[v][b][gd][xbin]=0;
                               }}}}


                     for (sline = 0; sline < num_scans; sline++)
                      {
                      if(l1cfile.format == FT_SPEXONE){
                          if( status = l1cstr.readl1b_spex_l1c(&l1cstr, &l1cfile,sline)>1){
                            printf("-E- %lu: status=2  error reading line#%d..from L1B file for sensor...............\n", sline,format.type);
                            exit(1);
                         }
                       }

                      else if(l1cfile.format == FT_HARP2){
                           if( status = l1cstr.readl1b_harp2_l1c(&l1cstr, &l1cfile,sline)>1){
                            printf("-E- %lu: status=2  error reading line#%d..from L1B file for sensor...............\n", sline,format.type);
                            exit(1);
                         }
                      }
                      else{  
                          if( status = l1cstr.readl1b_ocis_l1c(&l1cstr, &l1cfile,sline)>1){                         
                            printf("-E- %lu: status=2  error reading line#%d..from L1B file for sensor...............\n", sline,format.type);
                            exit(1);
                         }
                       }

                       cout<<"------------------------------------------------------------------------------------------------------------------"<<endl;
                       cout<<"sline............................................................................................#"<<sline+1<<endl;       
                       cout<<"------------------------------------------------------------------------------------------------------------------"<<endl;

                       ptl1c->binL1C_sbs_line3(swtd,ptl1c,&l1cstr,&l1cfile,&l1cinput,swtdid,fileid,nfiles_swt,binLt,bincount,binLt2,bincount2,sline,fi); 

   
                                                                         
                          if (sline  % 100 == 0)  cout <<"line #..."<<sline+1<<endl; 

                       }


                   if(l1cfile.format == FT_SPEXONE){
                       if(status=l1cstr.closel1b_spex_l1c(&l1cstr,&l1cfile)>1)
                       {
                      printf("-E- %d: status=2  error closing L1B file for sensor...............\n",format.type);
                        exit(1);
                       }                    
                   }
                   else if(l1cfile.format == FT_HARP2){
                      if(status=l1cstr.closel1b_harp2_l1c(&l1cstr,&l1cfile)>1)
                       {
                      printf("-E- %d: status=2  error closing L1B file for sensor...............\n",format.type);
                        exit(1);
                       }
                   }  
                   else{
                    if(status=l1cstr.closel1b_ocis_l1c(&l1cstr,&l1cfile)>1)
                       {
                      printf("-E- %d: status=2  error closing L1B file for sensor...............\n",format.type);
                        exit(1);
                       }
                     }

                 if(binLt!=nullptr)
                        delete [] (binLt);
                 if(bincount!=nullptr)
                      delete [] (bincount);
                     
                 if(binLt2!=nullptr)
                        delete [] (binLt2);
                 if(bincount2!=nullptr)
                      delete [] (bincount2);
                     }


                 if (l1cfile.lat_gd != nullptr)
                   delete[](l1cfile.lat_gd);
                 if (l1cfile.lon_gd != nullptr)
                    delete[](l1cfile.lon_gd);
                 if (l1cfile.alt_gd != nullptr)
                    delete[](l1cfile.alt_gd);
                 if (l1cfile.lat_asort != nullptr)
                     delete[](l1cfile.lat_asort);
                 if (l1cfile.index_xy != nullptr)
                     delete[](l1cfile.index_xy);
              
                 if (ptl1c != nullptr)
                   delete ptl1c;
                 if (fileid != nullptr)
                    delete [](fileid);
                 if (swtdid != nullptr)
                   delete [](swtdid);
                 if (odir != nullptr)
                   delete [](odir);
                 if (nfiles_swt != nullptr)
                    delete [](nfiles_swt);

                 auto stop1 = high_resolution_clock::now();
                 auto duration = duration_cast<microseconds>(stop1 - start);
                 etime = duration.count();
                 etime /= (1000000 * 60);
                 cout<<"done processing options #2..(Binning granule using saddleback search processing (LINE-BY-LINE)"<<etime<<"minutes..."<<endl;                   
                    }

 //LEGACY OPTION 3 (binning area-weighting)
        else if (l1cfile.l1c_pflag == 3)
            {
                int sfiles=0,fi=0,****bincount=nullptr,****bincount2=nullptr;
                float ****binLt=nullptr,****binLt2=nullptr;
                cout<<"------------------------------area-weighting method--------------------"<<endl;
                cout<<"*********************************legacy sensor ************************************"<<endl;
                cout << "number of files to be processed for creating the L1C grid..." << l1cfile.ifiles.size() << endl;

                nfiles = l1cfile.ifiles.size();
                 try{
                  int result;
                  result=l1cfile.selgran[0];
                  if(result>0){
                      while(l1cfile.selgran[sfiles]>0){
                       cout<<"selected granule #..........................."<<l1cfile.selgran[sfiles]<<endl;
                       sfiles++;
                     }}
                  else{
                  throw(result);
                  }}
                catch(int e){
                   if(e<0)
                    {
                     cout<<"all L1B granules will be processed...."<<nfiles<<endl;                     
                     sfiles=nfiles;
                     cout<<"the number of files must be <=10!!..."<<endl;

                     for(int j=0;j<sfiles;j++) l1cfile.selgran[j]=j+1; //the number of files must be <=10

                     cout<<"ALL GRANULES are processed TOO IF selgran input argument was not setup in CLI ..please use selgran=-1 for binning ALL granules or selgran>0 for speccific granules"<<endl;
                     cout<<"otherwise all granules will be processed by default!!"<<endl;
                     }
                    if(e==0)
                    {
                     cout<<"wrong granule #..."<<e<<"must be >0..."<<endl;
                     exit(1);
                    }
                  }

                cout<<"number of granules to be processed...."<<sfiles<<endl;

                l1cfile.nswath = 1;

                swtdid = (int16_t*)calloc(nfiles , sizeof(int16_t));
                odir = (int16_t*)calloc(l1cfile.nswath , sizeof(int16_t));
                fileid = (int16_t*)calloc(nfiles, sizeof(int16_t));
                nfiles_swt = (int16_t*)calloc(l1cfile.nswath , sizeof(int16_t));

                nfiles_swt[0]=nfiles;
                odir[0]=1;
                int swtd=1;
                swtdid[0]=1;

                l1cfile.nbinx=519;

                ptl1c->openL1Cgrid3(swtd,&l1cstr,&l1cfile,&l1cinput,swtdid,fileid,nfiles_swt);
                cout<<"OK after sorting latitudes..."<<endl; 

                for(int g=0;g<sfiles;g++){

                   if(l1cfile.format == FT_SPEXONE || l1cfile.format ==FT_HARP2){
               //intensity bands        
                      binLt = allocate4d_float(l1cfile.n_views,l1cfile.nband_view,l1cfile.num_gridlines,l1cfile.nbinx);
                      bincount = allocate4d_int(l1cfile.n_views,l1cfile.nband_view,l1cfile.num_gridlines,l1cfile.nbinx);
               //polarized bands
                      binLt2 = allocate4d_float(l1cfile.n_views,l1cfile.npol_band_view,l1cfile.num_gridlines,l1cfile.nbinx);
                      bincount2 = allocate4d_int(l1cfile.n_views,l1cfile.npol_band_view,l1cfile.num_gridlines,l1cfile.nbinx);
                      }
                   else{
                      binLt = allocate4d_float(l1cfile.n_views,l1cfile.nbands,l1cfile.num_gridlines,l1cfile.nbinx);
                      bincount = allocate4d_int(l1cfile.n_views,l1cfile.nbands,l1cfile.num_gridlines,l1cfile.nbinx);
                      }

                      for(size_t v=0;v<l1cfile.n_views;v++){
                           for(size_t b=0;b<l1cfile.nbands;b++){
                                for(int gd=0;gd<l1cfile.num_gridlines;gd++){
                                     for(int xbin=0;xbin<l1cfile.nbinx;xbin++){
                                      binLt[v][b][gd][xbin]=0;
                                      bincount[v][b][gd][xbin]=0;
                               }}}}
 

                   fi=l1cfile.selgran[g];
                   l1cfile.l1b_name=l1cfile.ifiles[fi-1];
                   str=l1cfile.ifiles[fi-1];
                   ptstr = str.c_str();

                   cout<<"binning granule....................................................................................."<<fi<<"filename......................................................#"<<l1cfile.l1b_name<<endl;

                   if(status = l1cstr.openl1b_ocis_l1c(&l1cstr, &l1cfile)>1){
                      printf("-E- %d: status=2  error opening..L1B file for sensor...............\n", format.type);
                        exit(1);
                    }

                   cout<<"sfiles.."<<sfiles<<"#gridlines.."<<l1cfile.num_gridlines<<"nbinx.."<<l1cfile.nbinx<<"num_scans.."<<l1cfile.nscan<<endl;
              
                   pix_size_u = allocate2d_float(l1cfile.nscan-1,l1cfile.npix-1);
                   pix_size_v = allocate2d_float(l1cfile.nscan-1,l1cfile.npix-1);

                   Ltfracsum = allocate2d_float(l1cfile.num_gridlines,l1cfile.nbinx);
                   areabinsum = allocate2d_float(l1cfile.num_gridlines,l1cfile.nbinx);
                   nobs_perbin = allocate2d_float(l1cfile.num_gridlines,l1cfile.nbinx);

                   for (sline = 0; sline <l1cfile.nscan-1; sline++)
                      {
                     if( status = l1cstr.readl1b_ocis_l1c(&l1cstr, &l1cfile,sline)>1)
                         {
                          printf("-E- %lu: status=2  error reading line#%d..from L1B file for sensor...............\n", sline,format.type);
                          exit(1);
                         }
                                 
               //        ptl1c->xy_pixsize_sf4(ptstr,&l1cstr,&l1cfile,&l1cinput,pix_size_u,pix_size_v,Ltfracsum,areabinsum,nobs_perbin,binLt,bincount,sline);
                         ptl1c->binL1C_wgranule_aw3(swtd,&l1cfile,&l1cinput,&l1cstr,Ltfracsum,areabinsum,nobs_perbin,sline);                     
                       
                  if(sline%10 == 0) cout<<"Processing scanline.#............................................................."<<sline+1<<endl;
                   }

                    if(status=l1cstr.closel1b_ocis_l1c(&l1cstr,&l1cfile)>1)
                       {
                      printf("-E- %d: status=2  error closing L1B file for sensor...............\n",format.type);
                        exit(1);
                       }

                    cout<<"producing nc files with area-weighted Lts................filename.=.."<<ptstr<<endl;
               //     ptl1c->savebinL1C_v3(swtd,&l1cinput,&l1cfile,&l1cstr,l1cfile.lat_gd,l1cfile.lon_gd,Ltfracsum,areabinsum,nobs_perbin);

                   
                   if(binLt!=nullptr)
                        delete [] (binLt);
                   if(bincount!=nullptr)
                      delete [] (bincount);
                   if(binLt2!=nullptr)
                        delete [] (binLt2);
                   if(bincount2!=nullptr)
                      delete [] (bincount2);                     
                   if (pix_size_u != nullptr)
                       delete [](pix_size_u);
                   if (pix_size_v != nullptr)
                     delete [](pix_size_v);
                   if (Ltfracsum != nullptr)
                       delete [](Ltfracsum);
                   if (areabinsum != nullptr)
                      delete [](areabinsum);
                   if (nobs_perbin != nullptr)
                     delete [](nobs_perbin);

                       } 
  
              if (l1cfile.lat_gd != nullptr)
                   delete[](l1cfile.lat_gd);
              if (l1cfile.lon_gd != nullptr)
                    delete[](l1cfile.lon_gd);
              if (l1cfile.alt_gd != nullptr)
                    delete[](l1cfile.alt_gd);
              if (l1cfile.lat_asort != nullptr)
                     delete[](l1cfile.lat_asort);
              if (l1cfile.index_xy != nullptr)
                     delete[](l1cfile.index_xy);

              if (ptl1c != nullptr)
                delete ptl1c;
              if (fileid != nullptr)
                delete [](fileid);
              if (swtdid != nullptr)
                delete [](swtdid);
              if (odir != nullptr)
                delete [](odir);
              if (nfiles_swt != nullptr)
                delete [](nfiles_swt);


            auto stop1 = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop1 - start);
            etime = duration.count();
            etime /= (1000000 * 60);
            cout<<"done processing options #3..Binning granule based on area-weighting in "<<etime<<"minutes..."<<endl;
            }        
     //----- READING L2 FILES ---READING WITH SDS --------------------------------------------------------------
     //**********************************************************************************************************-
            else if (l1cfile.l1c_pflag == 4){              
                float ****binprod=nullptr;
                int ****bincount=nullptr;
                int sfiles=0,fi=0;


                size_t *ybins=nullptr;
                size_t *xbins=nullptr;
                lat_gd = allocate2d_float(4000, 519);
                lon_gd = allocate2d_float(4000, 519);
            
                open_l1c(l1cinput.l1c_grid,ybins,xbins,lat_gd,lon_gd);
                cout<<"lat_gd[0][0]"<<lat_gd[0][0]<<"lon_gd[0][0]"<<lon_gd[0][0]<<"alt[0][0]"<<alt[0][0]<<endl;
                exit(1);

                cout<<"processing L2 file to L1C ----------"<<endl;
                cout<<"SDS binning products........"<<l1cinput.l2prod<<endl;

                cout<<"------------------------------SBS method--------------------"<<endl;
                cout<<"processing legacy sensor..."<<"l1cgen option #...."<<l1cfile.l1c_pflag<<endl;
                cout << "loading L1C grid and processing L2 file to L1C granule............................................................................." << endl;

                nfiles=l1cfile.ifiles.size();//listfile as txt in ifile

                try{
                  int result;
                  result=l1cfile.selgran[0];//example selgran=[1], selgran=[1,2,3]
                  if(result>0){
                      while(l1cfile.selgran[sfiles]>0){
                       cout<<"selected granule #..........................."<<l1cfile.selgran[sfiles]<<endl;
                       sfiles++;
                     }}
                  else{
                  throw(result);
                  }}
                catch(int e){
                   if(e<0)
                    {
                     cout<<"all L2 granules will be processed...."<<nfiles<<endl;
                     sfiles=nfiles;
                     }
                    if(e==0)
                    {
                     cout<<"wrong granule #..."<<e<<"must be >0..."<<endl;
                     exit(1);
                     }
                  }

                cout<<"number of granules to be processed...."<<sfiles<<endl;

                l1cfile.nswath = 1;

                swtdid = (int16_t*)calloc(nfiles , sizeof(int16_t));
                fileid = (int16_t*)calloc(nfiles, sizeof(int16_t));
                nfiles_swt = (int16_t*)calloc(l1cfile.nswath , sizeof(int16_t));

                nfiles_swt[0]=nfiles;
                swtdid[0]=1;

                ptl1c->openL1Cgrid(&l1cstr,&l1cfile,&l1cinput);                                                        
       //L2 granule loop
                binprod = allocate4d_float(l1cfile.n_views,l2str.nl2prod,l1cfile.num_gridlines,l1cfile.nbinx);
                bincount = allocate4d_int(l1cfile.n_views,l2str.nl2prod,l1cfile.num_gridlines,l1cfile.nbinx);
            
                for(int j=0;j<sfiles;j++){
                 l1cfile.l1b_name=l1cfile.ifiles[0];//first l2 file, fi=0
                 cout<<"binning filename......................................................................................#"<<l1cfile.l1b_name<<endl;

                 if(status = l2str.openl2_ocis_l1c(&l1cinput,&l2str,&l1cfile,fileid)>1){
                      printf("-E- %d: status=2  error opening..L2B file for sensor...............\n", l1cfile.format);
                        exit(1);
                    }
                  
                     for(size_t v=0;v<l1cfile.n_views;v++){
                      for(size_t p=0;p<l2str.nl2prod;p++){
                                for(int gd=0;gd<l1cfile.num_gridlines;gd++){
                                     for(int xbin=0;xbin<l1cfile.nbinx;xbin++){
                                      binprod[v][p][gd][xbin]=0.0;
                                      bincount[v][p][gd][xbin]=0.0;
                               }}}}

                    cout<<"number of views.."<<l1cfile.n_views<<"number of products.."<<l2str.nl2prod<<"#gridlines.."<<l1cfile.num_gridlines<<"number of xbins.."<<l1cfile.nbinx<<endl;

                    for (sline = 0; sline <l1cfile.nscan; sline++)
                      {
                      if( status = l2str.readl2_ocis_l1c(&l2str, &l1cfile, fileid, sline)>1)
                      {
                       printf("-E- %lu: status=2  error reading line#%d..from L2 file for sensor...............\n", sline,l1cfile.format);
                        exit(1);
                       }
                        ptl1c->binL1C_sbs_line_l2(ptl1c,&l2str,&l1cfile,&l1cinput,binprod,bincount,sline,fi);

                     if (sline  % 100 == 0)  cout <<"line #..."<<sline+1<<endl;
                       }

                    if(status=l2str.closel2_ocis_l1c(&l2str,&l1cfile)>1)
                       {
                        printf("-E- %d: status=2  error closing L2 file for sensor...............\n",l1cfile.format);
                        exit(1);
                       }
                   }

                if(binprod!=nullptr)
                        delete [] (binprod);
                if(bincount!=nullptr)
                      delete [] (bincount);
            
                if (l1cfile.lat_gd != nullptr)
                  delete[](l1cfile.lat_gd);
                if (l1cfile.lon_gd != nullptr)
                  delete[](l1cfile.lon_gd);   
                if (l1cfile.alt_gd != nullptr)
                  delete[](l1cfile.alt_gd);                
                if (l1cfile.lat_asort != nullptr)
                  delete[](l1cfile.lat_asort);
                if (l1cfile.index_xy != nullptr)
                  delete[](l1cfile.index_xy);

                if (ptl1c != nullptr)
                  delete ptl1c;
                if (fileid != nullptr)
                  delete [](fileid);
                if (swtdid != nullptr)
                  delete [](swtdid);
                if (nfiles_swt != nullptr)
                  delete [](nfiles_swt);

                auto stop1 = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop1 - start);
                etime = duration.count();
                etime /= (1000000 * 60);
                cout<<"done processing options #4..(Binning L2 file  using saddleback search processing "<<etime<<"minutes..."<<endl;
             } 
     //----- READING HKT FILES -----------------------------------------------------------------
     //**********************************************************************************************************-
            else if (l1cfile.l1c_pflag == 5){
                cout<<"Reading telemetry from L1A files..--> SOCEA -- L1C grid.."<<endl;
       
              if(status = ptl1c->open_l1atol1c3(&l1cinput,&l1cfile)>1){
                      printf("-E- %d: status=2  error opening..L1A file for sensor...............\n", format.type);
                        exit(1);
                    }

                 if (ptl1c != nullptr)
                    delete ptl1c;

                 auto stop1 = high_resolution_clock::now();
                 auto duration = duration_cast<microseconds>(stop1 - start);
                 etime = duration.count();
                 etime /= (1000000 * 60);
                 cout<<"done processing options #5..(creating SOCEA L1C common grid from telemetry "<<etime<<"minutes..."<<endl;
                    }

 
//----------------------------------------------------------------------------------------------------------
//FIXED-BEARING SWATH PROJECTION-----TESTED WITH OCI------------------------------------------------------- 
//---------------------------------------------------------------------------------------------------------           
            else if (l1cfile.l1c_pflag == 6)
              { 
                size_t sfiles=0;
                float *lati2=nullptr, *loni2=nullptr,*az_gd=nullptr;             
                double tcross1=-1;
                float loncross1=-1;

                float **lat_tot=nullptr, **lon_tot=nullptr;
                int16_t nfiles_swath=1;

                cout << "producing L1C grid...swath fixed bearing..................................................." << endl;

                nfiles = l1cfile.ifiles.size();
                num_scans = l1cfile.nscan;

                 try{
                  int result;
                  result=l1cfile.selgran[0];
                  if(result<0){
                cout << "number of files to be processed..." << l1cfile.ifiles.size() << endl;
                        while(l1cfile.selgran[sfiles]<0 && sfiles<l1cfile.ifiles.size()){
                             l1cfile.selgran[sfiles]=sfiles+1;
                             cout<<"selected granule #..........................."<<l1cfile.selgran[sfiles]<<endl;
                             sfiles++;                 
                        }}
                  else{
                  throw(result);
                  }}

                catch(int e){
                   if(e>=0)
                    {
                     cout<<"ERROR selgran>=0 ..all granules for the swath are used to produce the L1C grid, not specific ones.."<<"l1cfile.selgran[0]<=0...."<<e<<endl;
                     cout<<"DO NOT USE selgran option with l1c_pflag==1"<<endl;
                     exit(1);
                     }
                  }
                
                int16_t c = 0;

                l1cfile.nbinx = 519;

                for (int i = 0; i < nfiles; i++)
                {
                    str = l1cfile.ifiles[i];
                    ptstr = str.c_str();

                    if (status = ptl1c->ect_sf2(ptstr, &l1cinput,&l1cfile) > 1)
                    {
                        cout<<"error computing equatorial crossing time...."<<tcross1<<" for file%s.................\n"<<ptstr<<endl;;
                        exit(1);
                    }

                    if (l1cfile.eqt > 0)
                    {
                        tcross1=l1cfile.eqt;
                        loncross1=l1cfile.orbit_node_lon;
                        c++;
                    }
                } 
                 cout << "number of eq crossings per day..." << c << "tcross1..."<<tcross1<<"loncross1..."<<loncross1<<endl;
                l1cfile.nswath = c;

                l1cfile.eqt=tcross1;
                l1cfile.orbit_node_lon=loncross1;


                fileid = (int16_t*)calloc(nfiles, sizeof(int16_t));
                swtdid = (int16_t*)calloc(nfiles,sizeof(int16_t));
                ect_d = (double*)calloc(nfiles,sizeof(double));
                nfiles_swt = (int16_t*)calloc(l1cfile.nswath,sizeof(int16_t));
                tiod = (int16_t*)calloc(l1cfile.nswath,sizeof(int16_t));
                odir = (int16_t*)calloc(l1cfile.nswath,sizeof(int16_t));

                if (ect_d != nullptr)
                 delete [](ect_d);
                if (tiod != nullptr)
                 delete [](tiod);
                if (odir != nullptr)
                 delete [](odir);

                num_scans = l1cfile.nscan;
                num_pixels = l1cfile.npix;
                nfiles_swath=nfiles;

                lat_tot = allocate2d_float(num_scans * nfiles_swath, num_pixels);
                lon_tot = allocate2d_float(num_scans * nfiles_swath, num_pixels);     

                nfiles_swt[0]=nfiles;
                for (int i = 0; i < nfiles; i++){
                     swtdid[i]=1;
                     fileid[i]=i+1; 
                }


                cout<<"swath_scans.."<<num_scans * nfiles_swath<<"num pixels.."<<num_pixels<<endl;
                ptl1c->swath_latlon(1,&l1cfile,&l1cinput,swtdid,fileid,nfiles_swt,lat_tot,lon_tot);

                ptl1c->azmean_swt3(swtd,&l1cinput,&l1cfile,lat_tot,lon_tot);

                l1cfile.num_gridlines=4000;//assuming5.2. bin resolution for L1C grid

                lati2 = (float*)calloc(l1cfile.num_gridlines,sizeof(float));
                loni2 = (float*)calloc(l1cfile.num_gridlines,sizeof(float));
                az_gd = (float*)calloc(l1cfile.num_gridlines,sizeof(float));

                ptl1c->calc_biny_dist(1,&l1cfile,&l1cinput,lati2,loni2);

                num_gridlines = l1cfile.num_gridlines;

                lat_gd = allocate2d_float(num_gridlines, l1cfile.nbinx);
                lon_gd = allocate2d_float(num_gridlines, l1cfile.nbinx);
                      
                ptl1c->across_gridlines_l1c2(swtd, &l1cfile, &l1cinput, swtdid, fileid, nfiles_swt, lati2, loni2, lat_gd, lon_gd, az_gd);
                
                if (lati2 != nullptr)
                    delete [](lati2);
                if (loni2 != nullptr)
                    delete [](loni2);
                if (lat_tot != nullptr)
                    delete [](lat_tot);
                if (lon_tot != nullptr)
                    delete [](lon_tot);

                if (lat_gd != nullptr)
                    delete [](lat_gd);
                if (lon_gd != nullptr)
                    delete [](lon_gd);
                if (az_gd != nullptr)
                    delete [](az_gd);
                if (fileid != nullptr)
                delete [](fileid);
                if (swtdid != nullptr)
                delete [](swtdid);
                if (nfiles_swt != nullptr)
                delete [](nfiles_swt);
                if (ptl1c != nullptr)
                delete ptl1c;
 
                auto stop1 = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop1 - start);
                etime = duration.count();
                etime /= (1000000 * 60);
                cout<<"done processing options #7..(creating L1C grid using SWATH FIXED BEARING PROJ) in "<<etime<<"minutes..."<<endl;          
            }
     //**********************************************************************************************************************
    else if (l1cfile.l1c_pflag == 7)
            {                                                                         
               string senstr;
               string l1c_str=l1cinput.l1c_grid;
                 
                cout << "producing L1C grid from l2/l3/anc granules and CTH ......................................................" << endl;


                nfiles = l1cfile.ifiles.size();
                if(l1cinput.sensor==34){
                 senstr="SPEXone";
                 l1cfile.nbinx=25;
                       }
                else if (l1cinput.sensor==30){
                 senstr="OCI";
                 l1cfile.nbinx=519;
                    }
                else if (l1cinput.sensor==35){
                 senstr="HARP2";
                 l1cfile.nbinx=457;
                       }
                else if (l1cinput.sensor==28){
                    senstr="MISR";
                    l1cfile.nbinx=81;
                       }
                else{cout<<"sensor by default is OCI option 2....."<<endl;
                 senstr="OCI";
                 l1cfile.nbinx=519;
                    }

            
                 cout << "number of files to be processed..." << 1 << endl;
                 cout << "number of L1C granules for screening scantime..." << l1cinput.files_l1c.size() << endl;
                            
                 NcFile* nc_l12;
                 string l_12_str=l1cfile.ifiles[0];
                  const char *l_12=l_12_str.c_str();

                 try {
                  nc_l12 = new NcFile(l_12, NcFile::replace);
                    }
                 catch (NcException& e) {
                    e.what();
                    cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: "
                     + l_12_str << endl;
                    exit(1);
                    }

                 ptl1c->l1b_cloud_correct(&l1cinput,&l1cfile,nc_l12);

                 nc_l12->close();
              
                 
                if (ptl1c != nullptr)
                delete ptl1c;
 
                auto stop1 = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop1 - start);
                etime = duration.count();
                etime /= (1000000 * 60);
                cout<<"done processing options #8..(creating L1C grid USING CTH) in "<<etime<<"minutes..."<<endl;          
            }
    else if (l1cfile.l1c_pflag == 9)
            {                 
                string senstr;
                string l1c_str=l1cinput.l1c_grid;
                cout << "producing L1C grid from L1C granules  and at CTH ......................................................" << endl;


                nfiles = l1cfile.ifiles.size();
                if(l1cinput.sensor==34){
                 senstr="SPEXone";
                 l1cfile.nbinx=25;
                       }
                else if (l1cinput.sensor==30){
                 senstr="OCI";
                 l1cfile.nbinx=519;
                    }
                else if (l1cinput.sensor==35){
                 senstr="HARP2";
                 l1cfile.nbinx=457;
                       }
                else if (l1cinput.sensor==28){
                    senstr="MISR";
                    l1cfile.nbinx=81;
                       }
                else{cout<<"sensor by default is OCI option 2....."<<endl;
                 senstr="OCI";
                 l1cfile.nbinx=519;
                    }

            
                 cout << "number of L1B files to be processed..." << 1 << endl;
                 cout << "number of L1C granules for screening scantime..." << l1cinput.files_l1c.size() << endl;
                                 
                 ptl1c->l1c_cloud_correct(&l1cinput,&l1cfile);
              
              
                if (ptl1c != nullptr)
                delete ptl1c;
 
                auto stop1 = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop1 - start);
                etime = duration.count();
                etime /= (1000000 * 60);
                cout<<"done processing options #9..(creating L1C grid USING L1C and CTH) in "<<etime<<"minutes..."<<endl;          
            }

//end option 9      
      //**********************************************************************************************************************
    else if (l1cfile.l1c_pflag == 8)
        {
          cout << "searching row/col of L1C grid given a L1B granule and L1C grid inputs......................................................" << endl;
         //we are using here Dons L1B libraries for open, read and close, one line at the time!!
//---------------------------------------------------------
   //TEST LIBRARY -------------------------            
//---------------------------------------------------------             
              cout<<"init filehandle and add options.."<<endl; 
              filehandle_init(&l1file); 
              l1_input_init();
              clo_optionList_t* optionList = clo_createList();
              l1_add_options(optionList);
      //add extra options from CLI and coming in msl12_defaults.par
              clo_addOption(optionList, "suite", CLO_TYPE_STRING, "OC", "suite of OC products");
              clo_addOption(optionList, "aermodfile", CLO_TYPE_STRING, "Unspecified", "path dir for aerosol files anc data");       
              clo_addOption(optionList, "aer_wave_short", CLO_TYPE_INT, "765", "default shortest wavelength used for aerosol correction, epsilon");
              clo_addOption(optionList, "aer_wave_long", CLO_TYPE_INT, "865", "default longest wavelength used for aerosol correction, epsilon");
              clo_addOption(optionList, "mumm_alpha", CLO_TYPE_FLOAT, "1.72", "mumm alpha for AC Rudick turbid waters");
              clo_addOption(optionList, "mumm_gamma", CLO_TYPE_FLOAT, "1.0", "mumm gamma for AC Rudick turbid waters");
              clo_addOption(optionList, "mumm_epsilon", CLO_TYPE_FLOAT, "1.0", "mumm epsilon for AC Rudick turbid waters");
              clo_addOption(optionList, "chloc2_wave", CLO_TYPE_INT, "[-1,-1]", "sensor wavelengths for OC2 chlorophyll\n        algorithm");
              clo_addOption(optionList, "chloc2_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]", "coefficients for OC2\n        chlorophyll algorithm");
              clo_addOption(optionList, "chloc3_wave", CLO_TYPE_INT, "[-1,-1]", "sensor wavelengths for OC3 chlorophyll\n        algorithm");
              clo_addOption(optionList, "chloc3_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]", "coefficients for OC3\n        chlorophyll algorithm");
              clo_addOption(optionList, "chloc4_wave", CLO_TYPE_INT, "[-1,-1]", "sensor wavelengths for OC4 chlorophyll\n        algorithm");
              clo_addOption(optionList, "chloc4_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]", "coefficients for OC4\n        chlorophyll algorithm");
              clo_addOption(optionList, "kd2_wave", CLO_TYPE_INT, "[-1,-1]", "sensor wavelengths for polynomial Kd(490)\n        algorithm");
              clo_addOption(optionList, "kd2_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0,0.0]", "sensor wavelengths\n        for polynomial Kd(490) algorithm");
              clo_addOption(optionList, "coccolith", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]", "coefficients for coccolith  algorithm");
              clo_addOption(optionList, "qaa_wave", CLO_TYPE_INT, NULL,"sensor wavelengths for QAA algorithm");
              clo_addOption(optionList, "giop_wave", CLO_TYPE_FLOAT, "-1", "optimization comma-separated list, default is all visible bands (400-700nm)");


              l1_read_default_files(optionList, &l1file, ifile_char);
              l1_load_options(optionList, &l1file);
              clo_deleteList(optionList);
       
              cout<<"open l1.."<<endl;
              openl1(&l1file); 
              string l1c_str=l1cinput.l1c_grid;
              const char *l1c_grid=l1c_str.c_str();

              bin_str binl1c;
              binl1c.version=l1cfile.version;
              binl1c.history=l1cinput.history;

              check_swath_time(&l1file,l1c_grid,&binl1c);          
        
              //allocate mem for l1rec
                 if (alloc_l1(&l1file, &l1rec) == 0) {
                     cout<<"Unable to allocate L1 record....."<<endl;
                       exit(1);
               }
               
              int32_t spix=MAX(l1_input->spixl-1,0);
              int32_t epix=MAX(l1_input->epixl-1,l1file.npix-1);
              int32_t dpix=MAX(l1_input->dpixl,1);
              int32_t sscan=MAX(l1_input->sline-1,0);
              int32_t escan=MAX(l1_input->eline-1,l1file.nscan-1);
              int32_t dscan=MAX(l1_input->dline-1,1);


              l1file.spix=spix;
              l1file.epix=epix;

              int32_t npix=(epix-spix)/dpix+1;
              int32_t nscan=(escan-sscan)/dscan+1;

              l1file.npix=npix;

              cout<<"npix.."<<l1file.npix<<"npix.."<<npix<<"nscan.."<<l1file.nscan<<endl; 
              cout<<"spix.."<<spix<<"epix.."<<epix<<"sscan.."<<sscan<<"escan.."<<escan<<"dpix.."<<dpix<<"dscan.."<<dscan<<endl;

              short **gdindex=nullptr;

              cout<<"number of 'working' bands..."<<l1file.nbands<<endl;
              string senstr;
              char *ifile_char=strdup(l1cinput.l1c_grid);
              file_format format = getFormat(ifile_char);
              if(format.type==FT_HKT || format.type==FT_L1C) format.type=FT_OCIS; //if HKT then FT_OCIS----

              if(format.type==FT_SPEXONE){
                         senstr = "SPEXONE";
                     }
              else if(format.type==FT_HARP2){
                         senstr = "HARP2";                                                                                              
                     }
              else if(format.type==FT_OCIS || format.type==FT_OCIL1B){
                         senstr = "OCI";
                     }
              else//OCIS
                     {
                         senstr = "OCI";
                     }
           
              NcFile* nc_output;
              string file_str=l1c_str.substr(5,15);
              string l1c_full_str="PACE_"+senstr+"."+file_str+"Z.L1C.5.2km.nc";
              const char *l1c_full=l1c_full_str.c_str();

              cout<<"creating FULL L1C file---------------"<<l1c_full_str<<endl;
            
               try {
                  nc_output = new NcFile(l1c_full, NcFile::replace);
                    }
               catch (NcException& e) {
                 e.what();
                 cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: "
                 + l1c_full_str << endl;
                 exit(1);
              }

              meta_l1c_full(&l1file,&binl1c,l1cinput.l1c_grid,nc_output);
             
              gdindex = allocate2d_short(npix,2);              

                  
//--------SCANTIME-------------------------------------------------------------
             l1cfile.l1b_name=l1cfile.ifiles[0];
             double scantime=-1;

          //time_offset loop   
             if(status = l1cstr.openl1b_ocis_l1c(&l1cstr, &l1cfile)>1){
                      printf("-E- %d: status=2  error opening..L1B file for sensor...............\n", format.type);
                        exit(1);
                    }
             for (int sline = 0; sline <nscan; sline++)
                   {
                     if( status = l1cstr.readl1b_ocis_l1c(&l1cstr, &l1cfile,sline)>1)
                         {
                          printf("-E- %d: status=2  error reading line#%d..from L1B file for sensor...............\n", sline,format.type);
                          exit(1);
                         }
             
                      scantime=l1cstr.timepix[sline];
                      cout<<"scantime..."<<scantime<<endl;
            
                      readl1(&l1file,sline,&l1rec);
                      search_l1c(&l1file,&l1rec,&binl1c,gdindex);                
                      bintime_l1c(&l1file,&l1rec,&binl1c,gdindex,scantime,nc_output);
                    /*  
                     if(gdindex[635][0]<3710 || gdindex[635][1]<3710) {
                         cout<<"ERROR row/col in the south"<<endl;
                         cout<<"scanline.."<<sline+1<<"scantime..."<<scantime<<"lat.."<<l1rec.lat[635]<<endl;
                         exit(1);
                         }
                     */    


                  if(sline%10 == 0) cout<<"Processing scanline.#............................................................."<<sline+1<<endl;
 
                   }

//----------------------------------------------------------------------------

//--------NO TIME BINNED VARS-----------------------------------------------
              for(int sline=0;sline<nscan;sline++){                                    
                 readl1(&l1file,sline,&l1rec);
                 search_l1c(&l1file,&l1rec,&binl1c,gdindex);                                 
                 bin_l1c(&l1file,&l1rec,&binl1c,gdindex,nc_output,l1cstr.timepix[sline]);                             
               if((sline % 10)==0) cout<<"scanline#.."<<sline+1<<"inpix.."<<binl1c.inpix<<"outpix.."<<binl1c.outpix<<endl;            
              }//lines
      
              for(int sline=0;sline<nscan;sline++){                    
                 readl1(&l1file,sline,&l1rec);
                 search_l1c(&l1file,&l1rec,&binl1c,gdindex);
                 rmse_l1c_alt(&l1file,&binl1c,&l1rec,gdindex,l1cstr.timepix[sline]);
                 if((sline % 10)==0) cout<<"scanline#.."<<sline+1<<endl;   
              }      

              if(status=l1cstr.closel1b_ocis_l1c(&l1cstr,&l1cfile)>1)

                       {
                      printf("-E- %d: status=2  error closing L1B file for sensor...............\n",format.type);
                        exit(1);
                       }

//---------------------------------------------------------------------------     
              delete [] (gdindex);

                                                                               
              meta_l1c_bin(&l1file,&binl1c,nc_output); 
              meta_l1c_altvar(&binl1c,nc_output);

              nc_output->close();
              closel1(&l1file);
              binl1c.close_bin(&binl1c);        
                          
              delete(l1_input);


              auto stop1 = high_resolution_clock::now();
              auto duration = duration_cast<microseconds>(stop1 - start);
              etime = duration.count();
              etime /= (1000000 * 60);
              cout<<"#binned pixels.."<<binl1c.inpix<<"#pixels grid excluded.."<<binl1c.outpix<<"total # pixels.."<<nscan*npix<<endl;
              cout<<"done processing options #8..(using L1C functions from libl1 and l1 readers) in "<<etime<<"minutes..."<<endl;

//--------------------------------------------------------
//---------------------------------------------------------------
         }//end option 8
      } //end ft diff of invalid 

 
          
    else
    {
        printf("-E- %s: Error opening %s for reading unknown sensor...................\n", argv[0], l1cfile.l1b_name.c_str());
        return (1);
    }

    return (0);
}
