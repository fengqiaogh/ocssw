
//  Created by Martin Montes on 10/26/20.
//latest update 12/6/2022
//****************************************************8
#include <filehandle.h>
#include <l1.h>
#include "l1c.h"
#include "l1c_filehandle.h"
#include "l1c_input.h"
#include "l1c_str.h"
#include "l2_str.h"
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <list>
#include <stack>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <netcdf.h>
#include <nc4utils.h>
#include "hawkeye_methods.h"
#include "allocate4d.h"
#include <allocate3d.h>
#include <allocate2d.h>
#include <algorithm> 
#include <bits/stdc++.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <timeutils.h>
#include <genutils.h>
#include <stdint.h>
#include "l2prod.h"
#include "Strassen.h"
#include <libnav.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/Poisson_reconstruction_function.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>
#include <CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <boost/iterator/transform_iterator.hpp>

#include <cgal_interp.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/draw_triangulation_2.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/draw_triangulation_3.h>


#include <CGAL/Triangulation_data_structure_3.h>
#include <cassert>


#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/boost/graph/IO/OBJ.h>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>

#include<netcdf>
#include<l1c_latlongrid.h>


//TYPES
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3   Point;
typedef CGAL::Point_set_3<Point> Point_set;


typedef CGAL::Triangulation_2<K>                            Triangulation;
typedef Triangulation::Point                                Point2;

typedef CGAL::Delaunay_triangulation_3<K>                   DT3;
typedef CGAL::Creator_uniform_3<double,K::Point_3>          Creator;


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Sphere_3 Sphere_3;
typedef CGAL::Point_set_3<Point_3, Vector_3> Point_set;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Vertex_iterator        Vertex_iterator;

typedef std::pair<Point, Vector> Point_with_normal;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;

typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

typedef CGAL::Triangulation_3<K>      Triangulation3;
typedef Triangulation3::Finite_vertices_iterator Finite_vertices_iterator;
typedef Triangulation3::Finite_facets_iterator Finite_facets_iterator;

typedef CGAL::Triangulation_vertex_base_with_info_2<std::string, Gt> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Delaunay3;


typedef CGAL::Triangulation_data_structure_3<>      Tds2;
typedef Tds2::size_type                              size_type;
typedef Tds2::Cell_handle                            Cell_handle;
typedef Tds2::Vertex_handle                          Vertex_handle;


namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef K::FT                                                           FT;
typedef CGAL::Surface_mesh<Point_3>                                     Mesh;
typedef typename boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor;
typedef typename boost::graph_traits<Mesh>::face_descriptor             face_descriptor;
namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Mesh, FT>                                    Face_location;
typedef K::Ray_3                                                        Ray_3;


using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Projection_traits = CGAL::Projection_traits_xy_3<Kernel>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;
// Triangulated Irregular Network
using TIN = CGAL::Delaunay_triangulation_2<Projection_traits>;




static size_t INPIX,TOTPIX,OUTPIX;
static size_t num_scans,num_pixels,num_frames,nviews,num_records,nframes;
static size_t num_blue_bands, num_red_bands, num_SWIR_bands;
static int ncid_L1B,ncid_L2,ncid_L1A,ncid_L1C;


//Navigation data
static int navGrp, scGrp, geoGrp;
static int ovId, opId, otId, latId, lonId, altId,solzId,senzId,senaId,timeId;
static float** lat, ** lon;
static short** solz,**senz=nullptr,**sena=nullptr;;
static float degrad = M_PI / 180, Re = 6378.137;//Re earth radius in km at equator
//static float *orbp;
//static float *attang;
//static double *timeArr;


//L1C vars
static int16_t nswath, nswtfiles, ndswaths;


using namespace std;
using namespace boost::assign;
using namespace std::chrono;
using namespace GeographicLib;


using namespace netCDF;
using namespace netCDF::exceptions;



namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> Point_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::box<Point_t> Box_t;


namespace l1c {


    L1C::L1C() {
    }


    L1C::~L1C() {
    }


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


int32_t L1C::l1c_cloud_correct(L1C_input *l1cinput,l1c_filehandle *l1cfile){
   int dimid,status;
   string str,str2;
   double *ptime=nullptr,*ptime_l1c=nullptr;
   const char *ptstr_l1c;
   size_t num_ybin,num_xbin;
   int gfound=-1;
   float **latpix=nullptr,**lonpix=nullptr,**l1clat=nullptr,**l1clon=nullptr;
   short **senapix=nullptr,**senzpix=nullptr;
   float **cth_l1c=nullptr;
   float **lat_new=nullptr,**lon_new=nullptr;
   string senstr;
   float **posr=nullptr,**velr=nullptr;
   size_t n_orb_rec;
   float Rpole=6356;;
   float Xcorr,Ycorr,Zcorr;
   size_t att_len;
   float Rsurf,radius_ratio=Re/Rpole;
   float latcorr; 

   NcFile* nc_l12;
   string l12_str=l1cinput->files[0];
   const char *l_12=l12_str.c_str();

   try {
        nc_l12 = new NcFile(l_12, NcFile::read);
                    }
        catch (NcException& e) {
                    e.what();
                    cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: "
                     + l12_str << endl;
                    exit(1);
         }

//open L1B/L2 granule ------


   nc_l12->close();

/*

   status = nc_open(ptstr, NC_NOWRITE, &ncid_L2);
        if (status != NC_NOERR) {
                fprintf(stderr, "-E- Error failed nc_open.\n");
                        exit(EXIT_FAILURE);
            }
   cout<<"Opening L1B/L2 granule ..."<<l1cinput->files[0]<<"....for parallax-cloud corrections......."<<endl;


 //global attributes
   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "time_coverage_start", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* time_str = (char*)malloc(att_len + 1); // + 1 for trailing null 
        time_str[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "time_coverage_start", time_str);
        string tswt_ini(time_str);

   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "time_coverage_end", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* time_str2 = (char*)malloc(att_len + 1); // + 1 for trailing null
        time_str2[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "time_coverage_end", time_str2);
        string tswt_end(time_str2);

   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "startDirection", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* sdir = (char*)malloc(att_len + 1); // + 1 for trailing null
        sdir[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "startDirection", sdir);
        string start_dir(sdir);

   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "endDirection", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* edir = (char*)malloc(att_len + 1); // + 1 for trailing null
        edir[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "endDirection", edir);
        string end_dir(edir);
  
   //opening groups
   status = nc_inq_grp_ncid(ncid_L2, "scan_line_attributes", &scGrp);
        check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(scGrp, "time", &timeId);
        check_err(status, __LINE__, __FILE__);

   status = nc_inq_grp_ncid(ncid_L2, "geolocation_data", &geoGrp);
        check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(geoGrp, "latitude", &latId);
                check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(geoGrp, "longitude", &lonId);
                   check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(geoGrp, "sensor_zenith", &senzId);
                   check_err(status, __LINE__, __FILE__);                   
   status = nc_inq_varid(geoGrp, "sensor_azimuth", &senaId);
                   check_err(status, __LINE__, __FILE__);

   status = nc_inq_grp_ncid(ncid_L2, "navigation_data", &navGrp);
        check_err(status, __LINE__, __FILE__);
   ////open velocity vectors
   status = nc_inq_varid(navGrp, "orb_pos", &opId);
        check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(navGrp, "orb_vel", &ovId);
        check_err(status, __LINE__, __FILE__);
                   
   status = nc_inq_dimid(ncid_L2, "number_of_scans", &dimid);
        check_err(status, __LINE__, __FILE__);
   nc_inq_dimlen(ncid_L2, dimid, &num_scans);
   status = nc_inq_dimid(ncid_L2, "ccd_pixels", &dimid);
        check_err(status, __LINE__, __FILE__);
   nc_inq_dimlen(ncid_L2, dimid, &num_pixels);
   n_orb_rec=num_scans;

   posr = allocate2d_float(n_orb_rec, 3);
   velr = allocate2d_float(n_orb_rec, 3);
   ptime=(double*)calloc(num_scans,sizeof(double));     

   latpix = allocate2d_float(num_scans,num_pixels);
   lonpix = allocate2d_float(num_scans,num_pixels);
   senzpix = allocate2d_short(num_scans,num_pixels);
   senapix = allocate2d_short(num_scans,num_pixels);


   // get vars
   status = nc_get_var_float(navGrp, opId, &posr[0][0]);
        check_err(status, __LINE__, __FILE__);
   status = nc_get_var_float(navGrp, ovId, &velr[0][0]);
        check_err(status, __LINE__, __FILE__);
   status = nc_get_var_double(scGrp, timeId, &ptime[0]);
      check_err(status, __LINE__, __FILE__);
   status = nc_get_var_float(geoGrp, latId, &latpix[0][0]);
      check_err(status, __LINE__, __FILE__);
   status = nc_get_var_float(geoGrp, lonId, &lonpix[0][0]);
        check_err(status, __LINE__, __FILE__);
   status = nc_get_var_short(geoGrp, senzId, &senzpix[0][0]);
      check_err(status, __LINE__, __FILE__);
   status = nc_get_var_short(geoGrp, senaId, &senapix[0][0]);
        check_err(status, __LINE__, __FILE__);
          
//initial satellite position
  //Earth radius at pixel positions
//   Rsurf=Re/(sqrt(cos(glat*M_PI/180)*cos(glat*M_PI/180)+radius_ratio*radius_ratio * sin(glat*M_PI/180)*sin(glat*M_PI/180)));


// l1c_grid_str=1cinput->files_l1c[0];        

//--- OPEN L1C FILES ---
   for (unsigned int i = 0; i <l1cinput->files_l1c.size(); i++)  { 
    if(gfound<0){
          str=l1cinput->files_l1c[i];
       ptstr_l1c = str.c_str();
       status = nc_open(ptstr_l1c, NC_NOWRITE, &ncid_L1C);
        if (status != NC_NOERR) {
                fprintf(stderr, "-E- Error failed nc_open.\n");
                        exit(EXIT_FAILURE);
            }
       cout<<"Opening L1C .#.."<<i+1<<"...."<<ptstr_l1c<<"for parallax-cloud corrections......."<<endl;

  //dimensions
  //num_scans num_pixels
       status = nc_inq_grp_ncid(ncid_L1C, "geolocation_data", &geoGrp);
                check_err(status, __LINE__, __FILE__);

       status = nc_inq_dimid(ncid_L1C, "bins_along_track", &dimid);
        check_err(status, __LINE__, __FILE__);
       nc_inq_dimlen(ncid_L1C, dimid, &num_ybin);
       status = nc_inq_dimid(ncid_L1C, "bins_across_track", &dimid);
        check_err(status, __LINE__, __FILE__);
       nc_inq_dimlen(ncid_L1C, dimid, &num_xbin);


       ptime_l1c=(double*)calloc(num_ybin,sizeof(double));
       l1clat = allocate2d_float(num_ybin,num_xbin);
       l1clon = allocate2d_float(num_ybin,num_xbin);
       cth_l1c = allocate2d_float(num_ybin,num_xbin);
  //open groups and variables
       status = nc_inq_varid(geoGrp, "latitude", &latId);
                check_err(status, __LINE__, __FILE__);
       status = nc_inq_varid(geoGrp, "longitude", &lonId);
                   check_err(status, __LINE__, __FILE__);
       status = nc_inq_varid(geoGrp, "time_nadir", &timeId);
           
       status = nc_get_var_double(geoGrp, timeId, &ptime_l1c[0]);       
          check_err(status, __LINE__, __FILE__);
       status = nc_get_var_float(geoGrp, latId, &l1clat[0][0]);
          check_err(status, __LINE__, __FILE__);
       status = nc_get_var_float(geoGrp, lonId, &l1clon[0][0]);
        check_err(status, __LINE__, __FILE__);
  
       cout<<"num along bins L1C..."<<num_ybin<<"num_across bins L1C...."<<num_xbin<<"l1c time_ini.."<<ptime_l1c[0]<<"l1c time end.."<<ptime_l1c[num_ybin-1]<<endl;
       if ((status = nc_close(ncid_L1C)))
            check_err(status, __LINE__, __FILE__);


//mathing L1C granule with L1B list----
//-----------------------------------       
       for (unsigned int j=0;j<num_scans;j++){  
               if(ptime[j]>=ptime_l1c[0] && ptime[j]<=ptime_l1c[num_ybin-1]){
                  cout<<"L1C granule found!! #......"<<i+1<<"......l2 ptime[j].."<<ptime[j]<<"l1c time_ini.."<<ptime_l1c[0]<<"l1c time end.."<<ptime_l1c[num_ybin-1]<<endl;
                  gfound=1;
                  break;}                             
                }

        //interpolate CTH from L2 to l1c grid geograph coordinates
       if(gfound>0){
             lat_new = allocate2d_float(num_ybin,num_xbin);
             lon_new = allocate2d_float(num_ybin,num_xbin);

//          int nl2=num_scans*num_pixels;
//          int nl1c=num_ybin*num_xbin;

            cout<<"num_scans.."<<num_scans<<"num pixels.."<<num_pixels<<"nybins.."<<num_ybin<<"nxbin.."<<num_xbin<<endl;
       
*/

            
         //MAPPING L1C TO L2 ---
/*           int nnc_tag_l2 = cgal_nnc(nl2, &latpix[0][0], &lonpix[0][0],nl1c, &l1clat[0][0], &l1clon[0][0]);
           cout<<"nnc_tag_l2.."<<nnc_tag_l2<<endl;

               cgal_interp2(nl2,&latpix[0][0],&lonpix[0][0],&cth[0][0],nl1c, &cth_l1c[0][0], nnc_tag_l2);//this is using infor from cgal_nnc

             if(nnc_tag_l2>=0){
               lat_new = allocate2d_float(num_ybin,num_xbin);
               lon_new = allocate2d_float(num_ybin,num_xbin);  
                 // Define the dimensions,vars and attributes at the root level
               NY = num_ybin;
               NX = num_xbin;

               lat_out = allocate2d_float(NY, NX);
               lon_out = allocate2d_float(NY, NX);
               time_nad_out=(double*)calloc(NY,sizeof(double));//time of the day in seconds
                }
              else{cout<<"wrong nnc_tag_l2 <0.."<<nnc_tag_l2<<endl; exit(1);}
  */        

    /*
             for (unsigned int i = 0; i <num_ybin; i++)  {
                for (unsigned int j = 0; j <num_xbin; j++)  {
                  //   cout<<"ybin.."<<i+1<<"xbin..."<<j+1<<"cth_l1c.."<<cth_l1c[i][j]<<endl;                          
                 
              */      

   //compute new lat/lon based on CTH ---
   //spehrical coors at CTH
                  //  R=1000*cth_l1c[i][j]; //CTH in meters
               /*
                    R=1;
                    x_new=R*cos(l1clon[i][j]*M_PI/180)*cos(l1clat[i][j]*M_PI/180);
                    y_new=R*sin(l1clon[i][j]*M_PI/180)*cos(l1clat[i][j]*M_PI/180);
                    z_new=R*sin(l1clat[i][j]*M_PI/180);
    
  //geographic coordinates
                    rad_new=sqrt(x_new*x_new +y_new*y_new+z_new*z_new);
                    lat_new[i][j]=asin(z_new/rad_new)*180/M_PI;
                    lon_new[i][j]=atan2(y_new,x_new)*180/M_PI;
         //           cout<<"rad_new.."<<rad_new<<"lat_new..."<<lat_new<<"lon_new.."<<lon_new<<endl;         
                */

//               if(l1cinput->cloud_correct==2){
/*
                           Xp=Rsurf*cos(l1clat[i][j]*M_PI/180)*sin(l1clon[i][j]*M_PI/180);
                           Yp=Rsurf*sin(l1clat[i][j]*M_PI/180);
                           Zp=Rsurf*cos(l1clat[i][j]*M_PI/180)*cos(l1clon[i][j]*M_PI/180);

                           radius_ratio_cloud=((Re+cth[i][j])/(Rpole+cth[i][j]))*((Re+cth[i][j])/(Rpole+cth[i][j]));


                           Xdiff=Xs-Xp;
                           Ydiff=Ys-Yp;
                           Zdiff=Zs-Zp;

                           //local zenith angle

                           xfact=sqrt(Xdiff*Xdiff+Ydiff*Ydiff+Zdiff*Zdiff);
                           zen=(Xdiff*Xp+Ydiff*Yp+Zdiff*Zp)/(mean_radius*xfact);
                           zen=acos(zen);
                           zen=zen*180/M_PI;

       //                    equation to solve for the line of sight at height Z

                           e1 = Xdiff*Xdiff + radius_ratio_cloud*Ydiff*Ydiff + Zdiff*Zdiff;
                           e2 = 2.0 * (Xp*Xdiff + radius_ratio_cloud*Yp*Ydiff + Zp*Zdiff);
                           e3 = Xp*Xp + Zp*Zp + radius_ratio_cloud*Yp*Yp -(Re+cth[i][j])*(Re+cth[i][j]);

                           dcorr = (sqrt(e2*e2 - 4.0*e1*e3) - e2)/(2.0*e1);
                    //    corrected surface coordinates

                           Xcorr = Xp + dcorr*Xdiff;
                           Ycorr = Xp + dcorr*Ydiff;
                           Zcorr = Xp + dcorr*Zdiff;
*/

/*                    
 //     convert back to latitude and longitude
                             //Earth radius at pixel positions
                           Rsurf=Re/(sqrt(cos(l1clat[i][j]*M_PI/180)*cos(l1clat[i][j]*M_PI/180)+radius_ratio*radius_ratio * sin(l1clat[i][j]*M_PI/180)*sin(l1clat[i][j]*M_PI/180)));
                           Xcorr=(l1cinput->cloud_height+Rsurf)*cos(l1clat[i][j]*M_PI/180)*sin(l1clon[i][j]*M_PI/180);//lon speherical and geodetic should be equal
                           Ycorr=(l1cinput->cloud_height+Rsurf)*sin(l1clat[i][j]*M_PI/180);
                           Zcorr=(l1cinput->cloud_height+Rsurf)*cos(l1clat[i][j]*M_PI/180)*cos(l1clon[i][j]*M_PI/180);
                             
                           latcorr = atan(Ycorr/sqrt(Xcorr*Xcorr + Zcorr*Zcorr));
                           lat_new[i][j] = atan(tan(latcorr)/radius_ratio*radius_ratio) * 180.0/M_PI;
                           lon_new[i][j] = atan2(Xcorr,Zcorr) * 180.0/M_PI;
                        

                           lon_new[i][j] =atan2(Ycorr,Xcorr)*180/M_PI;
                        //   gnorm = sqrt(Xcorr*Xcorr+Ycorr*Ycorr+Zcorr*Zcorr);
//                           omf2p = (omf2*rem + gnorm - rem)/gnorm;
                       //    pxy = Xcorr*Xcorr+Ycorr*Ycorr;
                        //   temp = sqrt(Zcorr*Zcorr + omf2p*omf2p*pxy);
                        //   lat_new[i][j] =asin(Zcorr/temp)*180/M_PI;

                           cout<<"l1cinput->cloud_height.."<<l1cinput->cloud_height<<"latnew.."<<lat_new[i][j]<<"lon_new.."<<lon_new[i][j]<<"l1clat.."<<l1clat[i][j]<<"l1clon.."<<l1clon[i][j]<<endl;
                           
                        //   rad_new=sqrt(Xc*Xc +Yc*Yc+Zc*Zc); //OK FOR SPHERE NOT ELLIPS 
 
                            //altitude
                   //        clatg = cos(atan(omf2*tan(lat_new[i][j]*M_PI/180.)));
                       //    rl = Re*(1.-fe)/sqrt(1.-(2.-fe)*fe*clatg*clatg);
                   //        altitude = (gnorm - rl)*1000;//in meters

              //        }
       
                }
       
             }

          delete [] (lat_new);
          delete [] (lon_new);

           }//gfound>0

          delete [] (ptime_l1c);
          delete [] (l1clat);
          delete [] (l1clon);
          delete [] (cth_l1c);
        }//gfound<0
     }//l1c granules loop

*/     

//   if ((status = nc_close(ncid_L2)))
//        check_err(status, __LINE__, __FILE__);

return 0;
}


int32_t L1C::l1b_cloud_correct(L1C_input *l1cinput,l1c_filehandle *l1cfile, NcFile *nc_l12){
   int dimid,status,ncid_out,x_dimid, y_dimid, varid1, varid2, varid3;
   string str;
   int grp_coor;
   double *ptime=nullptr;
   const char* ptstr;
   float **latpix=nullptr,**lonpix=nullptr;
   short **senapix=nullptr,**senzpix=nullptr;
   float **cth=nullptr;
   float **lat_new=nullptr,**lon_new=nullptr;
   float fill_value=-999.;
   double fill_value2=-999.;
   int NVIEWS, NBANDS,v_dimid,b_dimid;
   int32_t NX, NY;
   const char* filename_lt;
   int NDIMS=2,NDIMS2=1;
   int dimids[NDIMS],*dimids2;
   string senstr;
   int shuffle, deflate, deflate_level;
   float **posr=nullptr,**velr=nullptr;
   size_t n_orb_rec;
   double oangle,G[3],glat,glon,gnorm,rem=6371,omf2,omf2p,pxy,temp;
   double rl2,pos_norm,clatg2,fe=1/298.257;
   double v1[3],v2[3],vecout[3],orbnorm[3],nvec;
   float *latsat=nullptr,*lonsat=nullptr;
   float dv,Hsat=676.5;//sensor height in km above ellipsoid
   size_t att_len;
  
   shuffle = NC_SHUFFLE;
   deflate = 1;
   deflate_level = 5;
      

   str=l1cinput->files[0];
   ptstr = str.c_str();

   status = nc_open(ptstr, NC_NOWRITE, &ncid_L2);
        if (status != NC_NOERR) {
                fprintf(stderr, "-E- Error failed nc_open.\n");
                        exit(EXIT_FAILURE);
            }
   cout<<"Opening L1B/L2 granule ..."<<ptstr<<"....for parallax-cloud corrections......."<<endl;


 //global attributes
   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "time_coverage_start", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* time_str = (char*)malloc(att_len + 1); // + 1 for trailing null 
        time_str[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "time_coverage_start", time_str);
        string tswt_ini(time_str);

   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "time_coverage_end", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* time_str2 = (char*)malloc(att_len + 1); // + 1 for trailing null
        time_str2[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "time_coverage_end", time_str2);
        string tswt_end(time_str2);

   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "startDirection", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* sdir = (char*)malloc(att_len + 1); // + 1 for trailing null
        sdir[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "startDirection", sdir);
        string start_dir(sdir);

   status = nc_inq_attlen(ncid_L2, NC_GLOBAL, "endDirection", &att_len);
            check_err(status, __LINE__, __FILE__);
   char* edir = (char*)malloc(att_len + 1); // + 1 for trailing null
        edir[att_len] = '\0';
   status = nc_get_att_text(ncid_L2, NC_GLOBAL, "endDirection", edir);
        string end_dir(edir);
  
   //opening groups
   status = nc_inq_grp_ncid(ncid_L2, "scan_line_attributes", &scGrp);
        check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(scGrp, "time", &timeId);
        check_err(status, __LINE__, __FILE__);

   status = nc_inq_grp_ncid(ncid_L2, "geolocation_data", &geoGrp);
        check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(geoGrp, "latitude", &latId);
                check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(geoGrp, "longitude", &lonId);
                   check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(geoGrp, "sensor_zenith", &senzId);
                   check_err(status, __LINE__, __FILE__);                   
   status = nc_inq_varid(geoGrp, "sensor_azimuth", &senaId);
                   check_err(status, __LINE__, __FILE__);

   status = nc_inq_grp_ncid(ncid_L2, "navigation_data", &navGrp);
        check_err(status, __LINE__, __FILE__);
   ////open velocity vectors
   status = nc_inq_varid(navGrp, "orb_pos", &opId);
        check_err(status, __LINE__, __FILE__);
   status = nc_inq_varid(navGrp, "orb_vel", &ovId);
        check_err(status, __LINE__, __FILE__);
                   
   status = nc_inq_dimid(ncid_L2, "number_of_scans", &dimid);
        check_err(status, __LINE__, __FILE__);
   nc_inq_dimlen(ncid_L2, dimid, &num_scans);
   status = nc_inq_dimid(ncid_L2, "ccd_pixels", &dimid);
        check_err(status, __LINE__, __FILE__);
   nc_inq_dimlen(ncid_L2, dimid, &num_pixels);
   n_orb_rec=num_scans;

   posr = allocate2d_float(n_orb_rec, 3);
   velr = allocate2d_float(n_orb_rec, 3);
   ptime=(double*)calloc(num_scans,sizeof(double));     

   latpix = allocate2d_float(num_scans,num_pixels);
   lonpix = allocate2d_float(num_scans,num_pixels);
   senzpix = allocate2d_short(num_scans,num_pixels);
   senapix = allocate2d_short(num_scans,num_pixels);

   cth = allocate2d_float(num_scans,num_pixels);
   latsat=(float*)calloc(num_scans,sizeof(float));
   lonsat=(float*)calloc(num_scans,sizeof(float));

   lat_new = allocate2d_float(num_scans,num_pixels);
   lon_new = allocate2d_float(num_scans,num_pixels);

   // get vars
   status = nc_get_var_float(navGrp, opId, &posr[0][0]);
        check_err(status, __LINE__, __FILE__);
   status = nc_get_var_float(navGrp, ovId, &velr[0][0]);
        check_err(status, __LINE__, __FILE__);
   status = nc_get_var_double(scGrp, timeId, &ptime[0]);
      check_err(status, __LINE__, __FILE__);
   status = nc_get_var_float(geoGrp, latId, &latpix[0][0]);
      check_err(status, __LINE__, __FILE__);
   status = nc_get_var_float(geoGrp, lonId, &lonpix[0][0]);
        check_err(status, __LINE__, __FILE__);
   status = nc_get_var_short(geoGrp, senzId, &senzpix[0][0]);
      check_err(status, __LINE__, __FILE__);
   status = nc_get_var_short(geoGrp, senaId, &senapix[0][0]);
        check_err(status, __LINE__, __FILE__);


//compute geodetic latitude/longitude for satellite sub-orbital locations ---
   //angle subsat track----
   omf2=(1-fe)*(1-fe);

   for(unsigned int i=0;i<num_scans;i++){
       pos_norm=sqrt(posr[i][0]*posr[i][0]+posr[i][1]*posr[i][1]+posr[i][2]*posr[i][2]);
       clatg2=sqrt(posr[i][0]*posr[i][0]+posr[i][1]*posr[i][1])/pos_norm;
                        rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));
               //ground pos
                        v1[0]=(posr[i][0])*rl2/pos_norm;
                        v1[1]=(posr[i][1])*rl2/pos_norm;
                        v1[2]=(posr[i][2])*rl2/pos_norm;
               //ground vel
                        v2[0]=(velr[i][0])*rl2/pos_norm;
                        v2[1]=(velr[i][1])*rl2/pos_norm;
                        v2[2]=(velr[i][2])*rl2/pos_norm;

                        cross_product_double(v1,v2,vecout);
                        nvec=cross_product_norm_double(v1,v2);//length of orb norm vect

                        orbnorm[0]=vecout[0]/nvec;
                        orbnorm[1]=vecout[1]/nvec;
                        orbnorm[2]=vecout[2]/nvec;
                         
                        oangle=0;

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
      
                        latsat[i]=glat;
                        lonsat[i]=glon;

                      //altitude
                    //    clatg = cos(atan(omf2*tan(glat*M_PI/180.)));
                    //    rl = Re*(1.-fe)/sqrt(1.-(2.-fe)*fe*clatg*clatg);
                    //    altitude = (gnorm - rl)*1000;//in meters
   }



   //fill up cth with random numbers ----
   for (unsigned int i = 0; i <num_scans; i++)  {
       for (unsigned int j = 0; j <num_pixels; j++)  {
//          cth[i][j]=(rand()%18);// 0-15 km random cth, result in meters  
          cth[i][j]=l1cinput->cloud_height;
       }}


  //cloud in cartesian coordinates--   
   for(unsigned int i=0;i<num_scans;i++){
        //satellite in cartesian coordinates--
    //   Xs=Hsat*cos(latsat[i]*M_PI/180)*sin(lonsat[i]*M_PI/180);
   //    Ys=Hsat*sin(lonsat[i]*M_PI/180);
  //     Zs=Hsat*cos(latsat[i]*M_PI/180)*cos(lonsat[i]*M_PI/180);

       for(unsigned int j=0;j<num_pixels;j++){
       //KOENIG VICENTE APPROACH
       /*    
          term1=1-(1-cos(2*latpix[i][j]*M_PI/180))/2;
          term2=Rratio*Rratio*(1-cos(2*latpix[i][j]*M_PI/180))/2;
          Rlocal=Re/sqrt(term1+term2);

          Xc=Rlocal*cos(latpix[i][j]*M_PI/180)*sin(lonpix[i][j]*M_PI/180);
          Yc=Rlocal*sin(lonpix[i][j]*M_PI/180);
          Zc=Rlocal*cos(latpix[i][j]*M_PI/180)*cos(lonpix[i][j]*M_PI/180);

          Rlocal_cloud=((Re+cth[i][j])*(Re+cth[i][j]))/((Rpole+cth[i][j])*(Rpole+cth[i][j]));

          Xdiff=Xs-Xc;
          Ydiff=Ys-Yc;
          Zdiff=Zs-Zc;

          //correction for the line of sight
          e1=Xdiff*Xdiff+Rlocal_cloud*Ydiff*Ydiff+Zdiff*Zdiff;
          e2=2*(Xc*Xdiff+Rlocal_cloud*Yc*Ydiff+Zc*Zdiff);
          e3=Xc*Xc+Rlocal_cloud*Yc*Yc+Zc*Zc-(Re+cth[i][j])*(Re+cth[i][j]);

          cgf=(sqrt(e2*e2-4*e1*e3)-e2)/(2*e1);

          Xcorr=Xc+cgf*Xdiff;
          Ycorr=Yc+cgf*Ydiff;
          Zcorr=Zc+cgf*Zdiff;

          term3=tan(atan2(Ycorr,sqrt(Xcorr*Xcorr+Zcorr*Zcorr)));
          term4=Rratio*Rratio;

          lat_new[i][j]=atan2(term3,term4);
          lon_new[i][j]=atan2(Xcorr,Zcorr);

          rad_new=sqrt(Xcorr*Xcorr +Ycorr*Ycorr+Zcorr*Zcorr);
          lat_new[i][j]=asin(Zcorr/rad_new)*180/M_PI;
          lon_new[i][j]=atan2(Ycorr,Xcorr)*180/M_PI;
*/



          //Displacement vector -- wang et al. 2011 ----
          temp=senapix[i][j]*0.01;

          if(senzpix[i][j]*0.01>=0 && temp>=-180 && temp<=180){

             dv=Hsat*cth[i][j]*tan(senzpix[i][j]*0.01*M_PI/180)/(Hsat-cth[i][j]);

       

             lat_new[i][j]=latpix[i][j]*M_PI/180-dv*cos(temp*M_PI/180+M_PI)/Re;
             lon_new[i][j]=lonpix[i][j]*M_PI/180-dv*sin((temp*M_PI/180+M_PI))/(Re*cos(lat_new[i][j]*M_PI/180));


             if(lon_new[i][j]*180/M_PI>180){
                 lon_new[i][j]-=2*M_PI;
                 }
             if(lon_new[i][j]*180/M_PI<-180){
                 lon_new[i][j]+=2*M_PI;
              }
            
           }

          lat_new[i][j]*=180/M_PI;
          lon_new[i][j]*=180/M_PI;                              

            }}

   cout<<"l2 time_ini.."<<ptime[0]<<"l2 time end.."<<ptime[num_scans-1]<<endl;   
   
   if ((status = nc_close(ncid_L2)))
        check_err(status, __LINE__, __FILE__);


    //writing new L1C file-------------------
    //--------------------------------------
          std::string timestr,missionstr,fname_out, pathstr,monstr, daystr, yearstr, secstr,mistr,hstr,prodstr, gdstr, swtstr, swtnum, extstr,ofilestr,dirstr,datetimestr1,datetimestr2,fdatetimestr1;
          pathstr = "out/";
          fdatetimestr1="2022";
          string GATT_VAL2,GATT_VAL32;

          if(l1cinput->sensor==34){
              senstr="SPEXone";
              GATT_VAL2="SPEXone";
              GATT_VAL32="3567";
              NVIEWS=5;
              NBANDS=400;
              }
          else if (l1cinput->sensor==30){
               senstr="OCI";
              GATT_VAL2="OCI";
              GATT_VAL32="635";
              NVIEWS=2;
              NBANDS=249;
                }
          else if (l1cinput->sensor==35){
               senstr="HARP2";
               GATT_VAL2="HARP2";
               GATT_VAL32="127";
               NVIEWS=90;
               NBANDS=1;
               }
          else{
              cout<<"sensor by default is OCI option 2....."<<endl;
              senstr="OCI";
              GATT_VAL2="OCI";
              GATT_VAL32="635";
              NVIEWS=2;
              NBANDS=249;
               }

          fname_out=pathstr+"PACE_"+senstr+"."+fdatetimestr1+".L1C.5.2km.nc";   
    
      
          string  ATT_NAME="Units", ATT_VAL="degrees",GATT_NAME1="title",GATT_VAL1="PACE OCI Level-1C Data",GATT_NAME2="instrument",GATT_NAME3="processing_version",GATT_VAL3="V1.0",GATT_NAME4="Conventions",GATT_VAL4="CF-1.6";
          string GATT_NAME5="institution",GATT_VAL5="NASA Goddard Space Flight Center, Ocean Biology Processing Group",GATT_NAME6="license",GATT_VAL6="http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
          string GATT_NAME7="naming_authority",GATT_VAL7="gov.nasa.gsfc.sci.oceancolor",GATT_NAME8="keywords_vocabulary",GATT_VAL8="NASA Global Change Master Directory (GCMD) Science Keywords";
          string GATT_NAME9="stdname_vocabulary",GATT_VAL9="NetCDF Climate and Forecast (CF) Metadata Convention",GATT_NAME10="creator_name",GATT_VAL10="NASA/GSFC",GATT_NAME11="creator_email",GATT_VAL11="data@oceancolor.gsfc.nasa.gov";
          string GATT_NAME12="creator_url",GATT_VAL12="http://oceancolor.gsfc.nasa.gov",GATT_NAME13="project",GATT_VAL13="PACE Project",GATT_NAME14="publisher_name",GATT_VAL14="NASA/GSFC";
          string GATT_NAME15="publisher_email",GATT_VAL15="data@oceancolor.gsfc.nasa.gov",GATT_NAME16="publisher_url",GATT_VAL16="http://oceancolor.gsfc.nasa.gov",GATT_NAME17="processing_level",GATT_VAL17="L1C";
          string GATT_NAME18="cdm_data_type",GATT_VAL18="swath",GATT_NAME19="orbit_number",GATT_VAL19="xxx",GATT_NAME20="history",GATT_VAL20="",GATT_NAME21="CDL_version_date",GATT_VAL21="2021-09-10",GATT_NAME22="product_name",GATT_VAL22=fname_out;
          string GATT_NAME23="startDirection",GATT_VAL23=start_dir,GATT_NAME24="endDirection",GATT_VAL24=end_dir,GATT_NAME25="time_coverage_start",GATT_VAL25=tswt_ini,GATT_NAME26="time_coverage_end",GATT_VAL26=tswt_end,GATT_NAME27="date_created",GATT_VAL27="2022-03-25T15:12:41Z",GATT_NAME28="sun_earth_distance",GATT_VAL28="0.990849042172323",GATT_NAME29="terrain_data_source",GATT_VAL29="",GATT_NAME30="spectral_response_function",GATT_VAL30="",GATT_NAME31="systematic_uncertainty_model",GATT_VAL31="",GATT_NAME32="nadir_bin",GATT_NAME33="bin_size_at_nadir",GATT_VAL33="5.2km2";
          string ATT_NAME2="long_name",ATT_NAME3="valid_min",ATT_NAME4="valid_max";;
          string ATT_VAL2="seconds of the day",ATT_VAL3="Latitudes of bin locations",ATT_VAL4="Longitudes of bin locations",ATT_VAL5="UTC time of the day for each gridline",ATT_VAL6="-90",ATT_VAL7="-180",ATT_VAL8="0";
          string ATT_VAL9="+90",ATT_VAL10="+180",ATT_VAL11="86400";
  
          filename_lt = fname_out.c_str();

          cout<<"creating file for L1C coor..."<<filename_lt<<endl;
          NY = num_scans;
          NX = num_pixels;

       
          if ((status = nc_create(filename_lt, NC_CLOBBER | NC_NETCDF4, &ncid_out)))
            check_err(status, __LINE__, __FILE__);


           //DEF DIMENSIONS
             if ((status = nc_def_dim(ncid_out, "bins_across_track", NX, &x_dimid)))
                 check_err(status, __LINE__, __FILE__);
             if ((status = nc_def_dim(ncid_out, "bins_along_track", NY, &y_dimid)))
                 check_err(status, __LINE__, __FILE__);
           //dims for output var
           dimids[0] = y_dimid;
           dimids[1] = x_dimid;
           NDIMS=2;
           dimids2=&y_dimid;

           if ((status = nc_def_dim(ncid_out, "number_of_views", NVIEWS, &v_dimid)))
                  check_err(status, __LINE__, __FILE__);
           if ((status = nc_def_dim(ncid_out, "intensity_bands_per_view", NBANDS, &b_dimid)))
                  check_err(status, __LINE__, __FILE__);
              //define attributes
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME1.c_str(), strlen(GATT_VAL1.c_str()),
                 GATT_VAL1.c_str()))
                 check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME2.c_str(), strlen(GATT_VAL2.c_str()),
                 GATT_VAL2.c_str()))
                 check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME3.c_str(), strlen(GATT_VAL3.c_str()),
                 GATT_VAL3.c_str()))
                check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME4.c_str(), strlen(GATT_VAL4.c_str()),
                 GATT_VAL4.c_str()))
                 check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME5.c_str(), strlen(GATT_VAL5.c_str()),
                 GATT_VAL5.c_str()))
                   check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME6.c_str(), strlen(GATT_VAL6.c_str()),
                 GATT_VAL6.c_str()))
                 check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME7.c_str(), strlen(GATT_VAL7.c_str()),
                GATT_VAL7.c_str()))
                check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME8.c_str(), strlen(GATT_VAL8.c_str()),
                  GATT_VAL8.c_str()))
                check_err(status, __LINE__, __FILE__);
           if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME9.c_str(), strlen(GATT_VAL9.c_str()),
                GATT_VAL9.c_str()))
               check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME10.c_str(), strlen(GATT_VAL10.c_str()),
                GATT_VAL10.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME11.c_str(), strlen(GATT_VAL11.c_str()),
                GATT_VAL11.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME12.c_str(), strlen(GATT_VAL12.c_str()),
                GATT_VAL12.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME13.c_str(), strlen(GATT_VAL13.c_str()),
                 GATT_VAL13.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME14.c_str(), strlen(GATT_VAL14.c_str()),
                   GATT_VAL14.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME15.c_str(), strlen(GATT_VAL15.c_str()),
                GATT_VAL15.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME16.c_str(), strlen(GATT_VAL16.c_str()),
                GATT_VAL16.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME17.c_str(), strlen(GATT_VAL17.c_str()),
                 GATT_VAL17.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME18.c_str(), strlen(GATT_VAL18.c_str()),
                 GATT_VAL18.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME19.c_str(), strlen(GATT_VAL19.c_str()),
                 GATT_VAL19.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME20.c_str(), strlen(GATT_VAL20.c_str()),
                 GATT_VAL20.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME21.c_str(), strlen(GATT_VAL21.c_str()),
                GATT_VAL21.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME22.c_str(), strlen(GATT_VAL22.c_str()),
                GATT_VAL22.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME23.c_str(), strlen(GATT_VAL23.c_str()),
                GATT_VAL23.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME24.c_str(), strlen(GATT_VAL24.c_str()),
                GATT_VAL24.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME25.c_str(), strlen(GATT_VAL25.c_str()),
                 GATT_VAL25.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME26.c_str(), strlen(GATT_VAL26.c_str()),
                GATT_VAL26.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME27.c_str(), strlen(GATT_VAL27.c_str()),
                 GATT_VAL27.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME28.c_str(), strlen(GATT_VAL28.c_str()),
                GATT_VAL28.c_str()))
               check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME29.c_str(), strlen(GATT_VAL29.c_str()),
                GATT_VAL29.c_str()))
               check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME30.c_str(), strlen(GATT_VAL30.c_str()),
                GATT_VAL30.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME31.c_str(), strlen(GATT_VAL31.c_str()),
                GATT_VAL31.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME32.c_str(), strlen(GATT_VAL32.c_str()),
                GATT_VAL32.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME33.c_str(), strlen(GATT_VAL33.c_str()),
                   GATT_VAL33.c_str()))
                 check_err(status, __LINE__, __FILE__);

           //define groups                  check_err(status, __LINE__, __FILE__);
             if ((status = nc_def_grp(ncid_out, "geolocation_data", &grp_coor)))  //netcdf-4
                check_err(status, __LINE__, __FILE__);
             //def var
             if ((status = nc_def_var(grp_coor, "latitude", NC_FLOAT, NDIMS,
                dimids, &varid1)))
                check_err(status, __LINE__, __FILE__);
             if ((status = nc_def_var(grp_coor, "longitude", NC_FLOAT, NDIMS,
                dimids, &varid2)))
                check_err(status, __LINE__, __FILE__);
             if ((status = nc_def_var(grp_coor, "time_nadir", NC_DOUBLE, NDIMS2,
                dimids2, &varid3)))
                check_err(status, __LINE__, __FILE__);

             if ((status = nc_def_var_deflate(grp_coor, varid1, shuffle, deflate,deflate_level)))
               check_err(status, __LINE__, __FILE__);
             if ((status = nc_def_var_deflate(grp_coor, varid2, shuffle, deflate,deflate_level)))
               check_err(status, __LINE__, __FILE__);
             if ((status = nc_def_var_deflate(grp_coor, varid3, shuffle, deflate,deflate_level)))
               check_err(status, __LINE__, __FILE__);

             if (status=nc_def_var_fill(grp_coor, varid1, 1, &fill_value))
                 check_err(status, __LINE__, __FILE__);
             if (status=nc_def_var_fill(grp_coor, varid2, 1, &fill_value))
                check_err(status, __LINE__, __FILE__);
             if (status=nc_def_var_fill(grp_coor, varid3, 1, &fill_value2))
                 check_err(status, __LINE__, __FILE__);

            //define attributes of vars in groups
             //long_name
             if (status=nc_put_att_text(grp_coor, varid1, ATT_NAME2.c_str(), strlen(ATT_VAL3.c_str()), ATT_VAL3.c_str()))
               check_err(status, __LINE__, __FILE__);
          
             if (status=nc_put_att_text(grp_coor, varid2, ATT_NAME2.c_str(), strlen(ATT_VAL4.c_str()), ATT_VAL4.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (status=nc_put_att_text(grp_coor, varid3, ATT_NAME2.c_str(), strlen(ATT_VAL5.c_str()), ATT_VAL5.c_str()))
                 check_err(status, __LINE__, __FILE__);
            //units
             if (status=nc_put_att_text(grp_coor, varid1, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
                 check_err(status, __LINE__, __FILE__);
           
             if (status=nc_put_att_text(grp_coor, varid2, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (status=nc_put_att_text(grp_coor, varid3, ATT_NAME.c_str(), strlen(ATT_VAL2.c_str()), ATT_VAL2.c_str()))
                 check_err(status, __LINE__, __FILE__);


             if (status=nc_put_att_float(grp_coor, varid1,"_FillValue",NC_FLOAT,1 ,&fill_value))
                 check_err(status, __LINE__, __FILE__);
             if (status=nc_put_att_float(grp_coor, varid2,"_FillValue",NC_FLOAT,1 ,&fill_value))
                 check_err(status, __LINE__, __FILE__);
             if (status=nc_put_att_double(grp_coor, varid3,"_FillValue",NC_DOUBLE,1 ,&fill_value2))
                 check_err(status, __LINE__, __FILE__);

         //valid_min
             if (status=nc_put_att_text(grp_coor, varid1, ATT_NAME3.c_str(), strlen(ATT_VAL6.c_str()), ATT_VAL6.c_str()))
               check_err(status, __LINE__, __FILE__);            
             if (status=nc_put_att_text(grp_coor, varid2, ATT_NAME3.c_str(), strlen(ATT_VAL7.c_str()), ATT_VAL7.c_str()))
                check_err(status, __LINE__, __FILE__);
             if (status=nc_put_att_text(grp_coor, varid3, ATT_NAME3.c_str(), strlen(ATT_VAL8.c_str()), ATT_VAL8.c_str()))
               check_err(status, __LINE__, __FILE__);

        //valid_max
             if (status=nc_put_att_text(grp_coor, varid1, ATT_NAME4.c_str(), strlen(ATT_VAL9.c_str()), ATT_VAL9.c_str()))
               check_err(status, __LINE__, __FILE__);
             if (status=nc_put_att_text(grp_coor, varid2, ATT_NAME4.c_str(), strlen(ATT_VAL10.c_str()), ATT_VAL10.c_str()))
                 check_err(status, __LINE__, __FILE__);
             if (status=nc_put_att_text(grp_coor, varid3, ATT_NAME4.c_str(), strlen(ATT_VAL11.c_str()), ATT_VAL11.c_str()))
                 check_err(status, __LINE__, __FILE__);

          //leave define mode-----------------------
             if ((status = nc_enddef(grp_coor))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);

          //writing the whole thing
             if ((status = nc_put_var_float(grp_coor, varid1, &lat_new[0][0])))
                 check_err(status, __LINE__, __FILE__);

             if ((status = nc_put_var_float(grp_coor, varid2, &lon_new[0][0])))
                 check_err(status, __LINE__, __FILE__);

             if ((status = nc_put_var_double(grp_coor, varid3, &ptime[0])))
                  check_err(status, __LINE__, __FILE__);

            //close file
             if ((status = nc_close(ncid_out)))
                 check_err(status, __LINE__, __FILE__);                                   

  delete [] (lat_new);
  delete [] (lon_new);

  delete [] (posr);
  delete [] (velr);
  delete [] (ptime);
  delete [] (latpix);
  delete [] (lonpix);
  delete [] (senapix);
  delete [] (senzpix);
  delete [] (cth);
  delete [] (latsat);
  delete [] (lonsat);

  return 0;

}

    
    int32_t L1C::azmean_swt3(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,float **lat_tot,float **lon_tot){
  float az,lat1_l,lat2,lon1_l,lon2,dlambda,az_r;
  float sum=0,tot=0;
  int32_t swath_scans=-1;
  size_t azc=0;
  int nadpix=-1;
  nadpix=l1cfile->nadpix;
  cout<<"in azmean_swt3......nadpix............"<<nadpix<<endl;


   swath_scans = l1cfile->swath_scans;

  
     for (int i = 0;i < swath_scans-1;i++) {
        if(lat_tot[i][nadpix]>=-10. && lat_tot[i][nadpix]<=10.){
            //assign central positions---
            lat1_l = lat_tot[i][nadpix] * M_PI / 180.;
            lat2 = lat_tot[i + 1][nadpix] * M_PI / 180.;
            lon1_l = lon_tot[i][nadpix] * M_PI / 180.;
            lon2 = lon_tot[i + 1][nadpix] * M_PI / 180.;
            dlambda = (lon2 - lon1_l);
            //bearing---
            //make sure there are not NAN values ----------
            az = atan2(sin(dlambda) * cos(lat2), cos(lat1_l) * sin(lat2) - sin(lat1_l) * cos(lat2) * cos(dlambda));//this is the satellite track bearing in radians


            if (az > M_PI || az < -M_PI) {
                cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "az in degrees.." << az * 180 / M_PI << endl;
                exit(1);
                    }
            if (isnan(az) < 1) {
                az_r = (az + M_PI / 2.);
                sum = az_r * 180. / M_PI;
                if (sum < 0.0001 && sum>0) sum = 0.0;
                tot += sum;
                azc++;
                 }
            else {
                sum = NAN;
                  }
  
              }

           }

        float mean_az_east = tot / (azc);

        cout << "mean_az_east.." << mean_az_east << "swath_scans.." <<swath_scans<< endl;

        l1cfile->mean_az_east = mean_az_east;

    return 0;
}


int32_t L1C::calc_biny_dist(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,float* lati2, float* loni2) {
        std::string str,pathstr,filestr,ofilestr;
        int32_t  num_gridlines = 0;//4000 gridlines at 5.2 km and 8000 at 2.6 km resolution
        float res = 0;
        double s12;
        float mean_az_east;
        float az_track;
        double newlat_b, newlon_b,newlat_f,newlon_f;
        float lat0=0.;
        double tcross1=-1.;
        float loncross1=-999.;
        int bina=0,binb=0;
        float oldlat,oldlon;
        vector<float> latvec_f, lonvec_f,latvec_b, lonvec_b;

        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());

        num_gridlines = l1cfile->num_gridlines;//first estimate based on resolution of 5.2 km 
        res = (l1cinput->grid_resolution);
        tcross1=l1cfile->eqt;
        loncross1=l1cfile->orbit_node_lon;


        cout<<"swath#.."<<swtd<<"tcross..."<<tcross1<<"loncross.."<<loncross1<<"num_gridlines..."<<num_gridlines<<endl;

        mean_az_east = l1cfile->mean_az_east;
        if(mean_az_east<-180){cout<<"ERROR mean_az_east should be >=-180..."<<mean_az_east<<endl; exit(1);}

        az_track = mean_az_east - 90;

        if(tcross1>0. && loncross1>=-180. && loncross1<=180.){
       
            //compute first pair of gridcenter points  bracketing equator...
            geod.Direct(lat0, loncross1, az_track, res * 1000/2, newlat_f, newlon_f);
            latvec_f.push_back(newlat_f);
            lonvec_f.push_back(newlon_f);
            geod.Direct(lat0, loncross1, az_track+180., res * 1000/2, newlat_b, newlon_b);
            latvec_b.push_back(newlat_b);
            lonvec_b.push_back(newlon_b);      

            while(bina<(num_gridlines)/2 && binb<(num_gridlines)/2){
//forward

             oldlat=newlat_f;
             oldlon=newlon_f;

             geod.Direct(oldlat,oldlon,az_track, res * 1000, newlat_f, newlon_f);
             bina++;
             latvec_f.push_back(newlat_f);
             lonvec_f.push_back(newlon_f);

//backward         
             oldlat=newlat_b;
             oldlon=newlon_b;

             geod.Direct(oldlat,oldlon,az_track+180, res * 1000, newlat_b, newlon_b);
             binb++;
             latvec_b.push_back(newlat_b);
             lonvec_b.push_back(newlon_b);             
                  }

        }
       
 
        //forward
        int k=0;
        for(int i=(num_gridlines)/2;i<num_gridlines;i++){
           lati2[i]=latvec_f[k];
           loni2[i]=lonvec_f[k];
           k++;
          }

        //backward
        k=0;
        for(int i=(num_gridlines)/2-1;i>=0;i--){
           lati2[i]=latvec_b[k];
           loni2[i]=lonvec_b[k];
           k++;
          }

        for(int i=0;i<num_gridlines-1;i++){    
           geod.Inverse(lati2[i],loni2[i],lati2[i+1],loni2[i+1], s12);     
        }

        cout<<"bina..."<<bina<<"binb..."<<binb<<"number of gridlines..."<<binb+bina<<endl;

        latvec_f.clear();
        lonvec_f.clear();
        latvec_b.clear();
        lonvec_b.clear();

        return 0;
    }


int32_t L1C::calc_biny(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv){
  int16_t ix=-1,ngridlines,gd=0,bina,binb;
  double tg,mot=0.,tini=0.,tend=0.;


  mot = ((l1cinput->grid_resolution) * 1000) /mgv;//in seconds
  l1cfile->mot=mot;


  if((l1cinput->grid_resolution) * 1000==5200)
   ngridlines=4000;
  else if((l1cinput->grid_resolution) * 1000==5200/2)
   ngridlines=8000;
  else{cout<<"error wrong number of gridlines at calc_biny...."<<endl; exit(1);}  

  tini=tswt[0];
  ix=norbs-1;
  tend=tswt[ix];


        int flag_time=-1;

    if (tcross > 0.0) {
          flag_time=0;
          cout<<"computing time series assuming mean swath velocity..for swath#."<<swt<<endl;

         tg = tcross - mot / 2;          
         gd=ngridlines/2-1;
         binb=1;

         while(binb<ngridlines/2 && tg>=tini+mot){
             binb++;
             tg -= mot;
             gd--;
         }

         tg = tcross + mot / 2;
         gd=ngridlines/2;
         bina=1;

         while(bina<ngridlines/2 && tg<=tend-mot){
             bina++; 
             tg += mot;
             gd++;        
         }

         cout<<"bina.."<<bina<<"binb.."<<binb<<endl;

         l1cfile->num_gridlines = binb + bina;
         l1cfile->binyb = binb; 
         l1cfile->binya = bina; 

            cout<<"number of L1C gridlines along-track..."<<l1cfile->num_gridlines<<"for swath #.."<<swt<<endl;
             }//end equat crossing
     else{
       cout<<"time series not possible for swath #.."<<swt<<"tcross<0...."<<endl;
       flag_time=1;
       }


   return flag_time;
}


int32_t L1C::swath_latlon(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,int16_t* swtd_id, int16_t* file_id, int16_t* nfiles_swt,float** lat_tot, float** lon_tot) {
        const char* ptstr;
        int dimid, status;
        std::string str,pathstr,filestr,ofilestr;
        int16_t fi = 0;
        int32_t  num_scans_tot = 0, swt_orb_rec = 0, num_orb_rec;
        int32_t n_files, iscan = 0;//4000 gridlines at 5.2 km and 8000 at 2.6 km resolution
     //   size_t first_line=5000,last_line=5000;
        int asc_mode=-1;

        asc_mode = l1cfile->orb_dir;
        nswath = l1cfile->nswath;
        num_pixels = l1cfile->npix;
        n_files = l1cfile->ifiles.size();
        int16_t nfiles_swath = nfiles_swt[swtd - 1];
        int nadpix = l1cfile->nadpix;

        cout<<"nfiles.."<<n_files<<"swtd_id.."<<swtd_id[0]<<"nadpix..."<<nadpix<<"asc_mode..."<<asc_mode<<endl;

        //opening swath files----
        for (int i = 0;i < n_files;i++) {
            if (swtd == swtd_id[i]) {

                fi = file_id[i];
                str = l1cfile->ifiles[fi - 1];
                ptstr = str.c_str();

                cout<<"fi.."<<fi<<"opening file..."<<ptstr<<endl;

                status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
                if (status != NC_NOERR) {
                      fprintf(stderr, "-E- failed nc_open\n");
                    exit(EXIT_FAILURE);
               }
                //Open dimensions
                // num_scans
                status = nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid);
                if (status != NC_NOERR) {
                    fprintf(stderr, "-E- Error reading number_of_scans.\n");
                    exit(EXIT_FAILURE);
                }
                nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

                status = nc_inq_dimid(ncid_L1B, "ccd_pixels", &dimid);
                if (status != NC_NOERR) {
                    fprintf(stderr, "-E- Error reading number_of_scans.\n");
                    exit(EXIT_FAILURE);
                }
                nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);

                lat = allocate2d_float(num_scans, num_pixels);
                lon = allocate2d_float(num_scans, num_pixels);

                //open geo data
                status = nc_inq_grp_ncid(ncid_L1B, "geolocation_data", &geoGrp);
                check_err(status, __LINE__, __FILE__);
                status = nc_inq_varid(geoGrp, "latitude", &latId);//scans x velements
                check_err(status, __LINE__, __FILE__);
                status = nc_get_var_float(geoGrp, latId, &lat[0][0]);
                check_err(status, __LINE__, __FILE__);

                status = nc_inq_varid(geoGrp, "longitude", &lonId);//scans x velements
                check_err(status, __LINE__, __FILE__);
                status = nc_get_var_float(geoGrp, lonId, &lon[0][0]);
                check_err(status, __LINE__, __FILE__);

                //concatenate 1-D time and 2-D geo vectors for the swath---
                //check orbit direction changes (ascending to descending) within the same file if the file is the lastof the swath

                num_scans_tot += num_scans;
                num_orb_rec = num_scans;
                swt_orb_rec += num_orb_rec;

                for (unsigned int j = 0;j < num_scans;j++) {
                     for (unsigned int k = 0;k < num_pixels;k++) {
                        lat_tot[j + iscan][k] = lat[j][k];
                        lon_tot[j + iscan][k] = lon[j][k];
                    }
                }

                //close file if not the last
                status = nc_close(ncid_L1B);
                check_err(status, __LINE__, __FILE__);

                delete[](lat);
                delete[](lon);

                iscan += num_scans;

/*
 //ascending
              for (size_t k =num_scans-1;k>0 ;k--) {
                    //checking partial granules--
                   if(asc_mode==1 && i==0 && lat[k-1][nadpix]>lat[k][nadpix] && first_line==5000){
                         first_line=k;
                         cout<<"first line.."<<first_line<<endl;
                         }
                    }
              for (size_t k = 0;k<num_scans-1 ;k++) {
                     if(asc_mode==1 && lat[k+1][nadpix]<lat[k][nadpix] && i==n_files-1 && last_line==5000){//assuming change of direction NH, last granule
                         last_line=k;
                         cout<<"last line.."<<last_line<<endl;
                         }
                     }
//descending
              for (size_t k =num_scans-1;k>0 ;k--) {                     //descending
                      if(asc_mode==0 && lat[k-1][nadpix]<lat[k][nadpix] && i==0 && first_line==5000){//assuming change of direction SH, first granule
                         first_line=k;
                         cout<<"first line.."<<first_line<<endl;
                         }
                   }

              for (size_t k = 0;k<num_scans-1 ;k++) {
                     if(asc_mode==0 && lat[k+1][nadpix]>lat[k][nadpix] && i==n_files-1 && last_line==5000){//assuming change of direction NH, last granule
                         last_line=k;
                         cout<<"last line.."<<last_line<<endl;
                         }
                   }
  */
            }//end if
        }//end for

       l1cfile->swath_scans=num_scans_tot;

       cout << "num_scans_tot.." << num_scans_tot << "should be equal to numfiles x numscans/file.." << num_scans * nfiles_swath << "and swath records.." << swt_orb_rec << endl;

        return 0;
    }


int32_t L1C::binL1C_sbs_line_l2(L1C *l1c,l2_str *l2str,l1c_filehandle *l1cfile,L1C_input *l1cinput,float ****binmean_prod,int ****bincount,size_t sline,int granid){
        std::string str,ofilestr;
        int status, ncid_out,p_dimid,v_dimid,x_dimid,y_dimid,varid1,varid2, varid_count,varid_prod,NDIMS,NDIMS3;
        int32_t bin_xpix, bin_ypix;
        int dimids[2],dimids3[4];
        const char* filename_lt;
        int32_t NY = -1, NX = -1;
        float****data_out2=nullptr,** lat_out=nullptr, ** lon_out=nullptr;
        int ****data_out=nullptr;
        int flag_inpix;//first index file id, second index line #
        bool boolbin1;
        char* ifile_char;
        std::string ifile_str;
        int16_t selyear = -1, selmon = -1, selday = -1;
        string gridname, azeast_name;
        int grp_obs,grp_coor;
        size_t  NVIEWS=-1,NL2PRODS=-1;
        int view=0;
        short gd_row, gd_col;
        std::string fname_out, pathstr, senstr, monstr, daystr, yearstr, prodstr, gdstr, swtstr, swtnum, extstr, granstr,timestr,azstr,missionstr;
        size_t  iprod;

        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());

        num_pixels = l1cfile->npix;
        num_scans=l1cfile->nscan;
        num_pixels=l1cfile->npix;

        ifile_str = l1cinput->files[0];
        ifile_char = &ifile_str[0];
        file_format format = getFormat(ifile_char);


                    //allocate mem for lat/lon/azeast pointers---
        //determine # of gd groups to be processede
        //gdlines_group=num_gridlines;//ONE BIG GROUP

        //--NORMALIZED  gridline longitude is -180 to 180 degrees before comparison with longitude extracted from image pixels---
        for (int i = 0; i < l1cfile->num_gridlines; i++) {
            for (int j = 0; j < l1cfile->nbinx; j++) {
                if (l1cfile->lon_gd[i][j] < -180.)l1cfile->lon_gd[i][j] += 360.;
                if (l1cfile->lon_gd[i][j] > 180.)l1cfile->lon_gd[i][j] -= 360.;
            }
            }

//********* LOOP FILES ************************************

                //radiance variables
                if(format.type==FT_SPEXONE){
                   NVIEWS= 5;//should be 5 but only provided 1
                }
                else if(format.type==FT_HARP){
                   NVIEWS= 1;//HARP-2 10 but 60 for 669 nm
                }
                else{//OCIS
                   NVIEWS= l1cfile->n_views;//for OCI
                   NL2PRODS=l2str->nl2prod;
                }



//BIG SCANLINE LOOP----------------------------------
                //********************************************************************************
                 //-----reading line by line-----------------------------------                  

   //*********** BIG LOOP M_PIXEL *****************************************************************
                     for (unsigned int pix = 0;pix < num_pixels;pix++) {
                            gd_row=0,gd_col=0;
                            flag_inpix = 0;
                            iprod=0;
                             //     cout<<"pix..."<<pix+1<<"latpix..."<<l1cstr->latpix[pix]<<"lonpix..."<<l1cstr->lonpix[pix]<<"sensor azimuth.."<<l1cstr->senazpix[pix]/100<<endl;

                            boolbin1 = sbs2_l1c(l1cinput, l1cfile->num_gridlines, l1cfile->nbinx, l1cfile->lat_asort, l1cfile->index_xy,l2str->latpix[pix], l2str->lonpix[pix], l1cfile->lon_gd, &gd_row, &gd_col);

                            bin_ypix = gd_row + 1;
                            bin_xpix = gd_col + 1;
                   //         cout<<"boolbin1.."<<boolbin1<<"bin_ypix............................."<<bin_ypix<<"bin_xpix........................."<<bin_xpix<<endl;

                            if (boolbin1 == 1) flag_inpix = 1;
                            //************ BINNING ***************************
                         //assign identified pixel to Lt and bin stat arrays-----


                            if (flag_inpix == 1 && bin_xpix >= 1 && bin_ypix >= 1 && bin_xpix <= l1cfile->nbinx && bin_ypix <= l1cfile->num_gridlines && l2str->latpix[pix]>=l1cinput->south && l2str->latpix[pix]<=l1cinput->north && l2str->lonpix[pix]>=l1cinput->west && l2str->lonpix[pix]<=l1cinput->east) {
                                INPIX++;
                              if(format.type==FT_SPEXONE){
                                  view=0;
                               /*    while (sb<NBANDS) {
                                      if (l1cstr->Lt[pix][ib] > 0.) {
                                              binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] + l1cstr->Lt[pix][ib];
                                              bincount[view][ib][bin_ypix - 1][bin_xpix - 1] += 1;
                                           }
                                    ib++;
                                    sb++;
                                     }//end while intensity bands
                                  */
                                 }

                               else if(format.type==FT_HARP){
                                  view=0;
                            /*      if (l1cstr->Lt[pix][ib] > 0.) {
                                              binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] + l1cstr->Lt[pix][ib];
                                              bincount[view][ib][bin_ypix - 1][bin_xpix - 1] += 1;
                                           }
                                    ib++;
                                    sb++;*/
                                 }
                               else{ //OCIS

                               if(l2str->tilt[sline]<=0) view=0; else view=1;
                                 while (iprod<NL2PRODS) {
                       //                  cout<<"binning OCIS L2 products........................................................................................................."<<endl;
                         //                cout<<"l2str->l2prod[iprod][pix]..."<<l2str->l2prod[iprod][pix]<<endl;
                                      if (l2str->l2prod[iprod][pix] > 0.) {
                                              cout<<"**********************************************************************************************************************************************************"<<endl;
                                              cout<<"l2str->l2prod[iprod][pix]..."<<l2str->l2prod[iprod][pix]*l2str->slopeprod[iprod]+l2str->offsetprod[iprod]<<"bin_ypix.."<<bin_ypix<<"bin_xpix.."<<bin_xpix<<endl;
                                              cout<<"pix..."<<pix+1<<"latpix..."<<l2str->latpix[pix]<<"lonpix..."<<l2str->lonpix[pix]<<endl;
                                               cout<<"**********************************************************************************************************************************************************"<<endl;

                                              binmean_prod[view][iprod][bin_ypix - 1][bin_xpix - 1] = binmean_prod[view][iprod][bin_ypix - 1][bin_xpix - 1] + l2str->l2prod[iprod][pix]*l2str->slopeprod[iprod]+l2str->offsetprod[iprod];
                                              bincount[view][iprod][bin_ypix - 1][bin_xpix - 1] += 1;
                                           }
                                    iprod++;

                                }//end products

                             }

                            }//end if INPIX==1



    //********** pixel AREA WEIGHTING ************************
                            if (flag_inpix == 0) OUTPIX++;

                            TOTPIX++;

                         }//end pixels loop

                     //END OF FILE--****** reset starting sline index if last line processed and within k box
                        //go for another granule if it is not the last one---------------

         //***** skip lines and files if sline is not within k group ********
              //increase k and screen files and lines again starting from the last binning


                //**********************************************************************************************************************************
                //writing each granule ---------------------------------
                //*****************************************************************************

                if (sline==num_scans-1) {
                    cout << "writing file #....................................................................................................................."  <<l1cfile->l1b_name<< endl;
                    //write mean Lt as nc file---
                   //**********************************************************************
                   //**********************************************************************8
                   //create Lt file
                    pathstr = "";
                    missionstr="PACE";
                    senstr = "OCIS_";
                    timestr="T00:00:00Z";
                     selday = l1cinput->selday;
                    selmon = l1cinput->selmon;
                    selyear = l1cinput->selyear;
                    monstr = std::to_string(selmon);
                    daystr = std::to_string(selday);
                    yearstr = std::to_string(selyear);
                    prodstr = "binL2_sbs_LINEBYLINE";
                    swtstr = "_swt";
                    granstr = "_" + std::to_string(granid);


                    extstr = ".nc";
                    string GATT_NAME1,GATT_VAL1,GATT_NAME2,GATT_VAL2;

                     if(format.type==FT_SPEXONE){
                         senstr = "SPEXONE";
                         GATT_NAME1="title",GATT_VAL1="PACE SPEXone Level-1C Data",GATT_NAME2="instrument",GATT_VAL2="SPEXone";
                     }
                     else//OCIS
                     {
                         senstr = "OCI";
                         GATT_NAME1="title",GATT_VAL1="PACE OCI Level-1C Data",GATT_NAME2="instrument",GATT_VAL2="OCI";
                     }

            //        fname_out = pathstr + missionstr + "_" + senstr + "." + yearstr + monstr + daystr + timestr + prodstr+granstr+extstr;
                    ofilestr=std::string(l1cinput->ofile);
                    fname_out = pathstr + ofilestr+"_"+ prodstr+granstr+extstr;

                    string  ATT_NAME="Units", ATT_VAL="degrees",GATT_NAME3="processing_version",GATT_VAL3="V1.0",GATT_NAME4="Conventions",GATT_VAL4="CF-1.6";
                    string GATT_NAME5="institution",GATT_VAL5="NASA Goddard Space Flight Center, Ocean Biology Processing Group",GATT_NAME6="license",GATT_VAL6="http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
                    string GATT_NAME7="naming_authority",GATT_VAL7="gov.nasa.gsfc.sci.oceancolor",GATT_NAME8="keywords_vocabulary",GATT_VAL8="NASA Global Change Master Directory (GCMD) Science Keywords";
                    string GATT_NAME9="stdname_vocabulary",GATT_VAL9="NetCDF Climate and Forecast (CF) Metadata Convention",GATT_NAME10="creator_name",GATT_VAL10="NASA/GSFC",GATT_NAME11="creator_email",GATT_VAL11="data@oceancolor.gsfc.nasa.gov";
                    string GATT_NAME12="creator_url",GATT_VAL12="http://oceancolor.gsfc.nasa.gov",GATT_NAME13="project",GATT_VAL13="PACE Project",GATT_NAME14="publisher_name",GATT_VAL14="NASA/GSFC";
                    string GATT_NAME15="publisher_email",GATT_VAL15="data@oceancolor.gsfc.nasa.gov",GATT_NAME16="publisher_url",GATT_VAL16="http://oceancolor.gsfc.nasa.gov",GATT_NAME17="processing_level",GATT_VAL17="L1C";
                    string GATT_NAME18="cdm_data_type",GATT_VAL18="swath",GATT_NAME19="orbit_number",GATT_VAL19="12345",GATT_NAME20="history",GATT_VAL20="",GATT_NAME21="CDL_version_date",GATT_VAL21="2021-09-10",GATT_NAME22="product_name",GATT_VAL22=fname_out;
                    string GATT_NAME23="startDirection",GATT_VAL23="Ascending",GATT_NAME24="endDirection",GATT_VAL24="Ascending",GATT_NAME25="time_coverage_start",GATT_VAL25=yearstr+"-"+monstr+"-"+daystr+"-"+timestr,GATT_NAME26="time_coverage_end",GATT_VAL26=yearstr+"-"+monstr+"-"+daystr+"-"+timestr,GATT_NAME27="date{_created",GATT_VAL27="2021-09-10T15:12:41Z",GATT_NAME28="sun_earth_distance",GATT_VAL28="0.990849042172323",GATT_NAME29="terrain_data_source",GATT_VAL29="",GATT_NAME30="spectral_response_function",GATT_VAL30="",GATT_NAME31="systematic_uncertainty_model",GATT_VAL31="",GATT_NAME32="nadir_bin",GATT_VAL32="12345",GATT_NAME33="bin_size_at_nadir",GATT_VAL33="5.2km2",GATT_NAME34="L2_products",GATT_VAL34=l1cinput->l2prod;

                    l1cfile->gridname = fname_out.c_str();
                    filename_lt = fname_out.c_str();
                    cout<<"filename_lt.."<<filename_lt<<endl;

                    if ((status = nc_create(filename_lt, NC_CLOBBER | NC_NETCDF4, &ncid_out)))
                        check_err(status, __LINE__, __FILE__);
                     //define dims
                     // Define the dimensions,vars and attributes at the root level
                    NDIMS=2;
                    NDIMS3=4;
                    NY = l1cfile->num_gridlines;
                    NX = l1cfile->nbinx;
                    //NVIEWS= 2;//for OCI
                   // NBANDS=249;//for OCI

/*
                    NY1 = l1cfile->NY1;//these are indexes not line #!!!!!!!!!!!!!!!!!
                    NY2 = l1cfile->NY2;

                    cout << "rowgrid ini.." << NY1 << "rowgrid end.." << NY2 << endl;
                    NY = NY2 - NY1 + 1;
                    cout << "NY.." << NY << "NX.." << NX << "NY1.."<<NY1<<"NY2.."<<NY2<<endl;
*/
  //DEF DIMENSIONS
                    if ((status = nc_def_dim(ncid_out, "bins_across_track", NX, &x_dimid)))
                        check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_dim(ncid_out, "bins_along_track", NY, &y_dimid)))
                           check_err(status, __LINE__, __FILE__);
                    //dims for output var
                    dimids[0] = y_dimid;
                    dimids[1] = x_dimid;

                    if ((status = nc_def_dim(ncid_out, "number_of_views", NVIEWS, &v_dimid)))
                       check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_dim(ncid_out, "l2_products", NL2PRODS, &p_dimid)))
                      check_err(status, __LINE__, __FILE__);

                    dimids3[0] = y_dimid;
                    dimids3[1] = x_dimid;
                    dimids3[2] = v_dimid;//
                    dimids3[3] = p_dimid;

                    //define attributes
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME1.c_str(), strlen(GATT_VAL1.c_str()),
                        GATT_VAL1.c_str()))
                       check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME2.c_str(), strlen(GATT_VAL2.c_str()),
                         GATT_VAL2.c_str()))
                       check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME3.c_str(), strlen(GATT_VAL3.c_str()),
                         GATT_VAL3.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME4.c_str(), strlen(GATT_VAL4.c_str()),
                            GATT_VAL4.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME5.c_str(), strlen(GATT_VAL5.c_str()),
                         GATT_VAL5.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME6.c_str(), strlen(GATT_VAL6.c_str()),
                          GATT_VAL6.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME7.c_str(), strlen(GATT_VAL7.c_str()),
                          GATT_VAL7.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME8.c_str(), strlen(GATT_VAL8.c_str()),
                          GATT_VAL8.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME9.c_str(), strlen(GATT_VAL9.c_str()),
                          GATT_VAL9.c_str()))
                         check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME10.c_str(), strlen(GATT_VAL10.c_str()),
                          GATT_VAL10.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME11.c_str(), strlen(GATT_VAL11.c_str()),
                           GATT_VAL11.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME12.c_str(), strlen(GATT_VAL12.c_str()),
                            GATT_VAL12.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME13.c_str(), strlen(GATT_VAL13.c_str()),
                             GATT_VAL13.c_str()))
                     check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME14.c_str(), strlen(GATT_VAL14.c_str()),
                             GATT_VAL14.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME15.c_str(), strlen(GATT_VAL15.c_str()),
                                 GATT_VAL15.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME16.c_str(), strlen(GATT_VAL16.c_str()),
                             GATT_VAL16.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME17.c_str(), strlen(GATT_VAL17.c_str()),
                             GATT_VAL17.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME18.c_str(), strlen(GATT_VAL18.c_str()),
                             GATT_VAL18.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME19.c_str(), strlen(GATT_VAL19.c_str()),
                             GATT_VAL19.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME20.c_str(), strlen(GATT_VAL20.c_str()),
                             GATT_VAL20.c_str()))
                              check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME21.c_str(), strlen(GATT_VAL21.c_str()),
                             GATT_VAL21.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME22.c_str(), strlen(GATT_VAL22.c_str()),
                             GATT_VAL22.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME23.c_str(), strlen(GATT_VAL23.c_str()),
                              GATT_VAL23.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME24.c_str(), strlen(GATT_VAL24.c_str()),
                              GATT_VAL24.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME25.c_str(), strlen(GATT_VAL25.c_str()),
                           GATT_VAL25.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME26.c_str(), strlen(GATT_VAL26.c_str()),
                           GATT_VAL26.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME27.c_str(), strlen(GATT_VAL27.c_str()),
                             GATT_VAL27.c_str()))
                              check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME28.c_str(), strlen(GATT_VAL28.c_str()),
                             GATT_VAL28.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME29.c_str(), strlen(GATT_VAL29.c_str()),
                             GATT_VAL29.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME30.c_str(), strlen(GATT_VAL30.c_str()),
                             GATT_VAL30.c_str()))
                              check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME31.c_str(), strlen(GATT_VAL31.c_str()),
                             GATT_VAL31.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME32.c_str(), strlen(GATT_VAL32.c_str()),
                              GATT_VAL32.c_str()))
                                check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME33.c_str(), strlen(GATT_VAL33.c_str()),
                                GATT_VAL33.c_str()))
                                check_err(status, __LINE__, __FILE__);
                     if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME34.c_str(), strlen(GATT_VAL34.c_str()),
                                GATT_VAL34.c_str()))
                                check_err(status, __LINE__, __FILE__);

                    //define groups                  check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_grp(ncid_out, "geolocation_data", &grp_coor)))  //netcdf-4
                         check_err(status, __LINE__, __FILE__);
                     //def var grp1
                    if ((status = nc_def_var(grp_coor, "latitude", NC_FLOAT, NDIMS,
                         dimids, &varid1)))
                          check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_var(grp_coor, "longitude", NC_FLOAT, NDIMS,
                            dimids, &varid2)))
                          check_err(status, __LINE__, __FILE__);
                    //leave define mode-----------------------
                    if ((status = nc_enddef(grp_coor))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);


                    if ((status = nc_def_grp(ncid_out, "observation_data", &grp_obs)))  //netcdf-4
                        check_err(status, __LINE__, __FILE__);
                    //def var grp2
                        //counts and Lt
                    if ((status = nc_def_var(grp_obs, "obs_per_view", NC_INT, NDIMS3,
                        dimids3, &varid_count)))
                        check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_var(grp_obs, "l2_product", NC_FLOAT, NDIMS3,
                        dimids3, &varid_prod)))
                        check_err(status, __LINE__, __FILE__);

                    //leave define mode-----------------------
                    if ((status = nc_enddef(grp_obs))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);


                    lat_out = allocate2d_float(NY, NX);
                    lon_out = allocate2d_float(NY, NX);

                    int c = 0;

                    for (int i = 0; i < NY; i++) {
                        for (int j = 0; j < NX; j++) {
                            lat_out[NY - 1 - c][NX - 1 - j] = l1cfile->lat_gd[i][j];
                            lon_out[NY - 1 - c][j] = l1cfile->lon_gd[i][j];
                        }
                        c++;
                    }

                    if ((status = nc_put_var_float(grp_coor, varid1, &lat_out[0][0])))
                        check_err(status, __LINE__, __FILE__);

                    if ((status = nc_put_var_float(grp_coor, varid2, &lon_out[0][0])))
                        check_err(status, __LINE__, __FILE__);


           //alloc mem for output variables obs_per_view and I as part of the observation_data group
                    data_out = allocate4d_int(NY, NX,NVIEWS,NL2PRODS);
                    data_out2 = allocate4d_float(NY, NX,NVIEWS,NL2PRODS);

            //BINCOUNT ARRAY

                   c = 0;
                   for (size_t v = 0; v < NVIEWS; v++) {
                      for (size_t p = 0;p <NL2PRODS;p++) {
                        for (int i = 0; i < NY; i++) {
                            for (int j = 0; j < NX; j++) {
                                    if(bincount[v][p][i][j]>0)
                                         data_out[NY - 1 - c][NX - 1 - j][v][p]+=bincount[v][p][i][j];
                                    else
                                         data_out[NY - 1 - c][NX - 1 - j][v][p]+=0;
                                }
                            }
                            c++;
                        }
                        c=0;
                      }


                        cout << "writing countbin for all L2 products ."<< endl;
                       if ((status = nc_put_var_int(grp_obs, varid_count, &data_out[0][0][0][0])))
                            check_err(status, __LINE__, __FILE__);

            // l2 products ARRAYS
                    c=0;
                    for (size_t v = 0;v <NVIEWS;v++) {
                      for (size_t p = 0; p <NL2PRODS; p++) {
                        for (int i = 0; i < NY; i++) {
                             for (int j = 0; j < NX; j++) {
                                if (bincount[v][p][i][j] > 0)
                                    data_out2[NY - 1 - c][NX - 1 - j][v][p] = binmean_prod[v][p][i][j] / bincount[v][p][i][j];
                                else         data_out2[NY - 1 - c][NX - 1 - j][v][p] = NAN;
                            }
                            c++;
                        }
                        c=0;
                       }//end views
                     }//end for bands

                         cout << "writing L2 products .."<< endl;
                        if ((status = nc_put_var_float(grp_obs, varid_prod, &data_out2[0][0][0][0])))
                            check_err(status, __LINE__, __FILE__);

                    delete [] (data_out);
                    delete [] (data_out2);
                    delete[](lat_out);
                    delete[](lon_out);

                    if ((status = nc_close(ncid_out)))
                        check_err(status, __LINE__, __FILE__);

                    lat_out = nullptr;;
                    lon_out = nullptr;
                    data_out = nullptr;
                    data_out2 = nullptr;
    
                    cout << "total inpix.." << INPIX << "total pix.." << TOTPIX << endl;
                    printf("*** SUCCESS writing bincount and binLt rasters as nc..!\n");
                    cout<<"LINE-BY-LINE sbs binning method....."<<endl;
                }//end writing file


return 0;
}



int32_t L1C::swtime_swt3(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv){
  int16_t ix=-1,ngridlines,bina,binb,gd;
  double mot=0.,tini=0.,tend=0.,tg;
  int flag_time=-1;

  mot = ((l1cinput->grid_resolution) * 1000) /mgv;//in seconds
  ngridlines=l1cfile->num_gridlines;
  binb=l1cfile->binyb;
  bina=l1cfile->binya;

  cout<<"bin b.."<<binb<<"bina.."<<bina<<"tcross..."<<tcross<<endl;

  tini=tswt[0];
  ix=norbs-1;
  tend=tswt[ix];

  if(ngridlines>0 && ngridlines<=4000 && (l1cinput->grid_resolution)*1000==5200){
      flag_time=0;}
  else flag_time=1;

  if(ngridlines>0 && ngridlines<=8000 && (l1cinput->grid_resolution)*1000==5200/2){
      flag_time=0;}
  else flag_time=1;


   tg = tcross - mot / 2;
   gd=binb-1;

   tmgv[gd]=tg;
  while(gd>0){
      gd--;
      tg-=mot;

      tmgv[gd]=tg;
       if(tmgv[gd]<tini-mot){cout<<"error at swtime_swt3.....tmgv[gd]<tini-mot.."<<endl; exit(1);}
       }

   tg = tcross +  mot / 2;
   gd=binb;

   tmgv[gd]=tg;
  while(gd<(binb+bina)-1){
      gd++;
      tg+=mot;

      tmgv[gd]=tg;

       if(tmgv[gd]>tend+mot){cout<<"error at swtime_swt3.....tmgv[gd]>tend+mot.."<<endl; exit(1);}
       }

   return flag_time;
}




int32_t L1C::time_swt2(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput, double* ect_swtd, int16_t* swtd_id, int16_t* file_id, int16_t* nfiles_swt,  float* mgv_swath, double* time_mgv) {
        float mot, ect = -1;//mean time
        const char* ptstr;
        int dimid, status;
        std::string str;
        int16_t fi = 0;
        size_t num_scans_tot = 0, swt_orb_rec = 0, num_orb_rec = 0;
        double  tini, tend;//in seconds
        int32_t n_files, ffile = 1;//4000 gridlines at 5.2 km and 8000 at 2.6 km resolution

        nswath = l1cfile->nswath;
        num_scans = l1cfile->nscan;
        n_files = l1cfile->ifiles.size();

        mot = ((l1cinput->grid_resolution) * 1000) / mgv_swath[swtd - 1];//delta-t in seconds
//        cout << "mean swath vel.." << mgv_swath[swtd - 1] << "[mot..." << mot << "for day swath #.." << swtd << endl;

        cout<<"computing mean velocity time intervals......................................................................................for swath #...."<<swtd<<endl;

        for (int i = 0;i < n_files;i++) {
            if (swtd == swtd_id[i]) {
                //open files with the same crossing swath and compute time series based on mean time        
                fi = file_id[i];
                str = l1cfile->ifiles[fi - 1];
                ptstr = str.c_str();
                cout << "Opening file.." << ptstr << " for crossing swath..." << swtd << "file_id.." << file_id[i] << "swt_id.." << swtd_id[i] << endl;
                status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);

                if (status != NC_NOERR) {
                    fprintf(stderr, "-E- Error failed nc_open.\n");
                    exit(EXIT_FAILURE);
           }

                //Open dimensions
                // num_scans
                status = nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid);
                if (status != NC_NOERR) {
                    fprintf(stderr, "-E- Error reading number_of_scans.\n");
                    exit(EXIT_FAILURE);
                }
                nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

                num_scans_tot += num_scans;
                num_orb_rec = num_scans;
                swt_orb_rec += num_orb_rec;

    //            cout << "num_scans_tot..." << num_scans_tot << "num of orb records for swath..." << swt_orb_rec << "so far file.." << ptstr << "swath #.." << swtd << endl;

                //defining swath pointers
                double* orb_time = new double[num_scans]();//with Fred file num_scans = sdim

                //open scan_time (in Fred L1B file is orb_time interpolated)
                status = nc_inq_grp_ncid(ncid_L1B, "scan_line_attributes", &scGrp);
                check_err(status, __LINE__, __FILE__);
                status = nc_inq_varid(scGrp, "time", &otId);
                check_err(status, __LINE__, __FILE__);
                status = nc_get_var_double(scGrp, otId, orb_time);
                check_err(status, __LINE__, __FILE__);

                if (ffile == 1)
                    tini = orb_time[0];

                if (ffile != 1 & swtd != swtd_id[i + 1])//this is the last file of the swath
                    tend = orb_time[num_scans - 1];

                if (ect_swtd[i] > 0.0)//only 1 file is carrying the crossing time the other bunch of that swath is -999
                    ect = ect_swtd[i];
                //close file if not the last   
                status = nc_close(ncid_L1B);
                check_err(status, __LINE__, __FILE__);

                delete[](orb_time);
                ffile++;
            }//end if
         }//end for  

        if (ect > 0.0) {
  //          cout << "tini..." << tini << "tend..." << tend << "for swath..#.." << swtd << "mot..." << mot << "ect..." << ect << endl;
            //before eq crosssing time bins
  //         double  dtime = (ect - mot / 2) - tini;

            double tg = ect - mot / 2;
            int binb = 0, bina = 0;

            while (tg > (tini - mot / 2)) {
                binb++; //bins below eq crossing bin which is eq cros time +-mot/2
                tg -= mot;
            }

            tg = ect + mot / 2;

            while (tg < (tend + mot / 2)) {
                bina++;//time bins aove eq crossing time
                tg += mot;
            }

            //second way of computing gridlines---
//            dtime = tend - tini;


            l1cfile->num_gridlines = binb + bina;

            //before--  
            for (int bin = 0;bin < binb;bin++) //300 seconds are 5 minutes or 1 file from Fred
                time_mgv[bin] = tini + bin * mot;

            int binc = binb-1;
      //      time_mgv[binc] = ect;//equat crossing time bin

            //after
            for (int bin = 1;bin < bina;bin++) //300 seconds are 5 minutes or 1 file from Fred
                time_mgv[binc + bin] = time_mgv[binc] + bin * mot;

//        for (int bin = 0;bin <bina+binb-1;bin++) cout<<time_mgv[bin]<<endl;



        }//ect>0 
        return 0;
    }


int32_t L1C::swtime_swt2_segment(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv){
  int16_t bina=0,binb=0,ngridlines,gd=0;
  double tg=0.,mot=0,t_start=-1,t_end=-1;

  t_start=tswt[0];
  t_end=tswt[norbs-1];

  mot = ((l1cinput->grid_resolution) * 1000) /mgv;//in seconds
 // asc_mode=l1cfile->orb_dir;
  ngridlines=l1cfile->num_gridlines;

        int flag_time=-1;

    if (tcross > 0.0) {
          flag_time=0;
          cout<<"computing time series assuming mean swath velocity..for swath#."<<swt<<endl;

       //backward section of swath--before equat   
         tg = tcross - mot / 2;
         gd=ngridlines/2-1;
         tmgv[gd]=tg;
         binb=1;

         while(tg>=t_start){
             binb++;
             tg -= mot;
             gd--;
             tmgv[gd]=tg;
         }
      //forward section of swath--after equator   
         tg = tcross + mot / 2;
         gd=ngridlines/2;
         tmgv[gd]=tg;
         bina=1;

         while(tg<=t_end){
             bina++;
             tg += mot;
             gd++;
             tmgv[gd]=tg;
         }

         cout<<"bina.."<<bina<<"binb.."<<binb<<endl;

         l1cfile->num_gridlines = binb + bina;

            cout<<"number of L1C gridlines along-track..."<<l1cfile->num_gridlines<<"for swath #.."<<swt<<"t_start.."<<t_start<<"t_end..."<<t_end<<endl;
             }//end equat crossing
     else{
       cout<<"time series not possible for swath #.."<<swt<<"tcross<0...."<<endl;
       flag_time=1;
       }


     
   return flag_time;
}



int32_t L1C::swtime_swt2(int swt,L1C_input *l1cinput,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double tcross,double mgv,double *tmgv){
  int16_t bina=0,binb=0,ngridlines,gd=0;
  double tg=0.,mot=0;

  mot = ((l1cinput->grid_resolution) * 1000) /mgv;//in seconds
 // asc_mode=l1cfile->orb_dir;
  ngridlines=l1cfile->num_gridlines;

        int flag_time=-1;

    if (tcross > 0.0) {
          flag_time=0;
          cout<<"computing time series assuming mean swath velocity..for swath#."<<swt<<endl;

         tg = tcross - mot / 2;          
         gd=ngridlines/2-1;
         tmgv[gd]=tg;
         binb=1;

         while(binb<ngridlines/2){
             binb++;
             tg -= mot;
             gd--;
             tmgv[gd]=tg;
         }
         tg = tcross + mot / 2;
         gd=ngridlines/2;
         tmgv[gd]=tg;
         bina=1;

         while(bina<ngridlines/2+1){
             bina++; 
             tg += mot;
             gd++;
             tmgv[gd]=tg;        
         }

         cout<<"bina.."<<bina<<"binb.."<<binb<<endl;

         l1cfile->num_gridlines = binb + bina;

            cout<<"number of L1C gridlines along-track..."<<l1cfile->num_gridlines<<"for swath #.."<<swt<<endl;
             }//end equat crossing
     else{
       cout<<"time series not possible for swath #.."<<swt<<"tcross<0...."<<endl;
       flag_time=1;
       }


   return flag_time;
}


int32_t L1C::ect_swt(int swt,l1c_filehandle *l1cfile,int32_t norbs,double *tswt,double *latswt,double *lonswt,float*tcross,float*loncross){
    //determine equatorial crossing time--------------
        double t1=-1,t2=-1;
        int asc_mode=-1;
        size_t tindex=-1;
        int scan_brack[2] = { -1,-1 };
        float lat1, lat2, lon1, lon2;//lat,lon defined globally       
        int flag_ect=-1;

        //determine orbit direction--asc/desc
           if (latswt[0] > latswt[1])
                asc_mode = 0;//descending
            else
                asc_mode = 1;//ascending

        cout<<"computing equatorial crossing time and longitude at crossing for a specific swath...#..."<<swt<<"orbit direction 0: descending, 1: ascending..."<<asc_mode<<endl;
           for (int i = 0;i < norbs - 1;i++) { //improve with binary search
            if (asc_mode == 1 && latswt[i] < 0. && latswt[i + 1]>0.) {
                scan_brack[0] = i + 1;
                scan_brack[1] = i + 2;
                lat1 = latswt[i];
                lat2 = latswt[i + 1];
                lon1 = lonswt[i];
                lon2 = lonswt[i + 1];
                tindex=i;
                i=norbs;
                break;
              }

            else if (asc_mode == 0 && latswt[i] > 0. && latswt[i + 1] < 0) { //descending orbit
                scan_brack[0] = i + 2;//negative lat first convention
                scan_brack[1] = i + 1;
                lat1 = latswt[i + 1];
                lat2 = latswt[i];
                lon1 = lonswt[i + 1];
                lon2 = lonswt[i];
                tindex=i;
                i=norbs;
                break;
              }
             } //end for


         //interpolate time--linear
        if (scan_brack[0] > -1) {//if file with equat crossing
            t1 = tswt[tindex];
            t2 = tswt[tindex + 1];
            cout<<"lat1.."<<lat1<<"lat2..."<<lat2<<"t1.."<<t1<<"t2..."<<t2<<endl;
            *tcross = t1 - lat1 * (t2 - t1) / (lat2 - lat1);//equat crossing time

            float dtcross = (*tcross - t1) / (t2 - t1);
            *loncross = lon1 + (lon2 - lon1) * dtcross;//latcross is zero

            l1cfile->eqt =*tcross;
            l1cfile->orbit_node_lon = *loncross;
            l1cfile->orb_dir = asc_mode;
            flag_ect=0;
        }
        else {
            l1cfile->eqt = -999.0;
            l1cfile->orbit_node_lon = -999.0;
            l1cfile->orb_dir = asc_mode;
            flag_ect=1;
        }


    return flag_ect;
    }



//open telemetry file parameters for creating L1C grid----
   int32_t L1C::open_l1atol1c3(L1C_input *l1cinput,l1c_filehandle *l1cfile){
   int32_t n_files;
   const char* ptstr;
   char *ifile_char;
   string str,ifile_str;
   int status=-1,status1=-1,status2=-1;
   int telGrp,navGrp,scpId,otimeId,olatId,olonId;
   int dimhp,dimtp,dimor;
   unsigned char  **hkpackets=nullptr;
   uint8_t *apids=nullptr;
   size_t number_hkpack,number_scpack,number_orecords,number_hkpack_tot,number_scpack_tot,number_orecords_tot,nr;
   uint8_t ubnum1,ubnum2;
   double *sec=nullptr,*tai58_sec=nullptr;
   int16_t *year=nullptr,*mon=nullptr,*day=nullptr,*hour=nullptr,*min=nullptr;      
   int16_t oneyear,onemon,oneday,onehour,onemin;
   double onesec;
   double *orb_time=nullptr,*orb_lat=nullptr,*orb_lon=nullptr;
   float **lat_gd=nullptr,**lon_gd=nullptr,**alt=nullptr;
   double  omeg = 7.29211585494e-5;
   short n_ephem;
   string temp_str,tai_str;
   double vxyz=0,sum1=0,sum2=0,mov1=0.,mgv1=0.,mov2=0,mgv2,*orb_vel=nullptr,*grn_vel=nullptr;
   double sum3=0;
   int *orb_dir=nullptr;
   size_t nswath;
   int32_t  ix1=-1,ix2=-1,ix3=-1,ix4=-1,ix5=-1,ix6=-1;
   double *tmgv1=nullptr,*tmgv2=nullptr,*tmgvf=nullptr,*tmgvf2=nullptr;   
   int32_t num_gridlines=-1,norbs=-1;
   double *tswt=nullptr,*latswt=nullptr,*lonswt=nullptr;
   double *tswt2=nullptr,*latswt2=nullptr,*lonswt2=nullptr;
   const char* outlist;
   string ofile_str,senstr;        
   float latmin,latmax;
   string y_swt,mo_swt,d_swt,h_swt,mi_swt,s_swt,tswt_ini,tswt_end;
   int32_t gd_per_gran=-1,numgran=-1;
   double deltasec=-1;
   string  s_swt2;
   int length;
   string tswt_ini_file;
   double v1[3],v2[3],vecout[3],orbnorm[3],nvec,vi,toff;
   double rl2,pos_norm,clatg2,fe=1/298.257;
   double oangle,G[3],glat,glon,gnorm,rem=6371,omf2,omf2p,pxy,temp;
   int day_mode=-1;//0 is nighttime 1 is dayttime
   int day1,year1;
   float hour1,lat1,lon1,sunz1,suna1;
   double tai58unix;
   double tswt_ini_sec,tswt_end_sec;
   double clatg,rl;         
   vector<size_t> ix;
   double tfile_ini,tfile_end;
   size_t att_len;
   int16_t syear, smon, sday,shour,smin,syear2,smon2,sday2,shour2,smin2;
   double secs,second,secs2,second2;
   int logoff=-1;    
   double oangle_nad;
   double tfile_ini_sec,tfile_end_sec;
  
    
         //TOTAL ARRAYS---ASSUMING 6000 rcords
         number_orecords_tot=6000;
         size_t cc=0,c=0;

         int16_t *year_tot=(int16_t*)calloc(number_orecords_tot,sizeof(int16_t));
         int16_t *day_tot=(int16_t*)calloc(number_orecords_tot,sizeof(int16_t));
         int16_t *hour_tot=(int16_t*)calloc(number_orecords_tot,sizeof(int16_t));
         double *orb_time_tot=(double*)calloc(number_orecords_tot,sizeof(double));
         double *orb_lat_tot=(double*)calloc(number_orecords_tot,sizeof(double));
         double *orb_lon_tot=(double*)calloc(number_orecords_tot,sizeof(double));
         int *orb_dir_tot = (int*)calloc(number_orecords_tot,sizeof(int));
         orb_array2* posr_tot = new orb_array2[number_orecords_tot]();
         orb_array2* velr_tot = new orb_array2[number_orecords_tot]();
         double *orb_vel_tot = (double*)calloc(number_orecords_tot,sizeof(double));
         double *grn_vel_tot = (double*)calloc(number_orecords_tot,sizeof(double));


//*************************************************************

         ifile_str = l1cinput->files[0];
         ifile_char = &ifile_str[0];
         ofile_str=ifile_str.substr(0,24);
         outlist="list_l1c_granules.txt";//DEFAULT OFILE
  
         if(l1cinput->outlist[0]=='\0'){ //only when no txt output list file is provided
            strcpy(l1cinput->outlist,outlist);
            cout<<"L1C granules written to DEFAULT file........"<<outlist<<endl;
         }
         else {cout<<"L1C granules written to file......................................................................................."<<l1cinput->outlist<<endl;
         }

         l1cfile->gransize=l1cinput->gransize;
   
         file_format format = getFormat(ifile_char);
         cout<<"format.type.."<<format.type<<endl;
       

         l1cfile->sensor=l1cinput->sensor;//SPEX 1, OCI 2 and HARP 3

         if(l1cinput->sensor==34){
             senstr="SPEXone";
             l1cfile->nbinx=25;
             l1cfile->n_views=10;
              }
        else if (l1cinput->sensor==30){
            senstr="OCI";
            l1cfile->nbinx=519;
            l1cfile->n_views=2;
             }
        else if (l1cinput->sensor==35){
            senstr="HARP2";
            l1cfile->nbinx=457;
            l1cfile->n_views=90;
            }
        else{cout<<"sensor by default is OCI option 2....."<<endl;
             senstr="OCI";
             l1cfile->nbinx=514;
             l1cfile->n_views=2;                   
          }

         cout<<"PACE sensor to be griddded....."<<senstr<<endl;

         n_files = l1cfile->ifiles.size();
         cout<<"number of files in the list...."<<n_files<<endl;


        int32_t nfiles=l1cfile->ifiles.size(); 
        for (int fi = 0; fi <nfiles; fi++)
           {                           
            str = l1cfile->ifiles[fi];
            ptstr = str.c_str();
            cout<<"*********************************************************************************************************"<<endl;
            cout<<"Opening L1A file ..."<<ptstr<<"for telemetry......."<<endl;
            cout<<"*********************************************************************************************************"<<endl;

            status = nc_open(ptstr, NC_NOWRITE, &ncid_L1A);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- Error failed nc_open.\n");
                        exit(EXIT_FAILURE);
            }
   
   
            status = nc_inq_attlen(ncid_L1A, NC_GLOBAL, "time_coverage_start", &att_len);
            check_err(status, __LINE__, __FILE__);
            // allocate required space before retrieving values 
            char* time_str = (char*)malloc(att_len + 1); // + 1 for trailing null 
         // get global attribute values 
            status = nc_get_att_text(ncid_L1A, NC_GLOBAL, "time_coverage_start", time_str);
            check_err(status, __LINE__, __FILE__);
            time_str[att_len] = '\0';
            tfile_ini = isodate2unix(time_str);

            status = nc_inq_attlen(ncid_L1A, NC_GLOBAL, "time_coverage_end", &att_len);
            check_err(status, __LINE__, __FILE__);
            // allocate required space before retrieving values
            char* time_str2 = (char*)malloc(att_len + 1); // + 1 for trailing null
         // get global attribute values
            status = nc_get_att_text(ncid_L1A, NC_GLOBAL, "time_coverage_end", time_str2);
            check_err(status, __LINE__, __FILE__);
            time_str2[att_len] = '\0';           
            tfile_end = isodate2unix(time_str2);
            
            //this is default start time, start_timeflag=1 or using start time of swath file HKT------
            unix2ymds(tfile_ini, &syear, &smon, &sday,&secs);
            unix2ymds(tfile_end, &syear2, &smon2, &sday2,&secs2);

            cout<<"secs elapsed.."<<secs<<"initial granule #..."<<round(secs/(l1cfile->gransize*60))<<endl;
            unix2ymdhms(tfile_ini, &syear,&smon,&sday,&shour,&smin,&second);
            cout<<"HKT file start time................."<<"year.."<<syear<<"month..."<<smon<<"day..."<<sday<<"hour.."<<shour<<"min.."<<smin<<"sec..."<<second<<endl;
            unix2ymdhms(tfile_end, &syear2,&smon2,&sday2,&shour2,&smin2,&second2);
            cout<<"HKT file end time................."<<"year.."<<syear2<<"month..."<<smon2<<"day..."<<sday2<<"hour.."<<shour2<<"min.."<<smin2<<"sec..."<<second2<<endl;

            if(l1cinput->start_timeflag==0){   //forcing to hms =0:0:0 
               secs=0.;               
            }
               
            tfile_ini_sec=ymds2unix(syear,smon,sday,secs);
            tfile_end_sec=tfile_ini_sec+49*60*3;//49 minutes per half-orbit, 3x cause 1.5 orbit, three swaths

            int16_t y_ini,mo_ini,d_ini,h_ini,mi_ini;
            double sec_ini; 
            unix2ymdhms(tfile_ini_sec,&y_ini,&mo_ini,&d_ini, &h_ini, &mi_ini, &sec_ini);
            cout<<"tfile_ini_sec.."<<"YEAR.."<<y_ini<<"MONTH.."<<mo_ini<<"DAY.."<<d_ini<<"HOUR.."<<h_ini<<"MIN.."<<mi_ini<<"SEC.."<<sec_ini<<endl;
            unix2ymdhms(tfile_end_sec,&y_ini,&mo_ini,&d_ini, &h_ini, &mi_ini, &sec_ini);
            cout<<"tfile_end_sec.."<<"YEAR.."<<y_ini<<"MONTH.."<<mo_ini<<"DAY.."<<d_ini<<"HOUR.."<<h_ini<<"MIN.."<<mi_ini<<"SEC.."<<sec_ini<<endl;

            l1cfile->tfile_ini_sec=tfile_ini_sec;
            l1cfile->tfile_end_sec=tfile_end_sec;
            
       //open dimensions
          ///number of packets
            status = nc_inq_dimid(ncid_L1A, "SC_hkt_pkts", &dimhp);//number of spacecraft packets
             check_err(status, __LINE__, __FILE__);
            status = nc_inq_dimid(ncid_L1A, "max_SC_packet", &dimtp);//size of packets
             check_err(status, __LINE__, __FILE__);
            status = nc_inq_dimid(ncid_L1A, "orb_records", &dimor);//number of orbital records for SC track lat/lon
            check_err(status, __LINE__, __FILE__);

            nc_inq_dimlen(ncid_L1A, dimhp, &number_hkpack);
            nc_inq_dimlen(ncid_L1A, dimtp, &number_scpack);
            nc_inq_dimlen(ncid_L1A, dimor, &number_orecords);

            //total counters---
            number_hkpack_tot+=number_hkpack;
            number_scpack_tot+=number_scpack;
            nr+=number_orecords;

            cout<<"number_hkpack_tot.."<<number_hkpack_tot<<"number_scpack_tot.."<<number_scpack_tot<<"number of orbit records total.REAL."<<nr<<endl;

  //allocat mem
            hkpackets=allocate2d_uchar(number_hkpack, number_scpack);//NUMBER of ephem elements x max SIZE of package
            apids=(uint8_t*)calloc(number_hkpack,sizeof(uint8_t));

            orb_time=(double*)calloc(number_orecords,sizeof(double));
            orb_lat=(double*)calloc(number_orecords,sizeof(double));
            orb_lon=(double*)calloc(number_orecords,sizeof(double));
            orb_dir = (int*)calloc(number_orecords,sizeof(int));

   

        //open groups   
         status = nc_inq_grp_ncid(ncid_L1A, "housekeeping_data", &telGrp);
                check_err(status, __LINE__, __FILE__);   
         status = nc_inq_grp_ncid(ncid_L1A, "navigation_data", &navGrp);
                check_err(status, __LINE__, __FILE__);
          
         //open vars ids
         status = nc_inq_varid(telGrp, "SC_HKT_packets", &scpId);
                check_err(status, __LINE__, __FILE__);
         status = nc_inq_varid(navGrp, "orb_time", &otimeId);//second of the day
                check_err(status, __LINE__, __FILE__);
         status = nc_inq_varid(navGrp, "orb_lat", &olatId);
                check_err(status, __LINE__, __FILE__);
         status = nc_inq_varid(navGrp, "orb_lon", &olonId);
                check_err(status, __LINE__, __FILE__);

         //get data       
         status = nc_get_var_ubyte(telGrp, scpId, &hkpackets[0][0]);
            check_err(status, __LINE__, __FILE__);      


            status = nc_get_var_double(navGrp, otimeId, &orb_time[0]);
             check_err(status, __LINE__, __FILE__);
            status = nc_get_var_double(navGrp, olatId, &orb_lat[0]);
              check_err(status, __LINE__, __FILE__);
            status = nc_get_var_double(navGrp, olonId, &orb_lon[0]);
             check_err(status, __LINE__, __FILE__);


         for(size_t hk=0;hk<number_hkpack;hk++){
            ubnum1=(uint8_t)hkpackets[hk][0]; //-48;//48 is 0 ascii code
            ubnum2=(uint8_t)hkpackets[hk][1];

            apids[hk] = (ubnum1%8)*256 + ubnum2;
            if(apids[hk]==128){
                ix.push_back(hk);//packet index where ephemr               
               } 
               }

         cout<<"#number of ephem elements...."<<ix.size()<<"for HKT file..."<<ptstr<<"#orecords..."<<number_orecords<<endl;
         n_ephem=ix.size();
         
         // Reverse the byte order  from small endian to large endian, small endian last byte is stored first
         // convert to double and later swap endian

      //allocate mem
            orb_vel = (double*)calloc(n_ephem,sizeof(double));
            grn_vel = (double*)calloc(n_ephem,sizeof(double));

         sec=(double*)calloc(n_ephem,sizeof(double));
         year=(int16_t*)calloc(n_ephem,sizeof(int16_t));
         mon=(int16_t*)calloc(n_ephem,sizeof(int16_t));
         day=(int16_t*)calloc(n_ephem,sizeof(int16_t));
         hour=(int16_t*)calloc(n_ephem,sizeof(int16_t));
         min=(int16_t*)calloc(n_ephem,sizeof(int16_t));

            tai58_sec=(double*)calloc(n_ephem,sizeof(double));


     //process all packets----
         double tai58;
         c=0;
         //TIME   

            while (c<ix.size()){
               double *tai_ptr = (double*) (hkpackets[ix[c]] + 16);//moving pointer to position 16
               swapc_bytes((char*) tai_ptr, 8, 1);//swap 8 bytes once or 1 16-23 pos, lets say we have 16-32 indexes and we want to swap 8 bytes at the time, some we need ntime=2
               tai58 = *tai_ptr;
               tai58unix=tai58_to_unix(tai58);

               unix2ymdhms(tai58unix,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);

               year[c]=oneyear;
               mon[c]=onemon;
               day[c]=oneday;
               hour[c]=onehour;
               min[c]=onemin;
               sec[c]=onesec;

               tai58_sec[c]=tai58unix;

              c++;               
                }


        //POS/VEL alloc mem-----
              orb_array2* posr = new orb_array2[n_ephem]();
              orb_array2* velr = new orb_array2[n_ephem]();

         double dp1,dp2,dp3;

            //computed orbital velocity from ephemeris
         c=0;
            while (c<ix.size()){
               double *posi = (double*) (hkpackets[ix[c]] + 120);
               double *veli = (double*) (hkpackets[ix[c]] + 144);
               double *ecmat = (double*) (hkpackets[ix[c]] + 176);

             swapc_bytes((char*) posi, 8, 3);
             swapc_bytes((char*) veli, 8, 3);
             swapc_bytes((char*) ecmat, 8, 9);
    
       //assuming ecmat index 1-3 first row, 4-6 second row and 7-9 third row
       //dot product
            //position          
             dp1=posi[0]*ecmat[0]+posi[1]*ecmat[1]+posi[2]*ecmat[2];
             dp2=posi[0]*ecmat[3]+posi[1]*ecmat[4]+posi[2]*ecmat[5];
             dp3=posi[0]*ecmat[6]+posi[1]*ecmat[7]+posi[2]*ecmat[8];

             posr[c][0]=dp1;
             posr[c][1]=dp2;
             posr[c][2]=dp3;
       
       //velocity
             dp1=veli[0]*ecmat[0]+veli[1]*ecmat[1]+veli[2]*ecmat[2];
             dp2=veli[0]*ecmat[3]+veli[1]*ecmat[4]+veli[2]*ecmat[5];
             dp3=veli[0]*ecmat[6]+veli[1]*ecmat[7]+veli[2]*ecmat[8];

             velr[c][0]=dp1;
             velr[c][1]=dp2;
             velr[c][2]=dp3;

             velr[c][0]=velr[c][0] + posr[c][1]*omeg;
             velr[c][1]=velr[c][1] - posr[c][0]*omeg;

         //computing orbital velocity
               vxyz = sqrt(velr[c][0] * velr[c][0] + velr[c][1] * velr[c][1] + velr[c][2] * velr[c][2]);//units m/s
               orb_vel[c]=vxyz;               

               pos_norm=sqrt(posr[c][0]*posr[c][0]+posr[c][1]*posr[c][1]+posr[c][2]*posr[c][2]);
               clatg2=sqrt(posr[c][0]*posr[c][0]+posr[c][1]*posr[c][1])/pos_norm;
               rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));
               grn_vel[c]=vxyz*rl2/pos_norm;

         c++;    
         }


         for(short n=0;n<n_ephem;n++){
             orb_time_tot[cc]=orb_time[n];
             orb_lon_tot[cc]=orb_lon[n];
             orb_lat_tot[cc]=orb_lat[n];
             orb_vel_tot[cc]=orb_vel[n];
             grn_vel_tot[cc]=grn_vel[n];
             velr_tot[cc][0]=velr[n][0];
             velr_tot[cc][1]=velr[n][1];
             velr_tot[cc][2]=velr[n][2];
             posr_tot[cc][0]=posr[n][0];
             posr_tot[cc][1]=posr[n][1];
             posr_tot[cc][2]=posr[n][2];
             year_tot[cc]=year[n];
             day_tot[cc]=day[n];
             hour_tot[cc]=hour[n];
             cc++;
         }


         if (hkpackets != nullptr)
               delete[](hkpackets);
         hkpackets= nullptr;

         if (apids != nullptr)
               delete[](apids);
         apids= nullptr;

         ix.clear();

         if (tai58_sec != nullptr)
               delete[](tai58_sec);
         tai58_sec= nullptr;

         if (sec != nullptr)
               delete[](sec);
         sec= nullptr;

         if (year != nullptr)
               delete[](year);
         year= nullptr;

         if (mon != nullptr)
               delete[](mon);
         mon= nullptr;

         if (day != nullptr)
               delete[](day);
         day= nullptr;

         if (hour != nullptr)
               delete[](hour);
         hour= nullptr;

         if (min != nullptr)
               delete[](min);
         min= nullptr;

         if (orb_lat != nullptr)
               delete[](orb_lat);
         orb_lat= nullptr;

         if (orb_lon != nullptr)
               delete[](orb_lon);
         orb_lon= nullptr;

         if (orb_time != nullptr)
               delete[](orb_time);
         orb_time= nullptr;

         if (orb_vel != nullptr)
               delete[](orb_vel);
         orb_vel= nullptr;

         if (grn_vel != nullptr)
               delete[](grn_vel);
          grn_vel= nullptr;

         if (orb_dir != nullptr)
               delete[](orb_dir);
         orb_dir= nullptr;

         if (posr != nullptr)
               delete[](posr);
            posr= nullptr;

         if (velr != nullptr)
               delete[](velr);
            velr= nullptr;
         if (senz != nullptr)
               delete[](senz);
         senz= nullptr;

         if (sena != nullptr)
               delete[](sena);
         sena= nullptr;

         if ((status = nc_close(ncid_L1A)))
            check_err(status, __LINE__, __FILE__);

     }//end big loop files HKT batch processing


      
//determine ini/end time indexes for asc/desc passes
             for(size_t k=0;k<number_orecords_tot-1;k++){
                 if(orb_lat_tot[k+1]>orb_lat_tot[k])
                     orb_dir_tot[k]=1;//ascending
                 else orb_dir_tot[k]=0;//descending       
                    }

//find time limits for half-orbits
            int swt=1;
            int16_t tlimits_swt[2][2]={{-1,-1},{-1,-1}};//2 swath x time limits ini/end
            tlimits_swt[swt-1][0]=0;//only the record indexes


            for(size_t k=0;k<number_orecords_tot-1;k++){
              if(orb_dir_tot[k]!=orb_dir_tot[k+1]){
                if(swt==1){           
                  tlimits_swt[swt-1][1]=k;//this time is a 'local' end cause may continue circling the earth after ascending or descending
                  swt++;
                  tlimits_swt[swt-1][0]=k+1;
                    }         
                else if(swt==2){//another partial swath
                  tlimits_swt[swt-1][1]=k;
                  tlimits_swt[swt-2][1]=number_orecords_tot-1;    
                  swt++;
                  }                       
                 }
         
                if(k==number_orecords_tot-2 && swt==2)  tlimits_swt[swt-1][1]=number_orecords_tot-1;
                 }     

            nswath=swt;

       cout<<"nswath...................."<<nswath<<"number_orecords_tot.."<<number_orecords_tot<<endl;


//split half-orbits----  group lat/lon time and velocity vector
            for (size_t sw=0;sw<2;sw++){
              if(nswath==2){
                 if(sw==0){ 
                    ix1=tlimits_swt[sw][0];
                    ix2=tlimits_swt[sw][1];
                    cout<<"nswath#2......swat#1...ix1"<<ix1<<"ix2.."<<ix2<<endl;
                   }
                 else if(sw==1){
                    ix3=tlimits_swt[sw][0];
                    ix4=tlimits_swt[sw][1];
                    cout<<"nswath2.....swat#2...ix3"<<ix3<<"ix4.."<<ix4<<endl;
                    }   
                 }
              else if(nswath==3){ //>2 we have 4 limits for swath 1
                 if(sw==0){ 
                     ix1=tlimits_swt[sw][0];
                     ix2=tlimits_swt[sw+1][0]-1;
                     ix3=tlimits_swt[sw+1][1]+1;
                     ix4=tlimits_swt[sw][1];
                      cout<<"nswath3....swat#1...ix1"<<ix1<<"ix2.."<<ix2<<"ix3.."<<ix3<<"ix4..."<<ix4<<endl;
                     }
                 else if(sw==1){
                    ix5=tlimits_swt[sw][0];
                    ix6=tlimits_swt[sw][1];
                    cout<<"nswath3....swat#2...ix5"<<ix5<<"ix6.."<<ix6<<endl;
                   }               
                  }
              else{
                cout<<"number of swaths is less than 2..!! exit..."<<endl;
                 exit(1);            
                 }
              }

          int kk;
          float tcross1=-999.,loncross1=-999.,tcross2=-999.,loncross2=-999.;


        //call fortran blocks  
          cdata_();

         
//--------------------------------------------------------------------------------------          
//nswath#2
//----case type: half orbits are consecutive          
//swath1
//--------------------------------------------------------------------------------------
          if(nswath==2){
            if(ix1>=0 && ix5<0){
                norbs=ix2-ix1+1;
                kk=ix1;       
                c=0;
                latmin=orb_lat_tot[ix1];
                latmax=orb_lat_tot[ix2];

 //compute solar zenith angle for checking day/night conditions for each swath ---
   //-------------------------------------------------------
                year1=year_tot[ix1];
                day1=day_tot[ix1];
                hour1=hour_tot[ix1];
                lat1=orb_lat_tot[ix1];
                lon1=orb_lon_tot[ix1];
                    
                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                 if(sunz1<=93) day_mode=1; else day_mode=0;//nighttime

                year1=year_tot[ix2];
                day1=day_tot[ix2];
                hour1=hour_tot[ix2];
                lat1=orb_lat_tot[ix2];
                lon1=orb_lon_tot[ix2];

                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                 if(sunz1<=93 && day_mode==1) day_mode=1; else day_mode=0;


             if(day_mode==1){      
//orbit direction 
                l1cfile->orb_dir=orb_dir_tot[ix1];

//ini and end UTC--this is swath orbital time, not interpolated time
                tswt=(double*)calloc(norbs,sizeof(double));
                latswt=(double*)calloc(norbs,sizeof(double));
                lonswt=(double*)calloc(norbs,sizeof(double));

                sum1=0,sum3=0,c=0;
                while(kk<=ix2){
                   latswt[c]=orb_lat_tot[kk];
                   lonswt[c]=orb_lon_tot[kk];
                   tswt[c]=orb_time_tot[kk];

                   if(orb_lat_tot[kk]<=latmin){ 
                       latmin=orb_lat_tot[kk];                          
                        }
                   if(orb_lat_tot[kk]>=latmax){ 
                       latmax=orb_lat_tot[kk];
                        }   

                   sum1+=orb_vel_tot[kk];
                   sum3+=grn_vel_tot[kk];
                   kk++;
                   c++;
                     }
                //computing mean swath velocity (m/s) --------------------------
                mov1 = (sum1 /norbs);//mean velocity in m/s 
                mgv1 = (sum3 /norbs)*1000;//in meters/s

                cout<<"mov1..."<<mov1<<"mvg1..."<<mgv1<<endl;
                
                status=ect_swt(1,l1cfile,norbs,tswt,latswt,lonswt,&tcross1,&loncross1);
                cout<<"nswath==2 --tcross equat crossing in (s)..swath#1..."<<tcross1<<"loncross1..."<<loncross1<<endl;
                cout<<"latmin.."<<latmin<<"latmax..."<<latmax<<endl;


                if (latswt != nullptr)
                     delete[](latswt);
                latswt= nullptr;
                if (lonswt != nullptr)
                     delete[](lonswt);
                lonswt= nullptr;
                
                l1cfile->num_gridlines=4000;                
                if(status==0){
                   tmgv1=(double*)calloc(l1cfile->num_gridlines,sizeof(double));
                   tmgvf=(double*)calloc(l1cfile->num_gridlines,sizeof(double));

                   swtime_swt2(1,l1cinput,l1cfile,norbs,tswt,tcross1,mgv1,tmgv1);                                                        
                   //ini/end times for swath
                   num_gridlines = l1cfile->num_gridlines;
                                                           
                   cout<<"number of across bins L1C grid...#.."<<l1cfile->nbinx<<"l1cfile->n_views..."<<l1cfile->n_views<<endl;

                   lat_gd = allocate2d_float(4000, l1cfile->nbinx);
                   lon_gd = allocate2d_float(4000, l1cfile->nbinx);
                   alt=allocate2d_float(4000, l1cfile->nbinx);             
                   

                   orb_array2* posgrid = new orb_array2[num_gridlines]();//these are doubles
                   orb_array2* velgrid = new orb_array2[num_gridlines]();

                   orb_array2* posgrid2 = new orb_array2[num_gridlines]();//these are doubles
                   orb_array2* velgrid2 = new orb_array2[num_gridlines]();

//interp 1-- always based on orbitals---                   
     //interpola pos/veloc as a function of time                  
                   orb_interp2(number_orecords_tot, num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgv1, posgrid,velgrid);
                 
                   for(int i=0;i<num_gridlines-1;i++){
                          pos_norm=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1]+posgrid[i][2]*posgrid[i][2]);
                          clatg2=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1])/pos_norm;
                          rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                          v2[0]=velgrid[i][0]*rl2/pos_norm;
                          v2[1]=velgrid[i][1]*rl2/pos_norm;
                          v2[2]=velgrid[i][2]*rl2/pos_norm;

                          vi=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])*1000;
                          toff=vi/mgv1;    
                          tmgvf[i]=tmgv1[i]+toff;
                          }

                   orb_interp2(number_orecords_tot, num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgvf, posgrid2,velgrid2);


                  //angle subsat track----
                  omf2=(1-fe)*(1-fe);

                  for(int i=0;i<num_gridlines-1;i++){ 
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
                        nvec=cross_product_norm_double(v1,v2);//length of orb norm vect

                        orbnorm[0]=vecout[0]/nvec;
                        orbnorm[1]=vecout[1]/nvec;
                        orbnorm[2]=vecout[2]/nvec;

                   for(int j=0;j<l1cfile->nbinx;j++){
                      pos_norm=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);//first gridline 
                      oangle=asin((j-(l1cfile->nbinx-1)/2)*5.2/pos_norm);


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
                                                                      
                      //altitude
                      clatg = cos(atan(omf2*tan(glat*M_PI/180.)));
                      rl = Re*(1.-fe)/sqrt(1.-(2.-fe)*fe*clatg*clatg);
                      alt[i][j] = (gnorm - rl)*1000;//in meters
                        
                   }
                  }
                  
                   if (posgrid != nullptr)
                     delete[](posgrid);
               
                   if (velgrid != nullptr)
                     delete[](velgrid);  
                   
                   if (posgrid2 != nullptr)
                     delete[](posgrid2);

                   if (velgrid2 != nullptr)
                     delete[](velgrid2);           

                  l1cfile->num_gridlines = l1cfile->num_gridlines-1;

                  if(l1cinput->grantype==1){                               
                        tswt_ini_sec=tfile_ini_sec+tmgvf[0];
                        unix2ymdhms(tswt_ini_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);                    
                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_ini=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        tswt_ini_file=y_swt+mo_swt+d_swt+"T"+h_swt+mi_swt+s_swt2.substr(0,2);
                        
                        tswt_end_sec=tfile_ini_sec+tmgvf[num_gridlines-2];
                        unix2ymdhms(tswt_end_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);

                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_end=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        l1cfile->tswt_ini=tswt_ini;
                        l1cfile->tswt_end=tswt_end;
                        l1cfile->tswt_ini_file=tswt_ini_file;     
                      create_SOCEA2(1,l1cinput,l1cfile,lat_gd,lon_gd,alt,tmgvf);//THIS IS SWATH PROCESSING
                     }
                  else if(l1cinput->grantype==0)
                  {
//--------------------------------------------------------------
//granule processing---------------------------------------------                 
//-----------------------------------------------------------
                   deltasec=tmgvf[num_gridlines-2]-tmgvf[0]+1;
                   cout<<"deltasec..swath."<<deltasec<<endl;

                   if(tswt != nullptr)
                     delete[](tswt);
                   tswt= nullptr;

                   numgran=144*2;
                   l1cfile->numgran=numgran;
                   gd_per_gran=round(num_gridlines/10);//10 granules per half orbit
                   l1cfile->gd_per_gran=gd_per_gran;

                   cout<<"estimated # of granules to be processed..."<<numgran<<"gd_per_gran..."<<gd_per_gran<<"#gridlines.."<<num_gridlines<<endl;

                   write_L1C_granule2(1,l1cfile,l1cinput,tmgvf,lat_gd, lon_gd,alt);
                   
                
                  }
                  else{
                   cout<<"ERROR selecting grantype, must be 0: granules or 1: swath........................."<<endl;
                   exit(1);
                  }
//-----------------------------------------------------------------
//-----------------------------------------------------------------   

                   if (lat_gd != nullptr)
                      delete [](lat_gd);
                   lat_gd = nullptr;
                   if (lon_gd != nullptr)
                       delete [](lon_gd);
                   lon_gd = nullptr;
                               
                   if (tmgv1 != nullptr)
                     delete[](tmgv1);
                   tmgv1= nullptr;
                   if (tmgvf != nullptr)
                     delete[](tmgvf);
                   tmgvf= nullptr;
                   if (alt != nullptr)
                     delete[](alt);
                   alt= nullptr;
                  }

               } //end day_mode
              else{
                   cout<<"ERROR swath #1 day_mode = "<<day_mode<<"nightime...continue to swath2 (nswath#2).............."<<endl;
                 }


//nswath=2                 
//-------------------------------------------------------------------------------------------------------         
//swath#2    
//-------------------------------------------------------------------------------------------------------    
//              if(l1cinput->swath_num==2){

                kk=ix3;
                c=0;
                norbs=ix4-ix3+1;
                latmin=orb_lat_tot[ix3];
                latmax=orb_lat_tot[ix4];

               //compute solar zenith angle for checking day/night conditions for each swath ---
   //-------------------------------------------------------
                year1=year_tot[ix3];
                day1=day_tot[ix3];
                hour1=hour_tot[ix3];
                lat1=orb_lat_tot[ix3];
                lon1=orb_lon_tot[ix3];

                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                 if(sunz1<=93) day_mode=1; else day_mode=0;

                year1=year_tot[ix4];
                day1=day_tot[ix4];
                hour1=hour_tot[ix4];
                lat1=orb_lat_tot[ix4];
                lon1=orb_lon_tot[ix4];

                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                 if(sunz1<=93 || day_mode==1) day_mode=1; else day_mode=0;

              if(day_mode==1){
                //orbit direction
                l1cfile->orb_dir=orb_dir_tot[ix3];

                            tswt2=(double*)calloc(norbs,sizeof(double));
                latswt2=(double*)calloc(norbs,sizeof(double));
                lonswt2=(double*)calloc(norbs,sizeof(double));

                sum2=0,sum3=0,c=0;
                while(kk<=ix4){
                   latswt2[c]=orb_lat_tot[kk];
                   lonswt2[c]=orb_lon_tot[kk];
                   tswt2[c]=orb_time_tot[kk];

                   if(orb_lat_tot[kk]<=latmin){
                       latmin=orb_lat_tot[kk];
                        }
                   if(orb_lat_tot[kk]>=latmax){
                       latmax=orb_lat_tot[kk];
                        }

                   sum2+=orb_vel_tot[kk];
                   sum3+=grn_vel_tot[kk];
                   kk++;
                   c++;
                     }
                mov2 = (sum2 /norbs);//mean velocity in m/s x 10^3
                mgv2 = (sum3 /norbs)*1000;//in meters/s  

                cout<<"mov2..."<<mov2<<"mvg2..."<<mgv2<<endl;

                status=ect_swt(2,l1cfile,norbs,tswt2,latswt2,lonswt2,&tcross2,&loncross2);
                cout<<"nswath==2 ---tcross equat crossing in (s)..swath#2..."<<tcross2<<"loncross2..."<<loncross2<<endl;
                cout<<"latmin.."<<latmin<<"latmax..."<<latmax<<endl;

                if (latswt2 != nullptr)
                     delete[](latswt2);
                latswt2= nullptr;
                if (lonswt2 != nullptr)
                     delete[](lonswt2);
                lonswt2= nullptr;               

                l1cfile->num_gridlines=4000;

                if(status==0){
                  tmgv2=(double*)calloc(l1cfile->num_gridlines,sizeof(double)); 
                  tmgvf2=(double*)calloc(l1cfile->num_gridlines,sizeof(double));

                  swtime_swt2(2,l1cinput,l1cfile,norbs,tswt2,tcross2,mgv2,tmgv2);                 
                  num_gridlines = l1cfile->num_gridlines;

//granule processing ------------------------                  
   
                  numgran=144;         
                  l1cfile->numgran=numgran;
                  gd_per_gran=round(num_gridlines/10);
                  l1cfile->gd_per_gran=gd_per_gran;

                  cout<<"# of granules to be processed..."<<numgran<<"gd_per_gran..."<<gd_per_gran<<endl;

                  if(tswt2 != nullptr)
                     delete[](tswt2);
                  tswt2= nullptr;

                  cout<<"number of across bins L1C grid...#.."<<l1cfile->nbinx<<endl;
                  lat_gd = allocate2d_float(4000, l1cfile->nbinx);
                  lon_gd = allocate2d_float(4000, l1cfile->nbinx);
                  alt=allocate2d_float(4000, l1cfile->nbinx);

                  orb_array2* posgrid = new orb_array2[num_gridlines]();//these are doubles
                  orb_array2* velgrid = new orb_array2[num_gridlines]();
                  orb_array2* posgrid2 = new orb_array2[num_gridlines]();//these are doubles
                  orb_array2* velgrid2 = new orb_array2[num_gridlines]();

                  orb_interp2(number_orecords_tot, num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgv2, posgrid,velgrid);
                 
                  for(int i=0;i<num_gridlines-1;i++){
                          pos_norm=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1]+posgrid[i][2]*posgrid[i][2]);
                          clatg2=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1])/pos_norm;
                          rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                          v2[0]=velgrid[i][0]*rl2/pos_norm;
                          v2[1]=velgrid[i][1]*rl2/pos_norm;
                          v2[2]=velgrid[i][2]*rl2/pos_norm;

                          vi=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])*1000;
                          toff=vi/mgv2;    
                          tmgvf2[i]=tmgv2[i]+toff;
                            }


                   orb_interp2(number_orecords_tot, num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgvf2, posgrid2,velgrid2);

                  //angle subsat track----
                  omf2=(1-fe)*(1-fe);

                  for(int i=0;i<num_gridlines-1;i++){ 
                        pos_norm=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1]+posgrid2[i][2]*posgrid2[i][2]);
                        clatg2=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1])/pos_norm;
                        rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                        v1[0]=(posgrid2[i][0])*rl2/pos_norm;
                        v1[1]=(posgrid2[i][1])*rl2/pos_norm;
                        v1[2]=(posgrid2[i][2])*rl2/pos_norm;

                        v2[0]=(velgrid2[i][0])*rl2/pos_norm;
                        v2[1]=(velgrid2[i][1])*rl2/pos_norm;
                        v2[2]=(velgrid2[i][2])*rl2/pos_norm;         

                        cross_product_double(v1,v2,vecout);
                        nvec=cross_product_norm_double(v1,v2);//length of orb norm vect

                        orbnorm[0]=vecout[0]/nvec;
                        orbnorm[1]=vecout[1]/nvec;
                        orbnorm[2]=vecout[2]/nvec;

                   for(int j=0;j<l1cfile->nbinx;j++){

                      pos_norm=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);//first gridline 

                      oangle=asin((j-(l1cfile->nbinx-1)/2)*5.2/pos_norm);

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
                                            
                      //altitude
                      clatg = cos(atan(omf2*tan(glat*M_PI/180.)));
                      rl = Re*(1.-fe)/sqrt(1.-(2.-fe)*fe*clatg*clatg);
                      alt[i][j] = (gnorm - rl)*1000;//in meters           
                   }
                  }
                  
                   if (posgrid != nullptr)
                     delete[](posgrid);
               
                   if (velgrid != nullptr)
                     delete[](velgrid);  
                   
                   if (posgrid2 != nullptr)
                     delete[](posgrid2);

                   if (velgrid2 != nullptr)
                     delete[](velgrid2);           

                  l1cfile->num_gridlines = l1cfile->num_gridlines-1;

                  if(l1cinput->grantype==1){
                        tswt_ini_sec=tfile_ini_sec+tmgvf2[0];
                        unix2ymdhms(tswt_ini_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);
                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_ini=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        tswt_ini_file=y_swt+mo_swt+d_swt+"T"+h_swt+mi_swt+s_swt2.substr(0,2);

                        tswt_end_sec=tfile_ini_sec+tmgvf2[num_gridlines-2];
                        unix2ymdhms(tswt_end_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);

                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_end=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        l1cfile->tswt_ini=tswt_ini;
                        l1cfile->tswt_end=tswt_end;
                        l1cfile->tswt_ini_file=tswt_ini_file;
                     create_SOCEA2(2,l1cinput,l1cfile,lat_gd,lon_gd,alt,tmgvf2);//THIS IS SWATH PROCESSING
                     }
                  else if(l1cinput->grantype==0)
                    {
//--------------------------------------------------------------
//granule processing---------------------------------------------                 
//---------------------------------------------------------------
                   deltasec=tmgvf2[num_gridlines-2]-tmgvf2[0]+1;
                   cout<<"deltasec..swath."<<deltasec<<endl;

                   numgran=144*2;
                   l1cfile->numgran=numgran;
                   gd_per_gran=round(num_gridlines/10);//10 granules per half orbit
                   l1cfile->gd_per_gran=gd_per_gran;

                   cout<<"estimated # of granules to be processed..."<<numgran<<"gd_per_gran..."<<gd_per_gran<<"#gridlines.."<<num_gridlines<<endl;

                   write_L1C_granule2(2,l1cfile,l1cinput,tmgvf2,lat_gd, lon_gd,alt);
                    }
                  else{
                   cout<<"ERROR selecting grantype, must be 0: granules or 1: swath........................."<<endl;
                   exit(1);
                  }
//--------------------------------------------------------------------
//--------------------------------------------------------------------
                  if (lat_gd != nullptr)
                    delete [](lat_gd);
                  lat_gd = nullptr;
                  if (lon_gd != nullptr)
                       delete [](lon_gd);
                  lon_gd = nullptr;
                  if (tmgv2 != nullptr)
                     delete[](tmgv2);
                  tmgv2= nullptr;
                  if (tmgvf2 != nullptr)
                     delete[](tmgvf2);
                  tmgvf2= nullptr;
                  if (alt != nullptr)
                     delete[](alt);
                  alt= nullptr;
                    }
                 else{
                     cout<<"ERROR swath #2 does not cross the equator..NO L1C grid for swath #2"<<endl;
                 }

             }
            else{ cout<<"day_mode = 0 nightime...no L1C grid produced--exit (swath#2 --- nswath#2).............."<<day_mode<<endl;
               exit(1);
            }

           }//end index x1 and x5 condition
                     
          }//end nswath=2

          
 

//-------------HALF-ORBITS IN TWO NON-CONSECUTIVE PORTIONS ------        
//nswath =3  
//swath1
//------------------------------------------------------------- 
          c=0;
          if(nswath==3){
            if(ix1>=0 && ix5>0){
                norbs=ix2-ix1+1;
                kk=ix1;
                latmin=orb_lat_tot[ix1];
                latmax=orb_lat_tot[ix2];


               //compute solar zenith angle for checking day/night conditions for each swath ---
   //-------------------------------------------------------
                year1=year_tot[ix1];
                day1=day_tot[ix1];
                hour1=hour_tot[ix1];
                lat1=orb_lat_tot[ix1];
                lon1=orb_lon_tot[ix1];

                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                if(sunz1<=93) day_mode=1; else day_mode=0;

                cout<<"sunz1.."<<sunz1<<"lat1.."<<lat1<<endl;

                year1=year_tot[ix2];
                day1=day_tot[ix2];
                hour1=hour_tot[ix2];
                lat1=orb_lat_tot[ix2];
                lon1=orb_lon_tot[ix2];

                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                if(sunz1<=93 || day_mode==1) day_mode=1; else day_mode=0;
            
                cout<<"sunz1.."<<sunz1<<"lat1.."<<lat1<<endl;
                cout<<"day_mode at nswath #3 SWATH1--segment#1."<<day_mode<<endl;
           
              if(day_mode==1){//daylight
                //orbit direction
                l1cfile->orb_dir=orb_dir_tot[ix1];

                tswt=(double*)calloc(norbs,sizeof(double));
                latswt=(double*)calloc(norbs,sizeof(double));
                lonswt=(double*)calloc(norbs,sizeof(double));

                sum1=0,sum3=0,c=0;                
                while(kk<=ix2){
                   latswt[c]=orb_lat_tot[kk];
                   lonswt[c]=orb_lon_tot[kk];
                   tswt[c]=orb_time_tot[kk];
                   sum1+=orb_vel_tot[kk];
                   sum3+=grn_vel_tot[kk];
                   kk++;
                   c++;
                     }
  //computing mean swath velocity (m/s) --------------------------
                mov1 = (sum1 /norbs);//mean velocity in m/s x 10^3
                mgv1 = (sum3 /norbs)*1000;
                cout<<"mov1..."<<mov1<<"mvg1..."<<mgv1<<endl;

                status1=ect_swt(1,l1cfile,norbs,tswt,latswt,lonswt,&tcross1,&loncross1);
                cout<<"nswath==3 --tcross equat crossing in (s)..swath#1.(segment #1).."<<tcross1<<"loncross1..."<<loncross1<<endl;


                c=0;                 
                if(tcross1<0.){

                  if(tswt != nullptr)
                     delete[](tswt);
                  tswt= nullptr;
                  if (latswt != nullptr)
                     delete[](latswt);
                  latswt= nullptr;
                  if (lonswt != nullptr)
                     delete[](lonswt);
                  lonswt= nullptr;  

                  norbs=ix4-ix3+1;  

                  tswt=(double*)calloc(norbs,sizeof(double));
                  latswt=(double*)calloc(norbs,sizeof(double));
                  lonswt=(double*)calloc(norbs,sizeof(double));

                  kk=ix3;
                  sum1=0.;
                  tcross1=-999,loncross1=-999.;


                  l1cfile->orb_dir=orb_dir_tot[ix3];

                   //compute solar zenith angle for checking day/night conditions for each swath ---
   //-------------------------------------------------------
                  year1=year_tot[ix3];
                  day1=day_tot[ix3];
                  hour1=hour_tot[ix3];
                  lat1=orb_lat_tot[ix3];
                  lon1=orb_lon_tot[ix3];

                  sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                  if(sunz1<=93) day_mode=1; else day_mode=0;

                  year1=year_tot[ix4];
                  day1=day_tot[ix4];
                  hour1=hour_tot[ix4];
                  lat1=orb_lat_tot[ix4];
                  lon1=orb_lon_tot[ix4];

                  sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                  if(sunz1<=93 || day_mode==1) day_mode=1; else day_mode=0;

                  cout<<"day_mode at nswath #3 SWATH1--segment#2."<<day_mode<<endl;

                            sum1=0,sum3=0,c=0;
                 while(kk<=ix4){
                      latswt[c]=orb_lat_tot[kk];
                      lonswt[c]=orb_lon_tot[kk];
                      tswt[c]=orb_time_tot[kk];        
                   sum1+=orb_vel_tot[kk];
                      sum3+=grn_vel_tot[kk];
                   kk++;
                      c++;
                     }
                //computing mean swath velocity (m/s) --------------------------
                  mov1 = (sum1 /norbs);//mean velocity in m/s x 10^3
                  mgv1 = (sum3 /norbs)*1000;
                  cout<<"mov1..."<<mov1<<"mvg1..."<<mgv1<<endl;

                  status2=ect_swt(1,l1cfile,norbs,tswt,latswt,lonswt,&tcross1,&loncross1);
                  cout<<"nswath==3 -swath1---tcross equat crossing in (s)..swath#1.(segment #2).."<<tcross1<<"loncross1..."<<loncross1<<endl;

                } //end segment2

                 l1cfile->num_gridlines=4000;


                 if(status1==0 || status2==0 && day_mode==1){                 
                   tmgv1=(double*)calloc(l1cfile->num_gridlines,sizeof(double));
                

                   cout<<"#gridlines..."<<l1cfile->num_gridlines<<"norbs.."<<norbs<<endl;

                   for(int i=0;i<l1cfile->num_gridlines;i++)   tmgv1[i]=-1;

                   swtime_swt2_segment(1,l1cinput,l1cfile,norbs,tswt,tcross1,mgv1,tmgv1);//number of gridlines may change here!!!!!!!!! not anymore 4000 and assymetric around the equator so bina and binb not the same!!!                               
                   num_gridlines = l1cfile->num_gridlines;

                   double *tmgv1_segment=(double*)calloc(l1cfile->num_gridlines,sizeof(double));
                   tmgvf=(double*)calloc(l1cfile->num_gridlines,sizeof(double));

                   int ss=0;
                   for(int i=0;i<4000;i++)
                       {
                        if(tmgv1[i]>=0){                              
                            tmgv1_segment[ss]=tmgv1[i];
                            ss++;
                           }
                        }

                    cout<<"#segment gridlines.ss counter.."<<ss<<"equal to num_gridlines..."<<num_gridlines<<endl;
               
                   if(tswt != nullptr)
                     delete[](tswt);
                   tswt= nullptr;
                   if (latswt != nullptr)
                     delete[](latswt);
                   latswt= nullptr;
                   if (lonswt != nullptr)
                     delete[](lonswt);
                   lonswt= nullptr;
              
                   cout<<"number of across bins L1C grid...#.."<<l1cfile->nbinx<<"l1cfile->n_views..."<<l1cfile->n_views<<endl;
                   lat_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                   lon_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                   alt=allocate2d_float(num_gridlines, l1cfile->nbinx);

                   orb_array2* posgrid = new orb_array2[num_gridlines]();//these are doubles
                   orb_array2* velgrid = new orb_array2[num_gridlines]();

                   orb_array2* posgrid2 = new orb_array2[num_gridlines]();//these are doubles
                   orb_array2* velgrid2 = new orb_array2[num_gridlines]();

                   orb_interp2(number_orecords,num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgv1_segment,posgrid,velgrid);
                                               
                   for(int i=0;i<num_gridlines-1;i++){
                          pos_norm=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1]+posgrid[i][2]*posgrid[i][2]);
                          clatg2=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1])/pos_norm;
                          rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                          v2[0]=velgrid[i][0]*rl2/pos_norm;
                          v2[1]=velgrid[i][1]*rl2/pos_norm;
                          v2[2]=velgrid[i][2]*rl2/pos_norm;

                          vi=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])*1000;
                          toff=vi/mgv1;    
                          tmgvf[i]=tmgv1_segment[i]+toff;
                            }



                   deltasec=tmgvf[num_gridlines-2]-tmgvf[0]+1;
                   cout<<"deltasec..swath."<<deltasec<<endl;

                   orb_interp2(number_orecords,num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgvf,posgrid2,velgrid2);
                               
                  //angle subsat track----
                   omf2=(1-fe)*(1-fe);

                   for(int i=0;i<num_gridlines-1;i++){ 
                        pos_norm=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1]+posgrid2[i][2]*posgrid2[i][2]);
                        clatg2=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1])/pos_norm;
                        rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                        v1[0]=(posgrid2[i][0])*rl2/pos_norm;
                        v1[1]=(posgrid2[i][1])*rl2/pos_norm;
                        v1[2]=(posgrid2[i][2])*rl2/pos_norm;

                        v2[0]=(velgrid2[i][0])*rl2/pos_norm;
                        v2[1]=(velgrid2[i][1])*rl2/pos_norm;
                        v2[2]=(velgrid2[i][2])*rl2/pos_norm;         

                        cross_product_double(v1,v2,vecout);
                        nvec=cross_product_norm_double(v1,v2);//length of orb norm vect

                        orbnorm[0]=vecout[0]/nvec;
                        orbnorm[1]=vecout[1]/nvec;
                        orbnorm[2]=vecout[2]/nvec;

                    for(int j=0;j<l1cfile->nbinx;j++){

                      pos_norm=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);//first gridline 

                      oangle=asin((j-(l1cfile->nbinx-1)/2)*5.2/pos_norm);

                      if(j==(l1cfile->nbinx-1)/2){
                          oangle_nad=asin((j-(l1cfile->nbinx-1)/2)*5.2/pos_norm);
                          //Geocentric vector
                          G[0]=v1[0]*cos(oangle_nad)-orbnorm[0]*pos_norm*sin(oangle_nad);
                          G[1]=v1[1]*cos(oangle_nad)-orbnorm[1]*pos_norm*sin(oangle_nad);
                          G[2]=v1[2]*cos(oangle_nad)-orbnorm[2]*pos_norm*sin(oangle_nad);                  
                        }

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
                                                                                      
                      //altitude
                      clatg = cos(atan(omf2*tan(glat*M_PI/180.)));
                      rl = Re*(1.-fe)/sqrt(1.-(2.-fe)*fe*clatg*clatg);
                      alt[i][j] = (gnorm - rl)*1000;//in meters                                 
                   }
                  }
                                         
                   if (posgrid != nullptr)
                     delete[](posgrid);
               
                   if (velgrid != nullptr)
                     delete[](velgrid);  
                   
                   if (posgrid2 != nullptr)
                     delete[](posgrid2);

                   if (velgrid2 != nullptr)
                     delete[](velgrid2);           

                   l1cfile->num_gridlines = l1cfile->num_gridlines-1;

                   if(l1cinput->grantype==1){
                        tswt_ini_sec=tfile_ini_sec+tmgvf[0];
                        unix2ymdhms(tswt_ini_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);
                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_ini=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        tswt_ini_file=y_swt+mo_swt+d_swt+"T"+h_swt+mi_swt+s_swt2.substr(0,2);

                        tswt_end_sec=tfile_ini_sec+tmgvf[num_gridlines-2];
                        unix2ymdhms(tswt_end_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);

                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_end=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        l1cfile->tswt_ini=tswt_ini;
                        l1cfile->tswt_end=tswt_end;
                        l1cfile->tswt_ini_file=tswt_ini_file;  
                     create_SOCEA2(1,l1cinput,l1cfile,lat_gd,lon_gd,alt,tmgvf);
                     }
                   else if(l1cinput->grantype==0)
                    {
//--------------------------------------------------------------
//granule processing---------------------------------------------                 
//---------------------------------------------------------------
                   numgran=144*2;
                   l1cfile->numgran=numgran;
                   gd_per_gran=round(num_gridlines/10);
                   l1cfile->gd_per_gran=gd_per_gran;

                   cout<<"estimated # of granules to be processed..."<<numgran<<"gd_per_gran..."<<gd_per_gran<<"#gridlines.."<<num_gridlines<<endl;

                   write_L1C_granule2(1,l1cfile,l1cinput,tmgvf,lat_gd, lon_gd,alt);

                    }
                   else{
                    cout<<"ERROR selecting grantype, must be 0: granules or 1: swath........................."<<endl;
                    exit(1);
                    }
    
//-----------------------------------------------------------------
//-----------------------------------------------------------------                                    
                   if (lat_gd != nullptr)
                      delete [](lat_gd);
                   lat_gd = nullptr;
                   if (lon_gd != nullptr)
                       delete [](lon_gd);
                   lon_gd = nullptr;
                   if (tmgv1 != nullptr)
                     delete[](tmgv1);
                   tmgv1= nullptr;
                   if (tmgv1_segment != nullptr)
                     delete[](tmgv1_segment);                    
                   if (tmgvf != nullptr)
                     delete[](tmgvf);
                   tmgvf= nullptr;
                   if (alt != nullptr)
                     delete[](alt);
                   alt= nullptr;
                                   
                   }//end tcross (all segments)               
                 else{
                        cout<<"ERROR swath #1 does not cross the equator..."<<endl;
                         cout<<"checking swath #2..."<<endl;
                   }   
                 } //end day_mode
            
                 else{
                   cout<<"day_mode==0 (nighttime)...."<<day_mode<<endl;
                   cout<<"checking swath #2..."<<endl;
                 }

               }//end indexes x1 and x5 condition          

//-----------------------------------------------------------------------------------
//nswath =3
//swath 2
//----------------------------------------------------------------------------------

                kk=ix5;
                c=0;
                norbs=ix6-ix5+1;
                latmin=orb_lat_tot[ix5];
                latmax=orb_lat_tot[ix6];


      //compute solar zenith angle for checking day/night conditions for each swath ---
   //-------------------------------------------------------
                year1=year_tot[ix5];
                day1=day_tot[ix5];
                hour1=hour_tot[ix5];
                lat1=orb_lat_tot[ix5];
                lon1=orb_lon_tot[ix5];

                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                if(sunz1<=93) day_mode=1; else day_mode=0;
                cout<<"sunz1.."<<sunz1<<"lat1.."<<lat1<<endl;

                year1=year_tot[ix6];
                day1=day_tot[ix6];
                hour1=hour_tot[ix6];
                lat1=orb_lat_tot[ix6];
                lon1=orb_lon_tot[ix6];

                sunangs_(&year1,&day1,&hour1,&lon1,&lat1,&sunz1,&suna1);
                if(sunz1<=93 || day_mode==1) day_mode=1; else day_mode=0;

                cout<<"day_mode at nswath #3 SWATH2."<<day_mode<<endl;
                cout<<"sunz1.."<<sunz1<<"lat1.."<<lat1<<endl;

               if(day_mode==1){
                //orbit direction
                l1cfile->orb_dir=orb_dir_tot[ix5];

                tswt2=(double*)calloc(norbs,sizeof(double));
                latswt2=(double*)calloc(norbs,sizeof(double));
                lonswt2=(double*)calloc(norbs,sizeof(double));

                sum2=0,sum3=0,c=0;
                while(kk<=ix6){
                   latswt2[c]=orb_lat_tot[kk];
                   lonswt2[c]=orb_lon_tot[kk];
                   tswt2[c]=orb_time_tot[kk];

                   if(orb_lat_tot[kk]<=latmin){
                       latmin=orb_lat_tot[kk];
                        }
                   if(orb_lat_tot[kk]>=latmax){
                       latmax=orb_lat_tot[kk];
                        }

                   sum2+=orb_vel_tot[kk];
                   sum3+=grn_vel_tot[kk];
                   kk++;
                   c++;
                     }
                mov2 = (sum2 /norbs);//mean velocity in m/s x 10^3
                mgv2 = (sum3 /norbs)*1000;//in meters/s  

                cout<<"mov2..."<<mov2<<"mvg2..."<<mgv2<<endl;

                status=ect_swt(2,l1cfile,norbs,tswt2,latswt2,lonswt2,&tcross2,&loncross2);
                cout<<"nswath==3 -swath2--tcross equat crossing in (s)..swath#2..."<<tcross2<<"loncross2..."<<loncross2<<endl;
                cout<<"latmin.."<<latmin<<"latmax..."<<latmax<<endl;

                if (latswt2 != nullptr)
                     delete[](latswt2);
                latswt2= nullptr;
                if (lonswt2 != nullptr)
                     delete[](lonswt2);
                lonswt2= nullptr;               

                l1cfile->num_gridlines=4000;

                if(status==0 && day_mode==1){
                  tmgv2=(double*)calloc(l1cfile->num_gridlines,sizeof(double)); 
                  tmgvf2=(double*)calloc(l1cfile->num_gridlines,sizeof(double));

                  swtime_swt2(2,l1cinput,l1cfile,norbs,tswt2,tcross2,mgv2,tmgv2);                 
                  num_gridlines = l1cfile->num_gridlines;

                  if(tswt2 != nullptr)
                     delete[](tswt2);
                  tswt2= nullptr;

                  cout<<"number of across bins L1C grid...#.."<<l1cfile->nbinx<<"l1cfile->n_views..."<<l1cfile->n_views<<endl;
                  lat_gd = allocate2d_float(4000, l1cfile->nbinx);
                  lon_gd = allocate2d_float(4000, l1cfile->nbinx);
                  alt=allocate2d_float(4000, l1cfile->nbinx);

                  orb_array2* posgrid = new orb_array2[num_gridlines]();//these are doubles
                  orb_array2* velgrid = new orb_array2[num_gridlines]();
                  orb_array2* posgrid2 = new orb_array2[num_gridlines]();//these are doubles
                  orb_array2* velgrid2 = new orb_array2[num_gridlines]();

                  orb_interp2(number_orecords_tot, num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgv2, posgrid,velgrid);
                 
                  for(int i=0;i<num_gridlines-1;i++){
                          pos_norm=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1]+posgrid[i][2]*posgrid[i][2]);
                          clatg2=sqrt(posgrid[i][0]*posgrid[i][0]+posgrid[i][1]*posgrid[i][1])/pos_norm;
                          rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                          v2[0]=velgrid[i][0]*rl2/pos_norm;
                          v2[1]=velgrid[i][1]*rl2/pos_norm;
                          v2[2]=velgrid[i][2]*rl2/pos_norm;

                          vi=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])*1000;
                          toff=vi/mgv2;    
                          tmgvf2[i]=tmgv2[i]+toff;
                            }

                  orb_interp2(number_orecords_tot, num_gridlines,orb_time_tot,posr_tot,velr_tot,tmgvf2, posgrid2,velgrid2);

                  //angle subsat track----
                  omf2=(1-fe)*(1-fe);

                  for(int i=0;i<num_gridlines-1;i++){ 
                        pos_norm=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1]+posgrid2[i][2]*posgrid2[i][2]);
                        clatg2=sqrt(posgrid2[i][0]*posgrid2[i][0]+posgrid2[i][1]*posgrid2[i][1])/pos_norm;
                        rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                        v1[0]=(posgrid2[i][0])*rl2/pos_norm;
                        v1[1]=(posgrid2[i][1])*rl2/pos_norm;
                        v1[2]=(posgrid2[i][2])*rl2/pos_norm;

                        v2[0]=(velgrid2[i][0])*rl2/pos_norm;
                        v2[1]=(velgrid2[i][1])*rl2/pos_norm;
                        v2[2]=(velgrid2[i][2])*rl2/pos_norm;         

                        cross_product_double(v1,v2,vecout);
                        nvec=cross_product_norm_double(v1,v2);//length of orb norm vect

                        orbnorm[0]=vecout[0]/nvec;
                        orbnorm[1]=vecout[1]/nvec;
                        orbnorm[2]=vecout[2]/nvec;

                   for(int j=0;j<l1cfile->nbinx;j++){

                      pos_norm=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);//first gridline 

                      oangle=asin((j-(l1cfile->nbinx-1)/2)*5.2/pos_norm);

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
                                                   
                      //altitude
                      clatg = cos(atan(omf2*tan(glat*M_PI/180.)));
                      rl = Re*(1.-fe)/sqrt(1.-(2.-fe)*fe*clatg*clatg);
                      alt[i][j] = (gnorm - rl)*1000;//in meters              
                   }
                  }
                  
                   if (posgrid != nullptr)
                     delete[](posgrid);
               
                   if (velgrid != nullptr)
                     delete[](velgrid);  
                   
                   if (posgrid2 != nullptr)
                     delete[](posgrid2);

                   if (velgrid2 != nullptr)
                     delete[](velgrid2);           

                   l1cfile->num_gridlines = l1cfile->num_gridlines-1;
                
                   if(l1cinput->grantype==1){
                        tswt_ini_sec=tfile_ini_sec+tmgvf2[0];
                        unix2ymdhms(tswt_ini_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);
                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_ini=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        tswt_ini_file=y_swt+mo_swt+d_swt+"T"+h_swt+mi_swt+s_swt2.substr(0,2);

                        tswt_end_sec=tfile_ini_sec+tmgvf2[num_gridlines-2];
                        unix2ymdhms(tswt_end_sec,&oneyear,&onemon,&oneday, &onehour, &onemin, &onesec);

                        y_swt=std::to_string(oneyear);
                        mo_swt=std::to_string(onemon);
                        d_swt=std::to_string(oneday);
                        h_swt=std::to_string(onehour);
                        mi_swt=std::to_string(onemin);
                        s_swt2=std::to_string(round(onesec));

                        length = (int) floor( log10 (onemon) ) + 1;
                        if(length==1) mo_swt="0"+mo_swt;
                        length = (int) floor( log10 (oneday) ) + 1;
                        if(length==1) d_swt="0"+d_swt;
                        if(onehour==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onehour+logoff)) + 1;
                        if(length==1) h_swt="0"+h_swt;
                        if(onemin==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (onemin+logoff)) + 1;
                        if(length==1) mi_swt="0"+mi_swt;
                        if(onesec==0) logoff=1; else logoff=0;
                        length = (int) floor( log10 (round(onesec)+logoff)) + 1;
                        if(length==1) s_swt2="0"+s_swt2;

                        tswt_end=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                        l1cfile->tswt_ini=tswt_ini;
                        l1cfile->tswt_end=tswt_end;
                        l1cfile->tswt_ini_file=tswt_ini_file;
                     create_SOCEA2(1,l1cinput,l1cfile,lat_gd,lon_gd,alt,tmgvf2);
                     }
                   else if(l1cinput->grantype==0)
                    {                  
                  //--------------------------------------------------------------
//granule processing---------------------------------------------                 
//---------------------------------------------------------------


                   deltasec=tmgvf2[num_gridlines-2]-tmgvf2[0]+1;
                   cout<<"deltasec..swath."<<deltasec<<endl;

                   numgran=144*2;

                   l1cfile->numgran=numgran;
                   gd_per_gran=round(num_gridlines/10);//10 granules per half orbit
                   l1cfile->gd_per_gran=gd_per_gran;

                   cout<<"estimated # of granules to be processed..."<<numgran<<"gd_per_gran..."<<gd_per_gran<<"#gridlines.."<<num_gridlines<<endl;

                   write_L1C_granule2(2,l1cfile,l1cinput,tmgvf2,lat_gd, lon_gd,alt);
               

                    }
                   else{
                    cout<<"ERROR selecting grantype, must be 0: granules or 1: swath........................."<<endl;
                    exit(1);
                    }


//---------------------------------------------------------------
//---------------------------------------------------------------
                  if (lat_gd != nullptr)
                    delete [](lat_gd);
                  lat_gd = nullptr;
                  if (lon_gd != nullptr)
                       delete [](lon_gd);
                  lon_gd = nullptr;
                  if (tmgv2 != nullptr)
                     delete[](tmgv2);
                  tmgv2= nullptr;
                  if (tmgvf2 != nullptr)
                     delete[](tmgvf2);
                  tmgvf2= nullptr;
                  if (alt != nullptr)
                     delete[](alt);
                  alt= nullptr;


           }
        else {
                     cout<<"ERROR swath #2 does not cross the equator..NO L1C grid for swath #2"<<endl;
              }

         }//end day_mode==1
         else{ cout<<"day_mode = 0 nightime...no L1C grid produced--exit (swath#2 ---nswath#3).............."<<day_mode<<endl;
               exit(1);
            }      

    }//end nswath =3

          
       delete [] (year_tot);
       delete [] (day_tot);
       delete [] (hour_tot);
       delete [] (orb_time_tot);
       delete [] (orb_lon_tot);
       delete [] (orb_lat_tot);
       delete [] (orb_dir_tot);
       delete [] (posr_tot);
       delete [] (velr_tot);
       delete [] (orb_vel_tot);
       delete [] (grn_vel_tot); 

         
    return 0;
    }




    bool L1C::binIntersectsPix4corn4_l1c2(l1c_filehandle* l1cfile, L1C_input* l1cinput, short row, short col, float** lat_gd, float** lon_gd, Polygon_t& pixelPoly, double areaFracBox[3][3], double areabinBox[3][3]) {
        bool result = false;
        double  ws_lat, wn_lat, es_lat, en_lat, ws_lon, wn_lon, es_lon, en_lon;
        double ws_lat2, wn_lat2, es_lat2, en_lat2, ws_lon2, wn_lon2, es_lon2, en_lon2;
        double ws_lat3, wn_lat3, es_lat3, en_lat3, ws_lon3, wn_lon3, es_lon3, en_lon3;
        double ws_lat4, wn_lat4, es_lat4, en_lat4, ws_lon4, wn_lon4, es_lon4, en_lon4;
        double ws_lat5, wn_lat5, es_lat5, en_lat5, ws_lon5, wn_lon5, es_lon5, en_lon5;
        double ws_lat6, wn_lat6, es_lat6, en_lat6, ws_lon6, wn_lon6, es_lon6, en_lon6;
        double ws_lat7, wn_lat7, es_lat7, en_lat7, ws_lon7, wn_lon7, es_lon7, en_lon7;
        double ws_lat8, wn_lat8, es_lat8, en_lat8, ws_lon8, wn_lon8, es_lon8, en_lon8;
        double ws_lat9, wn_lat9, es_lat9, en_lat9, ws_lon9, wn_lon9, es_lon9, en_lon9;
        bool ingeom = false;
        int pc = 1;
        double intersectArea = 0, binAreapix = 0, binAreagrid = 0.;
        std::string ws_lon_str, ws_lat_str, wn_lon_str, wn_lat_str, en_lon_str, en_lat_str, es_lon_str, es_lat_str, gridstr;
        std::string ws2_lon_str, ws2_lat_str, wn2_lon_str, wn2_lat_str, en2_lon_str, en2_lat_str, es2_lon_str, es2_lat_str;
        std::string ws3_lon_str, ws3_lat_str, wn3_lon_str, wn3_lat_str, en3_lon_str, en3_lat_str, es3_lon_str, es3_lat_str;
        std::string ws4_lon_str, ws4_lat_str, wn4_lon_str, wn4_lat_str, en4_lon_str, en4_lat_str, es4_lon_str, es4_lat_str;
        std::string ws5_lon_str, ws5_lat_str, wn5_lon_str, wn5_lat_str, en5_lon_str, en5_lat_str, es5_lon_str, es5_lat_str;
        std::string ws6_lon_str, ws6_lat_str, wn6_lon_str, wn6_lat_str, en6_lon_str, en6_lat_str, es6_lon_str, es6_lat_str;
        std::string ws7_lon_str, ws7_lat_str, wn7_lon_str, wn7_lat_str, en7_lon_str, en7_lat_str, es7_lon_str, es7_lat_str;
        std::string ws8_lon_str, ws8_lat_str, wn8_lon_str, wn8_lat_str, en8_lon_str, en8_lat_str, es8_lon_str, es8_lat_str;
        std::string ws9_lon_str, ws9_lat_str, wn9_lon_str, wn9_lat_str, en9_lon_str, en9_lat_str, es9_lon_str, es9_lat_str;
        float binres = 0, azgc,deltaphi, deltalam, azpixc;
        float aterm, bterm, dist_u, dist_v, theta, thetares,thetagrad,dist_corn2;
        double res1, res2, res3, res4, lat1, lon1, lat2, lon2;
        float az1, az2, az3, az4;
        bool infull = false;
        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());

        binres = (l1cinput->grid_resolution);//in km

//        cout << "binIntersectsPix4corn4_l1c...........................row......" << row << "col...." << col << endl;

        //*******************************************************************************************
        //center center bin
        //********************************************************************************************************
        ingeom = 0, binAreagrid = 0., binAreapix = 0, intersectArea = 0., pc = 1;
        deltaphi = (lat_gd[row][col + 1]-lat_gd[row][col]) * degrad;
        deltalam = (lon_gd[row][col + 1]-lon_gd[row][col] ) * degrad;

        //across-track grid size---HARVERSINE NOT VICENTY SO LESS ACCURATE!!!
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row][col] * degrad) * cos(lat_gd[row][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row + 1][col]-lat_gd[row][col]) * degrad;
        deltalam = (lon_gd[row + 1][col]-lon_gd[row][col]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row][col] * degrad) * cos(lat_gd[row + 1][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        lat1 = lat_gd[row][col];
        lon1 = lon_gd[row][col];
        lat2 = lat_gd[row][col + 1];
        lon2 = lon_gd[row][col + 1];
        geod.Inverse(lat1, lon1, lat2, lon2, res1);

        dist_corn2 = 0.5 * sqrt(binres * binres + binres * binres) * 1000;

        //azimuth grid cell calculation
        azgc=atan2(sin(deltalam)*cos(lat_gd[row+1][col]*degrad),cos(lat_gd[row][col]*degrad)*sin(lat_gd[row+1][col]*degrad)-sin(lat_gd[row][col]*degrad)*cos(lat_gd[row+1][col]*degrad)*cos(deltaphi));

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;
        azgc = azgc * 180 / M_PI;

     //   cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<"azgc.."<<azgc<<"theta..."<<theta*180 / M_PI<<endl;
       

        //nw corner 
        azpixc = (azgc - thetares) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "NW azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az1 = azpixc;
        geod.Direct(lat_gd[row][col], lon_gd[row][col], azpixc, dist_corn2, wn_lat, wn_lon);
        geod.Inverse(lat_gd[row][col], lon_gd[row][col], wn_lat, wn_lon, res1);
        //    cout<<"azgc.Center-Center."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //    cout<<"res1..nw."<<res1<<endl;

        //ne corner 
        azpixc = (azgc + thetares) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "NE azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az2 = azpixc;
        geod.Direct(lat_gd[row][col], lon_gd[row][col], azpixc, dist_corn2, en_lat, en_lon);
        geod.Inverse(lat_gd[row][col], lon_gd[row][col], en_lat, en_lon, res2);
        //  cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
      //sw corner 
        azpixc = (azgc - thetares - thetagrad) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "SW azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az3 = azpixc;
        geod.Direct(lat_gd[row][col], lon_gd[row][col], azpixc, dist_corn2, ws_lat, ws_lon);
        geod.Inverse(lat_gd[row][col], lon_gd[row][col], ws_lat, ws_lon, res3);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;

       //se corner 
        azpixc = (azgc + thetares + thetagrad) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "SE azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az4 = azpixc;
        geod.Direct(lat_gd[row][col], lon_gd[row][col], azpixc, dist_corn2, es_lat, es_lon);
        geod.Inverse(lat_gd[row][col], lon_gd[row][col], es_lat, es_lon, res4);
        //    cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

        Polygon_t gridPoly;
        ws_lon_str = std::to_string(ws_lon);
        ws_lat_str = std::to_string(ws_lat);
        wn_lon_str = std::to_string(wn_lon);
        wn_lat_str = std::to_string(wn_lat);
        en_lon_str = std::to_string(en_lon);
        en_lat_str = std::to_string(en_lat);
        es_lon_str = std::to_string(es_lon);
        es_lat_str = std::to_string(es_lat);
        gridstr = "POLYGON((" + ws_lon_str + " " + ws_lat_str + "," + wn_lon_str + " " + wn_lat_str + "," + en_lon_str + " " + en_lat_str + "," + es_lon_str + " " + es_lat_str + "," + ws_lon_str + " " + ws_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2
    /*
        cout<<"center center-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col]<<"lon_gd.."<<lon_gd[row][col]<<"wn_lat.."<<wn_lat<<"wn_lon.."<<wn_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col]<<"lon_gd.."<<lon_gd[row][col]<<"ws_lat.."<<ws_lat<<"ws_lon.."<<ws_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col]<<"lon_gd.."<<lon_gd[row][col]<<"en_lat.."<<en_lat<<"en_lon.."<<en_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col]<<"lon_gd.."<<lon_gd[row][col]<<"es_lat.."<<es_lat<<"es_lon.."<<es_lon<<endl;
    */
        areabinBox[1][1] = binAreagrid;

        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[1][1] = intersectArea / binAreagrid;

        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                std::deque<Polygon_t> output;//this is a list
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[1][1] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[1][1] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else

    //   cout<<"center center.."<<"areaFracBox[1][1].........................................."<<areaFracBox[1][1]<<endl;

    //***********************************************************************************************************
       //left center bin
        //********************************************************************************************************
        ingeom = 0, binAreagrid = 0., binAreapix = 0, intersectArea = 0., pc = 1;
        if (col >= 1) {

            deltaphi = (lat_gd[row][col]-lat_gd[row][col - 1]) * degrad;
            deltalam = (lon_gd[row][col]-lon_gd[row][col - 1]) * degrad;
            //across-track grid size
            aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row][col - 1] * degrad) * cos(lat_gd[row][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
            bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
            dist_u = Re * bterm;         //horizontal distance Harversine in km

            deltaphi = (lat_gd[row + 1][col - 1]-lat_gd[row][col - 1]) * degrad;
            deltalam = (lon_gd[row + 1][col - 1]-lon_gd[row][col - 1]) * degrad;
            //along-track distance
            aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row][col - 1] * degrad) * cos(lat_gd[row + 1][col - 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
            bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
            dist_v = Re * bterm;

            theta = atan2(dist_v, dist_u);
            thetagrad = 2 * theta * 180 / M_PI;
            thetares = (M_PI / 2 - theta) * 180 / M_PI;


            //nw corner
            azpixc = az1;
            geod.Direct(lat_gd[row][col - 1], lon_gd[row][col - 1], azpixc, dist_corn2, wn_lat2, wn_lon2);
            geod.Inverse(lat_gd[row][col - 1], lon_gd[row][col - 1], wn_lat2, wn_lon2, res1);
            //   cout<<"azgc..CENTER-LEFT"<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
            //   cout<<"res1..nw."<<res1<<endl;
           //ne corner 
            azpixc = az2;
            geod.Direct(lat_gd[row][col - 1], lon_gd[row][col - 1], azpixc, dist_corn2, en_lat2, en_lon2);
            geod.Inverse(lat_gd[row][col - 1], lon_gd[row][col - 1], en_lat2, en_lon2, res2);
            //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
           //sw corner 
            azpixc = az3;
            geod.Direct(lat_gd[row][col - 1], lon_gd[row][col - 1], azpixc, dist_corn2, ws_lat2, ws_lon2);
            geod.Inverse(lat_gd[row][col - 1], lon_gd[row][col - 1], ws_lat2, ws_lon2, res3);
            //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
           //se corner 
            azpixc = az4;
            geod.Direct(lat_gd[row][col - 1], lon_gd[row][col - 1], azpixc, dist_corn2, es_lat2, es_lon2);
            geod.Inverse(lat_gd[row][col - 1], lon_gd[row][col - 1], es_lat2, es_lon2, res4);
            //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

            Polygon_t gridPoly;
            ws2_lon_str = std::to_string(ws_lon2);
            ws2_lat_str = std::to_string(ws_lat2);
            wn2_lon_str = std::to_string(wn_lon2);
            wn2_lat_str = std::to_string(wn_lat2);
            en2_lon_str = std::to_string(en_lon2);
            en2_lat_str = std::to_string(en_lat2);
            es2_lon_str = std::to_string(es_lon2);
            es2_lat_str = std::to_string(es_lat2);
            gridstr = "POLYGON((" + ws2_lon_str + " " + ws2_lat_str + "," + wn2_lon_str + " " + wn2_lat_str + "," + wn_lon_str + " " + wn_lat_str + "," + ws_lon_str + " " + ws_lat_str + "," + ws2_lon_str + " " + ws2_lat_str + "))";
            bg::read_wkt(gridstr, gridPoly);

            ingeom = within(pixelPoly, gridPoly);
            binAreagrid = bg::area(gridPoly);
            binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2
        /*
            cout<<"center LEFT-----------------"<<endl;
            cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
            cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
            cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
            cout<<"lat_gd.."<<lat_gd[row][col-1]<<"lon_gd.."<<lon_gd[row][col-1]<<"wn_lat.."<<wn_lat2<<"wn_lon.."<<wn_lon2<<endl;
            cout<<"lat_gd.."<<lat_gd[row][col-1]<<"lon_gd.."<<lon_gd[row][col-1]<<"ws_lat.."<<ws_lat2<<"ws_lon.."<<ws_lon2<<endl;
            cout<<"lat_gd.."<<lat_gd[row][col-1]<<"lon_gd.."<<lon_gd[row][col-1]<<"en_lat.."<<wn_lat<<"en_lon.."<<wn_lon<<endl;
            cout<<"lat_gd.."<<lat_gd[row][col-1]<<"lon_gd.."<<lon_gd[row][col-1]<<"es_lat.."<<ws_lat<<"es_lon.."<<ws_lon<<endl;
        */

            areabinBox[1][0] = binAreagrid;

            if (ingeom > 0) { //pixel fully inside the grid cell
                infull = true;
                intersectArea = binAreapix;
                result = true;
                areaFracBox[1][0] = intersectArea / binAreagrid;
            }
            else {
                if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                    std::deque<Polygon_t> output;
                    if (bg::intersection(pixelPoly, gridPoly, output)) {
                        BOOST_FOREACH(Polygon_t const& p, output)
                        {
                            intersectArea = bg::area(p);
                            if (intersectArea > 0.0 && pc == 1) {
                                result = true;
                                areaFracBox[1][0] = intersectArea / binAreagrid;
                                if (intersectArea / binAreagrid < 0.000001) areaFracBox[1][0] = 0.0;
                                binAreapix = bg::area(pixelPoly);
                                pc++;
                            }
                        }
                    }
                }
            }//end else

        }//end col>=1 left center bin

       // cout<<"center LEFT.."<<"areaFracBox[1][0]........................................."<<areaFracBox[1][0]<<endl;


       //***************************************************************************
        //right center bin
       //**************************************************************************8
        ingeom = 0, binAreagrid = 0., binAreapix = 0, intersectArea = 0.;

        deltaphi = (lat_gd[row][col + 2]-lat_gd[row][col + 1]) * degrad;
        deltalam = (lon_gd[row][col + 2]-lon_gd[row][col + 1]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row][col + 1] * degrad) * cos(lat_gd[row][col + 2] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row][col + 1] - lat_gd[row + 1][col + 1]) * degrad;
        deltalam = (lon_gd[row][col + 1] - lon_gd[row + 1][col + 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row][col + 1] * degrad) * cos(lat_gd[row + 1][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(lat_gd[row][col + 1], lon_gd[row][col + 1], azpixc, dist_corn2, wn_lat3, wn_lon3);
        geod.Inverse(lat_gd[row][col + 1], lon_gd[row][col + 1], wn_lat3, wn_lon3, res1);
        //   cout<<"azgc..CENTER-RIGHT..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;

       //ne corner 
        azpixc = az2;
        geod.Direct(lat_gd[row][col + 1], lon_gd[row][col + 1], azpixc, dist_corn2, en_lat3, en_lon3);
        geod.Inverse(lat_gd[row][col + 1], lon_gd[row][col + 1], en_lat3, en_lon3, res2);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(lat_gd[row][col + 1], lon_gd[row][col + 1], azpixc, dist_corn2, ws_lat3, ws_lon3);
        geod.Inverse(lat_gd[row][col + 1], lon_gd[row][col + 1], ws_lat3, ws_lon3, res3);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(lat_gd[row][col + 1], lon_gd[row][col + 1], azpixc, dist_corn2, es_lat3, es_lon3);
        geod.Inverse(lat_gd[row][col + 1], lon_gd[row][col + 1], es_lat3, es_lon3, res4);
        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;
        ws3_lon_str = std::to_string(ws_lon3);
        ws3_lat_str = std::to_string(ws_lat3);
        wn3_lon_str = std::to_string(wn_lon3);
        wn3_lat_str = std::to_string(wn_lat3);
        en3_lon_str = std::to_string(en_lon3);
        en3_lat_str = std::to_string(en_lat3);
        es3_lon_str = std::to_string(es_lon3);
        es3_lat_str = std::to_string(es_lat3);
        gridstr = "POLYGON((" + es_lon_str + " " + es_lat_str + "," + en_lon_str + " " + en_lat_str + "," + en3_lon_str + " " + en3_lat_str + "," + es3_lon_str + " " + es3_lat_str + "," + es_lon_str + " " + es_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2

    /*
        cout<<"center RIGHT-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col+1]<<"lon_gd.."<<lon_gd[row][col+1]<<"wn_lat.."<<en_lat<<"wn_lon.."<<en_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col+1]<<"lon_gd.."<<lon_gd[row][col+1]<<"ws_lat.."<<es_lat<<"ws_lon.."<<es_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col+1]<<"lon_gd.."<<lon_gd[row][col+1]<<"en_lat.."<<en_lat3<<"en_lon.."<<en_lon3<<endl;
        cout<<"lat_gd.."<<lat_gd[row][col+1]<<"lon_gd.."<<lon_gd[row][col+1]<<"es_lat.."<<es_lat3<<"es_lon.."<<es_lon3<<endl;
    */

        areabinBox[1][2] = binAreagrid;
        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[1][2] = intersectArea / binAreagrid;
        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0.) {
                std::deque<Polygon_t> output;
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[1][2] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[1][2] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else

    //  cout<<"center RIGHT.."<<"areaFracBox[1][2]........................................."<<areaFracBox[1][2]<<endl;

    //******************************************************************************************************
     //upper center bin
    //*****************************************************************************************************
        ingeom = 0, binAreagrid = 0., binAreapix = 0, intersectArea = 0., pc = 1;

        deltaphi = (lat_gd[row + 1][col] - lat_gd[row + 1][col + 1]) * degrad;
        deltalam = (lon_gd[row + 1][col] - lon_gd[row + 1][col + 1]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row + 1][col] * degrad) * cos(lat_gd[row + 1][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row + 2][col] - lat_gd[row + 1][col]) * degrad;
        deltalam = (lon_gd[row + 2][col] - lon_gd[row + 1][col]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row + 1][col] * degrad) * cos(lat_gd[row + 2][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(lat_gd[row + 1][col], lon_gd[row + 1][col], azpixc, dist_corn2, wn_lat4, wn_lon4);
        geod.Inverse(lat_gd[row + 1][col], lon_gd[row + 1][col], wn_lat4, wn_lon4, res1);
        //   cout<<"azgc..CENTER-UPPER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;
       //ne corner 
        azpixc = az2;
        geod.Direct(lat_gd[row + 1][col], lon_gd[row + 1][col], azpixc, dist_corn2, en_lat4, en_lon4);
        geod.Inverse(lat_gd[row + 1][col], lon_gd[row + 1][col], en_lat4, en_lon4, res2);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(lat_gd[row + 1][col], lon_gd[row + 1][col], azpixc, dist_corn2, ws_lat4, ws_lon4);
        geod.Inverse(lat_gd[row + 1][col], lon_gd[row + 1][col], ws_lat4, ws_lon4, res3);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;

       //se corner 
        azpixc = az4;
        geod.Direct(lat_gd[row + 1][col], lon_gd[row + 1][col], azpixc, dist_corn2, es_lat4, es_lon4);
        geod.Inverse(lat_gd[row + 1][col], lon_gd[row + 1][col], es_lat4, es_lon4, res4);
        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

        ws4_lon_str = std::to_string(ws_lon4);
        ws4_lat_str = std::to_string(ws_lat4);
        wn4_lon_str = std::to_string(wn_lon4);
        wn4_lat_str = std::to_string(wn_lat4);
        en4_lon_str = std::to_string(en_lon4);
        en4_lat_str = std::to_string(en_lat4);
        es4_lon_str = std::to_string(es_lon4);
        es4_lat_str = std::to_string(es_lat4);
        gridstr = "POLYGON((" + wn_lon_str + " " + wn_lat_str + "," + wn4_lon_str + " " + wn4_lat_str + "," + en4_lon_str + " " + en4_lat_str + "," + en_lon_str + " " + en_lat_str + "," + wn_lon_str + " " + wn_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2

    /*
        cout<<"center UPPER-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col]<<"lon_gd.."<<lon_gd[row+1][col]<<"wn_lat.."<<wn_lat4<<"wn_lon.."<<wn_lon4<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col]<<"lon_gd.."<<lon_gd[row+1][col]<<"ws_lat.."<<wn_lat<<"ws_lon.."<<wn_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col]<<"lon_gd.."<<lon_gd[row+1][col]<<"en_lat.."<<en_lat4<<"en_lon.."<<en_lon4<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col]<<"lon_gd.."<<lon_gd[row+1][col]<<"es_lat.."<<en_lat<<"es_lon.."<<en_lon<<endl;
    */

        areabinBox[2][1] = binAreagrid;
        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[2][1] = intersectArea / binAreagrid;
            //        cout<<"pixel fully inside gricell [2][1].."<<endl;
        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                std::deque<Polygon_t> output;
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[2][1] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[2][1] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else


    //  cout<<"CENTER UPPER.."<<"areaFracBox[2][1]........................................."<<areaFracBox[2][1]<<endl;


    //********************************************************
    //upper left  bin
    //*******************************************************
        ingeom = 0, binAreagrid = 0., binAreapix = 0, intersectArea = 0., pc = 1;

        deltaphi = (lat_gd[row + 1][col - 1] - lat_gd[row + 1][col]) * degrad;
        deltalam = (lon_gd[row + 1][col - 1] - lon_gd[row + 1][col]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row + 1][col - 1] * degrad) * cos(lat_gd[row + 1][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row + 2][col - 1] - lat_gd[row + 1][col - 1]) * degrad;
        deltalam = (lon_gd[row + 2][col - 1] - lon_gd[row + 1][col - 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row + 1][col - 1] * degrad) * cos(lat_gd[row + 2][col - 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], azpixc, dist_corn2, wn_lat5, wn_lon5);
        geod.Inverse(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], wn_lat5, wn_lon5, res1);
        //   cout<<"azgc..LEFT-UPPER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;
       //ne corner 
        azpixc = az2;
        geod.Direct(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], azpixc, dist_corn2, en_lat5, en_lon5);
        geod.Inverse(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], en_lat5, en_lon5, res2);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], azpixc, dist_corn2, ws_lat5, ws_lon5);
        geod.Inverse(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], ws_lat5, ws_lon5, res3);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], azpixc, dist_corn2, es_lat5, es_lon5);
        geod.Inverse(lat_gd[row + 1][col - 1], lon_gd[row + 1][col - 1], es_lat5, es_lon5, res4);
        //  cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;
        ws5_lon_str = std::to_string(ws_lon5);
        ws5_lat_str = std::to_string(ws_lat5);
        wn5_lon_str = std::to_string(wn_lon5);
        wn5_lat_str = std::to_string(wn_lat5);
        en5_lon_str = std::to_string(en_lon5);
        en5_lat_str = std::to_string(en_lat5);
        es5_lon_str = std::to_string(es_lon5);
        es5_lat_str = std::to_string(es_lat5);
        gridstr = "POLYGON((" + wn2_lon_str + " " + wn2_lat_str + "," + wn5_lon_str + " " + wn5_lat_str + "," + wn4_lon_str + " " + wn4_lat_str + "," + wn_lon_str + " " + wn_lat_str + "," + wn2_lon_str + " " + wn2_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2 
    /*
        cout<<"LEFT UPPER-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col-1]<<"lon_gd.."<<lon_gd[row+1][col-1]<<"wn_lat.."<<wn_lat5<<"wn_lon.."<<wn_lon5<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col-1]<<"lon_gd.."<<lon_gd[row+1][col-1]<<"ws_lat.."<<wn_lat2<<"ws_lon.."<<wn_lon2<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col-1]<<"lon_gd.."<<lon_gd[row+1][col-1]<<"en_lat.."<<wn_lat4<<"en_lon.."<<wn_lon4<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col-1]<<"lon_gd.."<<lon_gd[row+1][col-1]<<"es_lat.."<<wn_lat<<"es_lon.."<<wn_lon<<endl;
    */

        areabinBox[2][0] = binAreagrid;
        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[2][0] = intersectArea / binAreagrid;
            //        cout<<"pixel fully inside gricell [2][0].."<<endl;
        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                std::deque<Polygon_t> output;
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[2][0] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[2][0] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else

    //  cout<<"LEFT UPPER.."<<"areaFracBox[2][0]........................................."<<areaFracBox[2][0]<<endl;

    //**************************************************************************************
    //upper right bin
    //**************************************************************************************
        ingeom = 0, binAreagrid = 0., binAreapix = 0, intersectArea = 0., pc = 1;

        deltaphi = (lat_gd[row + 1][col + 1] - lat_gd[row + 1][col + 2]) * degrad;
        deltalam = (lon_gd[row + 1][col + 1] - lon_gd[row + 1][col + 2]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row + 1][col + 1] * degrad) * cos(lat_gd[row + 1][col + 2] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row + 2][col + 1] - lat_gd[row + 1][col + 1]) * degrad;
        deltalam = (lon_gd[row + 2][col + 1] - lon_gd[row + 1][col + 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row + 1][col + 1] * degrad) * cos(lat_gd[row + 2][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], azpixc, dist_corn2, wn_lat6, wn_lon6);
        geod.Inverse(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], wn_lat6, wn_lon6, res1);
        //  cout<<"azgc..RIGHT-UPPER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
       //   cout<<"res1..nw."<<res1<<endl;
      //ne corner 
        azpixc = az2;
        geod.Direct(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], azpixc, dist_corn2, en_lat6, en_lon6);
        geod.Inverse(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], en_lat6, en_lon6, res2);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], azpixc, dist_corn2, ws_lat6, ws_lon6);
        geod.Inverse(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], ws_lat6, ws_lon6, res3);
        //  cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
      //se corner 
        azpixc = az4;
        geod.Direct(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], azpixc, dist_corn2, es_lat6, es_lon6);
        geod.Inverse(lat_gd[row + 1][col + 1], lon_gd[row + 1][col + 1], es_lat6, es_lon6, res4);
        //  cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

        ws6_lon_str = std::to_string(ws_lon6);
        ws6_lat_str = std::to_string(ws_lat6);
        wn6_lon_str = std::to_string(wn_lon6);
        wn6_lat_str = std::to_string(wn_lat6);
        en6_lon_str = std::to_string(en_lon6);
        en6_lat_str = std::to_string(en_lat6);
        es6_lon_str = std::to_string(es_lon6);
        es6_lat_str = std::to_string(es_lat6);
        gridstr = "POLYGON((" + en_lon_str + " " + en_lat_str + "," + en4_lon_str + " " + en4_lat_str + "," + en6_lon_str + " " + en6_lat_str + "," + en3_lon_str + " " + en3_lat_str + "," + en_lon_str + " " + en_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2
    /*
        cout<<"RIGHT UPPER-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col+1]<<"lon_gd.."<<lon_gd[row+1][col+1]<<"wn_lat.."<<en_lat4<<"wn_lon.."<<en_lon4<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col+1]<<"lon_gd.."<<lon_gd[row+1][col+1]<<"ws_lat.."<<en_lat<<"ws_lon.."<<en_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col+1]<<"lon_gd.."<<lon_gd[row+1][col+1]<<"en_lat.."<<en_lat6<<"en_lon.."<<en_lon6<<endl;
        cout<<"lat_gd.."<<lat_gd[row+1][col+1]<<"lon_gd.."<<lon_gd[row+1][col+1]<<"es_lat.."<<en_lat3<<"es_lon.."<<en_lon3<<endl;
    */

        areabinBox[2][2] = binAreagrid;
        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[2][2] = intersectArea / binAreagrid;
            //        cout<<"pixel fully inside gricell [2][2].."<<endl;
        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                std::deque<Polygon_t> output;
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[2][2] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[2][2] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else

    //  cout<<"RIGHT UPPER.."<<"areaFracBox[2][0]........................................."<<areaFracBox[2][2]<<endl;

    //************************************************************
    //lower center  bin
    //*************************************************************

        ingeom = 0, binAreagrid = 0., binAreapix = 0, intersectArea = 0., pc = 1;

        deltaphi = (lat_gd[row - 1][col] - lat_gd[row - 1][col + 1]) * degrad;
        deltalam = (lon_gd[row - 1][col] - lon_gd[row - 1][col + 1]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row - 1][col] * degrad) * cos(lat_gd[row - 1][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row - 1][col] - lat_gd[row][col]) * degrad;
        deltalam = (lon_gd[row - 1][col] - lon_gd[row][col]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row - 1][col] * degrad) * cos(lat_gd[row][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        /*    azpixc=(azgc-thetares)*degrad;//pixel bearing
            if(azpixc>M_PI | azpixc<-M_PI){
                    cout<<"problem with BEARING in across-gridline method...az<-180 or >180...."<<"NW azpixc in degrees.."<<azpixc*180/M_PI<<endl;
                     }
            azpixc=azpixc*180/M_PI;
        */
        azpixc = az1;
        geod.Direct(lat_gd[row - 1][col], lon_gd[row - 1][col], azpixc, dist_corn2, wn_lat7, wn_lon7);
        geod.Inverse(lat_gd[row - 1][col], lon_gd[row - 1][col], wn_lat7, wn_lon7, res1);
        //   cout<<"azgc..LOWER CENTER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;

       //ne corner 
        azpixc = az2;
        geod.Direct(lat_gd[row - 1][col], lon_gd[row - 1][col], azpixc, dist_corn2, en_lat7, en_lon7);
        geod.Inverse(lat_gd[row - 1][col], lon_gd[row - 1][col], en_lat7, en_lon7, res2);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(lat_gd[row - 1][col], lon_gd[row - 1][col], azpixc, dist_corn2, ws_lat7, ws_lon7);
        geod.Inverse(lat_gd[row - 1][col], lon_gd[row - 1][col], ws_lat7, ws_lon7, res3);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;

       //se corner 
        azpixc = az4;
        geod.Direct(lat_gd[row - 1][col], lon_gd[row - 1][col], azpixc, dist_corn2, es_lat7, es_lon7);
        geod.Inverse(lat_gd[row - 1][col], lon_gd[row - 1][col], es_lat7, es_lon7, res4);
        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

        ws7_lon_str = std::to_string(ws_lon7);
        ws7_lat_str = std::to_string(ws_lat7);
        wn7_lon_str = std::to_string(wn_lon7);
        wn7_lat_str = std::to_string(wn_lat7);
        en7_lon_str = std::to_string(en_lon7);
        en7_lat_str = std::to_string(en_lat7);
        es7_lon_str = std::to_string(es_lon7);
        es7_lat_str = std::to_string(es_lat7);
        gridstr = "POLYGON((" + ws7_lon_str + " " + ws7_lat_str + "," + ws_lon_str + " " + ws_lat_str + "," + es_lon_str + " " + es_lat_str + "," + es7_lon_str + " " + es7_lat_str + "," + ws7_lon_str + " " + ws7_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2 
    /*
        cout<<"LOWER CENTER-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col]<<"lon_gd.."<<lon_gd[row-1][col]<<"wn_lat.."<<ws_lat<<"wn_lon.."<<ws_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col]<<"lon_gd.."<<lon_gd[row-1][col]<<"ws_lat.."<<ws_lat7<<"ws_lon.."<<ws_lon7<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col]<<"lon_gd.."<<lon_gd[row-1][col]<<"en_lat.."<<es_lat<<"en_lon.."<<es_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col]<<"lon_gd.."<<lon_gd[row-1][col]<<"es_lat.."<<es_lat7<<"es_lon.."<<es_lon7<<endl;
    */
        areabinBox[0][1] = binAreagrid;
        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[0][1] = intersectArea / binAreagrid;
            //        cout<<"pixel fully inside gricell [0][1].."<<endl;
        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                std::deque<Polygon_t> output;
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[0][1] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[0][1] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else


    //  cout<<"CENTER LOWER.."<<"areaFracBox[0][1]........................................."<<areaFracBox[0][1]<<endl;


    //******************************************************************
    //lower left bin
    //****************************************************************
        ingeom = 0, binAreagrid = 0, binAreapix = 0, intersectArea = 0., pc = 1;

        deltaphi = (lat_gd[row - 1][col - 1] - lat_gd[row - 1][col]) * degrad;
        deltalam = (lon_gd[row - 1][col - 1] - lon_gd[row - 1][col]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row - 1][col - 1] * degrad) * cos(lat_gd[row - 1][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row - 1][col - 1] - lat_gd[row][col - 1]) * degrad;
        deltalam = (lon_gd[row - 1][col - 1] - lon_gd[row][col - 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row - 1][col - 1] * degrad) * cos(lat_gd[row][col - 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], azpixc, dist_corn2, wn_lat8, wn_lon8);
        geod.Inverse(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], wn_lat8, wn_lon8, res1);
        //   cout<<"azgc..LOWER LEFT..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;

       //ne corner 
        azpixc = az2;
        geod.Direct(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], azpixc, dist_corn2, en_lat8, en_lon8);
        geod.Inverse(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], en_lat8, en_lon8, res2);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], azpixc, dist_corn2, ws_lat8, ws_lon8);
        geod.Inverse(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], ws_lat8, ws_lon8, res3);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;

       //se corner 
        azpixc = az4;
        geod.Direct(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], azpixc, dist_corn2, es_lat8, es_lon8);
        geod.Inverse(lat_gd[row - 1][col - 1], lon_gd[row - 1][col - 1], es_lat8, es_lon8, res4);
        //  cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;
        ws8_lon_str = std::to_string(ws_lon8);
        ws8_lat_str = std::to_string(ws_lat8);
        wn8_lon_str = std::to_string(wn_lon8);
        wn8_lat_str = std::to_string(wn_lat8);
        en8_lon_str = std::to_string(en_lon8);
        en8_lat_str = std::to_string(en_lat8);
        es8_lon_str = std::to_string(es_lon8);
        es8_lat_str = std::to_string(es_lat8);
        gridstr = "POLYGON((" + ws8_lon_str + " " + ws8_lat_str + "," + ws2_lon_str + " " + ws2_lat_str + "," + ws_lon_str + " " + ws_lat_str + "," + ws7_lon_str + " " + ws7_lat_str + "," + ws8_lon_str + " " + ws8_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2 
    /*
        cout<<"LOWER LEFT-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col-1]<<"lon_gd.."<<lon_gd[row-1][col-1]<<"wn_lat.."<<ws_lat2<<"wn_lon.."<<ws_lon2<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col-1]<<"lon_gd.."<<lon_gd[row-1][col-1]<<"ws_lat.."<<ws_lat8<<"ws_lon.."<<ws_lon8<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col-1]<<"lon_gd.."<<lon_gd[row-1][col-1]<<"en_lat.."<<ws_lat<<"en_lon.."<<ws_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col-1]<<"lon_gd.."<<lon_gd[row-1][col-1]<<"es_lat.."<<ws_lat7<<"es_lon.."<<ws_lon7<<endl;
    */

        areabinBox[0][0] = binAreagrid;
        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[0][0] = intersectArea / binAreagrid;
        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                std::deque<Polygon_t> output;
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[0][0] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[0][0] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else


    //  cout<<"LEFT LOWER.."<<"areaFracBox[0][0]........................................."<<areaFracBox[0][0]<<endl;

    //******************************************************************
    //lower right bin
    //****************************************************************
        ingeom = 0, binAreagrid = 0, binAreapix = 0, intersectArea = 0., pc = 1;

        deltaphi = (lat_gd[row - 1][col + 1] - lat_gd[row - 1][col + 2]) * degrad;
        deltalam = (lon_gd[row - 1][col + 1] - lon_gd[row - 1][col + 2]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row - 1][col + 1] * degrad) * cos(lat_gd[row - 1][col + 2] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (lat_gd[row - 1][col + 1] - lat_gd[row][col + 1]) * degrad;
        deltalam = (lon_gd[row - 1][col + 1] - lon_gd[row][col + 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(lat_gd[row - 1][col + 1] * degrad) * cos(lat_gd[row][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], azpixc, dist_corn2, wn_lat9, wn_lon9);
        geod.Inverse(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], wn_lat9, wn_lon9, res1);
        //  cout<<"azgc..LOWER RIGHT..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //  cout<<"res1..nw."<<res1<<endl;

      //ne corner 
        azpixc = az2;
        geod.Direct(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], azpixc, dist_corn2, en_lat9, en_lon9);
        geod.Inverse(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], en_lat9, en_lon9, res2);
        //  cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
      //sw corner 
        azpixc = az3;
        geod.Direct(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], azpixc, dist_corn2, ws_lat9, ws_lon9);
        geod.Inverse(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], ws_lat9, ws_lon9, res3);
        //  cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;

      //se corner 
        azpixc = az4;
        geod.Direct(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], azpixc, dist_corn2, es_lat9, es_lon9);
        geod.Inverse(lat_gd[row - 1][col + 1], lon_gd[row - 1][col + 1], es_lat9, es_lon9, res4);
        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

        ws9_lon_str = std::to_string(ws_lon9);
        ws9_lat_str = std::to_string(ws_lat9);
        wn9_lon_str = std::to_string(wn_lon9);
        wn9_lat_str = std::to_string(wn_lat9);
        en9_lon_str = std::to_string(en_lon9);
        en9_lat_str = std::to_string(en_lat9);
        es9_lon_str = std::to_string(es_lon9);
        es9_lat_str = std::to_string(es_lat9);
        gridstr = "POLYGON((" + es7_lon_str + " " + es7_lat_str + "," + es_lon_str + " " + es_lat_str + "," + es3_lon_str + " " + es3_lat_str + "," + es9_lon_str + " " + es9_lat_str + "," + es7_lon_str + " " + es7_lat_str + "))";
        bg::read_wkt(gridstr, gridPoly);

        ingeom = within(pixelPoly, gridPoly);
        binAreagrid = bg::area(gridPoly);
        binAreapix = bg::area(pixelPoly);//area in m2 --Andoyer method --- DIVIDE BY 10^6 FOR KM2 

    /*
        cout<<"LOWER RIGHT-----------------"<<endl;
        cout<<"gd_row.."<<row<<"gd_col.."<<col<<endl;
        cout<<"dist_u.."<<dist_u<<"dist_v.."<<dist_v<<endl;
        cout<<"ingeom.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix.."<<binAreapix<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col+1]<<"lon_gd.."<<lon_gd[row-1][col+1]<<"wn_lat.."<<es_lat<<"wn_lon.."<<es_lon<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col+1]<<"lon_gd.."<<lon_gd[row-1][col+1]<<"ws_lat.."<<es_lat7<<"ws_lon.."<<es_lon7<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col+1]<<"lon_gd.."<<lon_gd[row-1][col+1]<<"en_lat.."<<es_lat3<<"en_lon.."<<es_lon3<<endl;
        cout<<"lat_gd.."<<lat_gd[row-1][col+1]<<"lon_gd.."<<lon_gd[row-1][col+1]<<"es_lat.."<<es_lat9<<"es_lon.."<<es_lon9<<endl;
    */

        areabinBox[0][2] = binAreagrid;
        if (ingeom > 0) { //pixel fully inside the grid cell
            infull = true;
            intersectArea = binAreapix;
            result = true;
            areaFracBox[0][2] = intersectArea / binAreagrid;
        }
        else {
            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                std::deque<Polygon_t> output;
                if (bg::intersection(pixelPoly, gridPoly, output)) {
                    BOOST_FOREACH(Polygon_t const& p, output)
                    {
                        intersectArea = bg::area(p);
                        if (intersectArea > 0.0 && pc == 1) {
                            result = true;
                            areaFracBox[0][2] = intersectArea / binAreagrid;
                            if (intersectArea / binAreagrid < 0.000001) areaFracBox[0][2] = 0.0;
                            pc++;
                        }
                    }
                }
            }
        }//end else

    //  cout<<"RIGHT LOWER.."<<"areaFracBox[0][2]........................................."<<areaFracBox[0][2]<<endl;



        if (infull == true) {
     //       cout << "pixelbox fully inside gridbox..........................................." << endl;
            //    exit(1);
        }

        return result;

    }


int L1C::create_SOCEA2(int swtd,L1C_input* l1cinput, l1c_filehandle* l1cfile,float** lat_gd, float **lon_gd,float **altitude,double *time_nad){
   string tswt_ini,tswt_end,tswt_ini_file,senstr,date_created;
   int asc_mode=-1;
   std::string ifile_str;
   std::string y_create,m_create,d_create,t_create,binstr,titlestr,extstr=".nc";
   char *ifile_char;
   string dirstr,prodstr,verstr;
   const char* filename_lt;
   int Nwest=-1, Neast=-1,Ngring=-1,midix=-1,dp=-1;   
   int p,ix=-1;
   float latemp=-1,lontemp1=-1,lontemp2=-1,dlat_gd=-1,dlon_gd=-1,dlat20=-1,dlon20=-1,lon360=-1;
   int re=-1,rw=-1;


   ifile_str = l1cinput->files[0];
   ifile_char = &ifile_str[0];
   file_format format = getFormat(ifile_char);

     if(format.type==FT_SPEXONE){
             senstr="SPEXone";
             l1cfile->nbinx=25;
             binstr="12";
             midix=12;
             titlestr="PACE SPEXone Level-1C Data";
         //    NVIEWS=l1cfile->n_views;
         //    NBANDS=400;
              }
        else if(format.type==FT_OCIS){
             senstr="OCI";
             l1cfile->nbinx=519;
             binstr="259";
             midix=259;
             titlestr="PACE OCI Level-1C Data";
         //    NVIEWS=2;
         //    NBANDS=249;
             }
        else if(format.type==FT_HARP2){
             senstr="HARP2";
             l1cfile->nbinx=457;
             binstr="228";
             midix=228;
             titlestr="PACE HARP2 Level-1C Data";
        //     NVIEWS=l1cfile->n_views;
        //     NBANDS=4;
            }
        else{cout<<"sensor by default is OCI option 2....."<<endl;
             senstr="OCI";
             l1cfile->nbinx=519;
             binstr="259";
             midix=259;
             titlestr="PACE OCI Level-1C Data";
         //    NVIEWS=2;
         //    NBANDS=249;
          }
   
   //time nadir in seconds

   tswt_ini_file=l1cfile->tswt_ini_file;
   prodstr="PACE."+tswt_ini_file+".L1C"+extstr;    

   l1cfile->gridname = prodstr.c_str();

   filename_lt = prodstr.c_str();
   char *gridchar=strdup(filename_lt);
   string l1c_str=filename_lt;

   NcFile* nc_output;
          try {
              nc_output = new NcFile(filename_lt, NcFile::replace);
              }
          catch (NcException& e) {
              e.what();
              cerr << "l1cgen l1c_pflag= 5 : producing L1C grid: "
               + l1c_str << endl;
              exit(1);
              }

   meta_l1c_grid(gridchar,l1cfile->num_gridlines,nc_output);
    //gobal attrs--
   asc_mode=l1cfile->orb_dir;
     if(asc_mode==1) dirstr="Ascending";
     else if(asc_mode==0) dirstr="Descending";
   tswt_ini=l1cfile->tswt_ini;
   tswt_end=l1cfile->tswt_end;

   verstr=l1cfile->version;
   verstr="V"+verstr.substr(0,4);
   nc_output->putAtt("processing_version",verstr);
   nc_output->putAtt("history",l1cinput->history);
   nc_output->putAtt("product_name",prodstr);
   nc_output->putAtt("startDirection",dirstr);
   nc_output->putAtt("endDirection",dirstr);
   nc_output->putAtt("time_coverage_start",tswt_ini);
   nc_output->putAtt("time_coverage_end",tswt_end);
   //vars
   NcGroup ba_grp=nc_output->getGroup("bin_attributes");
   NcVar v1=ba_grp.getVar("nadir_view_time");
   v1.putVar(&time_nad[0]);
//   v1=ba_grp.getVar("view_time_offsets");
//   v1.putVar(&time_off[0][0][0]);
   NcGroup geo_grp=nc_output->getGroup("geolocation_data");
   v1=geo_grp.getVar("latitude");
   v1.putVar(&lat_gd[0][0]);
   v1=geo_grp.getVar("longitude");
   v1.putVar(&lon_gd[0][0]);
   v1=geo_grp.getVar("altitude");
   v1.putVar(&altitude[0][0]);

   //GRING-------
   //determine the number of GCpoint indexes
   //default 6 for the swath sides + number of coordinates every 20 degrees latitude

//ascending pass
   if(l1cfile->orb_dir==1)
   {
        
   Nwest=round((lat_gd[l1cfile->num_gridlines-1][0]-lat_gd[0][0])/20);
   Neast=round((lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1]-lat_gd[0][l1cfile->nbinx-1])/20);
   }
   else //descending
   {
   Neast=(round(lat_gd[0][0]-lat_gd[l1cfile->num_gridlines-1][0])/20);
   Nwest=round((lat_gd[0][l1cfile->nbinx-1]-lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1])/20);
   }

//first NGring estimate----   
   Ngring=Nwest+Neast+6;


   if(Ngring>0)
   {
      float *latarr=(float*)calloc(Ngring,sizeof(float));
      float *lonarr=(float*)calloc(Ngring,sizeof(float));
      int *narr=(int*)calloc(Ngring,sizeof(int));
      int *p_west=(int*)calloc(Ngring,sizeof(int));
      int *p_east=(int*)calloc(Ngring,sizeof(int));

      //corners--counterclockwise and ascending pass
     if(l1cfile->orb_dir==1){
         if(Ngring==6)
         {
          latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
          latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];
          latarr[3]=lat_gd[0][0];
          latarr[4]=lat_gd[0][midix];
          latarr[5]=lat_gd[0][l1cfile->nbinx-1];

          lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
          lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];
          lonarr[3]=lon_gd[0][0];
          lonarr[4]=lon_gd[0][midix];
          lonarr[5]=lon_gd[0][l1cfile->nbinx-1];
         }
         else
         {
         latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
         latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];

         lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
         lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];
    
         latemp=latarr[2];     
         latemp-=20;         
         p=1;
         rw=0;
         //west side
         while(latemp>lat_gd[0][0])
         {
         latarr[2+p]=latemp;
         cout<<"p west--"<<p<<"point# = "<<3+p<<"lat20 = "<<latemp<<endl;
         p_west[rw]=3+p;
         latemp-=20;
         p++;
         rw++;
         }   

         p--;
       
         latarr[2+p+1]=lat_gd[0][0];
         latarr[2+p+2]=lat_gd[0][midix];
         latarr[2+p+3]=lat_gd[0][l1cfile->nbinx-1];

         lonarr[2+p+1]=lon_gd[0][0];
         lonarr[2+p+2]=lon_gd[0][midix];
         lonarr[2+p+3]=lon_gd[0][l1cfile->nbinx-1];



         latemp=latarr[2+p+3];
         latemp+=20;
         p++;
       
         int c=1;
         re=0;
         //east side
         while(latemp<lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1])
         {
         latarr[5+p]=latemp;
         cout<<"p east--"<<c<<"point# = "<<6+p<<"lat20 = "<<latemp<<endl;       
         p_east[re]=6+p;
         latemp+=20;
         p++;
         c++;
         re++;
         }

         p--;
      
         cout<<"mid points west.."<<rw<<"mid points east.."<<re<<endl;
         cout<<"# points in GRING = "<<6+p<<"# mid points west--"<<rw<<"# mid points east--"<<re<<endl;
//west
         for(int i=0;i<rw;i++)
         {
             ix=p_west[i]-1;
             cout<<"lat:"<<latarr[ix]<<"pwest--"<<p_west[i]<<endl;
          for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {
              if(latarr[ix]>lat_gd[row][0] && latarr[ix]<=lat_gd[row+1][0])
              {
        //       cout<<"mid point west# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;            

               if(lon_gd[row][0]<0.) lontemp1=lon_gd[row][0]+360;else lontemp1=lon_gd[row][0];
               if(lon_gd[row+1][0]<0.) lontemp2=lon_gd[row+1][0]+360;else lontemp2=lon_gd[row+1][0];

               dlat_gd=abs(lat_gd[row+1][0]-lat_gd[row][0]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row+1][0]);
               dlon20=dlat20*dlon_gd/dlat_gd;
   //            if(lontemp1>lontemp2) lon360=lontemp1-dlon20;else lon360=lontemp1+dlon20;
               if(lontemp1>lontemp2) lon360=lontemp2+dlon20;else lon360=lontemp2-dlon20;
               if(lon360>180) lon360=lon360-360.; 
               lonarr[ix]=lon360;
       //        cout<<"lon_gd row+1.."<<lon_gd[row+1][0]<<"lon_gd row.."<<lon_gd[row][0]<<"lon cring.."<<lonarr[ix]<<endl;
                break;
              }
          }
         }
//east
         for(int i=0;i<re;i++)
         {
             ix=p_east[i]-1;
             cout<<"lat:"<<latarr[ix]<<"peast--"<<p_east[i]<<endl;
             for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {         
              if(latarr[ix]>lat_gd[row][l1cfile->nbinx-1] && latarr[ix]<=lat_gd[row+1][l1cfile->nbinx-1])
              {
     //          cout<<"mid point east# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;

               if(lon_gd[row][l1cfile->nbinx-1]<0.) lontemp1=lon_gd[row][l1cfile->nbinx-1]+360;else lontemp1=lon_gd[row][l1cfile->nbinx-1];
               if(lon_gd[row+1][l1cfile->nbinx-1]<0.) lontemp2=lon_gd[row+1][l1cfile->nbinx-1]+360;else lontemp2=lon_gd[row+1][l1cfile->nbinx-1];

               dlat_gd=abs(lat_gd[row+1][l1cfile->nbinx-1]-lat_gd[row][l1cfile->nbinx-1]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row][l1cfile->nbinx-1]);
               dlon20=dlat20*dlon_gd/dlat_gd;
               if(lontemp1>lontemp2) lon360=lontemp1-dlon20;else lon360=lontemp1+dlon20;
               if(lon360>180) lon360=lon360-360.;
               lonarr[ix]=lon360;
       //        cout<<"lon_gd row+1.."<<lon_gd[row+1][l1cfile->nbinx-1]<<"lon_gd row.."<<lon_gd[row][l1cfile->nbinx-1]<<"lon cring.."<<lonarr[ix]<<endl;

                 break;
              }
          }           
         }

         }
     }
     else //descending orb
     {
     if(Ngring==6)
         {
          latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
          latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];
          latarr[3]=lat_gd[0][0];
          latarr[4]=lat_gd[0][midix];
          latarr[5]=lat_gd[0][l1cfile->nbinx-1];

          lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
          lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];
          lonarr[3]=lon_gd[0][0];
          lonarr[4]=lon_gd[0][midix];
          lonarr[5]=lon_gd[0][l1cfile->nbinx-1];
         }
         else
         {
         latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
         latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];

         lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
         lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];

         latemp=latarr[2];     
         latemp+=20;         
         p=1;
         int rw=0;
         //west side
         while(latemp<lat_gd[0][0])
         {
         latarr[2+p]=latemp;
         cout<<"p west--"<<p<<"point# = "<<3+p<<"lat20 = "<<latemp<<endl;
         p_west[rw]=3+p;
         latemp+=20;
         p++;
         rw++;
         }   

         p--;
       
         latarr[2+p+1]=lat_gd[0][0];
         latarr[2+p+2]=lat_gd[0][midix];
         latarr[2+p+3]=lat_gd[0][l1cfile->nbinx-1];

         lonarr[2+p+1]=lon_gd[0][0];
         lonarr[2+p+2]=lon_gd[0][midix];
         lonarr[2+p+3]=lon_gd[0][l1cfile->nbinx-1];
         latemp=latarr[2+p+3];
         latemp-=20;
         p++;
       
         int c=1;
         int re=0;
         //east side
         while(latemp>lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1])
         {
         latarr[5+p]=latemp;
         cout<<"p east--"<<c<<"point# = "<<6+p<<"lat20 = "<<latemp<<endl;       
         p_east[re]=6+p;
         latemp-=20;
         p++;
         c++;
         re++;
         }

         p--;
      
         cout<<"mid points west.."<<rw<<"mid points east.."<<re<<endl;
         cout<<"# points in GRING = "<<6+p<<"# mid points west--"<<rw<<"# mid points east--"<<re<<endl;

//west side
         for(int i=0;i<rw;i++)
         {
             ix=p_west[i]-1;
             cout<<"lat:"<<latarr[ix]<<"pwest--"<<p_west[i]<<endl;
          for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {
              if(latarr[ix]<=lat_gd[row][0] && latarr[ix]>lat_gd[row+1][0])
              {
            //      cout<<"mid point west# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;            
               if(lon_gd[row][0]<0.) lontemp1=lon_gd[row][0]+360;else lontemp1=lon_gd[row][0];
               if(lon_gd[row+1][0]<0.) lontemp2=lon_gd[row+1][0]+360;else lontemp2=lon_gd[row+1][0];

               dlat_gd=abs(lat_gd[row][0]-lat_gd[row+1][0]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row+1][0]);
               dlon20=dlat20*dlon_gd/dlat_gd;

               if(lontemp1>lontemp2) lon360=lontemp2+dlon20;else lon360=lontemp2-dlon20;
               if(lon360>180) lon360=lon360-360.; 
               lonarr[ix]=lon360;
      //         cout<<"lon_gd row+1.."<<lon_gd[row+1][0]<<"lon_gd row.."<<lon_gd[row][0]<<"lon cring.."<<lonarr[ix]<<endl;
                break;
              }
          }
         }
//east side
         for(int i=0;i<re;i++)
         {
             ix=p_east[i]-1;
             cout<<"lat:"<<latarr[ix]<<"peast--"<<p_east[i]<<endl;
             for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {         
              if(latarr[ix]<=lat_gd[row][l1cfile->nbinx-1] && latarr[ix]>lat_gd[row+1][l1cfile->nbinx-1])
              {
             //     cout<<"mid point east# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;
               if(lon_gd[row][l1cfile->nbinx-1]<0.) lontemp1=lon_gd[row][l1cfile->nbinx-1]+360;else lontemp1=lon_gd[row][l1cfile->nbinx-1];
               if(lon_gd[row+1][l1cfile->nbinx-1]<0.) lontemp2=lon_gd[row+1][l1cfile->nbinx-1]+360;else lontemp2=lon_gd[row+1][l1cfile->nbinx-1];

               dlat_gd=abs(lat_gd[row][l1cfile->nbinx-1]-lat_gd[row+1][l1cfile->nbinx-1]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row][l1cfile->nbinx-1]);
               dlon20=dlat20*dlon_gd/dlat_gd;
               if(lontemp1>lontemp2) lon360=lontemp1-dlon20;else lon360=lontemp1+dlon20;
               if(lon360>180) lon360=lon360-360.;
               lonarr[ix]=lon360;
      //         cout<<"lon_gd row+1.."<<lon_gd[row+1][l1cfile->nbinx-1]<<"lon_gd row.."<<lon_gd[row][l1cfile->nbinx-1]<<"lon cring.."<<lonarr[ix]<<endl;

                 break;
              }
          }           
         }

         }//else GRING>6

     }//end descending
    
     //sequence- GRING
     if(Ngring==6) dp=6;else dp=rw+re+6;
      for(int s=0;s<dp;s++){        
         narr[s]=s+1; 
      }  

      geo_grp.putAtt("GRingPointLatitude",ncFloat,dp,latarr);
      geo_grp.putAtt("GRingPointLongitude",ncFloat,dp,lonarr);
      geo_grp.putAtt("GRingPointSequenceNo",ncInt,dp,narr);

      delete [] (latarr);
      delete [] (lonarr);
      delete [] (narr);
      delete [] (p_west);
      delete [] (p_east);
   }
   else
   {
   cout<<"ERROR EXTRACTING GRING coordinates!!-----"<<endl;
   exit(1);
   }



   nc_output->close();
return 0;    
}    


//derived from search_rc_l1c2
//same as search_SOCEA but faster with richard's addition
int L1C::search_rc_l1c5(L1C_input* l1cinput, l1c_filehandle* l1cfile, l1c_str *l1cstr,short **gdindex) {
   int32_t num_gridlines, nbinx,ic;
    short gdrow=-1,gdcol=-1;
    float pvec[3],gnvm,db,pdotgn_first,pdotgn_last;
    float **gnvec=nullptr,**cnvec=nullptr,**dc=nullptr;;
    int flag_out=-1;
    size_t pix;
    int32_t i;
    size_t inpix=0,outpix=0;
    float dotprod,dotprod2,dcm;
    
        flag_out=-1;
        num_gridlines=l1cfile->num_gridlines;
        nbinx=l1cfile->nbinx;
        num_pixels = l1cfile->npix;
        db = (l1cinput->grid_resolution)/6371/2; //Half of bin size in radians

        gnvec=allocate2d_float(3,num_gridlines);
        cnvec=allocate2d_float(3,num_gridlines);
        dc=allocate2d_float(3,num_gridlines);

        ic = nbinx/2;

        for(i=0;i<num_gridlines;i++){  
          //normal vectors for L1C rows            
             if(l1cfile->lat_gd[i][nbinx-1]>90 || l1cfile->lat_gd[i][nbinx-1]<-90 || l1cfile->lon_gd[i][nbinx-1]<-180 || l1cfile->lon_gd[i][nbinx-1]>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}
             if(l1cfile->lat_gd[i][0]>90 || l1cfile->lat_gd[i][0]<-90 || l1cfile->lon_gd[i][0]<-180 || l1cfile->lon_gd[i][0]>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}
         //compute normal vector for each rowgrid----
              gnvec[0][i] = sin(l1cfile->lon_gd[i][nbinx-1]*M_PI/180)*cos(l1cfile->lat_gd[i][nbinx-1]*M_PI/180)*sin(l1cfile->lat_gd[i][0]*M_PI/180) - sin(l1cfile->lat_gd[i][nbinx-1]*M_PI/180)*sin(l1cfile->lon_gd[i][0]*M_PI/180)*cos(l1cfile->lat_gd[i][0]*M_PI/180);
              gnvec[1][i] = sin(l1cfile->lat_gd[i][nbinx-1]*M_PI/180)*cos(l1cfile->lon_gd[i][0]*M_PI/180)*cos(l1cfile->lat_gd[i][0]*M_PI/180) - cos(l1cfile->lon_gd[i][nbinx-1]*M_PI/180)*cos(l1cfile->lat_gd[i][nbinx-1]*M_PI/180)*sin(l1cfile->lat_gd[i][0]*M_PI/180);
              gnvec[2][i] = cos(l1cfile->lon_gd[i][nbinx-1]*M_PI/180)*cos(l1cfile->lat_gd[i][nbinx-1]*M_PI/180)*sin(l1cfile->lon_gd[i][0]*M_PI/180)*cos(l1cfile->lat_gd[i][0]*M_PI/180) - sin(l1cfile->lon_gd[i][nbinx-1]*M_PI/180)*cos(l1cfile->lat_gd[i][nbinx-1]*M_PI/180)*cos(l1cfile->lon_gd[i][0]*M_PI/180)*cos(l1cfile->lat_gd[i][0]*M_PI/180);

              gnvm=sqrt(gnvec[0][i]*gnvec[0][i]+gnvec[1][i]*gnvec[1][i]+gnvec[2][i]*gnvec[2][i]);
              if(isnan(gnvm)==1){cout<<"NAN value for gnvm.."<<endl; exit(1);}
              if(gnvm==0){ cout<<"ERROR gnvm == 0--- WE CANT NORMALIZE..."<<endl;exit(1);}
              //normalization
              gnvec[0][i] = gnvec[0][i]/gnvm;
              gnvec[1][i] = gnvec[1][i]/gnvm;
              gnvec[2][i] = gnvec[2][i]/gnvm;

         // Compute normals to center columns
              cnvec[0][i] = sin(l1cfile->lon_gd[i][ic]*M_PI/180)*cos(l1cfile->lat_gd[i][ic]*M_PI/180)*gnvec[2][i] - sin(l1cfile->lat_gd[i][ic]*M_PI/180)*gnvec[1][i];
              cnvec[1][i] = sin(l1cfile->lat_gd[i][ic]*M_PI/180)*gnvec[0][i] - cos(l1cfile->lon_gd[i][ic]*M_PI/180)*cos(l1cfile->lat_gd[i][ic]*M_PI/180)*gnvec[2][i];
              cnvec[2][i] = cos(l1cfile->lon_gd[i][ic]*M_PI/180)*cos(l1cfile->lat_gd[i][ic]*M_PI/180)*gnvec[1][i] - sin(l1cfile->lon_gd[i][ic]*M_PI/180)*cos(l1cfile->lat_gd[i][ic]*M_PI/180)*gnvec[0][i];

            //; Compute grid row nadir resolution
              dc[0][i] = cos(l1cfile->lon_gd[i][ic+1]*M_PI/180)*cos(l1cfile->lat_gd[i][ic+1]*M_PI/180)- cos(l1cfile->lon_gd[i][ic]*M_PI/180)*cos(l1cfile->lat_gd[i][ic]*M_PI/180);
              dc[1][i] = sin(l1cfile->lon_gd[i][ic+1]*M_PI/180)*cos(l1cfile->lat_gd[i][ic+1]*M_PI/180)- sin(l1cfile->lon_gd[i][ic]*M_PI/180)*cos(l1cfile->lat_gd[i][ic]*M_PI/180);
              dc[2][i] = sin(l1cfile->lat_gd[i][ic+1]*M_PI/180)- sin(l1cfile->lat_gd[i][ic]*M_PI/180);

           }   
        

      //dot product
        for(pix=0;pix<num_pixels;pix++){
          for(i=0;i<num_gridlines;i++){       
             if(l1cstr->latpix[pix]>90 || l1cstr->latpix[pix]<-90 || l1cstr->lonpix[pix]<-180 || l1cstr->lonpix[pix]>180){cout<<"latpix lonpix out of the boundaries.."<<endl; exit(1);} 

            pvec[0]=cos(l1cstr->lonpix[pix]*M_PI/180)*cos(l1cstr->latpix[pix]*M_PI/180);
            pvec[1]=sin(l1cstr->lonpix[pix]*M_PI/180)*cos(l1cstr->latpix[pix]*M_PI/180);
            pvec[2]=sin(l1cstr->latpix[pix]*M_PI/180);

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
           gdindex[pix][0]=gdrow;
           gdindex[pix][1]=gdcol;

           gdrow=-1;
           gdcol=-1;
            } //end pixels


   //bmcm only checked one pixel!!I
     if(dc!=nullptr)
        delete [] (dc);
     if(gnvec!=nullptr)
        delete [] (gnvec);
     if(cnvec!=nullptr)
        delete [] (cnvec);


//     cout<<"#of pixels binned....:"<<inpix<<"pixels outside the grid..."<<outpix<<endl;

    if(flag_out==0) 
       return 0;
    else
       return 1;
}




//fred2 based on orbit vectors  pixel by pixel search L1C grid algo
bool L1C::search_rc_l1c3(L1C_input* l1cinput, l1c_filehandle* l1cfile,int32_t ydim,int32_t xdim,float latpix,float lonpix,float **lat_gd,float **lon_gd,short *erow,short *ecol) {          
    float latgd,longd,latgd2,longd2;
    int32_t num_gridlines, nbinx;
    short gdrow=-1,gdcol=-1;
    float pvec[3],gnvm,pdotgn,db,pdotgn_first=0.,pdotgn_last=0.;
    float gvec[3];
    int flag_out=-1;
    float bcm,bcm_min=10;
    int i,j;
    size_t inpix=0,outpix=0;
    float v1,v2,v3,gvec11,gvec12,gvec13,gvec21,gvec22,gvec23;

        flag_out=-1;
        num_gridlines=ydim;
        nbinx=xdim;
        num_pixels = l1cfile->npix;
        db = (l1cinput->grid_resolution)/6371/2; //Half of bin size in radians

       
         pvec[0]=cos(lonpix*M_PI/180)*cos(latpix*M_PI/180);
         pvec[1]=sin(lonpix*M_PI/180)*cos(latpix*M_PI/180);
         pvec[2]=sin(latpix*M_PI/180);
                  
        for(i=0;i<num_gridlines;i++){            
              latgd=lat_gd[i][0];
              longd=lon_gd[i][0];
              gvec11=cos(longd*M_PI/180)*cos(latgd*M_PI/180);
              gvec12=sin(longd*M_PI/180)*cos(latgd*M_PI/180);
              gvec13=sin(latgd*M_PI/180);

              if(latgd>90 || latgd<-90 || longd<-180 || longd>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}

              latgd=lat_gd[i][nbinx-1];
              longd=lon_gd[i][nbinx-1];
              gvec21=cos(longd*M_PI/180)*cos(latgd*M_PI/180);
              gvec22=sin(longd*M_PI/180)*cos(latgd*M_PI/180);
              gvec23=sin(latgd*M_PI/180);             
           
              if(latgd>90 || latgd<-90 || longd<-180 || longd>180){cout<<"lat lon out of the boundaries.."<<endl; exit(1);}


              v1 = gvec22*gvec13 - gvec23*gvec12;
              v2 = gvec23*gvec11 - gvec21*gvec13;
              v3 = gvec21*gvec12 - gvec22*gvec11;

              if(isnan(v1)==1 || isnan(v2)==1 || isnan(v3)==1 || isinf(v1)==1 || isinf(v2)==1 || isinf(v3)==1){cout<<"NAN value for gnvm.."<<endl; exit(1);}

              gnvm=sqrt(v1*v1+v2*v2+v3*v3);

              if(isnan(gnvm)==1){cout<<"NAN value for gnvm.."<<endl; exit(1);}
             
              if(gnvm==0){ cout<<"ERROR gnvm == 0--- WE CANT NORMALIZE..."<<endl;exit(1);}
              //normalization                   
              v1 = v1/gnvm;
              v2 = v2/gnvm;
              v3 = v3/gnvm;
                        
           //dot product between gridrow normal and pixel vector
              pdotgn=v1*pvec[0]+v2*pvec[1]+v3*pvec[2];

         //     cout<<"pdotgn..."<<pdotgn<<endl;
      
             if(pdotgn<=db && pdotgn>-db && gdrow<0){                 
                 gdrow=i;
                 *erow=gdrow;      
        //check column index 
              if(gdrow>=0 && gdrow<=num_gridlines-1){
                 for(j=0;j<nbinx;j++){                
                    latgd2=lat_gd[gdrow][j];
                    longd2=lon_gd[gdrow][j];
                    gvec[0]=cos(longd2*M_PI/180)*cos(latgd2*M_PI/180);
                    gvec[1]=sin(longd2*M_PI/180)*cos(latgd2*M_PI/180);
                    gvec[2]=sin(latgd2*M_PI/180);
                    bcm=sqrt((pvec[0]-gvec[0])*(pvec[0]-gvec[0])+(pvec[1]-gvec[1])*(pvec[1]-gvec[1])+(pvec[2]-gvec[2])*(pvec[2]-gvec[2]));
                   if(bcm<bcm_min){ 
                       bcm_min=bcm;
                       gdcol=j;
                       *ecol=gdcol;
                      }
                  }

                }
               }
                      
              if(i==0){
                 pdotgn_first=pdotgn;//first gridline
                 }
              if(i==num_gridlines-1){
                  pdotgn_last=pdotgn;//last line
                 }
      
                 
             }//end gridlines

            //flag outside L1C grid pixels
              if(pdotgn_first>db && pdotgn_last<-db && gdrow<0 && gdcol<0){
                 flag_out=1;           
                 outpix++; 
           //      cout<<"************************* pixel out of the L1C grid ****************************"<<endl;
           //      cout<<"inside search row index.."<<i<<"col index...."<<j<<endl;
              }  
              else if(gdrow>=0 && gdcol>=0) {
                  flag_out=0;               
              //    cout<<"pix..inside the grid."<<pix+1<<"row index.."<<gdrow<<"gdcol.."<<gdcol<<endl;
                   inpix++;                   
                  if(gdcol==-1){
                     cout<<"gdcol==-1..."<<"bcm_min.."<<bcm_min<<endl;
                     }
                }
              else {//whole line outside, big dotproduct
               outpix++;
           //    cout<<"************************* pixel out of the L1C grid ****************************"<<endl;
               flag_out=-1;
              }

    if(flag_out==0) 
       return true;
    else
       return false;


    }



bool L1C::sbs2_l1c4(L1C_input* l1cinput, int32_t ydim, int32_t xdim, float** alat, short** alat_index, float latpix, float lonpix, float **lat_gd,float** lon_gd, short* rowindex, short* colindex)
    {
         //mat is is the sorted matrix in ascending order

        short i = ydim - 1, j = 0, erow, ecol;  //set indexes for bottom left element 
        float binres, diflat,diflat2;
        float latpos;
        double dtogd,dtogd2;
        float lat1,lat2,lon1,lon2,az1,dlambda,hyp;
        float co2,ca2,h2,g_bearing;


        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());

        binres = (l1cinput->grid_resolution) * 1000; //L1C grid resolution in meters
        hyp=sqrt(0.5*binres*0.5*binres+0.5*binres*0.5*binres);

        while ((i - 1) >= 0 && j < xdim)
        {
            latpos = latpix + 90;//convert -90,+90 to lat 0-180


            if (latpos >= alat[i - 1][j] && latpos <= alat[i][j]) {


                diflat = abs(latpos - alat[i][j]);
                diflat2 = abs(latpos - alat[i - 1][j]);
                //check lat distance to gridpoint-
                dtogd = diflat * 111.321 * 1000;
                dtogd2 = diflat2 * 111.321 * 1000;


                if (dtogd < dtogd2 && dtogd <= binres/2) {
                    erow = alat_index[i][j];
                        }
                else if (dtogd2 < dtogd && dtogd2 <= binres/2) {
                    erow = alat_index[i - 1][j];
                        }

                lat1 = (alat[i-1][j]-90) * M_PI / 180.;
                lat2 = (alat[i][j]-90) * M_PI / 180.;
                lon1 = lon_gd[i-1][j] * M_PI / 180.;
                lon2 = lon_gd[i][j] * M_PI / 180.;
                dlambda = (lon2 - lon1);

                g_bearing = (float)atan2(sin(dlambda) * cos(lat2), cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlambda));

                geod.Inverse(latpix,lonpix,alat[erow][j]-90,lon_gd[erow][j],dtogd);

 //compute azimuth for grid cells bracketing

                lat1 = latpix * M_PI / 180.;
                lat2 = (alat[erow][j]-90) * M_PI / 180.;
                lon1 = lonpix* M_PI / 180.;
                lon2 = lon_gd[erow][j] * M_PI / 180.;
                dlambda = (lon2 - lon1);
                az1 = (float)atan2(sin(dlambda) * cos(lat2), cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlambda));
                geod.Inverse(latpix,lonpix,alat[erow][j]-90,lon_gd[erow][j],dtogd);


                az1+=g_bearing;
      //CHECK QUADRANTS
         //1
                if(az1<=0.5*M_PI && az1>=0){
                 ca2=hyp*cos(az1);
                 co2=hyp*sin(az1);
                 h2=sqrt(co2*co2+ca2*ca2);
                 }
       //2   
                else if(az1>0.5*M_PI){
                 az1=az1-0.5*M_PI;
                 ca2=-hyp*sin(az1);
                 co2=hyp*cos(az1);
                 h2=sqrt(co2*co2+ca2*ca2);
                 }
                //3
                else if(az1<-0.5*M_PI){
                 az1=0.5*M_PI-(az1+M_PI);
                 ca2=-hyp*sin(az1);//
                 co2=-hyp*cos(az1);//
                 h2=sqrt(co2*co2+ca2*ca2);
                 }
         //4
                else if(az1>=-0.5*M_PI && az1<0){
                 ca2=hyp*cos(az1);//
                 co2=-hyp*sin(az1);//
                 h2=sqrt(co2*co2+ca2*ca2);
                  }


                if (dtogd < h2) {

                    ecol = j;
                    *rowindex = erow;
                    *colindex = ecol;

                    return true;
                }
                else {
                  j++;

                   }

            }

            if (alat[i][j] > latpos && alat[i - 1][j] > latpos) {
                i--;

            }
            if (i == 0) return false;

            if (alat[i][j] < latpos && alat[i - 1][j] < latpos) {
                j++;

                if (j == xdim) return false;
            }
            if (i == 0 || j == xdim) return false;

        }

        return false;
    }



    //saddleback search algo improved 
    //this approach is actually based on a radius search distance rather than a polygon pixel shape. This must be improved 
    bool L1C::sbs2_l1c2(L1C_input* l1cinput, int32_t ydim, int32_t xdim, float** alat, short** alat_index, float latpix, float lonpix, float** lon_gd, short* rowindex, short* colindex)
    {
        //pixval can be latpix or lonpix
        //flag_coor 0 latitude, 1 longitude
         //mat is is the sorted matrix in ascending order
         //mat_index IS THE ORIGINAL INDEX of lat or lon values L1C gridpoint 
         //this INDEXING MAY CHANGE BETWEEN COLUMNS OR ACROSS GRIDLINES!!! 

        short i = ydim - 1, j = 0, erow, ecol;  //set indexes for bottom left element 
        float diflon, diflat, diflat2, dtogd, dtogd2, binres, londegkm;
        float latpos, lonpos, lonpos_gd;


        binres = (l1cinput->grid_resolution) * 1000; //L1C grid resolution in meters


        while ((i - 1) >= 0 && j < xdim)
        {
            latpos = latpix + 90;//convert -90,+90 to lat 0-180

            if (latpos >= alat[i - 1][j] && latpos <= alat[i][j]) {
                     //check coor difference with respect to center gridpoint--
                diflat = abs(latpos - alat[i][j]);
                diflat2 = abs(latpos - alat[i - 1][j]);
                //check lat distance to gridpoint-
                dtogd = diflat * 111.321 * 1000;
                dtogd2 = diflat2 * 111.321 * 1000;

                if (dtogd < dtogd2 && dtogd <= binres / 2) {
                    erow = alat_index[i][j];                
                        }
                else if (dtogd2 < dtogd && dtogd2 <= binres / 2) {
                    erow = alat_index[i - 1][j];                
                        }                       
                //
        //         if(erow>=0 && erow<ydim){
                   if (lonpix < 0.) lonpos = lonpix + 360;else lonpos = lonpix;//longitude between 0 and 360
                   if (lon_gd[erow][j] < 0.) lonpos_gd = lon_gd[erow][j] + 360;else lonpos_gd = lon_gd[erow][j];
       //             }


      //          cout<<"lonpos..pix."<<lonpos<<"<lonpos_gd[erow][j]..."<<lonpos_gd<<endl;

                diflon = abs(lonpos - lonpos_gd);
                londegkm = cos(latpix * degrad) * 111.321;
                dtogd = diflon * londegkm * 1000;

                if (dtogd <= binres / 2) {
                    ecol = j;
                    *rowindex = erow;
                    *colindex = ecol;
/*                    cout<<"rowindex.."<<erow<<"colindex.."<<ecol<<endl;
                    cout<<"lonpos..pix."<<lonpos<<"<lonpos_gd[erow][j]..."<<lonpos_gd<<endl;
                    cout<<"row#............................"<<erow<<"latpos..pix."<<latpos<<"alat[i - 1][j]..."<<alat[i - 1][j]<<"index.."<<i<<"col.."<<j+1<<endl;
                    exit(1);
*/
                    return true;
                }
                else {
                  j++;//go to next col

                   }
             dtogd=binres;
            } //end if lat bracket  

            if (alat[i][j] > latpos && alat[i - 1][j] > latpos) { //pixval < mat
                i--;

            }
            if (i == 0) return false;

            if (alat[i][j] < latpos && alat[i - 1][j] < latpos) {
                j++;

                if (j == xdim) return false;
            }
            if (i == 0 || j == xdim) return false;

            //the sorting is done over the rows for each column so row order is actually changing!!

        }//end while

        return false;
    }


    //saddleback search algo improved 
    //this approach is actually based on a radius search distance rather than a polygon pixel shape. This must be improved 
    bool L1C::sbs2_l1c(L1C_input* l1cinput, int32_t ydim, int32_t xdim, float** alat, short** alat_index, float latpix, float lonpix, float** lon_gd, short* rowindex, short* colindex)
    {
        //pixval can be latpix or lonpix
        //flag_coor 0 latitude, 1 longitude
         //mat is is the sorted matrix in ascending order
         //mat_index IS THE ORIGINAL INDEX of lat or lon values L1C gridpoint 
         //this INDEXING MAY CHANGE BETWEEN COLUMNS OR ACROSS GRIDLINES!!! 

        short i = ydim - 1, j = 0, erow, ecol;  //set indexes for bottom left element 
        float diflon, diflat, diflat2, dtogd, dtogd2, binres, londegkm;
        float latpos, lonpos, lonpos_gd;


        binres = (l1cinput->grid_resolution) * 1000; //L1C grid resolution in meters

        while ((i - 1) >= 0 && j < xdim)
        {
            latpos = latpix + 90;//convert -90,+90 to lat 0-180

            if (latpos >= alat[i - 1][j] && latpos <= alat[i][j]) {
                     //check coor difference with respect to center gridpoint--
                diflat = abs(latpos - alat[i][j]);
                diflat2 = abs(latpos - alat[i - 1][j]);
                //check lat distance to gridpoint-
                dtogd = diflat * 111.321 * 1000;
                dtogd2 = diflat2 * 111.321 * 1000;

                if (dtogd < dtogd2 && dtogd <= binres / 2) {
                    erow = alat_index[i][j];                
                        }
                else if (dtogd2 < dtogd && dtogd2 <= binres / 2) {
                    erow = alat_index[i - 1][j];                
                        }        

                //check longitude estimate
         //       if(erow>=0 && erow<ydim){
                   if (lonpix < 0.) lonpos = lonpix + 180;else lonpos = lonpix;//longitude between 0 and 360
                   if (lon_gd[erow][j] < 0.) lonpos_gd = lon_gd[erow][j] + 180;else lonpos_gd = lon_gd[erow][j];
        //         }

  

                diflon = abs(lonpos - lonpos_gd);
                londegkm = cos(latpix * degrad) * 111.321;
                dtogd = diflon * londegkm * 1000;

     //           cout<<"dtogd..."<<dtogd<<endl;

                if (dtogd <= binres / 2) {
                    ecol = j;
                    *rowindex = erow;
                    *colindex = ecol;

                    return true;
                }
                else j++;//go to next colI

     //          if ((i - 1) < 0 || j >= xdim || j<0) return false;

            } //end if lat bracket  


      
             if (alat[i][j] > latpos && alat[i - 1][j] > latpos) { //pixval < mat
                i--;
               }
     
            if (i == 0) return false;

    
               if (alat[i][j] < latpos && alat[i - 1][j] < latpos) {
                j++;

                if (j == xdim) return false;
              }
   

            if (i == 0 || j == xdim) return false;

            //the sorting is done over the rows for each column so row order is actually changing!!

        }//end while


        
        return false;
    }


//version after Fred proj
  int32_t L1C::openL1Cgrid3(int swtd,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t* swtd_id,int16_t* file_id, int16_t* nfiles_swt){
        size_t ybins,xbins;
        int32_t num_gridlines, nbinx;//number of gridlines to be processed     
        std::string str;
        const char* ptstr;
        int status,ncid_grid,x_dimid,y_dimid;   
        int32_t NY = -1, NX = -1;
        char* ifile_char;
        std::string ifile_str;
        string gridname, azeast_name;
        int geoGrp;    
        char tmp[256];
        string tswt_ini_file;         
        size_t att_len; 
        std::string fname_out, pathstr, senstr, monstr, daystr, yearstr, prodstr, gdstr, swtstr, swtnum, extstr, granstr,timestr,azstr,missionstr,ofilestr;
        int proj_type=-1;

         //open file
        ifile_str = l1cinput->files[0];
        ifile_char = &ifile_str[0];
        file_format format = getFormat(ifile_char);   

         
        num_pixels = l1cfile->npix;
        num_scans=l1cfile->nscan;
        l1cfile->sensor=l1cinput->sensor;//SPEX 1, OCI 2 and HARP 3

        if(l1cinput->sensor==34 || format.type==FT_SPEXONE){
             senstr="SPEXone";
             l1cfile->nbinx=25;
              }
        else if (l1cinput->sensor==30 || format.type==FT_OCIS){
            senstr="OCI";
            l1cfile->nbinx=519;
             }
        else if (l1cinput->sensor==35 || format.type==FT_HARP2){
            senstr="HARP2";
            l1cfile->nbinx=457;
            }
        else{cout<<"sensor by default is OCI option 2....."<<endl;
            senstr="OCI";
             l1cfile->nbinx=519;
          }


    //    num_frames=l1cfile->nframes;

        cout<<"Opening L1C grid........................................................................."<<endl;

        getcwd(tmp, 256);
        pathstr=tmp;

//open lat/lon L1C grid --------------------------------------
        if(format.type == FT_OCIS){
           int ix=ifile_str.find("PACE");
           tswt_ini_file=ifile_str.substr(ix+13,15);

           l1cfile->tswt_ini_file=tswt_ini_file;
           }
        else if(format.type == FT_HARP2){
          int ix=ifile_str.find("PACE");
           tswt_ini_file=ifile_str.substr(ix+52,15);
           l1cfile->tswt_ini_file=tswt_ini_file;
        }
        else if(format.type == FT_SPEXONE){
           tswt_ini_file="18991231T112700";
           l1cfile->tswt_ini_file=tswt_ini_file;
        }

        proj_type=l1cfile->proj_type;
        if(proj_type==0){//socea
        ofilestr="PACE_"+senstr+"."+tswt_ini_file+".L1C.5.2km.nc";

        gridname = pathstr + "/out/" +ofilestr;
        }
        if(proj_type==1){ //fixed bearing projection ----
          gridname="/accounts/mamontes/images/OCIS/sean/out/FB_L1Cgrid.nc";
          }
        
        cout<<"gridname.."<<gridname<<endl;
        ptstr = gridname.c_str();
   
        status = nc_open(ptstr, NC_NOWRITE, &ncid_grid);
          if (status != NC_NOERR) {
              fprintf(stderr, "nc_open error.\n");
               exit(EXIT_FAILURE);
          }

        status = nc_inq_attlen(ncid_grid, NC_GLOBAL, "instrument", &att_len);
            check_err(status, __LINE__, __FILE__);
        char* in_str = (char*)malloc(att_len + 1); // + 1 for trailing null 
        in_str[att_len] = '\0';
        status = nc_get_att_text(ncid_grid, NC_GLOBAL, "instrument", in_str);
        string senstr2(in_str);  
        status = nc_inq_attlen(ncid_grid, NC_GLOBAL, "product_name", &att_len);
            check_err(status, __LINE__, __FILE__);
        char* fn_str = (char*)malloc(att_len + 1); // + 1 for trailing null 
        fn_str[att_len] = '\0';
        status = nc_get_att_text(ncid_grid, NC_GLOBAL, "product_name", fn_str);
        string tswt_ini_file2(fn_str);

        status = nc_inq_attlen(ncid_grid, NC_GLOBAL, "time_coverage_start", &att_len);
            check_err(status, __LINE__, __FILE__);
        char* time_str = (char*)malloc(att_len + 1); // + 1 for trailing null 
        time_str[att_len] = '\0';
        status = nc_get_att_text(ncid_grid, NC_GLOBAL, "time_coverage_start", time_str);
        string tswt_ini(time_str);

        status = nc_inq_attlen(ncid_grid, NC_GLOBAL, "time_coverage_end", &att_len);
            check_err(status, __LINE__, __FILE__);
        char* time_str2 = (char*)malloc(att_len + 1); // + 1 for trailing null
        time_str2[att_len] = '\0';
        status = nc_get_att_text(ncid_grid, NC_GLOBAL, "time_coverage_end", time_str2);
        string tswt_end(time_str2);

        status = nc_inq_attlen(ncid_grid, NC_GLOBAL, "startDirection", &att_len);
            check_err(status, __LINE__, __FILE__);
        char* sdir = (char*)malloc(att_len + 1); // + 1 for trailing null
        sdir[att_len] = '\0';
        status = nc_get_att_text(ncid_grid, NC_GLOBAL, "startDirection", sdir);
        string start_dir(sdir);

        status = nc_inq_attlen(ncid_grid, NC_GLOBAL, "endDirection", &att_len);
            check_err(status, __LINE__, __FILE__);
        char* edir = (char*)malloc(att_len + 1); // + 1 for trailing null
        edir[att_len] = '\0';
        status = nc_get_att_text(ncid_grid, NC_GLOBAL, "endDirection", edir);
        string end_dir(edir);

        status = nc_inq_attlen(ncid_grid, NC_GLOBAL, "bin_size_at_nadir", &att_len);
            check_err(status, __LINE__, __FILE__);
        char* binsize= (char*)malloc(att_len + 1); // + 1 for trailing null
        binsize[att_len] = '\0';
        status = nc_get_att_text(ncid_grid, NC_GLOBAL, "bin_size_at_nadir", binsize);
        string binstr(binsize);


        l1cfile->tswt_ini=tswt_ini;
    //    l1cfile->tswt_ini_file=tswt_ini_file2;
        l1cfile->tswt_end=tswt_end;
        l1cfile->sens_type=senstr2;
        l1cfile->start_dir=start_dir;
        l1cfile->end_dir=end_dir;
        l1cfile->binstr=binstr;
       

        delete [](time_str);
        delete [](time_str2);

       //groups
        status = nc_inq_grp_ncid(ncid_grid, "geolocation_data", &geoGrp);
        check_err(status, __LINE__, __FILE__);
       
        status = nc_inq_dimid(ncid_grid, "bins_along_track", &y_dimid);
            check_err(status, __LINE__, __FILE__);
        nc_inq_dimlen(ncid_grid, y_dimid, &ybins);
        status = nc_inq_dimid(ncid_grid, "bins_across_track", &x_dimid);
            check_err(status, __LINE__, __FILE__);

        nc_inq_dimlen(ncid_grid, x_dimid, &xbins);
        nc_inq_dimlen(ncid_grid, y_dimid, &ybins);

         //************************************************************
        //mem allocation for geolocation pointers---
        num_gridlines=ybins;
        nbinx=xbins;

        cout<<"num_gridlines..."<<num_gridlines<<"nbinx....."<<nbinx<<endl;

        l1cfile->lat_gd = allocate2d_float(num_gridlines, nbinx);
        l1cfile->lon_gd = allocate2d_float(num_gridlines, nbinx);
        l1cfile->alt_gd = allocate2d_float(num_gridlines, nbinx);
        l1cfile->index_xy = allocate2d_short(num_gridlines, nbinx);
        l1cfile->lat_asort = allocate2d_float(num_gridlines, nbinx);


        status = nc_inq_varid(geoGrp, "latitude", &latId);//scans x velements
                check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGrp, "longitude", &lonId);//scans x velements
                check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGrp, "altitude", &altId);//scans x velements
                check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(geoGrp, latId, &l1cfile->lat_gd[0][0]);
                check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(geoGrp, lonId, &l1cfile->lon_gd[0][0]);
                check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(geoGrp, altId, &l1cfile->alt_gd[0][0]);
                check_err(status, __LINE__, __FILE__);
//close file
        if ((status = nc_close(ncid_grid)))
           check_err(status, __LINE__, __FILE__);

       for(int i=0;i<num_gridlines;i++){
            for(int i=0;i<nbinx;i++){
            }

       }



        //recompute number of gridlines based on lat limits -80 to 80 degrees  
       //asc orbit goes from negative to positive latitude....
        NY = num_gridlines;
        NX = nbinx;
        l1cfile->num_gridlines=NY;
        l1cfile->nbinx=NX;


        //*********** SORTING LAT/LON ASCENDING AND TRACK INDEXES ****************************
        //***********************************************************************************

        //Create index matrix for i and j of each lat and lon L1C grid (-80 to 80)
        //num_gridlines x nbinx

        //define 1-D vector to sort 1 col L1C grid at the time---
        vector<pair<float, int> > vp;

        //fill vector with lat_gd column--
        cout << "*********** SORTING LAT/LON GRIDPOINTS AND TRACKING INDEXES ********************" << endl;
        

        for (int col = 0;col < NX;col++) {
            for (int row = 0;row < NY;row++) {
                vp.push_back(make_pair(l1cfile->lat_gd[row][col] + 90., row));
            }
            //sort 
            stable_sort(vp.begin(), vp.end());

            for (unsigned int i = 0; i < vp.size(); i++) {
                l1cfile->lat_asort[i][col] = vp[i].first;
                l1cfile->index_xy[i][col] = vp[i].second;

                 }
            vp.clear();
        }
       
        cout<<"ok after sorting.............................."<<endl;

     
   return 0;
   }



//OPEN L1C grid
   int32_t L1C::openL1Cgrid(l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput){ 
   const char* ptstr;   
   std::string ifile_str,pathstr;
   string gridname;
  
   cout<<"Opening L1C grid........................................................................."<<endl;

   pathstr="";
   gridname = pathstr+l1cinput->l1c_grid;  
   cout<<"gridname.."<<gridname<<endl;


   ptstr=gridname.c_str();

   string l1c_str=l1cinput->l1c_grid;

   NcFile* nc_l1cgrid;
   try {
        nc_l1cgrid = new NcFile(ptstr, NcFile::read);
        }
   catch (NcException& e) {
          e.what();
          cerr << "l1cgen l1c_pflag= 4:: Failure reading L1C grid: "
          + l1c_str << endl;
          exit(1);
   }
  
   NcDim yd=nc_l1cgrid->getDim("bins_along_track");
   l1cfile->num_gridlines=yd.getSize();
   NcDim xd=nc_l1cgrid->getDim("bins_across_track");
   l1cfile->nbinx=xd.getSize();

   l1cfile->lat_gd = allocate2d_float(l1cfile->num_gridlines, l1cfile->nbinx);
   l1cfile->lon_gd = allocate2d_float(l1cfile->num_gridlines, l1cfile->nbinx);
   l1cfile->alt_gd = allocate2d_float(l1cfile->num_gridlines, l1cfile->nbinx);
   l1cfile->index_xy = allocate2d_short(l1cfile->num_gridlines, l1cfile->nbinx);
   l1cfile->lat_asort = allocate2d_float(l1cfile->num_gridlines, l1cfile->nbinx);

   cout<<"num_gridlines"<<l1cfile->num_gridlines<<"nbinx"<<l1cfile->nbinx<<endl;
   NcGroup geo_grp=nc_l1cgrid->getGroup("geolocation_data");
   NcVar v1=geo_grp.getVar("latitude");
   v1.getVar(&l1cfile->lat_gd[0][0]);
   v1=geo_grp.getVar("longitude");
   v1.getVar(&l1cfile->lon_gd[0][0]);
   v1=geo_grp.getVar("altitude");
   v1.getVar(&l1cfile->alt_gd[0][0]);

//close file
   nc_l1cgrid->close();

        //*********** SORTING LAT/LON ASCENDING AND TRACK INDEXES ****************************
        //***********************************************************************************

        //Create index matrix for i and j of each lat and lon L1C grid (-80 to 80)
        //num_gridlines x nbinx

        //define 1-D vector to sort 1 col L1C grid at the time---
   vector<pair<float, int> > vp;

        //fill vector with lat_gd column--
   cout << "*********** SORTING LAT/LON GRIDPOINTS AND TRACKING INDEXES ********************" << endl;
        

   for (int col = 0;col <l1cfile->nbinx;col++)
   {
            for (int row = 0;row < l1cfile->num_gridlines;row++)
            {
                vp.push_back(make_pair(l1cfile->lat_gd[row][col] + 90., row));
            }
            //sort 
            stable_sort(vp.begin(), vp.end());

            for (unsigned int i = 0; i < vp.size(); i++)
            {
                l1cfile->lat_asort[i][col] = vp[i].first;
                l1cfile->index_xy[i][col] = vp[i].second;
            }
            vp.clear();
   }
       
   cout<<"ok after sorting.............................."<<endl;

   return 0;
   }



//same as :binL1C_sbs_line but with Fred, no need of az east for ebaring proj
//this version also works with SPEX and HARP2----
//
int32_t L1C::binL1C_sbs_line3(int swtd,L1C *l1c,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,int16_t* swtd_id,int16_t* file_id, int16_t* nfiles_swt,float ****binmean_Lt,int ****bincount,float ****binmean_Lt_pol,int ****bincount_pol,size_t sline,int granid){ 

        int32_t num_gridlines, nbinx;//number of gridlines to be processed
        std::string str,ofilestr;
        int status, ncid_out,b_dimid,b_dimid_pol,v_dimid,x_dimid,y_dimid,varid1,varid2,varid_v1,varid_v2,varid_v3,varid_count,varid_lt,varid_lt_pol,NDIMS,NDIMS2,NDIMS3,NDIMS4;
        int32_t bin_xpix, bin_ypix;
        int dimids[2],dimids2[3],dimids3[4],*dimids4,dimids5[2],dimids6[4],dimids7[2];
        const char* filename_lt;
        int32_t NY = -1, NX = -1, NY1 = -1, NY2 = -1;
        float****data_out2=nullptr,****data_out2_pol=nullptr,** lat_out=nullptr, ** lon_out=nullptr,*view_out1=nullptr,**view_out2=nullptr,**view_out2_pol=nullptr;;
        int ***data_out=nullptr,***data_out_pol=nullptr;;
        int flag_inpix;//first index file id, second index line #
        bool boolbin1;
        char* ifile_char;
        std::string ifile_str;    
        string gridname;
        int grp_obs,grp_coor,grp_views;
        size_t  sb=0,bb = 0, rb = 0, swb = 0,ib=0, NVIEWS=-1,NBANDS=-1,NBANDS_POL=-1;
        int view=0;
        short gd_row, gd_col;
        std::string fname_out, pathstr, senstr,produstr;
        std::string  tswt_ini, tswt_end, tswt_ini_file2,senstr2,start_dir,end_dir,binstr,nadpix_str;
        float fill_value=-999.;
        int shuffle, deflate, deflate_level;
        shuffle = NC_SHUFFLE;
        deflate = 1;
        deflate_level = 5;      

        short **gdindex=nullptr;
         

        //nadir pixel ofr OCI        
        int nadpix = l1cfile->nadpix;
        nadpix_str=to_string(nadpix);        
  //      cout<<"nadpix index......"<<nadpix<<endl;
  
        
   //metdatada------
        tswt_ini=l1cfile->tswt_ini;
        tswt_ini_file2=l1cfile->tswt_ini_file;
        tswt_end=l1cfile->tswt_end;
        senstr2=l1cfile->sens_type;
        start_dir=l1cfile->start_dir;
        end_dir=l1cfile->end_dir;
        binstr=l1cfile->binstr;    

        num_gridlines = l1cfile->num_gridlines;
        nbinx = l1cfile->nbinx;
        num_pixels = l1cfile->npix;
        num_scans=l1cfile->nscan;

        ifile_str = l1cinput->files[0];
        ifile_char = &ifile_str[0];
        file_format format = getFormat(ifile_char);

        if(l1cinput->sort_method==0){//orbital Fred binning sorting method
      //sorting pixels into L1C grid
           gdindex = allocate2d_short(num_pixels,2);

    //Fred SOCEA        
           search_rc_l1c5(l1cinput,l1cfile,l1cstr,gdindex);


           }
    //    num_frames=l1cfile->nframes;
       
                    //allocate mem for lat/lon/azeast pointers---
        //determine # of gd groups to be processede
        //gdlines_group=num_gridlines;//ONE BIG GROUP

        //--NORMALIZED  gridline longitude is -180 to 180 degrees before comparison with longitude extracted from image pixels---
        for (int i = 0; i < num_gridlines; i++) {
            for (int j = 0; j < nbinx; j++) {
                if (l1cfile->lon_gd[i][j] < -180.)l1cfile->lon_gd[i][j] += 360.;
                if (l1cfile->lon_gd[i][j] > 180.)l1cfile->lon_gd[i][j] -= 360.;

               }
            }


//********* LOOP FILES ************************************

                //radiance variables
                if(format.type==FT_SPEXONE){
                   NVIEWS=l1cfile->n_views;//should be 5 but only provided 1
                   NBANDS=l1cfile->nband_view;
                   NBANDS_POL=l1cfile->npol_band_view;
                   cout << "#total number of Intensity bands per view.." << NBANDS << endl;
                   cout << "#total number of Polarization bands per view.." << NBANDS_POL << endl;
                }
                else if(format.type==FT_HARP){
                   NVIEWS= 1;//HARP-2 10 but 60 for 669 nm
                   NBANDS=1;//
                   NBANDS_POL=1;
                   cout << "#total number of Intensity bands per view.." << NBANDS << endl;
                   cout << "#total number of Polarization bands per view.." << NBANDS_POL << endl;
                }
                else if(format.type==FT_HARP2){
                   NVIEWS=l1cfile->n_views;//HARP-2 10 but 60 for 669 nm
                   NBANDS=4;//
                   NBANDS_POL=4;
                   cout << "#total number of Intensity bands per view.." << NBANDS << endl;
                   cout << "#total number of Polarization bands per view.." << NBANDS_POL << endl;
                }
                else{//OCIS
                   NVIEWS= l1cfile->n_views;//for OCI
                   NBANDS=l1cfile->nbands;//for OCI
                   NBANDS_POL=0;
                   num_blue_bands=l1cfile->nband_blue;
                   num_red_bands=l1cfile->nband_red;
                   num_SWIR_bands=l1cfile->nband_swir;
                }

   //*********** BIG LOOP M_PIXEL *****************************************************************
                     for (unsigned int pix = 0;pix < num_pixels;pix++) {
                            sb=0,bb=0,rb=0,swb=0,ib=0;
                            gd_row=-1,gd_col=-1;
                            flag_inpix = 0;


                          if(format.type==FT_HARP2){ //with only view index 4 and blue group band              
                            boolbin1=search_rc_l1c3(l1cinput,l1cfile,num_gridlines, nbinx,l1cstr->latpix_3d[4][pix],l1cstr->lonpix_3d[4][pix],l1cfile->lat_gd,l1cfile->lon_gd,&gd_row,&gd_col);



                          }
                          else{//OCIS
                         
                              //SBS method
                               if(l1cinput->sort_method==1){                                  
                           //       boolbin1 = sbs2_l1c2(l1cinput, num_gridlines, nbinx, l1cfile->lat_asort, l1cfile->index_xy, l1cstr->latpix[pix], l1cstr->lonpix[pix], l1cfile->lon_gd, &gd_row, &gd_col);
                                   boolbin1 = sbs2_l1c4(l1cinput, num_gridlines, nbinx, l1cfile->lat_asort, l1cfile->index_xy, l1cstr->latpix[pix], l1cstr->lonpix[pix],l1cfile->lat_gd, l1cfile->lon_gd, &gd_row, &gd_col);
                      
                                  }
                     
                          //     FRED2 method
                               
                 //              boolbin1=search_rc_l1c3(l1cinput,l1cfile,num_gridlines, nbinx,l1cstr->latpix[pix],l1cstr->lonpix[pix],l1cfile->lat_gd,l1cfile->lon_gd,&gd_row,&gd_col);
                               if(l1cinput->sort_method==0){                         
                                    gd_row=gdindex[pix][0];
                                    gd_col=gdindex[pix][1];
                                  }

            //                   cout<<"boolbin1.."<<boolbin1<<"pix..."<<pix<<"gdrow.."<<gd_row<<"gdcol.."<<gd_col<<endl;                                                      
                          }


                            bin_ypix = gd_row + 1;
                            bin_xpix = gd_col + 1;
                //            cout<<"boolbin1.."<<boolbin1<<"bin_ypix.."<<bin_ypix<<"bin_xpix.."<<bin_xpix<<endl;


                            if (gd_row>=0 && gd_col>=0)
                                boolbin1=1;
                            else boolbin1=0;

                            if (boolbin1 == 1) flag_inpix = 1;
                            //************ BINNING ***************************
                         //assign identified pixel to Lt and bin stat arrays-----
                         //
                         float lontemp,latemp;
                         if(format.type==FT_HARP2){
                             latemp=l1cstr->latpix_3d[4][pix];
                             lontemp=l1cstr->lonpix_3d[4][pix];
                            }
                         else{
                            latemp=l1cstr->latpix[pix];
                            lontemp=l1cstr->lonpix[pix]; 
                            }

        
                            if (flag_inpix == 1 && bin_xpix >= 1 && bin_ypix >= 1 && bin_xpix <= nbinx && bin_ypix <= num_gridlines && latemp>=l1cinput->south && latemp<=l1cinput->north && lontemp>=l1cinput->west && lontemp<=l1cinput->east) {

                                INPIX++;
                              if(format.type==FT_SPEXONE){
                                  view=l1cstr->viewport[pix]%8;
                                   while (ib<NBANDS) {
                                      if (l1cstr->I[pix][ib] > 0.) {
                                              binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] + l1cstr->I[pix][ib];
                                              bincount[view][ib][bin_ypix - 1][bin_xpix - 1] += 1;

                                       //      cout<<"I[pix][ib]..."<<l1cstr->I[pix][ib]<<"bin_ypix.."<<bin_ypix<<"bin_xpix.."<<bin_xpix<<endl;
                                       //      cout<<"pix..."<<pix+1<<"latpix..."<<l1cstr->latpix[pix]<<"lonpix..."<<l1cstr->lonpix[pix]<<"sensor azimuth pixel.."<<l1cstr->senazpix[pix]<<endl;
                                           }
                                        ib++;                                    
                                     }//end while intensity bands

                                   ib=0;
                                   while (ib<NBANDS_POL) {
                                      if (l1cstr->I_polsample[pix][ib] > 0.) {
                                              binmean_Lt_pol[view][ib][bin_ypix - 1][bin_xpix - 1] = binmean_Lt_pol[view][ib][bin_ypix - 1][bin_xpix - 1] + l1cstr->I_polsample[pix][ib];
                                              bincount_pol[view][ib][bin_ypix - 1][bin_xpix - 1] += 1;

                                       //      cout<<"I[pix][ib]..."<<l1cstr->I[pix][ib]<<"bin_ypix.."<<bin_ypix<<"bin_xpix.."<<bin_xpix<<endl;
                                       //      cout<<"pix..."<<pix+1<<"latpix..."<<l1cstr->latpix[pix]<<"lonpix..."<<l1cstr->lonpix[pix]<<"sensor azimuth pixel.."<<l1cstr->senazpix[pix]<<endl;
                                           }
                                        ib++;
                                    }
                                 }
                              
                               else if(format.type==FT_HARP){
                                  view=0,sb=0;
                                  if (l1cstr->Lt[pix][ib] > 0.) {
                                              binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][ib][bin_ypix - 1][bin_xpix - 1] + l1cstr->Lt[pix][ib];
                                              bincount[view][ib][bin_ypix - 1][bin_xpix - 1] += 1;
                                           }
                                    ib++;
                                    sb++;
                                 }
                               else if(format.type==FT_HARP2){
                                  view=4,sb=0;
                                  if (l1cstr->I[4][pix] > 0.) {
                                              binmean_Lt[view][0][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][0][bin_ypix - 1][bin_xpix - 1] + l1cstr->I[4][pix];
                                              bincount[view][0][bin_ypix - 1][bin_xpix - 1] += 1;
                                           }
                                    ib++;
                                    sb++;
                                 }
                               else{ //OCIS
        //BLUE---                                                 
                               if(l1cstr->senazpix[nadpix]/100<0) view=0; else view=1;
                                while (sb<num_blue_bands) {
                                      if (l1cstr->Lt_blue[bb][pix] > 0.) {
                                        //      cout<<"Lt_pix[pix][ib]..."<<l1cstr->Lt_blue[bb][pix]<<"bin_ypix.."<<bin_ypix<<"bin_xpix.."<<bin_xpix<<endl;
                                        //      cout<<"pix..."<<pix+1<<"latpix..."<<l1cstr->latpix[pix]<<"lonpix..."<<l1cstr->lonpix[pix]<<"sensor azimuth pixel.."<<l1cstr->senazpix[pix]/100<<endl;

                                              binmean_Lt[view][sb][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][sb][bin_ypix - 1][bin_xpix - 1] + l1cstr->Lt_blue[bb][pix];
                                              bincount[view][sb][bin_ypix - 1][bin_xpix - 1] += 1;
                                           }
                                    bb++;
                                    sb++;
                                }//end while blue band band

        //RED---
                                while (sb <num_blue_bands+num_red_bands) {
                                    if (l1cstr->Lt_red[rb][pix] > 0.) {
                                        binmean_Lt[view][sb][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][sb][bin_ypix - 1][bin_xpix - 1] + l1cstr->Lt_red[rb][pix];
                                        bincount[view][sb][bin_ypix - 1][bin_xpix - 1] += 1;
                                    }
                                    rb++;
                                    sb++;
                                }//end red bands

       //SWIR----
                                while (sb <num_blue_bands + num_red_bands + num_SWIR_bands) {
                                    if (l1cstr->Lt_SWIR[swb][pix] > 0.) {
                                        binmean_Lt[view][sb][bin_ypix - 1][bin_xpix - 1] = binmean_Lt[view][sb][bin_ypix - 1][bin_xpix - 1] + l1cstr->Lt_SWIR[swb][pix];
                                        bincount[view][sb][bin_ypix - 1][bin_xpix - 1] += 1;
                                    }
                                    swb++;
                                    sb++;
                                }//end swir bands
                             }

                            }//end if INPIX==1



    //********** pixel AREA WEIGHTING ************************
                            if (flag_inpix == 0) OUTPIX++;

                            TOTPIX++;

                         }//end pixels loop

                  
                        //END OF FILE--****** reset starting sline index if last line processed and within k box
                        //go for another granule if it is not the last one---------------
              
         //***** skip lines and files if sline is not within k group ********
              //increase k and screen files and lines again starting from the last binning


                //**********************************************************************************************************************************
                //writing each granule ---------------------------------
                //*****************************************************************************

                if (sline==num_scans-1) {
                    cout << "writing file #....................................................................................................................."  <<l1cfile->l1b_name<< endl;

                    //********* time metadata **************************************************************
      //              if(i==0) start_time=stime[first_line];
      //              if(i==n_files-1) end_time=stime[last_line];

//ini time swath
               //     unix2ymds(start_time, &syear, &smon, &sday,&secs);
  //                  double tfile_ini_sec=ymds2unix(syear,smon,sday,secs);  
//                    unix2ymdhms(tfile_ini_sec,&syear,&smon,&sday, &shour, &smin, &seconds); 
//                    
           /*         unix2ymdhms(start_time, &syear2,&smon2,&sday2,&shour,&smin,&seconds);

                    y_swt=std::to_string(syear);
                    mo_swt=std::to_string(smon);
                    d_swt=std::to_string(sday);
                    h_swt=std::to_string(shour);
                    mi_swt=std::to_string(smin);
                    s_swt2=std::to_string(round(seconds));

                    length = (int) floor( log10 (smon) ) + 1;
                    if(length==1) mo_swt="0"+mo_swt;
                    length = (int) floor( log10 (sday) ) + 1;
                    if(length==1) d_swt="0"+d_swt;

                    if(shour==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (shour+logoff)) + 1;
                    if(length==1) h_swt="0"+h_swt;

                    if(smin==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (smin+logoff)) + 1;
                    if(length==1) mi_swt="0"+mi_swt;

                    if(seconds==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (round(seconds)+logoff)) + 1;
                    if(length==1) s_swt2="0"+s_swt2;

                    tswt_ini=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                    tswt_ini_file=y_swt+mo_swt+d_swt+"T"+h_swt+mi_swt+s_swt2.substr(0,2);

                    l1cfile->tswt_ini=tswt_ini;
                    l1cfile->tswt_ini_file=tswt_ini_file;


//end time swath
              
                    unix2ymdhms(end_time, &syear2,&smon2,&sday2,&shour,&smin,&seconds);

                    y_swt=std::to_string(syear);
                    mo_swt=std::to_string(smon);
                    d_swt=std::to_string(sday);
                    h_swt=std::to_string(shour);
                    mi_swt=std::to_string(smin);
                    s_swt2=std::to_string(round(seconds));

                    length = (int) floor( log10 (smon) ) + 1;
                    if(length==1) mo_swt="0"+mo_swt;
                    length = (int) floor( log10 (sday) ) + 1;
                    if(length==1) d_swt="0"+d_swt;

                    if(shour==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (shour+logoff)) + 1;
                    if(length==1) h_swt="0"+h_swt;

                    if(smin==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (smin+logoff)) + 1;
                    if(length==1) mi_swt="0"+mi_swt;

                    if(seconds==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (round(seconds)+logoff)) + 1;
                    if(length==1) s_swt2="0"+s_swt2;

                    tswt_end=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);

                    l1cfile->tswt_end=tswt_end;
                 
*/

//*********************************************************************************************************
                    //write mean Lt as nc file---
                   //**********************************************************************
                   //**********************************************************************8
                   //create Lt file

                    string GATT_NAME1,GATT_VAL1,GATT_NAME2,GATT_VAL2;

                     if(format.type==FT_SPEXONE){
                         senstr = "SPEXONE";
                         GATT_NAME1="title",GATT_VAL1="PACE SPEXone Level-1C Data",GATT_NAME2="instrument",GATT_VAL2=senstr2;
                     }
                     else if(format.type==FT_HARP2){
                         senstr = "HARP2";
                         GATT_NAME1="title",GATT_VAL1="PACE HARP2 Level-1C Data",GATT_NAME2="instrument",GATT_VAL2=senstr2;
                     }
                     else//OCIS
                     {
                         senstr = "OCI";
                         GATT_NAME1="title",GATT_VAL1="PACE OCI Level-1C Data",GATT_NAME2="instrument",GATT_VAL2=senstr2;
                     }

//                    ofilestr=std::string(l1cinput->ofile);

                    pathstr="out/";
                    fname_out =pathstr+"PACE_"+senstr+"."+tswt_ini_file2+"ZZ"+".L1C.5.2km.nc";
                    produstr="PACE_"+senstr+"."+tswt_ini_file2+".L1C.5.2km.nc";

                    string  ATT_NAME="Units", ATT_VAL="degrees",GATT_NAME3="processing_version",GATT_VAL3="V1.0",GATT_NAME4="Conventions",GATT_VAL4="CF-1.6";
                    string GATT_NAME5="institution",GATT_VAL5="NASA Goddard Space Flight Center, Ocean Biology Processing Group",GATT_NAME6="license",GATT_VAL6="http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
                    string GATT_NAME7="naming_authority",GATT_VAL7="gov.nasa.gsfc.sci.oceancolor",GATT_NAME8="keywords_vocabulary",GATT_VAL8="NASA Global Change Master Directory (GCMD) Science Keywords";
                    string GATT_NAME9="stdname_vocabulary",GATT_VAL9="NetCDF Climate and Forecast (CF) Metadata Convention",GATT_NAME10="creator_name",GATT_VAL10="NASA/GSFC",GATT_NAME11="creator_email",GATT_VAL11="data@oceancolor.gsfc.nasa.gov";
                    string GATT_NAME12="creator_url",GATT_VAL12="http://oceancolor.gsfc.nasa.gov",GATT_NAME13="project",GATT_VAL13="PACE Project",GATT_NAME14="publisher_name",GATT_VAL14="NASA/GSFC";
                    string GATT_NAME15="publisher_email",GATT_VAL15="data@oceancolor.gsfc.nasa.gov",GATT_NAME16="publisher_url",GATT_VAL16="http://oceancolor.gsfc.nasa.gov",GATT_NAME17="processing_level",GATT_VAL17="L1C";
                    string GATT_NAME18="cdm_data_type",GATT_VAL18="swath",GATT_NAME19="orbit_number",GATT_VAL19="12345",GATT_NAME20="history",GATT_VAL20="",GATT_NAME21="CDL_version_date",GATT_VAL21="2021-09-10",GATT_NAME22="product_name",GATT_VAL22=produstr;
                    string GATT_NAME23="startDirection",GATT_VAL23=start_dir,GATT_NAME24="endDirection",GATT_VAL24=end_dir,GATT_NAME25="time_coverage_start",GATT_VAL25=tswt_ini,GATT_NAME26="time_coverage_end",GATT_VAL26=tswt_end,GATT_NAME27="date{_created",GATT_VAL27="2021-09-10T15:12:41Z",GATT_NAME28="sun_earth_distance",GATT_VAL28="0.990849042172323",GATT_NAME29="terrain_data_source",GATT_VAL29="",GATT_NAME30="spectral_response_function",GATT_VAL30="",GATT_NAME31="systematic_uncertainty_model",GATT_VAL31="",GATT_NAME32="nadir_bin",GATT_VAL32=nadpix_str,GATT_NAME33="bin_size_at_nadir",GATT_VAL33=binstr;

                    l1cfile->gridname = fname_out.c_str();
                    filename_lt = fname_out.c_str();
                    cout<<"filename_lt.."<<filename_lt<<endl;               

                    if ((status = nc_create(filename_lt, NC_CLOBBER | NC_NETCDF4, &ncid_out)))
                        check_err(status, __LINE__, __FILE__);

                     //define dims
                     // Define the dimensions,vars and attributes at the root level
                    NDIMS=2;
                    NDIMS2=3;
                    NDIMS3=4;
                    NDIMS4=1;
                    NY = num_gridlines;
                    NX = nbinx;
                    //NVIEWS= 2;//for OCI
                   // NBANDS=249;//for OCI
                    l1cfile->NY1=0;//these are indexes not line #!!!!!!!!!!!!!!!!!
                    l1cfile->NY2=num_gridlines-1; 

                    NY1 = l1cfile->NY1;//these are indexes not line #!!!!!!!!!!!!!!!!!
                    NY2 = l1cfile->NY2;
                    cout << "rowgrid ini.." << NY1 << "rowgrid end.." << NY2 << endl;
                    NY = NY2 - NY1 + 1;
                    cout << "NY.." << NY << "NX.." << NX << "NY1.."<<NY1<<"NY2.."<<NY2<<endl;

  //DEF DIMENSIONS
                    if ((status = nc_def_dim(ncid_out, "bins_across_track", NX, &x_dimid)))
                        check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_dim(ncid_out, "bins_along_track", NY, &y_dimid)))
                           check_err(status, __LINE__, __FILE__);
                    //dims for output var
                    dimids[0] = y_dimid;
                    dimids[1] = x_dimid;

                    if ((status = nc_def_dim(ncid_out, "number_of_views", NVIEWS, &v_dimid)))
                       check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_dim(ncid_out, "intensity_bands_per_view", NBANDS, &b_dimid)))
                      check_err(status, __LINE__, __FILE__);

                    if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                       if ((status = nc_def_dim(ncid_out, "polarization_bands_per_view", NBANDS_POL, &b_dimid_pol)))
                         check_err(status, __LINE__, __FILE__);
                      }

                    dimids2[0] = y_dimid;
                    dimids2[1] = x_dimid;
                    dimids2[2] = v_dimid;//#o

                    dimids3[0] = y_dimid;
                    dimids3[1] = x_dimid;
                    dimids3[2] = v_dimid;//
                    dimids3[3] = b_dimid;

                    dimids4=&v_dimid;

                    dimids5[0]=v_dimid;
                    dimids5[1]=b_dimid;


                    if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                       dimids6[0] = y_dimid;
                       dimids6[1] = x_dimid;
                       dimids6[2] = v_dimid;//
                       dimids6[3] = b_dimid_pol;

                       dimids7[0]=v_dimid;
                       dimids7[1]=b_dimid_pol;
                    }

                                                     
                    //define attributes
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME1.c_str(), strlen(GATT_VAL1.c_str()),
                        GATT_VAL1.c_str()))
                       check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME2.c_str(), strlen(GATT_VAL2.c_str()),
                         GATT_VAL2.c_str()))
                       check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME3.c_str(), strlen(GATT_VAL3.c_str()),
                         GATT_VAL3.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME4.c_str(), strlen(GATT_VAL4.c_str()),
                            GATT_VAL4.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME5.c_str(), strlen(GATT_VAL5.c_str()),
                         GATT_VAL5.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME6.c_str(), strlen(GATT_VAL6.c_str()),
                          GATT_VAL6.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME7.c_str(), strlen(GATT_VAL7.c_str()),
                          GATT_VAL7.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME8.c_str(), strlen(GATT_VAL8.c_str()),
                          GATT_VAL8.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME9.c_str(), strlen(GATT_VAL9.c_str()),
                          GATT_VAL9.c_str()))
                         check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME10.c_str(), strlen(GATT_VAL10.c_str()),
                          GATT_VAL10.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME11.c_str(), strlen(GATT_VAL11.c_str()),
                           GATT_VAL11.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME12.c_str(), strlen(GATT_VAL12.c_str()),
                            GATT_VAL12.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME13.c_str(), strlen(GATT_VAL13.c_str()),
                             GATT_VAL13.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME14.c_str(), strlen(GATT_VAL14.c_str()),
                             GATT_VAL14.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME15.c_str(), strlen(GATT_VAL15.c_str()),
                                 GATT_VAL15.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME16.c_str(), strlen(GATT_VAL16.c_str()),
                             GATT_VAL16.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME17.c_str(), strlen(GATT_VAL17.c_str()),
                             GATT_VAL17.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME18.c_str(), strlen(GATT_VAL18.c_str()),
                             GATT_VAL18.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME19.c_str(), strlen(GATT_VAL19.c_str()),
                             GATT_VAL19.c_str()))
                          check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME20.c_str(), strlen(GATT_VAL20.c_str()),
                             GATT_VAL20.c_str()))
                              check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME21.c_str(), strlen(GATT_VAL21.c_str()),
                             GATT_VAL21.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME22.c_str(), strlen(GATT_VAL22.c_str()),
                             GATT_VAL22.c_str()))
                           check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME23.c_str(), strlen(GATT_VAL23.c_str()),
                              GATT_VAL23.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME24.c_str(), strlen(GATT_VAL24.c_str()),
                              GATT_VAL24.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME25.c_str(), strlen(GATT_VAL25.c_str()),
                           GATT_VAL25.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME26.c_str(), strlen(GATT_VAL26.c_str()),
                           GATT_VAL26.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME27.c_str(), strlen(GATT_VAL27.c_str()),
                             GATT_VAL27.c_str()))
                              check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME28.c_str(), strlen(GATT_VAL28.c_str()),
                             GATT_VAL28.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME29.c_str(), strlen(GATT_VAL29.c_str()),
                             GATT_VAL29.c_str()))
                             check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME30.c_str(), strlen(GATT_VAL30.c_str()),
                             GATT_VAL30.c_str()))
                              check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME31.c_str(), strlen(GATT_VAL31.c_str()),
                             GATT_VAL31.c_str()))
                            check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME32.c_str(), strlen(GATT_VAL32.c_str()),
                              GATT_VAL32.c_str()))
                                check_err(status, __LINE__, __FILE__);
                    if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME33.c_str(), strlen(GATT_VAL33.c_str()),
                                GATT_VAL33.c_str()))
                                check_err(status, __LINE__, __FILE__);

                    //define groups                
                    if ((status = nc_def_grp(ncid_out, "geolocation_data", &grp_coor)))  //netcdf-4
                         check_err(status, __LINE__, __FILE__);
                     //def var grp1
                    if ((status = nc_def_var(grp_coor, "latitude", NC_FLOAT, NDIMS,
                         dimids, &varid1)))
                          check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_var(grp_coor, "longitude", NC_FLOAT, NDIMS,
                            dimids, &varid2)))
                          check_err(status, __LINE__, __FILE__);

                    if ((status = nc_def_var_deflate(grp_coor, varid1, shuffle, deflate,deflate_level)))
                      check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_var_deflate(grp_coor, varid2, shuffle, deflate,deflate_level)))
                      check_err(status, __LINE__, __FILE__);

                    if ((status=nc_def_var_fill(grp_coor, varid1, 1, &fill_value)))
                           check_err(status, __LINE__, __FILE__);
                    if ((status=nc_def_var_fill(grp_coor, varid2, 1, &fill_value)))
                           check_err(status, __LINE__, __FILE__);


                    //leave define mode-----------------------
                    if ((status = nc_enddef(grp_coor))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);


                    if ((status = nc_def_grp(ncid_out, "observation_data", &grp_obs)))  //netcdf-4
                        check_err(status, __LINE__, __FILE__);
                    //def var grp2
                        //counts and Lt
                    if ((status = nc_def_var(grp_obs, "obs_per_view", NC_INT, NDIMS2,
                        dimids2, &varid_count)))
                        check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_var(grp_obs, "I", NC_FLOAT, NDIMS3,
                        dimids3, &varid_lt)))
                        check_err(status, __LINE__, __FILE__);

                    if ((status = nc_def_var_deflate(grp_obs, varid_count, shuffle, deflate,deflate_level)))
                      check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_var_deflate(grp_obs, varid_lt, shuffle, deflate,deflate_level)))
                      check_err(status, __LINE__, __FILE__);

                    if ((status=nc_def_var_fill(grp_obs, varid_count, 1, &fill_value)))
                           check_err(status, __LINE__, __FILE__);
                    if ((status=nc_def_var_fill(grp_obs, varid_lt, 1, &fill_value)))
                           check_err(status, __LINE__, __FILE__);


                    if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                      if ((status = nc_def_var(grp_obs, "I_polsample", NC_FLOAT, NDIMS3,
                        dimids6, &varid_lt_pol)))
                        check_err(status, __LINE__, __FILE__);
                     
                      if ((status = nc_def_var_deflate(grp_obs, varid_lt_pol, shuffle, deflate,deflate_level)))
                         check_err(status, __LINE__, __FILE__);
                      if ((status=nc_def_var_fill(grp_obs, varid_lt_pol, 1, &fill_value)))
                           check_err(status, __LINE__, __FILE__);
                    }


                    //leave define mode-----------------------
                    if ((status = nc_enddef(grp_obs))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);

                     //define groups
                    if ((status = nc_def_grp(ncid_out, "sensor_views_bands", &grp_views)))  //netcdf-4
                         check_err(status, __LINE__, __FILE__);
                     //def var grp1
                    if ((status = nc_def_var(grp_views, "view_angles", NC_FLOAT, NDIMS4,
                         dimids4, &varid_v1)))
                          check_err(status, __LINE__, __FILE__);
                    if ((status = nc_def_var(grp_views, "intensity_wavelengths", NC_FLOAT, NDIMS,
                            dimids5, &varid_v2)))
                          check_err(status, __LINE__, __FILE__);
                    if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                      if ((status = nc_def_var(grp_views, "polarization_wavelengths", NC_FLOAT, NDIMS,
                            dimids7, &varid_v3)))
                          check_err(status, __LINE__, __FILE__);
                     }

                    //leave define mode-----------------------
                    if ((status = nc_enddef(grp_views))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);


      //sensor view group
                    view_out1 = (float*)calloc(NVIEWS,sizeof(float));
                    view_out2 = allocate2d_float(NVIEWS,NBANDS);

                    if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                      view_out2_pol = allocate2d_float(NVIEWS,NBANDS_POL);
                         }

           //assign views
                   if(format.type==FT_SPEXONE){  
                        view_out1[0]=-58;
                        view_out1[1]=-20;
                        view_out1[2]=0;
                        view_out1[3]=20;
                        view_out1[4]=58;

                        for (size_t i = 0; i < NVIEWS; i++) {
                            for (size_t j = 0; j < NBANDS; j++) {
                                view_out2[i][j]=l1cstr->I_lambdas[i][j];                                         
                                }}
                        for (size_t i = 0; i < NVIEWS; i++) {
                            for (size_t j = 0; j < NBANDS_POL; j++) {
                                view_out2_pol[i][j]=l1cstr->pol_lambdas[i][j];
                                }}
                            
                     }
                   else{
                       view_out1[0]=-20;//OCI
                       view_out1[1]=20;

                        //assign 249 wavelengths coming from L1B file (actually they are 239 as may 17, 2022 as overlapping are removed

                       size_t sb=0;

                         for (size_t i = 0; i < NVIEWS; i++) {
                               for (size_t j = 0; j < num_blue_bands; j++) {
                                view_out2[i][sb]=l1cstr->blue_lambdas[j];
                  //              cout<<"sb..."<<sb+1<<"bluelam.."<<l1cstr->blue_lambdas[j]<<endl;
                                sb++;
                                }
                         for (size_t j = 0; j < num_red_bands; j++) {
                                view_out2[i][sb]=l1cstr->red_lambdas[j];
                    //            cout<<"sb..."<<sb+1<<"redlam.."<<l1cstr->red_lambdas[j]<<endl;
                                sb++;
                                }
                         for (size_t j = 0; j < num_SWIR_bands; j++) {
                                 view_out2[i][sb]=l1cstr->SWIR_lambdas[j];
                      //           cout<<"sb..."<<sb+1<<"swirlam.."<<l1cstr->SWIR_lambdas[j]<<endl;
                                 sb++;
                               }
                            sb=0;
                            }
                     }

   //PLACE VARS                                 
                    if ((status = nc_put_var_float(grp_views, varid_v1, &view_out1[0])))//view angles
                        check_err(status, __LINE__, __FILE__);
                    if ((status = nc_put_var_float(grp_views, varid_v2, &view_out2[0][0])))//wavelengths per view
                        check_err(status, __LINE__, __FILE__);

                    if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                             if ((status = nc_put_var_float(grp_views, varid_v3, &view_out2_pol[0][0])))//wavelengths per view
                                   check_err(status, __LINE__, __FILE__);
                            }

         //create the arrays for geoloc
                    lat_out = allocate2d_float(NY, NX);
                    lon_out = allocate2d_float(NY, NX);

                    int c = 0;

                    for (int i = NY1; i < NY2 + 1; i++) {
                        for (int j = 0; j < NX; j++) {
                            lat_out[i][j] = l1cfile->lat_gd[i][j];
                            lon_out[i][j] = l1cfile->lon_gd[i][j];
                        }
                        c++;
                    }

                    if ((status = nc_put_var_float(grp_coor, varid1, &lat_out[0][0])))
                        check_err(status, __LINE__, __FILE__);

                    if ((status = nc_put_var_float(grp_coor, varid2, &lon_out[0][0])))
                        check_err(status, __LINE__, __FILE__);


           //alloc mem for output variables obs_per_view and I as part of the observation_data group
                    data_out = allocate3d_int(NY, NX,NVIEWS);
                    data_out2 = allocate4d_float(NY, NX,NVIEWS,NBANDS);

                    if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                          data_out_pol = allocate3d_int(NY, NX,NVIEWS);
                          data_out2_pol = allocate4d_float(NY, NX,NVIEWS,NBANDS_POL);
                         }

            // ARRAY INIT
               for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                            for (int j = 0; j < NX; j++) {
                                 data_out[i][j][v]=0;
                            }}}

               for (size_t b = 0;b <NBANDS;b++) {
                    for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                            for (int j = 0; j < NX; j++) {
                                 data_out2[i][j][v][b]=0;
                            }}}}


                // ARRAY INIT POL
               if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                  for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                            for (int j = 0; j < NX; j++) {
                                 data_out_pol[i][j][v]=0;
                            }}}

                  for (size_t b = 0;b <NBANDS_POL;b++) {
                    for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                            for (int j = 0; j < NX; j++) {
                                 data_out2_pol[i][j][v][b]=0;
                            }}}}
                    }


//--------- intensity arrays ----------------------------------------------------------------
          //count  arrays---
                   
                   for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                            for (int j = 0; j < NX; j++) {
                                for (sb = 0;sb < NBANDS;sb++) {
                                    if(bincount[v][sb][i][j]>0)
                                         data_out[i][j][v]+=bincount[v][sb][i][j];
                                    else
                                         data_out[i][j][v]+=0;
                                }
                            }
                           
                        }
                        
                      }

                        cout << "writing countbin for all band ."<< endl;
                       if ((status = nc_put_var_int(grp_obs, varid_count, &data_out[0][0][0])))
                            check_err(status, __LINE__, __FILE__);

            // LT ARRAYS--
                    
                    for (size_t b = 0;b <NBANDS;b++) {
                      for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                             for (int j = 0; j < NX; j++) {
                                if (bincount[v][b][i][j] > 0)
                                    data_out2[i][j][v][b] = binmean_Lt[v][b][i][j] / bincount[v][b][i][j];
                                else         data_out2[i][j][v][b] = NAN;
                            }
                            
                        }
                        
                       }//end views
                     }//end for bands

                         cout << "writing Lt for all bands.."<< endl;
                        if ((status = nc_put_var_float(grp_obs, varid_lt, &data_out2[0][0][0][0])))
                            check_err(status, __LINE__, __FILE__);

//polarization arrays---count and Lts-----------------------------------------------------
          //count arrays---
                 if(format.type==FT_SPEXONE || format.type==FT_HARP2){
                   for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                            for (int j = 0; j < NX; j++) {
                                for (sb = 0;sb < NBANDS_POL;sb++) {
                                    if(bincount_pol[v][sb][i][j]>0)
                                         data_out_pol[i][j][v]+=bincount_pol[v][sb][i][j];
                                    else
                                         data_out_pol[i][j][v]+=0;
                                }
                            }
                           
                        }
                        
                      }


                        cout << "writing countbin for all band ."<< endl;
                        if ((status = nc_put_var_int(grp_obs, varid_count, &data_out_pol[0][0][0])))
                            check_err(status, __LINE__, __FILE__);

            // LT ARRAYS--
                    
                    for (size_t b = 0;b <NBANDS_POL;b++) {
                      for (size_t v = 0; v < NVIEWS; v++) {
                        for (int i = NY1; i < NY2 + 1; i++) {
                             for (int j = 0; j < NX; j++) {
                                if (bincount_pol[v][b][i][j] > 0)
                                    data_out2_pol[i][j][v][b] = binmean_Lt_pol[v][b][i][j] / bincount_pol[v][b][i][j];
                                else         data_out2_pol[i][j][v][b] = NAN;
                            }
                           
                        }
                        
                       }//end views
                     }//end for bands


                        cout << "writing Lt for all bands.."<< endl;
                        if ((status = nc_put_var_float(grp_obs, varid_lt_pol, &data_out2_pol[0][0][0][0])))
                            check_err(status, __LINE__, __FILE__);
                   }    

//polarization arrays---count and Lts-----------------------------------------------------
                  if (data_out != nullptr){
                    delete [] (data_out);}
                  if (data_out2 != nullptr){
                    delete [] (data_out2);}
                  if (data_out_pol != nullptr){
                    delete [] (data_out_pol);}
                  if (data_out2_pol != nullptr){
                    delete [] (data_out2_pol);}
                  if (lat_out != nullptr){
                    delete[](lat_out);}
                  if (lon_out != nullptr){
                    delete[](lon_out);}
                  if (view_out1 != nullptr){
                    delete[](view_out1);}
                  if (view_out2 != nullptr){
                    delete[](view_out2);}
                  if (view_out2_pol != nullptr){
                    delete[](view_out2_pol);}

                    if ((status = nc_close(ncid_out)))
                        check_err(status, __LINE__, __FILE__);

                    lat_out = nullptr;;
                    lon_out = nullptr;
                    data_out = nullptr;
                    data_out2 = nullptr;
                    view_out1 = nullptr;
                    view_out2 = nullptr;                 
                    view_out2_pol = nullptr;

                    cout << "number of files per swath.." << nfiles_swt[swtd - 1] << "for swath.." << swtd << endl;
                    cout << "total inpix.." << INPIX << "total pix.." << TOTPIX << endl;
                    cout<<"*** SUCCESS writing bincount and binLt rasters as nc..!.."<<endl;

                }//end writing file


    if(gdindex!=nullptr)
      delete [] (gdindex);

    cout<<"number of binned pixels..."<<INPIX<<endl;            
    return 0;
      }



int32_t L1C::pix_corners4_l1c(l1c_filehandle* l1cfile, L1C_input* l1cinput, float dist_u, float dist_v, float azpix, int32_t scanline, int32_t pix, float pixlat, float pixlon, float pixLt, float** lat_asort, short** index_xy, float** lat_gd, float** lon_gd, double areaFracBox[3][3], float** Ltfracsum, float** areafracsum, float** nobs_perbin) { //compute pixel corners for a specific line,pixel and file

        double lat_cnw, lat_cne, lat_csw, lat_cse, lon_cnw, lon_cne, lon_csw, lon_cse;
        float  azpixc;
        short gd_row, gd_col;
        int32_t num_gridlines, nbinx;
        double areabinBox[3][3];
        float theta, thetares, dist_corn, thetagrad;
        double res1, res2, res3, res4;//pixel resolution in meters

        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
        Polygon_t pixelPoly;

        num_gridlines = l1cfile->num_gridlines;
        nbinx = l1cfile->nbinx;


        //******* M_PIXEL CORNERS *****************************
        //***************************************************
        dist_corn = 0.5 * (sqrt(dist_u * dist_u + dist_v * dist_v)) * 1000;//distane center to corner in meters
\
        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = (azpix - thetares) * degrad;//pixel bearing
        if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in nw corner..azpixc" << azpixc << endl;exit(1); }
        azpixc = azpixc * 180 / M_PI;

      //  cout << "azpixc nw...." << azpixc << endl;
        geod.Direct(pixlat, pixlon, azpixc, dist_corn, lat_cnw, lon_cnw);
        geod.Inverse(pixlat, pixlon, lat_cnw, lon_cnw, res1);
        // lat_cne=asin(sin(pixlat*degrad)*cos(ad)+cos(pixlat*degrad)*sin(ad)*cos(azpixc*degrad))*180./M_PI;
        // lon_cne=pixlon+atan2(sin(azpixc*degrad)*sin(ad)*cos(pixlat*degrad),cos(ad)-sin(pixlat*degrad)*sin(lat_cne*degrad))*180./M_PI;


        //ne corner (sw-180)
        azpixc = (azpix + thetares) * degrad;//pixel bearing
        if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in ne corner..azpixc" << azpixc << endl;exit(1); }
        // azpixc=(azpix+45)*degrad;//pixel bearing
        /* if(azpixc>M_PI | azpixc<-M_PI){
              }
         if(azpixc<0.){
                  azpixc=360+azpixc*180/M_PI;
                     }
         else azpixc=azpixc*180/M_PI;
        */

        azpixc = azpixc * 180 / M_PI;
      //  cout << "azpixc ne...." << azpixc << endl;

        geod.Direct(pixlat, pixlon, azpixc, dist_corn, lat_cne, lon_cne);
        geod.Inverse(pixlat, pixlon, lat_cne, lon_cne, res2);

        //sw corner
        // azpixc=(azpix-45-90)*degrad;//pixel bearing
        azpixc = (azpix - thetares - thetagrad) * degrad;//pixel bearing
        if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in sw corner..azpixc" << azpixc << endl;exit(1); }
        /* if(azpixc>M_PI | azpixc<-M_PI){
              //      cout<<"problem with BEARING in across-gridline method...az<-180 or >180...."<<"SW azpixc in degrees.."<<azpixc*180/M_PI<<endl;
                     }
         if(azpixc<0.){
                  azpixc=360+azpixc*180/M_PI;
                     }
         else azpixc=azpixc*180/M_PI;
         */

        azpixc = azpixc * 180 / M_PI;
      //  cout << "azpixc sw....." << azpixc << endl;
      geod.Direct(pixlat, pixlon, azpixc, dist_corn, lat_csw, lon_csw);
        geod.Inverse(pixlat, pixlon, lat_csw, lon_csw, res3);
        //lat_csw=asin(sin(pixlat*degrad)*cos(ad)+cos(pixlat*degrad)*sin(ad)*cos(azpixc*degrad))*180./M_PI;
        // lon_csw=pixlon+atan2(sin(azpixc*degrad)*sin(ad)*cos(pixlat*degrad),cos(ad)-sin(pixlat*degrad)*sin(lat_csw*degrad))*180./M_PI;

        //se corner (nw-180)
        azpixc = (azpix + thetares + thetagrad) * degrad;//pixel bearingi
        if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in se corner..azpixc" << azpixc << endl;exit(1); }
        // azpixc=(azpix+45+90)*degrad;//pixel bearing
        /* if(azpixc>M_PI | azpixc<-M_PI){
              //      cout<<"problem with BEARING in across-gridline method...az<-180 or >180...."<<"SE azpixc in degrees.."<<azpixc*180/M_PI<<endl;
                     }
         if(azpixc<0.){
                  azpixc=360+azpixc*180/M_PI;
                     }
         else azpixc=azpixc*180/M_PI;
        */
        azpixc = azpixc * 180 / M_PI;

      //  cout << "azpix se...." << azpixc << endl;
        geod.Direct(pixlat, pixlon, azpixc, dist_corn, lat_cse, lon_cse);
        geod.Inverse(pixlat, pixlon, lat_cse, lon_cse, res4);
        // lat_cse=asin(sin(pixlat*degrad)*cos(ad)+cos(pixlat*degrad)*sin(ad)*cos(azpixc*degrad))*180./M_PI;
        // lon_cse=pixlon+atan2(sin(azpixc*degrad)*sin(ad)*cos(pixlat*degrad),cos(ad)-sin(pixlat*degrad)*sin(lat_cse*degrad))*180./M_PI;


       if(res1>4000 || res2>4000 || res3>4000 || res4>4000){
           cout<<"ERROR...half diagonal of pixel cant be larger than 3675 m for a 7.8 x 1.6 km  pixel..."<<endl;

           cout<<"pixlat.."<<pixlat<<"pixlon.."<<pixlon<<endl;
           cout << "scaline.."<<scanline+1<<"pix.." << pix + 1 << "azpix.." << azpix << "theta.." << thetagrad << "thetares.." << thetares << endl;
           cout << "dist_u.." << dist_u << "dist_v.." << dist_v << endl;
           cout << "res1." << res1 << "res2." << res2 << "res3.." << res3 << "res4.." << res4 << endl;
           exit(1);}
       else{
       /*    
        cout<<"pixlat.."<<pixlat<<"pixlon.."<<pixlon<<endl;
        cout << "scaline.."<<scanline+1<<"pix.." << pix + 1 << "azpix.." << azpix << "theta.." << thetagrad << "thetares.." << thetares << endl;
        cout << "dist_u.." << dist_u << "dist_v.." << dist_v << endl;
        cout << "res1." << res1 << "res2." << res2 << "res3.." << res3 << "res4.." << res4 << endl;
        cout << "pix sw corner..lat." << lat_csw << "sw corner lon.." << lon_csw << endl;
        cout << "pix nw corner..lat." << lat_cnw << "nw corner lon.." << lon_cnw << endl;
        cout << "pix ne corner..lat." << lat_cne << "ne corner lon.." << lon_cne << endl;
        cout << "pix se corner..lat." << lat_cse << "se corner lon.." << lon_cse << endl;
        */
        }



        std::string ws_lon_str, ws_lat_str, wn_lon_str, wn_lat_str, en_lon_str, en_lat_str, es_lon_str, es_lat_str, pixstr;
        ws_lon_str = std::to_string(lon_csw);
        ws_lat_str = std::to_string(lat_csw);
        wn_lon_str = std::to_string(lon_cnw);
        wn_lat_str = std::to_string(lat_cnw);
        en_lon_str = std::to_string(lon_cne);
        en_lat_str = std::to_string(lat_cne);
        es_lon_str = std::to_string(lon_cse);
        es_lat_str = std::to_string(lat_cse);


        pixstr = "POLYGON((" + ws_lon_str + " " + ws_lat_str + "," + wn_lon_str + " " + wn_lat_str + "," + en_lon_str + " " + en_lat_str + "," + es_lon_str + " " + es_lat_str + "," + ws_lon_str + " " + ws_lat_str + "))";
        bg::read_wkt(pixstr, pixelPoly);

        if (lon_csw > lon_cse || lon_cnw > lon_cne || lat_csw > lat_cnw || lat_cse > lat_cne) {
    //        cout << "lon_csw>lon_cse | lon_cnw>lon_cne | lat_csw>lat_cnw | lat_cse>lat_cne//" << lon_csw << lon_cse << lon_cnw << lon_cne << lat_csw << lat_cnw << lat_cse << lat_cne << endl;
    //        cout << "corners are switched for pixel due to the change of longitude from +180 to -180...." << pix + 1 << endl;
        }

        // exit(1);
        //*** L1C grid cell *******************************************************************************

        bool boolbin1 = 0;

        boolbin1 = sbs2_l1c(l1cinput, num_gridlines, nbinx, lat_asort, index_xy, pixlat, pixlon, lon_gd, &gd_row, &gd_col);
        //Ini areaBinbox
        for (int ar = 0;ar < 3;ar++) {
            for (int ac = 0;ac < 3;ac++) {
                areabinBox[ar][ac] = 0.;
                areaFracBox[ar][ac] = 0.;
                ac++;
            }
            ar++;
        }

        if (boolbin1 == 1 && gd_row >= 1 && gd_col >= 1 && gd_row < num_gridlines - 1 && gd_col < nbinx - 1) {

     //       cout<<"scanline.."<<scanline+1<<"pix..."<<pix+1<<"gd row.."<<gd_row<<"gd_col.."<<gd_col<<endl;

            binIntersectsPix4corn4_l1c2(l1cfile, l1cinput, gd_row, gd_col, lat_gd, lon_gd, pixelPoly, areaFracBox, areabinBox);

            int ar = 0, ac = 0;
            for (int gr = gd_row - 1;gr < gd_row + 2;gr++) {
                for (int gc = gd_col - 1;gc < gd_col + 2;gc++) {
                    if (pixLt > 0.0) {
                        Ltfracsum[gr][gc] += pixLt * areaFracBox[ar][ac];
                        areafracsum[gr][gc] += areaFracBox[ar][ac];
                        nobs_perbin[gr][gc] += 1;

         //           cout<<"Ltfracsum[gr][gc].."<<Ltfracsum[gr][gc]<<"areafracsum[gr][gc].."<<areafracsum[gr][gc]<<"nobs_perbin[gr][gc].."<<nobs_perbin[gr][gc]<<endl;
                    }
                    else {
                        Ltfracsum[gr][gc] += 0.0;
                        areafracsum[gr][gc] += 0.0;
                        nobs_perbin[gr][gc] += 0;
                    }

                    ac++;
                }
                ac = 0;
                ar++;
            }

        }
        else {
            //  cout<<"pix out of the L1Cgrid.........................."<<"scanline.."<<scanline+1<<"pix#..."<<pix+1<<"pixlat.."<<pixlat<<"pixlon..."<<pixlon<<endl;
        }


        return 0;
    }



    int32_t L1C::binL1C_wgranule_aw3(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput, l1c_str* l1cstr,float** Ltfracsum, float** areabinsum, float** nobs_perbin, size_t sline) {     
       int32_t num_gridlines, nbinx, n_files, inpix = 0, outpix = 0, totpix = 0;
       float minv = 0., maxv = 0.; 
       int16_t fi = 0;
       string str, pathstr, senstr, monstr, daystr, yearstr, prodstr, swtstr, granstr, swtnum, extstr, fname_out;
       float Re = 6378;
       float latmin_swt = 100, latmax_swt = 0, lonmin_swt = 0, lonmax_swt = 0;
       int32_t NY = -1,NY1 = -1,NY2 = -1;
       float aterm = 0.;
       bool boolbin1=0;
       int num_pixels;
       string gridname, azeast_name;
       float azpix, azpixc;
       double** latcornBox = nullptr, ** loncornBox = nullptr;
       float deltaphi, deltalam, bterm, dist_u, dist_v;
       double lat_cnw, lat_cne, lat_csw, lat_cse, lon_cnw, lon_cne, lon_csw, lon_cse;
       float theta, thetares, dist_corn, thetagrad;
       double intersectArea = 0, binAreapix = 0, binAreagrid = 0.;
       bool ingeom = false;
       int pc;
       short  gd_row = 0, gd_col = 0;
       int bb = 0;
       int pix;
       int flag_inpix;
       string ws_lon_str, ws_lat_str, wn_lon_str, wn_lat_str, en_lon_str, en_lat_str, es_lon_str, es_lat_str, pixstr;
   
        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
                  
     
        num_gridlines = l1cfile->num_gridlines;
        num_scans = l1cfile->nscan;
        nbinx = l1cfile->nbinx;
        num_pixels = l1cfile->npix;
        n_files = l1cfile->ifiles.size();
        size_t sfiles = 0;

        latcornBox = allocate2d_double(4, 9);//4 corners, 9 grid cells
        loncornBox = allocate2d_double(4, 9);


        if (l1cfile->selgran[0] < 0)//process all granules
        {
            for (unsigned int j = 0;j < sfiles;j++) l1cfile->selgran[j] = j + 1;
        }
        else {
            while (l1cfile->selgran[sfiles] >0) {
                cout << "selected granule #..........................." << l1cfile->selgran[sfiles] << endl;
                sfiles++;
            }
        }

        //determine crossing time of the swath
        cout << "nfiles.for the swath." << n_files << endl;

        //determine # of gd groups to be processede  
        //gdlines_group=num_gridlines;//ONE BIG GROUP
        //--NORMALIZED  gridline longitude is -180 to 180 degrees before comparison with longitude extracted from image pixels---
        for (int i = 0; i < num_gridlines; i++) {
            for (int j = 0; j < nbinx; j++) {
                if (l1cfile->lon_gd[i][j] < -180.)l1cfile->lon_gd[i][j] += 360.;
                if (l1cfile->lon_gd[i][j] > 180.)l1cfile->lon_gd[i][j] -= 360.;
            }        
         }  
        //get min/max lat/lon of the swath
        // Assume first element as maximum and minimum
        maxv = l1cfile->lat_gd[0][0];
        minv = l1cfile->lat_gd[0][0];

        //Find maximum and minimum in all array elements.
        for (int i = 0; i < num_gridlines; i++) {
            for (int j = 0; j < nbinx; j++) {
                if (l1cfile->lat_gd[i][j] > maxv)
                    maxv = l1cfile->lat_gd[i][j];
                if (l1cfile->lat_gd[i][j] < minv)
                    minv = l1cfile->lat_gd[i][j];
            }        
         }      
        latmin_swt = minv;
        latmax_swt = maxv;
        maxv = l1cfile->lon_gd[0][0];
        minv = l1cfile->lon_gd[0][0];
        //Find maximum and minimum in all array elements.
        for (int i = 0; i < num_gridlines; i++) {
            for (int j = 0; j < nbinx; j++) {
                if (l1cfile->lon_gd[i][j] > maxv)
                    maxv = l1cfile->lon_gd[i][j];
                if (l1cfile->lon_gd[i][j] < minv)
                    minv = l1cfile->lon_gd[i][j];
            }
          }
        lonmin_swt = minv;
        lonmax_swt = maxv;

        //check lat bound min/max of each gd group---
        //
         //recompute number of gridlines based on lat limits -80 to 80 degrees  
        //asc orbit goes from negative to positive latitude....
        NY = num_gridlines;
 

        // cout<<"gridlines limited to -80 to 80 degrees of latitutde.."<<"gdline ini.."<<NY1<<"gdline end.."<<NY2<<endl;
        NY1=0;
        l1cfile->NY1 = NY1;
        NY2=NY-1;
        l1cfile->NY2 = NY2;

  //      cout << "NY.." << NY << "NY1.." << NY1 << "NY2.." << NY2 << endl;

        num_gridlines = NY2 - NY1 + 1;//this is the # gridlines after -80 to 80 lat constraint---

        num_blue_bands = l1cfile->nband_blue;//this includes uv + visible bands
        num_red_bands = l1cfile->nband_red;
        num_SWIR_bands = l1cfile->nband_swir;

        //*********** SORTING LAT/LON ASCENDING AND TRACK INDEXES ****************************
        //***********************************************************************************

        //Create index matrix for i and j of each lat and lon L1C grid (-80 to 80)
        //num_gridlines x nbinx
  
 
        //  lon_asort=allocate2d_float(num_gridlines,nbinx);
             //fill indexe


        fi = l1cfile->selgran[0];
        str = l1cfile->ifiles[fi - 1];


        for (pix = 0;pix < num_pixels;pix++) {
            if (l1cstr->latpix[pix] < latmin_swt) latmin_swt = l1cstr->latpix[pix];
            if (l1cstr->latpix[pix] > latmax_swt) latmax_swt = l1cstr->latpix[pix];
            if (l1cstr->lonpix[pix] < lonmin_swt) lonmin_swt = l1cstr->lonpix[pix];
            if (l1cstr->lonpix[pix] > lonmax_swt) lonmax_swt = l1cstr->lonpix[pix];
        }


             cout<<"Processing line #..."<<sline+1<<endl;

//*********** BIG LOOP M_PIXEL *****************************************************************
        for (pix = 0;pix < num_pixels;pix++) {

            bb = 0, gd_row = 0, gd_col = 0;
            flag_inpix = 0;

            latcornBox = allocate2d_double(4, 9);//4 corners, 9 grid cells
            loncornBox = allocate2d_double(4, 9);

            for (int i = 0;i < 4;i++) {
                for (int j = 0;j < 9;j++) {
                    latcornBox[i][j] = 0.;
                    loncornBox[i][j] = 0.;
                }            
         }

//            cout<<"#gridlines.."<<num_gridlines<<"nbinx.."<<nbinx<<"pix......"<<pix+1<< "latpix..." << l1cstr->latpix[pix] << "lonpix..." << l1cstr->lonpix[pix] << endl;

            boolbin1 = sbs2_l1c(l1cinput, num_gridlines, nbinx, l1cfile->lat_asort, l1cfile->index_xy, l1cstr->latpix[pix], l1cstr->lonpix[pix], l1cfile->lon_gd, &gd_row, &gd_col);

            //           if (sline % 1 == 0){     cout<<"LINE NUMBER------------------------------------------------------------------------------------------------------------------------#"<<sline+1<<"pix #........"<<pix+1<<endl;
 
            if (boolbin1 == 1){
                flag_inpix = 1;
             //    cout << "boolbin1..." << boolbin1 <<"pix..."<<pix+1<< "row.." << gd_row << "col..." << gd_col << "latpix..." << l1cstr->latpix[pix] << "lonpix..." << l1cstr->lonpix[pix] << endl;
                 }



            //************ BINNING ***************************
         //assign identified pixel to Lt and bin stat arrays-----         


            if (flag_inpix == 1 && gd_row >= 2 && gd_col >= 2 && gd_col < nbinx - 2 && gd_row < num_gridlines - 2) {
                Polygon_t pixelPoly;
                Polygon_t gridPoly;

                inpix++;

            //across-track distance    
                deltaphi = (l1cstr->latpix[pix] - l1cstr->latpix[pix + 1]) * degrad;
                deltalam = (l1cstr->lonpix[pix] - l1cstr->lonpix[pix + 1]) * degrad;
                aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cstr->latpix[pix] * degrad) * cos(l1cstr->latpix[pix + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
                bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
                dist_u = Re * bterm;         //horizontal distance Harversine in km
            //along-track distance
                deltaphi = (l1cstr->latpix[pix] - l1cstr->latpix2[pix]) * degrad;
                deltalam = (l1cstr->lonpix[pix] - l1cstr->lonpix2[pix]) * degrad;
                aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cstr->latpix[pix] * degrad) * cos(l1cstr->latpix2[pix] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
                bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
                dist_v = Re * bterm;

//                cout<<"deltaphi.."<<deltaphi<<"deltalam.."<<deltalam<<"cos(l1cstr->latpix2[pix] * degrad).."<<cos(l1cstr->latpix2[pix] * degrad)<<"cos(l1cstr->latpix[pix] * degrad..."<<cos(l1cstr->latpix[pix] * degrad)<<endl;
//                cout<<"l1cstr->latpix[pix].."<<l1cstr->latpix[pix]<<"l1cstr->latpix[pix+1].."<<l1cstr->latpix[pix+1]<<"l1cstr->lonpix[pix].."<<l1cstr->lonpix[pix]<<"l1cstr->lonpix[pix+1].."<<l1cstr->lonpix[pix+1]<<endl;

                //pixel azimuth
                azpix = atan2(sin(deltalam) * cos(l1cstr->latpix2[pix] * degrad), cos(l1cstr->latpix[pix] * degrad) * sin(l1cstr->latpix2[pix] * degrad) - sin(l1cstr->latpix[pix] * degrad) * cos(l1cstr->latpix2[pix] * degrad) * cos(deltaphi));

                if (azpix > M_PI || azpix < -M_PI) {
                    cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "azpixc in degrees.." << azpix * 180 / M_PI << endl;
                    exit(1);
                }

                azpix = azpix * 180. / M_PI;

                //      cout<<"line.."<<sline+1<<"pix #:.."<<pix+1<<"pix size.U in km."<<dist_u<<"pixsize.V in km."<<dist_v<<"azpix before  pix_corners4_l1c....in degrees."<<azpix<<endl;
               //       cout<<"FOUND grid cell!!-----row.."<<gd_row<<"col.."<<gd_col<<"latpix.."<<l1cstr->latpix[pix]<<"lonpix.."<<l1cstr->lonpix[pix]<<"latpix2.."<<l1cstr->latpix2[pix]<<"lonpix2.."<<l1cstr->lonpix2[pix]<<endl;                             

                //******* M_PIXEL CORNERS ***********************************************************************************************
                dist_corn = 0.5 * (sqrt(dist_u * dist_u + dist_v * dist_v)) * 1000;//distane center to corner in meters
                theta = atan2(dist_v, dist_u);
                thetagrad = 2 * theta * 180 / M_PI;
                thetares = (M_PI / 2 - theta) * 180 / M_PI;
                //nw corner 
                azpixc = (azpix - thetares) * degrad;//pixel bearing
                if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in nw corner..azpixc" << azpixc << endl;exit(1); }
                azpixc = azpixc * 180 / M_PI;
                //      cout<<"azpixc nw...."<<azpixc<<endl;
                geod.Direct(l1cstr->latpix[pix], l1cstr->lonpix[pix], azpixc, dist_corn, lat_cnw, lon_cnw);
                //ne corner (sw-180)
                azpixc = (azpix + thetares) * degrad;//pixel bearing
                if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in ne corner..azpixc" << azpixc << endl;exit(1); }
                azpixc = azpixc * 180 / M_PI;
                //      cout<<"azpixc ne...."<<azpixc<<endl;
                geod.Direct(l1cstr->latpix[pix], l1cstr->lonpix[pix], azpixc, dist_corn, lat_cne, lon_cne);
                //sw corner               
                azpixc = (azpix - thetares - thetagrad) * degrad;//pixel bearing
                if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in sw corner..azpixc" << azpixc << endl;exit(1); }
                azpixc = azpixc * 180 / M_PI;
                //       cout<<"azpixc sw....."<<azpixc<<endl;
                geod.Direct(l1cstr->latpix[pix], l1cstr->lonpix[pix], azpixc, dist_corn, lat_csw, lon_csw);
                //se corner (nw-180)
                azpixc = (azpix + thetares + thetagrad) * degrad;//pixel bearingi
                if (azpixc<-M_PI || azpixc>M_PI) { cout << "error in se corner..azpixc" << azpixc << endl;exit(1); }
            //    azpixc = (azpix + 45 + 90) * degrad;//pixel bearing
                azpixc = azpixc * 180 / M_PI;
                //       cout<<"azpix se...."<<azpixc<<endl;
                geod.Direct(l1cstr->latpix[pix], l1cstr->lonpix[pix], azpixc, dist_corn, lat_cse, lon_cse);



                //                            cout<<"sline................................................................................................................"<<sline+1<<endl;
                                       //     cout<<"sline.."<<sline+1<<"dist_corn.."<<dist_corn<<"lat_cnw.."<<lat_cnw<<"lon_cnw.."<<lon_cnw<<"lat_cne.."<<lat_cne<<"lon_cne.."<<lon_cne<<"lat_csw.."<<lat_csw<<"lon_csw.."<<lon_csw<<"lat_cse.."<<lat_cse<<"lon_cse.."<<lon_cse<<endl;

                                     //       ws_lon_str=std::to_string(lon_csw);
                                     //       ws_lat_str=std::to_string(lat_csw);
                                     //       wn_lon_str=std::to_string(lon_cnw);
                                     //       wn_lat_str=std::to_string(lat_cnw);
                                     //       en_lon_str=std::to_string(lon_cne);
                                     //       en_lat_str=std::to_string(lat_cne);
                                     //       es_lon_str=std::to_string(lon_cse);
                                     //       es_lat_str=std::to_string(lat_cse);

                                     //       pixstr="POLYGON(("+ws_lon_str+" "+ws_lat_str+","+wn_lon_str+" "+wn_lat_str+","+en_lon_str+" "+en_lat_str+","+es_lon_str+" "+es_lat_str+","+ws_lon_str+" "+ws_lat_str+"))";
                                     //       bg::read_wkt(pixstr,pixelPoly);

                bg::append(pixelPoly.outer(), Point_t(lon_csw, lat_csw));
                bg::append(pixelPoly.outer(), Point_t(lon_cnw, lat_cnw));
                bg::append(pixelPoly.outer(), Point_t(lon_cne, lat_cne));
                bg::append(pixelPoly.outer(), Point_t(lon_cse, lat_cse));
                bg::append(pixelPoly.outer(), Point_t(lon_csw, lat_csw));



                // make sure the polygon is defined properly
                bg::correct(pixelPoly);

                if (lon_csw > lon_cse || lon_cnw > lon_cne || lat_csw > lat_cnw || lat_cse > lat_cne) {
                    cout << "lon_csw>lon_cse | lon_cnw>lon_cne | lat_csw>lat_cnw | lat_cse>lat_cne//" << lon_csw << lon_cse << lon_cnw << lon_cne << lat_csw << lat_cnw << lat_cse << lat_cne << endl;
                    cout << "corners are switched for pixel.." << pix + 1 << endl;
                }


                lon_csw = 0, lat_csw = 0, lon_cnw = 0, lat_cnw = 0, lon_cne = 0, lat_cne = 0, lon_cse = 0, lat_cse = 0;
                //    cout<<"azpixc nw...."<<azpixc<<"lat.."<<l1cstr->latpix[pix]<<"lon.."<<l1cstr->lonpix[pix]<<"dist_corn.."<<dist_corn<<"lat_cnw.."<<lat_cnw<<"lon_cnw.."<<lon_cnw<<endl;


                if (gd_row >= 1 && gd_col >= 1 && gd_row < num_gridlines - 2 && gd_col < nbinx - 2) {

                    //           cout<<"sline.."<<sline<<"gd_row.."<<gd_row<<"gd_col.."<<gd_col<<"bb.."<<bb<<endl;
                    gwindowTopix_l1c2(l1cfile, l1cinput, gd_row, gd_col,latcornBox, loncornBox);

 //                  cout<<"ok after gwindowTopix.."<<"pix..."<<pix+1<<endl;

                    for (int i = 0;i < 9;i++) { //for each intersection
              //          corn1_lon_str=std::to_string(loncornBox[0][i]);
              //          corn1_lat_str=std::to_string(latcornBox[0][i]);
              //          corn2_lon_str=std::to_string(loncornBox[1][i]);
              //          corn2_lat_str=std::to_string(latcornBox[1][i]);
              //          corn3_lon_str=std::to_string(loncornBox[2][i]);
              //          corn3_lat_str=std::to_string(latcornBox[2][i]);
             //           corn4_lon_str=std::to_string(loncornBox[3][i]);
              //          corn4_lat_str=std::to_string(latcornBox[3][i]);
             //           gridstr="POLYGON(("+corn1_lon_str+" "+corn1_lat_str+","+corn2_lon_str+" "+corn2_lat_str+","+corn3_lon_str+" "+corn3_lat_str+","+corn4_lon_str+" "+corn4_lat_str+","+corn1_lon_str+" "+corn1_lat_str+"))";
              //          bg::read_wkt(gridstr,gridPoly);


                        bg::append(gridPoly.outer(), Point_t(loncornBox[0][i], latcornBox[0][i]));
                        bg::append(gridPoly.outer(), Point_t(loncornBox[1][i], latcornBox[1][i]));
                        bg::append(gridPoly.outer(), Point_t(loncornBox[2][i], latcornBox[2][i]));
                        bg::append(gridPoly.outer(), Point_t(loncornBox[3][i], latcornBox[3][i]));
                        bg::append(gridPoly.outer(), Point_t(loncornBox[0][i], latcornBox[0][i]));
                        bg::correct(gridPoly);

                        ingeom = within(pixelPoly, gridPoly);
                        binAreagrid = bg::area(gridPoly);
                        binAreapix = bg::area(pixelPoly);//         

             //           if (sline % 10 == 0 & binAreagrid>10000000)            cout<<"sline.."<<sline+1<<"pix.."<<pix+1<<"gridcell#.."<<i<<"fully inside????.."<<ingeom<<"binAreagrid..."<<binAreagrid<<"binAreapix..."<<binAreapix<<endl;        
                  //compute intersections-----
//                   
                        intersectArea = 0., pc = 1;
                        if (ingeom > 0) { //pixel fully inside the grid cell                        
                            intersectArea = binAreapix;

                            if (l1cstr->Lt_blue[bb][pix] > 0.) {
                                Ltfracsum[gd_row][gd_col] += l1cstr->Lt_blue[bb][pix] * (intersectArea / binAreagrid);
                                areabinsum[gd_row][gd_col] += (intersectArea / binAreagrid);
                                nobs_perbin[gd_row][gd_col] += 1;
                            }
                        }
                        else {
                            if (!bg::disjoint(pixelPoly, gridPoly) && binAreagrid > 0) {
                                std::deque<Polygon_t> output;//this is a list
                                if (bg::intersection(pixelPoly, gridPoly, output)) {
                                    BOOST_FOREACH(Polygon_t const& p, output)
                                    {
                                        intersectArea = bg::area(p);
                                        if (intersectArea > 0.0 && pc == 1) {

                                            if (l1cstr->Lt_blue[bb][pix] > 0.) {
                                                Ltfracsum[gd_row][gd_col] += l1cstr->Lt_blue[bb][pix] * (intersectArea / binAreagrid);
                                                areabinsum[gd_row][gd_col] += (intersectArea / binAreagrid);
                                                nobs_perbin[gd_row][gd_col] += 1;
                                            }


                                            if (intersectArea / binAreagrid < 0.000001) {
                                                Ltfracsum[gd_row][gd_col] += 0.;
                                                areabinsum[gd_row][gd_col] += 0.;
                                                nobs_perbin[gd_row][gd_col] += 0.;
                                            }
                                            pc++;
                                        }
                                    }
                                }
                            }
                        }//end else

//                                       cout<<"gridcell#.."<<i<<"fully inside????.."<<ingeom<<"areafrac..."<<intersectArea / binAreagrid<<"Ltfrac..."<<Ltfracsum[gd_row][gd_col]<<"areabinsum.."<<areabinsum[gd_row][gd_col]<<"nobsperbin.."<<nobs_perbin[gd_row][gd_col]<<endl;
               //              cout<<"sline........"<<sline+1<<"pix#.........................."<<pix+1<<"gridcell#.."<<i<<endl;
                    }//end for intersections


                }//end if gdrow and gdcol are ok





/*
                         while (bb<num_blue_bands){
                          if(l1cstr->Lt_blue[bb][pix]>0.){
                            binmean_Lt[bb][bin_ypix-1][bin_xpix-1]=binmean_Lt[bb][bin_ypix-1][bin_xpix-1]+l1cstr->Lt_blue[bb][pix];
                            bincount_Lt[bb][bin_ypix-1][bin_xpix-1]+=1;
                           }
                        bb++;
                        }//end while blue band band
//RED---
                         while (rb<num_red_bands){
                          if(l1cstr->Lt_red[rb][pix]>0.){
                            binmean_Lt2[rb][bin_ypix-1][bin_xpix-1]=binmean_Lt2[rb][bin_ypix-1][bin_xpix-1]+l1cstr->Lt_red[rb][pix];
                            bincount_Lt2[rb][bin_ypix-1][bin_xpix-1]+=1;
                           }
                         rb++;
                         }//end red bands

//SWIR----
                        while (swb<num_SWIR_bands){
                          if(l1cstr->Lt_SWIR[swb][pix]>0.){
                            binmean_Lt3[swb][bin_ypix-1][bin_xpix-1]=binmean_Lt3[swb][bin_ypix-1][bin_xpix-1]+l1cstr->Lt_SWIR[swb][pix];
                            bincount_Lt3[swb][bin_ypix-1][bin_xpix-1]+=1;
                           }
                         swb++;
                         }//end swir bands

*/

            }//end if flag_inpix



//********** pixel AREA WEIGHTING ************************
            if (flag_inpix == 0) outpix++;

            totpix++;

            delete[](latcornBox);
            delete[](loncornBox);
            latcornBox = nullptr;
            loncornBox = nullptr;
        }//end pixels loop


        cout << "total inpix.." << inpix << "total pix.." << totpix << endl;

 
        return 0;
    }




    int32_t L1C::gwindowTopix_l1c2(l1c_filehandle* l1cfile, L1C_input* l1cinput, short row, short col, double** latcornBox, double** loncornBox) {
        double  ws_lat, wn_lat, es_lat, en_lat, ws_lon, wn_lon, es_lon, en_lon;
        double ws_lat2, wn_lat2, es_lat2, en_lat2, ws_lon2, wn_lon2, es_lon2, en_lon2;
        double ws_lat3, wn_lat3, es_lat3, en_lat3, ws_lon3, wn_lon3, es_lon3, en_lon3;
        double ws_lat4, wn_lat4, es_lat4, en_lat4, ws_lon4, wn_lon4, es_lon4, en_lon4;
        double ws_lat5, wn_lat5, es_lat5, en_lat5, ws_lon5, wn_lon5, es_lon5, en_lon5;
        double ws_lat6, wn_lat6, es_lat6, en_lat6, ws_lon6, wn_lon6, es_lon6, en_lon6;
        double ws_lat7, wn_lat7, es_lat7, en_lat7, ws_lon7, wn_lon7, es_lon7, en_lon7;
        double ws_lat8, wn_lat8, es_lat8, en_lat8, ws_lon8, wn_lon8, es_lon8, en_lon8;
        double ws_lat9, wn_lat9, es_lat9, en_lat9, ws_lon9, wn_lon9, es_lon9, en_lon9;
        //   shape->rowcol2bounds(row, col, n, s, e, w);
        float binres = 0, azgc, deltaphi, deltalam, azpixc;
    
        float aterm, bterm, dist_u, dist_v, theta, thetares, thetagrad, dist_corn2;
        float az1, az2, az3, az4;


        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());


        binres = (l1cinput->grid_resolution);//in km

    //    cout<<"gwindowTopix_l1c...........................gd_row......"<<row<<"gd_col...."<<col<<endl;

    //*******************************************************************************************
    //center center bin
    //********************************************************************************************************

        deltaphi = (l1cfile->lat_gd[row][col] - l1cfile->lat_gd[row][col + 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row][col] - l1cfile->lon_gd[row][col + 1]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row][col] * degrad) * cos(l1cfile->lat_gd[row][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row][col] - l1cfile->lat_gd[row + 1][col]) * degrad;
        deltalam = (l1cfile->lon_gd[row][col] - l1cfile->lon_gd[row + 1][col]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row][col] * degrad) * cos(l1cfile->lat_gd[row + 1][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        //    ad=0.5*(sqrt(dist_u*dist_u+dist_v*dist_v)/Re);   
    
        dist_corn2 = 0.5 * sqrt(binres * binres + binres * binres) * 1000;

     //azimuth grid cell calculation
        azgc=atan2(sin(deltalam)*cos(l1cfile->lat_gd[row+1][col]*degrad),cos(l1cfile->lat_gd[row][col]*degrad)*sin(l1cfile->lat_gd[row+1][col]*degrad)-sin(l1cfile->lat_gd[row][col]*degrad)*cos(l1cfile->lat_gd[row+1][col]*degrad)*cos(deltaphi));

        if (azgc > M_PI || azgc < -M_PI){cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "azpixc in degrees.." << azgc * 180 / M_PI << endl;exit(1);}

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;
        azgc = azgc * 180 / M_PI;

        //nw corner 
        azpixc = (azgc - thetares) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "NW azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az1 = azpixc;
        geod.Direct(l1cfile->lat_gd[row][col], l1cfile->lon_gd[row][col], azpixc, dist_corn2, wn_lat, wn_lon);

        //    cout<<"azgc.Center-Center."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //    cout<<"res1..nw."<<res1<<endl;

        //ne corner 
        azpixc = (azgc + thetares) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "NE azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az2 = azpixc;
        geod.Direct(l1cfile->lat_gd[row][col], l1cfile->lon_gd[row][col], azpixc, dist_corn2, en_lat, en_lon);

        //  cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
      //sw corner 
        azpixc = (azgc - thetares - thetagrad) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "SW azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az3 = azpixc;
        geod.Direct(l1cfile->lat_gd[row][col], l1cfile->lon_gd[row][col], azpixc, dist_corn2, ws_lat, ws_lon);

        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;

       //se corner 
        azpixc = (azgc + thetares + thetagrad) * degrad;//pixel bearing
        if (azpixc > M_PI || azpixc < -M_PI) {
            cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "SE azpixc in degrees.." << azpixc * 180 / M_PI << endl;
            exit(1);
        }

        azpixc = azpixc * 180 / M_PI;
        az4 = azpixc;
        geod.Direct(l1cfile->lat_gd[row][col], l1cfile->lon_gd[row][col], azpixc, dist_corn2, es_lat, es_lon);
        //   geod.Inverse(lat_gd[row][col],lon_gd[row][col],es_lat,es_lon,res4);

        latcornBox[0][0] = ws_lat;
        latcornBox[1][0] = wn_lat;
        latcornBox[2][0] = en_lat;
        latcornBox[3][0] = es_lat;
        loncornBox[0][0] = ws_lon;
        loncornBox[1][0] = wn_lon;
        loncornBox[2][0] = en_lon;
        loncornBox[3][0] = es_lon;

        //***********************************************************************************************************
           //left center bin
            //********************************************************************************************************

        deltaphi = (l1cfile->lat_gd[row][col - 1] - l1cfile->lat_gd[row][col]) * degrad;
        deltalam = (l1cfile->lon_gd[row][col - 1] - l1cfile->lon_gd[row][col]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row][col - 1] * degrad) * cos(l1cfile->lat_gd[row][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row][col - 1] - l1cfile->lat_gd[row + 1][col - 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row][col - 1] - l1cfile->lon_gd[row + 1][col - 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row][col - 1] * degrad) * cos(l1cfile->lat_gd[row + 1][col - 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;


        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row][col - 1], l1cfile->lon_gd[row][col - 1], azpixc, dist_corn2, wn_lat2, wn_lon2);

        //   cout<<"azgc..CENTER-LEFT"<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;
       //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row][col - 1], l1cfile->lon_gd[row][col - 1], azpixc, dist_corn2, en_lat2, en_lon2);

        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row][col - 1], l1cfile->lon_gd[row][col - 1], azpixc, dist_corn2, ws_lat2, ws_lon2);

        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row][col - 1], l1cfile->lon_gd[row][col - 1], azpixc, dist_corn2, es_lat2, es_lon2);

        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;


        latcornBox[0][1] = ws_lat2;
        latcornBox[1][1] = wn_lat2;
        latcornBox[2][1] = wn_lat;
        latcornBox[3][1] = ws_lat;
        loncornBox[0][1] = ws_lon2;
        loncornBox[1][1] = wn_lon2;
        loncornBox[2][1] = wn_lon;
        loncornBox[3][1] = ws_lon;


        //******************************************************************************************************************************
        //right center bin 
        //*******************************************************************************************************************************
        //***************************************************************************
         //right center bin
        //**************************************************************************8


        deltaphi = (l1cfile->lat_gd[row][col + 1] - l1cfile->lat_gd[row][col + 2]) * degrad;
        deltalam = (l1cfile->lon_gd[row][col + 1] - l1cfile->lon_gd[row][col + 2]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row][col + 1] * degrad) * cos(l1cfile->lat_gd[row][col + 2] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row][col + 1] - l1cfile->lat_gd[row + 1][col + 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row][col + 1] - l1cfile->lon_gd[row + 1][col + 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row][col + 1] * degrad) * cos(l1cfile->lat_gd[row + 1][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row][col + 1], l1cfile->lon_gd[row][col + 1], azpixc, dist_corn2, wn_lat3, wn_lon3);
        //   cout<<"azgc..CENTER-RIGHT..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;

       //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row][col + 1], l1cfile->lon_gd[row][col + 1], azpixc, dist_corn2, en_lat3, en_lon3);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row][col + 1], l1cfile->lon_gd[row][col + 1], azpixc, dist_corn2, ws_lat3, ws_lon3);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row][col + 1], l1cfile->lon_gd[row][col + 1], azpixc, dist_corn2, es_lat3, es_lon3);
        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

        latcornBox[0][2] = es_lat;
        latcornBox[1][2] = en_lat;
        latcornBox[2][2] = en_lat3;
        latcornBox[3][2] = es_lat3;
        loncornBox[0][2] = es_lon;
        loncornBox[1][2] = en_lon;
        loncornBox[2][2] = en_lon3;
        loncornBox[3][2] = es_lon3;

        //******************************************************************************************************
         //upper center bin
        //*****************************************************************************************************


        deltaphi = (l1cfile->lat_gd[row + 1][col] - l1cfile->lat_gd[row + 1][col + 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row + 1][col] - l1cfile->lon_gd[row + 1][col + 1]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row + 1][col] * degrad) * cos(l1cfile->lat_gd[row + 1][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row + 2][col] - l1cfile->lat_gd[row + 1][col]) * degrad;
        deltalam = (l1cfile->lon_gd[row + 2][col] - l1cfile->lon_gd[row + 1][col]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row + 1][col] * degrad) * cos(l1cfile->lat_gd[row + 2][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;
        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row + 1][col], l1cfile->lon_gd[row + 1][col], azpixc, dist_corn2, wn_lat4, wn_lon4);
        //   cout<<"azgc..CENTER-UPPER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;
       //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row + 1][col], l1cfile->lon_gd[row + 1][col], azpixc, dist_corn2, en_lat4, en_lon4);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row + 1][col], l1cfile->lon_gd[row + 1][col], azpixc, dist_corn2, ws_lat4, ws_lon4);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row + 1][col], l1cfile->lon_gd[row + 1][col], azpixc, dist_corn2, es_lat4, es_lon4);

        latcornBox[0][3] = wn_lat;
        latcornBox[1][3] = wn_lat4;
        latcornBox[2][3] = en_lat4;
        latcornBox[3][3] = en_lat;
        loncornBox[0][3] = wn_lon;
        loncornBox[1][3] = wn_lon4;
        loncornBox[2][3] = en_lon4;
        loncornBox[3][3] = en_lon;

        //********************************************************
        //upper left  bin
        //*******************************************************
     

        deltaphi = (l1cfile->lat_gd[row + 1][col - 1] - l1cfile->lat_gd[row + 1][col]) * degrad;
        deltalam = (l1cfile->lon_gd[row + 1][col - 1] - l1cfile->lon_gd[row + 1][col]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row + 1][col - 1] * degrad) * cos(l1cfile->lat_gd[row + 1][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row + 2][col - 1] - l1cfile->lat_gd[row + 1][col - 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row + 2][col - 1] - l1cfile->lon_gd[row + 1][col - 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row + 1][col - 1] * degrad) * cos(l1cfile->lat_gd[row + 2][col - 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row + 1][col - 1], l1cfile->lon_gd[row + 1][col - 1], azpixc, dist_corn2, wn_lat5, wn_lon5);
        //   cout<<"azgc..LEFT-UPPER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;
       //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row + 1][col - 1], l1cfile->lon_gd[row + 1][col - 1], azpixc, dist_corn2, en_lat5, en_lon5);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row + 1][col - 1], l1cfile->lon_gd[row + 1][col - 1], azpixc, dist_corn2, ws_lat5, ws_lon5);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row + 1][col - 1], l1cfile->lon_gd[row + 1][col - 1], azpixc, dist_corn2, es_lat5, es_lon5);
        //  cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;

        latcornBox[0][4] = wn_lat2;
        latcornBox[1][4] = wn_lat5;
        latcornBox[2][4] = wn_lat4;
        latcornBox[3][4] = wn_lat;
        loncornBox[0][4] = wn_lon2;
        loncornBox[1][4] = wn_lon5;
        loncornBox[2][4] = wn_lon4;
        loncornBox[3][4] = wn_lon;

        //**************************************************************************************
        //upper right bin
        //**************************************************************************************
     
        deltaphi = (l1cfile->lat_gd[row + 1][col + 1] - l1cfile->lat_gd[row + 1][col + 2]) * degrad;
        deltalam = (l1cfile->lon_gd[row + 1][col + 1] - l1cfile->lon_gd[row + 1][col + 2]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row + 1][col + 1] * degrad) * cos(l1cfile->lat_gd[row + 1][col + 2] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row + 2][col + 1] - l1cfile->lat_gd[row + 1][col + 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row + 2][col + 1] - l1cfile->lon_gd[row + 1][col + 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row + 1][col + 1] * degrad) * cos(l1cfile->lat_gd[row + 2][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row + 1][col + 1], l1cfile->lon_gd[row + 1][col + 1], azpixc, dist_corn2, wn_lat6, wn_lon6);

        //  cout<<"azgc..RIGHT-UPPER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
       //   cout<<"res1..nw."<<res1<<endl;
      //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row + 1][col + 1], l1cfile->lon_gd[row + 1][col + 1], azpixc, dist_corn2, en_lat6, en_lon6);

        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row + 1][col + 1], l1cfile->lon_gd[row + 1][col + 1], azpixc, dist_corn2, ws_lat6, ws_lon6);

        //  cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
      //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row + 1][col + 1], l1cfile->lon_gd[row + 1][col + 1], azpixc, dist_corn2, es_lat6, es_lon6);

        //  cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;
        latcornBox[0][5] = en_lat;
        latcornBox[1][5] = en_lat4;
        latcornBox[2][5] = en_lat6;
        latcornBox[3][5] = en_lat3;
        loncornBox[0][5] = en_lon;
        loncornBox[1][5] = en_lon4;
        loncornBox[2][5] = en_lon6;
        loncornBox[3][5] = en_lon3;

        //************************************************************
        //lower center  bin
        //*************************************************************

     

        deltaphi = (l1cfile->lat_gd[row - 1][col] - l1cfile->lat_gd[row - 1][col + 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row - 1][col] - l1cfile->lon_gd[row - 1][col + 1]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row - 1][col] * degrad) * cos(l1cfile->lat_gd[row - 1][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row - 1][col] - l1cfile->lat_gd[row][col]) * degrad;
        deltalam = (l1cfile->lon_gd[row - 1][col] - l1cfile->lon_gd[row][col]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row - 1][col] * degrad) * cos(l1cfile->lat_gd[row][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row - 1][col], l1cfile->lon_gd[row - 1][col], azpixc, dist_corn2, wn_lat7, wn_lon7);
        //   cout<<"azgc..LOWER CENTER..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;
       //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row - 1][col], l1cfile->lon_gd[row - 1][col], azpixc, dist_corn2, en_lat7, en_lon7);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row - 1][col], l1cfile->lon_gd[row - 1][col], azpixc, dist_corn2, ws_lat7, ws_lon7);

        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row - 1][col], l1cfile->lon_gd[row - 1][col], azpixc, dist_corn2, es_lat7, es_lon7);

        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;
        latcornBox[0][6] = ws_lat7;
        latcornBox[1][6] = ws_lat;
        latcornBox[2][6] = es_lat;
        latcornBox[3][6] = es_lat7;
        loncornBox[0][6] = ws_lon7;
        loncornBox[1][6] = ws_lon;
        loncornBox[2][6] = es_lon;
        loncornBox[3][6] = es_lon7;

        //******************************************************************
        //lower left bin
        //****************************************************************
    
        deltaphi = (l1cfile->lat_gd[row - 1][col - 1] - l1cfile->lat_gd[row - 1][col]) * degrad;
        deltalam = (l1cfile->lon_gd[row - 1][col - 1] - l1cfile->lon_gd[row - 1][col]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row - 1][col - 1] * degrad) * cos(l1cfile->lat_gd[row - 1][col] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row - 1][col - 1] - l1cfile->lat_gd[row][col - 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row - 1][col - 1] - l1cfile->lon_gd[row][col - 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row - 1][col - 1] * degrad) * cos(l1cfile->lat_gd[row][col - 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row - 1][col - 1], l1cfile->lon_gd[row - 1][col - 1], azpixc, dist_corn2, wn_lat8, wn_lon8);
        //   cout<<"azgc..LOWER LEFT..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //   cout<<"res1..nw."<<res1<<endl;
       //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row - 1][col - 1], l1cfile->lon_gd[row - 1][col - 1], azpixc, dist_corn2, en_lat8, en_lon8);
        //   cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
       //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row - 1][col - 1], l1cfile->lon_gd[row - 1][col - 1], azpixc, dist_corn2, ws_lat8, ws_lon8);
        //   cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
       //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row - 1][col - 1], l1cfile->lon_gd[row - 1][col - 1], azpixc, dist_corn2, es_lat8, es_lon8);

        //  cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;
        latcornBox[0][7] = ws_lat8;
        latcornBox[1][7] = ws_lat2;
        latcornBox[2][7] = ws_lat;
        latcornBox[3][7] = ws_lat7;
        loncornBox[0][7] = ws_lon8;
        loncornBox[1][7] = ws_lon2;
        loncornBox[2][7] = ws_lon;
        loncornBox[3][7] = ws_lon7;

        //******************************************************************
        //lower right bin
        //****************************************************************

        deltaphi = (l1cfile->lat_gd[row - 1][col + 1] - l1cfile->lat_gd[row - 1][col + 2]) * degrad;
        deltalam = (l1cfile->lon_gd[row - 1][col + 1] - l1cfile->lon_gd[row - 1][col + 2]) * degrad;
        //across-track grid size
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row - 1][col + 1] * degrad) * cos(l1cfile->lat_gd[row - 1][col + 2] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_u = Re * bterm;         //horizontal distance Harversine in km

        deltaphi = (l1cfile->lat_gd[row - 1][col + 1] - l1cfile->lat_gd[row][col + 1]) * degrad;
        deltalam = (l1cfile->lon_gd[row - 1][col + 1] - l1cfile->lon_gd[row][col + 1]) * degrad;
        //along-track distance
        aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cfile->lat_gd[row - 1][col + 1] * degrad) * cos(l1cfile->lat_gd[row][col + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
        bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
        dist_v = Re * bterm;

        theta = atan2(dist_v, dist_u);
        thetagrad = 2 * theta * 180 / M_PI;
        thetares = (M_PI / 2 - theta) * 180 / M_PI;

        //nw corner
        azpixc = az1;
        geod.Direct(l1cfile->lat_gd[row - 1][col + 1], l1cfile->lon_gd[row - 1][col + 1], azpixc, dist_corn2, wn_lat9, wn_lon9);
        //  cout<<"azgc..LOWER RIGHT..."<<azgc<<"azpixc.."<<azpixc<<"dist_corn.."<<dist_corn<<"dist_corn2.."<<dist_corn2<<endl;
        //  cout<<"res1..nw."<<res1<<endl;
      //ne corner 
        azpixc = az2;
        geod.Direct(l1cfile->lat_gd[row - 1][col + 1], l1cfile->lon_gd[row - 1][col + 1], azpixc, dist_corn2, en_lat9, en_lon9);
        //  cout<<"azpixc.."<<azpixc<<"res2 ne.."<<res2<<endl;
      //sw corner 
        azpixc = az3;
        geod.Direct(l1cfile->lat_gd[row - 1][col + 1], l1cfile->lon_gd[row - 1][col + 1], azpixc, dist_corn2, ws_lat9, ws_lon9);
        //  cout<<"azpixc.."<<azpixc<<"res3..sw.."<<res3<<endl;
      //se corner 
        azpixc = az4;
        geod.Direct(l1cfile->lat_gd[row - 1][col + 1], l1cfile->lon_gd[row - 1][col + 1], azpixc, dist_corn2, es_lat9, es_lon9);
        //   cout<<"azpixc.."<<azpixc<<"res4...se.."<<res4<<endl;
        latcornBox[0][8] = es_lat7;
        latcornBox[1][8] = es_lat;
        latcornBox[2][8] = es_lat3;
        latcornBox[3][8] = es_lat9;
        loncornBox[0][8] = es_lon7;
        loncornBox[1][8] = es_lon;
        loncornBox[2][8] = es_lon3;
        loncornBox[3][8] = es_lon9;

        return 0;
    }





//line by line version of xy_pixsize_sf----------- unlike sf3
  int32_t L1C::xy_pixsize_sf4(const char*ptstr,l1c_str *l1cstr,l1c_filehandle *l1cfile,L1C_input *l1cinput,float **pix_size_u,float **pix_size_v,float **Ltfracsum,float **areabinsum,float **nobs_perbin,float ****binLt,int ****bincount,size_t sline){
        float Re = 6378;
        float deltaphi = 0., deltalam = 0., aterm = 0., bterm, dist_u = 0., dist_v = 0., azpix = 0.;
        double areaFracBox[3][3];


        for (int ar = 0;ar < 3;ar++) {
            for (int ac = 0;ac < 3;ac++) {
                areaFracBox[ar][ac] = 0.;
                ac++;
            }
            ar++;
        }
 
        num_pixels = l1cfile->npix;
 
            if (sline % 100 == 1)  cout << "weight-binning sline..." << sline + 1 << endl;

            for (unsigned int pix = 0;pix < num_pixels - 1;pix++) {

                //across-track distance
                deltaphi = (l1cstr->latpix[pix] - l1cstr->latpix[pix + 1]) * degrad;
                deltalam = (l1cstr->lonpix[pix] - l1cstr->lonpix[pix + 1]) * degrad;
                aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cstr->latpix[pix] * degrad) * cos(l1cstr->latpix[pix + 1] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
                bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
                dist_u = Re * bterm;         //horizontal distance Harversine in km
            //along-track distance
                deltaphi = (l1cstr->latpix[pix] - l1cstr->latpix2[pix]) * degrad;
                deltalam = (l1cstr->lonpix[pix] - l1cstr->lonpix2[pix]) * degrad;
                aterm = sin(deltaphi / 2) * sin(deltaphi / 2) + cos(l1cstr->latpix[pix] * degrad) * cos(l1cstr->latpix2[pix] * degrad) * sin(deltalam / 2) * sin(deltalam / 2);
                bterm = 2 * atan2(sqrt(aterm), sqrt(1 - aterm));
                dist_v = Re * bterm;

                pix_size_u[sline][pix] = dist_u;
                pix_size_v[sline][pix] = dist_v;

                //pixel azimuth
                azpix = atan2(sin(deltalam) * cos(l1cstr->latpix2[pix] * degrad), cos(l1cstr->latpix[pix] * degrad) * sin(l1cstr->latpix2[pix] * degrad) - sin(l1cstr->latpix[pix] * degrad) * cos(l1cstr->latpix2[pix] * degrad) * cos(deltaphi));
                if (azpix > M_PI || azpix < -M_PI) {
                    cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "azpixc in degrees.." << azpix * 180 / M_PI << endl;
                    cout << "line.." << sline + 1 << "pix #:.." << pix + 1 << "pix size.U in km." << pix_size_u[sline][pix] << "pixsize.V in km." << pix_size_v[sline][pix] << "azpix before corners4.in degrees." << azpix << endl;
                    exit(1);
                     }


                azpix = azpix * 180. / M_PI;


                pix_corners4_l1c(l1cfile, l1cinput, dist_u, dist_v, azpix, sline, pix, l1cstr->latpix[pix], l1cstr->lonpix[pix], l1cstr->Lt_blue[0][pix], l1cfile->lat_asort, l1cfile->index_xy, l1cfile->lat_gd, l1cfile->lon_gd, areaFracBox, Ltfracsum, areabinsum, nobs_perbin);                                                

            }//end pixels

        return 0;
    }



int32_t L1C::write_L1C_granule2(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput,double *tmgv,float** lat_gd, float** lon_gd,float **alt_gd){
   int32_t num_gridlines, nbinx, NY1 = -1, NY2 = -1;
   int32_t NX, NY;
   const char* filename_lt;
   float **lat_out=nullptr,**lon_out=nullptr,**alt_out=nullptr;
   double *time_nad_out=nullptr;
   int asc_mode=-1;
   string senstr;
   double tg_ini,tg_end;
   int numgran;
   double tgridline,tfile_ini_sec,tfile_end_sec;
   int16_t *granid=nullptr;
   int32_t gransize;
   short **gdindex=nullptr;
   double **gdtime=nullptr;
   int16_t y_ini,mo_ini,d_ini,h_ini,mi_ini,y_end,mo_end,d_end,h_end,mi_end;
   double sec_ini,sec_end;
   string tswt_ini,tswt_end,tswt_ini_file;
   int logoff=-1;
   int16_t syear, smon, sday,syear2,smon2,sday2;
   double secs,secs2;
   double tgran_ini,tgran_ini_sec,tgran_end,tgran_end_sec;
   string gfull;
   int16_t gtime;
   int Nwest=-1, Neast=-1,Ngring=-1,midix=-1,dp=-1;
   int p,ix=-1;
   float latemp=-1,lontemp1=-1,lontemp2=-1,dlat_gd=-1,dlon_gd=-1,dlat20=-1,dlon20=-1,lon360=-1;
   int re=-1,rw=-1;   

        if(l1cinput->sensor==34){
             senstr="SPEXone";
             l1cfile->nbinx=25;
              }
        else if (l1cinput->sensor==30){
            senstr="OCI";
            l1cfile->nbinx=519;
             }
        else if (l1cinput->sensor==35){
            senstr="HARP2";
            l1cfile->nbinx=457;
            }
        else{cout<<"sensor by default is OCI option 2....."<<endl;
            senstr="OCI";
             l1cfile->nbinx=519;      
          }

        asc_mode=l1cfile->orb_dir;
        gransize=l1cfile->gransize;
      
        num_gridlines = l1cfile->num_gridlines;
        cout<<"swath#....................................................................."<<swtd<<"asc_mode..."<<asc_mode<<endl;

        granid=(int16_t*)calloc(num_gridlines,sizeof(int16_t));

        numgran=l1cfile->numgran;

        nbinx = l1cfile->nbinx;

        gdtime = allocate2d_double(numgran,2);
        gdindex= allocate2d_short(numgran,2);

        tfile_ini_sec=l1cfile->tfile_ini_sec;
        tfile_end_sec=l1cfile->tfile_end_sec;

        tg_ini=tfile_ini_sec;//always unix time
        tg_end=tg_ini+gransize*60;//always in seconds

        //first granule----
        if(l1cinput->start_time[0]!='\0' && l1cinput->end_time[0]!='\0' && l1cinput->grantype==0){
         cout<<"Processing L1C granules between ........................"<<l1cinput->start_time<<"and ........................."<<l1cinput->end_time<<endl;  
         tgran_ini = isodate2unix(l1cinput->start_time);
         tgran_end = isodate2unix(l1cinput->end_time);
         unix2ymds(tgran_ini, &syear, &smon, &sday,&secs);
         unix2ymds(tgran_end, &syear2, &smon2, &sday2,&secs2);
         tgran_ini_sec=ymds2unix(syear,smon,sday,secs);
         tg_ini=tgran_ini_sec;
         tg_end=tg_ini+gransize*60;
         tgran_end_sec=ymds2unix(syear2,smon2,sday2,secs2); 
         tfile_end_sec=tgran_end_sec+gransize*60;

         unix2ymdhms(tgran_ini_sec,&y_ini,&mo_ini,&d_ini, &h_ini, &mi_ini, &sec_ini);
         cout<<"tgran_ini_sec.."<<"YEAR.."<<y_ini<<"MONTH.."<<mo_ini<<"DAY.."<<d_ini<<"HOUR.."<<h_ini<<"MIN.."<<mi_ini<<"SEC.."<<sec_ini<<endl;
         unix2ymdhms(tgran_end_sec,&y_ini,&mo_ini,&d_ini, &h_ini, &mi_ini, &sec_ini);

         if(tgran_end_sec>tfile_end_sec){ cout<<"ERROR--tgran_end_sec>tfile_end_sec...wrong command line (gran_end_time).gran_end_time beyond swath time limits.."<<endl; exit(1);}

         cout<<"tgran_end_sec.."<<"YEAR.."<<y_ini<<"MONTH.."<<mo_ini<<"DAY.."<<d_ini<<"HOUR.."<<h_ini<<"MIN.."<<mi_ini<<"SEC.."<<sec_ini<<endl;
         unix2ymdhms(tg_ini,&y_ini,&mo_ini,&d_ini, &h_ini, &mi_ini, &sec_ini);
         cout<<"tg_ini.."<<"YEAR.."<<y_ini<<"MONTH.."<<mo_ini<<"DAY.."<<d_ini<<"HOUR.."<<h_ini<<"MIN.."<<mi_ini<<"SEC.."<<sec_ini<<endl;
         unix2ymdhms(tg_end,&y_ini,&mo_ini,&d_ini, &h_ini, &mi_ini, &sec_ini);
         cout<<"tg_end.."<<"YEAR.."<<y_ini<<"MONTH.."<<mo_ini<<"DAY.."<<d_ini<<"HOUR.."<<h_ini<<"MIN.."<<mi_ini<<"SEC.."<<sec_ini<<endl;
         }
        else{
         cout<<"ERROR Processing L1C granules between initial and final time --"<<"start_time.."<<l1cinput->start_time<<"end_time.."<<l1cinput->end_time<<"grantype..."<<l1cinput->grantype<<endl;        
        }

     
        for(int i=0;i<num_gridlines;i++){   
          granid[i]=-1;  
          }

        for(int i=0;i<numgran;i++){
             for(size_t j=0;j<2;j++){
            gdtime[i][j]=-1;
            gdindex[i][j]=-1;
        }}

   
       int neg=0,gmin=-1,c=0;
       for(int gran=0;gran<numgran;gran++){
        for(int i=0;i<num_gridlines;i++){
            tgridline=tfile_ini_sec+tmgv[i]; //seconds of the day are added to second since unix time reference
   
          if(tg_end>tgran_end_sec) tg_end=tgran_end_sec;  
          if(tgridline>=tg_ini && tgridline<tg_end){             
            if(l1cinput->start_time[0]!='\0' || l1cinput->end_time[0]!='\0')
             {
              if(tg_ini>=tgran_ini_sec && tg_end<=tgran_end_sec+gransize*60)
              {                      
               if(gmin<0){//first index
                  gdtime[c][0]=tg_ini;
                  if(i==num_gridlines-1) gdtime[c][1]=tgridline;else  gdtime[c][1]=tg_end;
                  gdindex[c][0]=i;
                  gmin=1;
               }
                  gdindex[c][1]=i;
                  if(i==num_gridlines-1) gdtime[c][1]=tgridline;else  gdtime[c][1]=tg_end; 
             }//gran time                   
            }//command line            
            else{
   
               if(gmin<0){//first index
                  gdtime[c][0]=tg_ini;
                  if(i==num_gridlines-1) gdtime[c][1]=tgridline;else  gdtime[c][1]=tg_end;
                  gdindex[c][0]=i;
                  gmin=1;
                         }
                  gdindex[c][1]=i;
                  if(i==num_gridlines-1) gdtime[c][1]=tgridline;else  gdtime[c][1]=tg_end;
            }
          }


              if(tmgv[i]<0 && gran==0){//previous day data
                 neg++;
              }
        }//gridlines
         //new granule----
           tg_ini=tg_end;
           tg_end=tg_ini+gransize*60;
           gmin=-1;
           c++;       

           if(tg_ini>tgridline){
              cout<<"gridlines with negative time..."<<neg<<endl;
              break;
           }
  //granule selection constraints----------------
           if(l1cinput->start_time[0]!='\0' || l1cinput->end_time[0]!='\0'){
               if(tg_ini>=tg_end ||  tg_ini>=tgran_end_sec || tg_end>(tgran_end_sec+gransize*60)){                                  
                 break;
               }                               
           }    

       }//granules

  std::string timestr,missionstr,fname_out, fname_out_nopath,pathstr,monstr, daystr, yearstr, secstr,mistr,hstr,prodstr, gdstr, swtstr, swtnum, extstr,ofilestr,dirstr,datetimestr1,datetimestr2,fdatetimestr1,date_created;
  std::string y_create,m_create,d_create,t_create;

  pathstr = "";
  missionstr="PACE";
  extstr = ".nc";
//write filenames to list----
  string outxt(l1cinput->outlist);
  outxt=pathstr+outxt;

  std::ofstream outf;

  outf.open(outxt, std::ofstream::out | std::ofstream::app);

  if (outf)
  {
  cout<<"writing L1C granules to outfile..."<<outxt<<endl; 
  }
  else
    {
        std::cerr <<"output file.."<<outxt<<" could not be opened for writing!\n";
        return 1;
    }


// Current date/time based on current system
  time_t now = time(0);   
  // Convert now to tm struct for UTC
  tm* gmtm = gmtime(&now);
  if (gmtm != NULL) {
//     cout << "The UTC date and time is: " << asctime(gmtm) << endl;
  }
  else {
    cerr << "Failed to get the UTC date and time" << endl;
    return EXIT_FAILURE;
  }

  y_create=std::to_string(1900 + gmtm->tm_year);
  m_create=std::to_string(1 + gmtm->tm_mon);
  d_create=std::to_string(gmtm->tm_mday);
  t_create=std::to_string(gmtm->tm_hour) +":"+std::to_string(gmtm->tm_min)+":"+std::to_string(gmtm->tm_sec);

  date_created=y_create+"-"+m_create+"-"+d_create+"T"+t_create;

   int16_t ngridlines,totlines=0;
//write files with granid>0 ---------------------
    for(int gran=0;gran<numgran;gran++){
       if(gdtime[gran][0]>0){
           cout<<"gran #..."<<gran+1<<endl;         

        tg_ini=gdtime[gran][0];
        tg_end=gdtime[gran][1];

        unix2ymdhms(tg_ini,&y_ini,&mo_ini,&d_ini, &h_ini, &mi_ini, &sec_ini);
        unix2ymdhms(tg_end,&y_end,&mo_end,&d_end, &h_end, &mi_end, &sec_end);

        if(mi_end*60==0)
        {
          gtime=((60*60+round(sec_end))-mi_ini*60-round(sec_ini))/60;
        }
        else
        {
        gtime=((mi_end*60+round(sec_end))-mi_ini*60-round(sec_ini))/60;
        }
    
        if(gtime>gransize)
        {
         cout<<"gtime = "<<gtime<<" is greater than granule size = "<<gransize<<endl;
         exit(1);
        }
        if((gtime % gransize)==0) gfull="1";else gfull="0";
        cout<<"gtime.."<<gtime<<"gfull.."<<gfull<<endl;

        //# gridlines per granule
        //gridline index 0-num_gridlines
        ngridlines=gdindex[gran][1]-gdindex[gran][0]+1;
        totlines+=ngridlines;
        cout<<"ngridlines..."<<ngridlines<<"tot gridlines..."<<totlines<<endl;
        

        if(asc_mode==1) dirstr="Ascending";
        else if(asc_mode==0) dirstr="Descending";

//time coverage start----
        secstr = std::to_string(sec_ini);
        mistr = std::to_string(mi_ini);
        hstr = std::to_string(h_ini);
        daystr = std::to_string(d_ini);
        monstr = std::to_string(mo_ini);    
        yearstr = std::to_string(y_ini);

        int length = (int) floor( log10 (mo_ini) ) + 1;
        if(length==1) monstr="0"+monstr;
        length = (int) floor( log10 (d_ini) ) + 1;
        if(length==1) daystr="0"+daystr;

        if(h_ini==0) logoff=1; else logoff=0;
        length = (int) floor( log10 (h_ini+logoff)) + 1;
        if(length==1) hstr="0"+hstr;
        if(mi_ini==0) logoff=1; else logoff=0;
        length = (int) floor( log10 (mi_ini+logoff)) + 1;
        if(length==1) mistr="0"+mistr;
        if(sec_ini==0) logoff=1; else logoff=0;
        length = (int) floor( log10 (round(sec_ini+logoff))) + 1;
        if(length==1) secstr="0"+secstr;


        fdatetimestr1=yearstr+monstr+daystr+"T"+hstr+mistr+secstr.substr(0,2);
        datetimestr1=yearstr+"-"+monstr+"-"+daystr+"T"+hstr+":"+mistr+":"+secstr.substr(0,2);
//time coevrage end----
//
        secstr = std::to_string(sec_end);
        mistr = std::to_string(mi_end);
        hstr = std::to_string(h_end);
        daystr = std::to_string(d_end);
        monstr = std::to_string(mo_end);
        yearstr = std::to_string(y_end);

        length = (int) floor( log10 (mo_end) ) + 1;
        if(length==1) monstr="0"+monstr;
        length = (int) floor( log10 (d_end) ) + 1;
        if(length==1) daystr="0"+daystr;

         if(h_end==0) logoff=1; else logoff=0;
        length = (int) floor( log10 (h_end+logoff)) + 1;
        if(length==1) hstr="0"+hstr;
        if(mi_end==0) logoff=1; else logoff=0;
        length = (int) floor( log10 (mi_end+logoff)) + 1;
        if(length==1) mistr="0"+mistr;
        if(sec_end==0) logoff=1; else logoff=0;
        length = (int) floor( log10 (round(sec_end+logoff))) + 1;
        if(length==1) secstr="0"+secstr;

        datetimestr2=yearstr+"-"+monstr+"-"+daystr+"T"+hstr+":"+mistr+":"+secstr.substr(0,2);
      
        fname_out=pathstr+"PACE."+fdatetimestr1+".L1C"+extstr;
        fname_out_nopath="PACE."+fdatetimestr1+".L1C"+extstr;
   
        cout<<"granule filename.."<<fname_out<<endl;
             
        //output file with list of granules                   
        outf <<fname_out_nopath<<","<<datetimestr1<<","<<datetimestr2<<","<<gfull<<"\n";
        //---------------------------------------------------------------------------------
  
        l1cfile->gridname = fname_out.c_str();
        filename_lt = fname_out.c_str();
        char *gridchar=strdup(filename_lt);
        string l1c_str=filename_lt;

        NcFile* nc_output;
          try {
              nc_output = new NcFile(filename_lt, NcFile::replace);
              }
          catch (NcException& e) {
              e.what();
              cerr << "l1cgen l1c_pflag= 5 : producing L1C grid: "
               + l1c_str << endl;
              exit(1);
              }

        meta_l1c_grid(gridchar,ngridlines,nc_output);

        string verstr=l1cfile->version;
        verstr="V"+verstr.substr(0,4);
        nc_output->putAtt("processing_version",verstr);
        nc_output->putAtt("history",l1cinput->history);
        nc_output->putAtt("product_name",fname_out_nopath);
        nc_output->putAtt("startDirection",dirstr);
        nc_output->putAtt("endDirection",dirstr);
        nc_output->putAtt("time_coverage_start",datetimestr1);
        nc_output->putAtt("time_coverage_end",datetimestr2);

        NY = ngridlines;
        NX = nbinx;
        //subset lat and lon arrays
        NY1=gdindex[gran][0];
        NY2=gdindex[gran][1];

        lat_out = allocate2d_float(NY, NX);
        lon_out = allocate2d_float(NY, NX);
        alt_out = allocate2d_float(NY, NX);
        time_nad_out=(double*)calloc(NY,sizeof(double));//time of the day in seconds
   
       int cc=0;
        for (int i = NY1; i < NY2 + 1; i++) {
            for (int j = 0; j < NX; j++) {
                lat_out[cc][j] = lat_gd[i][j];
                lon_out[cc][j] = lon_gd[i][j];
                alt_out[cc][j] = alt_gd[i][j];
                time_nad_out[cc]=tmgv[i];  
            }
            cc++;
        }

        //vars
        NcGroup ba_grp=nc_output->getGroup("bin_attributes");
        NcVar v1=ba_grp.getVar("nadir_view_time");
        v1.putVar(&time_nad_out[0]);
        NcGroup geo_grp=nc_output->getGroup("geolocation_data");
        v1=geo_grp.getVar("latitude");
        v1.putVar(&lat_out[0][0]);
        v1=geo_grp.getVar("longitude");
        v1.putVar(&lon_out[0][0]);
        v1=geo_grp.getVar("altitude");
        v1.putVar(&alt_out[0][0]);

   //GRING-------
   //determine the number of GCpoint indexes
   //default 6 for the swath sides + number of coordinates every 20 degrees latitude

//ascending pass
   if(l1cfile->orb_dir==1)
   {
        
   Nwest=round((lat_gd[l1cfile->num_gridlines-1][0]-lat_gd[0][0])/20);
   Neast=round((lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1]-lat_gd[0][l1cfile->nbinx-1])/20);
   }
   else //descending
   {
   Neast=(round(lat_gd[0][0]-lat_gd[l1cfile->num_gridlines-1][0])/20);
   Nwest=round((lat_gd[0][l1cfile->nbinx-1]-lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1])/20);
   }

//first NGring estimate----   
   Ngring=Nwest+Neast+6;


   if(Ngring>0)
   {
      float *latarr=(float*)calloc(Ngring,sizeof(float));
      float *lonarr=(float*)calloc(Ngring,sizeof(float));
      int *narr=(int*)calloc(Ngring,sizeof(int));
      int *p_west=(int*)calloc(Ngring,sizeof(int));
      int *p_east=(int*)calloc(Ngring,sizeof(int));

      //corners--counterclockwise and ascending pass
     if(l1cfile->orb_dir==1){
         if(Ngring==6)
         {
          latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
          latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];
          latarr[3]=lat_gd[0][0];
          latarr[4]=lat_gd[0][midix];
          latarr[5]=lat_gd[0][l1cfile->nbinx-1];

          lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
          lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];
          lonarr[3]=lon_gd[0][0];
          lonarr[4]=lon_gd[0][midix];
          lonarr[5]=lon_gd[0][l1cfile->nbinx-1];
         }
         else
         {
         latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
         latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];

         lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
         lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];
    
         latemp=latarr[2];     
         latemp-=20;         
         p=1;
         rw=0;
         //west side
         while(latemp>lat_gd[0][0])
         {
         latarr[2+p]=latemp;
     //    cout<<"p west--"<<p<<"point# = "<<3+p<<"lat20 = "<<latemp<<endl;
         p_west[rw]=3+p;
         latemp-=20;
         p++;
         rw++;
         }   

         p--;
       
         latarr[2+p+1]=lat_gd[0][0];
         latarr[2+p+2]=lat_gd[0][midix];
         latarr[2+p+3]=lat_gd[0][l1cfile->nbinx-1];

         lonarr[2+p+1]=lon_gd[0][0];
         lonarr[2+p+2]=lon_gd[0][midix];
         lonarr[2+p+3]=lon_gd[0][l1cfile->nbinx-1];



         latemp=latarr[2+p+3];
         latemp+=20;
         p++;
       
         int c=1;
         re=0;
         //east side
         while(latemp<lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1])
         {
         latarr[5+p]=latemp;
      //   cout<<"p east--"<<c<<"point# = "<<6+p<<"lat20 = "<<latemp<<endl;       
         p_east[re]=6+p;
         latemp+=20;
         p++;
         c++;
         re++;
         }

         p--;
      
      //   cout<<"mid points west.."<<rw<<"mid points east.."<<re<<endl;
     //    cout<<"# points in GRING = "<<6+p<<"# mid points west--"<<rw<<"# mid points east--"<<re<<endl;
//west
         for(int i=0;i<rw;i++)
         {
             ix=p_west[i]-1;
   //          cout<<"lat:"<<latarr[ix]<<"pwest--"<<p_west[i]<<endl;
          for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {
              if(latarr[ix]>lat_gd[row][0] && latarr[ix]<=lat_gd[row+1][0])
              {
        //       cout<<"mid point west# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;            

               if(lon_gd[row][0]<0.) lontemp1=lon_gd[row][0]+360;else lontemp1=lon_gd[row][0];
               if(lon_gd[row+1][0]<0.) lontemp2=lon_gd[row+1][0]+360;else lontemp2=lon_gd[row+1][0];

               dlat_gd=abs(lat_gd[row+1][0]-lat_gd[row][0]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row+1][0]);
               dlon20=dlat20*dlon_gd/dlat_gd;
   //            if(lontemp1>lontemp2) lon360=lontemp1-dlon20;else lon360=lontemp1+dlon20;
               if(lontemp1>lontemp2) lon360=lontemp2+dlon20;else lon360=lontemp2-dlon20;
               if(lon360>180) lon360=lon360-360.; 
               lonarr[ix]=lon360;
       //        cout<<"lon_gd row+1.."<<lon_gd[row+1][0]<<"lon_gd row.."<<lon_gd[row][0]<<"lon cring.."<<lonarr[ix]<<endl;
                break;
              }
          }
         }
//east
         for(int i=0;i<re;i++)
         {
             ix=p_east[i]-1;
        //     cout<<"lat:"<<latarr[ix]<<"peast--"<<p_east[i]<<endl;
             for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {         
              if(latarr[ix]>lat_gd[row][l1cfile->nbinx-1] && latarr[ix]<=lat_gd[row+1][l1cfile->nbinx-1])
              {
     //          cout<<"mid point east# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;

               if(lon_gd[row][l1cfile->nbinx-1]<0.) lontemp1=lon_gd[row][l1cfile->nbinx-1]+360;else lontemp1=lon_gd[row][l1cfile->nbinx-1];
               if(lon_gd[row+1][l1cfile->nbinx-1]<0.) lontemp2=lon_gd[row+1][l1cfile->nbinx-1]+360;else lontemp2=lon_gd[row+1][l1cfile->nbinx-1];

               dlat_gd=abs(lat_gd[row+1][l1cfile->nbinx-1]-lat_gd[row][l1cfile->nbinx-1]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row][l1cfile->nbinx-1]);
               dlon20=dlat20*dlon_gd/dlat_gd;
               if(lontemp1>lontemp2) lon360=lontemp1-dlon20;else lon360=lontemp1+dlon20;
               if(lon360>180) lon360=lon360-360.;
               lonarr[ix]=lon360;
       //        cout<<"lon_gd row+1.."<<lon_gd[row+1][l1cfile->nbinx-1]<<"lon_gd row.."<<lon_gd[row][l1cfile->nbinx-1]<<"lon cring.."<<lonarr[ix]<<endl;

                 break;
              }
          }           
         }

         }
     }
     else //descending orb
     {
     if(Ngring==6)
         {
          latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
          latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];
          latarr[3]=lat_gd[0][0];
          latarr[4]=lat_gd[0][midix];
          latarr[5]=lat_gd[0][l1cfile->nbinx-1];

          lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
          lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
          lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];
          lonarr[3]=lon_gd[0][0];
          lonarr[4]=lon_gd[0][midix];
          lonarr[5]=lon_gd[0][l1cfile->nbinx-1];
         }
         else
         {
         latarr[0]=lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         latarr[1]=lat_gd[l1cfile->num_gridlines-1][midix];
         latarr[2]=lat_gd[l1cfile->num_gridlines-1][0];

         lonarr[0]=lon_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1];
         lonarr[1]=lon_gd[l1cfile->num_gridlines-1][midix];
         lonarr[2]=lon_gd[l1cfile->num_gridlines-1][0];

         latemp=latarr[2];     
         latemp+=20;         
         p=1;
         int rw=0;
         //west side
         while(latemp<lat_gd[0][0])
         {
         latarr[2+p]=latemp;
       //  cout<<"p west--"<<p<<"point# = "<<3+p<<"lat20 = "<<latemp<<endl;
         p_west[rw]=3+p;
         latemp+=20;
         p++;
         rw++;
         }   

         p--;
       
         latarr[2+p+1]=lat_gd[0][0];
         latarr[2+p+2]=lat_gd[0][midix];
         latarr[2+p+3]=lat_gd[0][l1cfile->nbinx-1];

         lonarr[2+p+1]=lon_gd[0][0];
         lonarr[2+p+2]=lon_gd[0][midix];
         lonarr[2+p+3]=lon_gd[0][l1cfile->nbinx-1];
         latemp=latarr[2+p+3];
         latemp-=20;
         p++;
       
         int c=1;
         int re=0;
         //east side
         while(latemp>lat_gd[l1cfile->num_gridlines-1][l1cfile->nbinx-1])
         {
         latarr[5+p]=latemp;
         //cout<<"p east--"<<c<<"point# = "<<6+p<<"lat20 = "<<latemp<<endl;       
         p_east[re]=6+p;
         latemp-=20;
         p++;
         c++;
         re++;
         }

         p--;
      
         //cout<<"mid points west.."<<rw<<"mid points east.."<<re<<endl;
        // cout<<"# points in GRING = "<<6+p<<"# mid points west--"<<rw<<"# mid points east--"<<re<<endl;

//west side
         for(int i=0;i<rw;i++)
         {
             ix=p_west[i]-1;
       //      cout<<"lat:"<<latarr[ix]<<"pwest--"<<p_west[i]<<endl;
          for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {
              if(latarr[ix]<=lat_gd[row][0] && latarr[ix]>lat_gd[row+1][0])
              {
 
               if(lon_gd[row][0]<0.) lontemp1=lon_gd[row][0]+360;else lontemp1=lon_gd[row][0];
               if(lon_gd[row+1][0]<0.) lontemp2=lon_gd[row+1][0]+360;else lontemp2=lon_gd[row+1][0];

               dlat_gd=abs(lat_gd[row][0]-lat_gd[row+1][0]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row+1][0]);
               dlon20=dlat20*dlon_gd/dlat_gd;

               if(lontemp1>lontemp2) lon360=lontemp2+dlon20;else lon360=lontemp2-dlon20;
               if(lon360>180) lon360=lon360-360.; 
               lonarr[ix]=lon360;
    
                break;
              }
          }
         }
//east side
         for(int i=0;i<re;i++)
         {
             ix=p_east[i]-1;
       //      cout<<"lat:"<<latarr[ix]<<"peast--"<<p_east[i]<<endl;
             for(int row=0;row<l1cfile->num_gridlines-1;row++)
          {         
              if(latarr[ix]<=lat_gd[row][l1cfile->nbinx-1] && latarr[ix]>lat_gd[row+1][l1cfile->nbinx-1])
              {
   
               if(lon_gd[row][l1cfile->nbinx-1]<0.) lontemp1=lon_gd[row][l1cfile->nbinx-1]+360;else lontemp1=lon_gd[row][l1cfile->nbinx-1];
               if(lon_gd[row+1][l1cfile->nbinx-1]<0.) lontemp2=lon_gd[row+1][l1cfile->nbinx-1]+360;else lontemp2=lon_gd[row+1][l1cfile->nbinx-1];

               dlat_gd=abs(lat_gd[row][l1cfile->nbinx-1]-lat_gd[row+1][l1cfile->nbinx-1]);
               dlon_gd=abs(lontemp1-lontemp2);
               dlat20=abs(latarr[ix]-lat_gd[row][l1cfile->nbinx-1]);
               dlon20=dlat20*dlon_gd/dlat_gd;
               if(lontemp1>lontemp2) lon360=lontemp1-dlon20;else lon360=lontemp1+dlon20;
               if(lon360>180) lon360=lon360-360.;
               lonarr[ix]=lon360;
  

                 break;
              }
          }           
         }

         }//else GRING>6

     }//end descending
    
     //sequence- GRING
     if(Ngring==6) dp=6;else dp=rw+re+6;
      for(int s=0;s<dp;s++){        
         narr[s]=s+1; 
      }  

      geo_grp.putAtt("GRingPointLatitude",ncFloat,dp,latarr);
      geo_grp.putAtt("GRingPointLongitude",ncFloat,dp,lonarr);
      geo_grp.putAtt("GRingPointSequenceNo",ncInt,dp,narr);

      delete [] (latarr);
      delete [] (lonarr);
      delete [] (narr);
      delete [] (p_west);
      delete [] (p_east);
   }
   else
   {
   cout<<"ERROR EXTRACTING GRING coordinates!!-----"<<endl;
   exit(1);
   }


        

        nc_output->close();

        if(lat_out!=nullptr)
         delete [] (lat_out); 
        if(lon_out!=nullptr)
         delete [] (lon_out);
        if(alt_out!=nullptr)
         delete [] (alt_out);
        if(time_nad_out!=nullptr)
         delete [] (time_nad_out);
       }//gdtime >0
      }//end granule loop

     if(gdtime!=nullptr)
         delete[] (gdtime);
     gdtime=nullptr;

      if(gdindex!=nullptr)
         delete[] (gdindex);
     gdindex=nullptr;

     if(granid!=nullptr)
         delete[] (granid);
     granid=nullptr;

outf.close();



return 0;
}    



//same as l1c but files are saved as nc
 int32_t L1C::across_gridlines_l1c2(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput, int16_t* swtd_id, int16_t* file_id, int16_t* nfiles_swt, float* lati2, float* loni2, float** lat_gd, float** lon_gd, float* az_east) {
        int32_t num_gridlines, nbinx, NY1 = -1, NY2 = -1, c1 = 1;
        float lat1_l, lat1_r, lon1_l, lon1_r,az_r, az_l, h,  *lati_l, *loni_l, *lati_r, *loni_r, htot = 0.;
        int  ncid_out,x_dimid, y_dimid, varid1, varid2, status;       
        int NDIMS=2;
        int dimids[NDIMS];
        int32_t NX, NY;
        int16_t selyear = -1, selmon = -1, selday = -1;
        double  latih, lonih, latih2, lonih2, res1, res2;
        float gres;
        int fnan1=0, fnan2=0, fnan3=0, fnan4=0;
        float latfirst = -999;
        int badline = 0;
        int grp_coor;
        const char* filename_lt;
        float **lat_out=nullptr,**lon_out=nullptr;
//        float **lt_out=nullptr;
        int NVIEWS, NBANDS,v_dimid,b_dimid, asc_mode=-1;
        vector <int32_t> cutoff_r,cutoff_l;
        int32_t ix=0,ix2=0;


        Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());

        //compute azimuth angle for each along-track position (nadir pixel or centered swath pixel) of the swath
        //azimuth angle in rads
        //lon+540/360-180 is a normalization factor for -180 to 180 degrees  

        int nadpix = l1cfile->nadpix;
        cout<<"nadpix..."<<nadpix<<endl;

        num_gridlines = l1cfile->num_gridlines;
        num_pixels = l1cfile->npix;
        nbinx = l1cfile->nbinx;
        lati_l = (float*)calloc(nbinx / 2 , sizeof(float));
        lati_r = (float*)calloc(nbinx / 2 , sizeof(float));
        loni_l = (float*)calloc(nbinx / 2 , sizeof(float));
        loni_r = (float*)calloc(nbinx / 2 , sizeof(float));

      //big loop

        int gd = 0;
        asc_mode=l1cfile->orb_dir;


        //compute mean az_east---
  /*      
        for (int i = 0;i < num_gridlines - 2;i++) {

        //    cout<<"gd#.."<<i+1<<"lati2[i].."<<lati2[i]<<endl;

            //assign central positions---
            lat1_l = lati2[i] * M_PI / 180.;
            lat1_r = lati2[i] * M_PI / 180.;
            lat2 = lati2[i + 1] * M_PI / 180.;
            lon1_l = loni2[i] * M_PI / 180.;
            lon1_r = loni2[i] * M_PI / 180.;
            //lon2 = loni2[i + 1] * M_PI / 180.;

 //           dlambda = (lon2 - lon1_l);
            //bearing---
            //make sure there are not NAN values ----------
 //           az = atan2(sin(dlambda) * cos(lat2), cos(lat1_l) * sin(lat2) - sin(lat1_l) * cos(lat2) * cos(dlambda));

            //     cout<<lat2<<"..."<<lat1_l<<".."<<dlambda<<az*180/M_PI<<endl;
         
            if (az > M_PI || az < -M_PI) {
                cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "az in degrees.." << az * 180 / M_PI << endl;
                exit(1);
            }
            
            if (isnan(az) < 1) {
                //       cout<<"az.."<<az<<"azcount"<<azc<<endl;
                az_r = (az + M_PI / 2.);
                az_east[i] = az_r * 180. / M_PI;
                if (az_east[i] < 0.0001 && az_east[i]>0) az_east[i] = 0.0;
                sum = az_east[i];
                tot += sum;
                azc++;
            }
            else {
                az_east[i] = NAN;
            }
        }
        */
//        float mean_az_east = tot / (azc);

        float mean_az_east=l1cfile->mean_az_east;

        cout << "mean_az_east.." << mean_az_east << "numgridlines.." << num_gridlines <<"nbinx/2..."<<nbinx/2<<endl;
//        exit(1);



        for (int j = 0;j < num_gridlines;j++) {
            az_east[j] = mean_az_east;
        }

        l1cfile->mean_az_east = mean_az_east;

        for (int i = 0;i < num_gridlines;i++) {

            ix=0;
            ix2=0;

            if (i % 1 == 0) {
                //assign central positions---
                lat1_l = lati2[i] * M_PI / 180.;
                lat1_r = lati2[i] * M_PI / 180.;
   //             lat2 = lati2[i + 1] * M_PI / 180.;
                lon1_l = loni2[i] * M_PI / 180.;
                lon1_r = loni2[i] * M_PI / 180.;
     //           lon2 = loni2[i + 1] * M_PI / 180.;

    //            dlambda = (lon2 - lon1_l);


                gres = (l1cinput->grid_resolution) * 1000; //grid resolution in meters

               //bearing---
               //this is Harverside--law f cosines is more accurate---
    //            az = atan2(sin(dlambda) * cos(lat2), cos(lat1_l) * sin(lat2) - sin(lat1_l) * cos(lat2) * cos(dlambda));
     /*           if (az > M_PI || az < -M_PI) {
                    cout << "problem with BEARING in across-gridline method...az<-180 or >180...." << "az in degrees.." << az * 180 / M_PI << endl;
                    exit(1);
                }
*/


                for (int k = 0;k < (nbinx-1) / 2;k++) {
                    h = (l1cinput->grid_resolution);//distance between along positions in km, must be 5.2 k m for L1C
                    if (k == 0) h /= 2;//if nadpix
                    htot += h;

                    //right pixels---
                    az_r = mean_az_east * degrad;

                    //left pixels---

                    az_l = az_r - M_PI;
                    //     az_l=(az-M_PI/2.);//this should be equal to az_r+M_PI=az_l

                   //bearing should be always positive (clockwise) before pixel location calculation based on Harverside---     
                    if (az_r < 0.) {
                        az_r = 360 + az_r * 180 / M_PI;
                        az_r = az_r * M_PI / 180;
                    }

                    if (az_l < 0.) {
                        az_l = 360 + az_l * 180 / M_PI;
                        az_l = az_l * M_PI / 180;
                    }
          
                    geod.Direct(lat1_r * 180 / M_PI, lon1_r * 180 / M_PI, az_r * 180 / M_PI, gres, latih, lonih);
                    geod.Inverse(lat1_r * 180 / M_PI, lon1_r * 180 / M_PI, latih, lonih, res1);
                    fnan1 = isnan(latih);
                    fnan2 = isnan(lonih);                    
//------------------------------------------------------------------
                    if(latih<lat1_r* 180 / M_PI && asc_mode==1){//ascending
         //              cout<<"latih..."<<latih<<"<lat1_r...."<<lat1_r<<"at scanline#.."<<i+1<<"at pixel# from nadir pos.."<<k+1<<endl;                                 
                       ix++;
                       if(ix==1)  cutoff_r.push_back(i);
                     }
         
                         
                    if(latih>lat1_r* 180 / M_PI && asc_mode==0){//descending
       //                cout<<"latih..."<<latih<<">lat1_r...."<<lat1_r<<"at scanline#.."<<i+1<<"at pixel# from nadir pos.."<<k+1<<endl;
                       ix++;
                       if(ix==1)  cutoff_r.push_back(i);
                       }
                      
//-----------------------------------------------------------------                    

                    lati_r[k] = latih;
                    loni_r[k] = lonih;
                

                    lat1_r = lati_r[k] * M_PI / 180;
                    lon1_r = loni_r[k] * M_PI / 180;

                    geod.Direct(lat1_l * 180 / M_PI, lon1_l * 180 / M_PI, az_l * 180 / M_PI, gres, latih2, lonih2);
                    geod.Inverse(lat1_l * 180 / M_PI, lon1_l * 180 / M_PI, latih2, lonih2, res2);
                    fnan3 = isnan(latih2);
                    fnan4 = isnan(lonih2);

//----------------------------------------------------------------
                    if(latih2>lat1_l* 180 / M_PI && asc_mode==1){
           //            cout<<"latih2..."<<latih2<<">lat1_l...."<<lat1_l<<"at scanline#.."<<i+1<<"at pixel.from nadir pos."<<k+1<<endl;                                                    
                       ix2++;            
                       if(ix2==1)  cutoff_l.push_back(i);
                       }
                    if(latih2<lat1_l* 180 / M_PI && asc_mode==0){
             //          cout<<"latih2..."<<latih2<<"<lat1_l...."<<lat1_l<<"at scanline#.."<<i+1<<"at pixel.from nadir pos."<<k+1<<endl;                                                              
                       ix2++;
                       if(ix2==1)  cutoff_l.push_back(i);
                          }
           /*            
                    if(asc_mode==1 && lon1_l* 180 / M_PI<0 && lonih2<0)
                        cout<<"lon - sign left..at line.."<<i<<"binx.."<<k+1<<endl;
                     if(asc_mode==1 && lon1_l* 180 / M_PI>0 && lonih2>0)
                        cout<<"lon + sign left..at line.."<<i<<"binx.."<<k+1<<endl;
             */           
//----------------------------------------------------------------

                    lati_l[k] = latih2;
                    loni_l[k] = lonih2;
                    //     cout<<"leftside.."<<latih2<<"..."<<lonih2<<"res2"<<res2<<endl;
                    //
                    if(latih>90 || latih<-90 || latih2>90 || latih2<-90 || lati2[i]>90 || lati2[i]<-90 ){
                        cout<<"WRONG latitude beyon -90 or 90 latitude degrees............."<<"lati2[i]....."<<lati2[i]<<"latih.."<<latih<<"latih2.."<<latih2<<endl;
                        exit(1);
                    }



                    lat1_l = lati_l[k] * M_PI / 180;
                    lon1_l = loni_l[k] * M_PI / 180;

                      //right swath side-------------
                    if (asc_mode==1 && lati2[i]>lati2[i+1] || asc_mode==0 && lati2[i]<lati2[i+1] || fnan1 > 0 || fnan2 > 0 || res1 > (gres + 6) || res1 < (gres - 6) || lati2[i + 1] < lati2[i] || lati2[i] < latfirst) { //nan values
                        badline++;
                        if (badline == 1) {
                            latfirst = lati2[i];
                        }

                 //       cout << "bad interpolation - missing orbital data at near nadir grid cell..NAN right side of line.." << i + 1 << "corresponding to lat..."<<latfirst <<"lati2[i].."<<lati2[i]<<endl;
                 //       cout<<"gd.."<<gd+1<<"XBIN.."<<k+1<<"lati_l[k]"<<lati_l[k]<<"loni_l[k]"<<loni_l[k]<<"lati_r[k]"<<lati_r[k]<<"loni_r[k]"<<loni_r[k]<<endl;
                        lati2[i] = NAN;
                        loni2[i] = NAN;

                        for (int k = 0;k < (nbinx-1) / 2;k++) {
                            lat_gd[gd][k + (nbinx-1)/ 2] = NAN;
                            lon_gd[gd][k + (nbinx-1) / 2] = NAN;
                        }
                        //       cout<<"nan value.right.."<<endl;
                    }
                    else {
                        lat_gd[gd][k + (nbinx-1) / 2] = latih;
                        lon_gd[gd][k + (nbinx-1) / 2] = lonih;                    
                }

                    //left side
                    if (asc_mode==1 && lati2[i]>lati2[i+1] || asc_mode==0 && lati2[i]<lati2[i+1] || fnan3 > 0 || fnan4 > 0 || res2 > (gres + 6) || res2 < (gres - 6) || lati2[i + 1] < lati2[i] || lati2[i] < latfirst) {
                        //     if(fnan3>0 | fnan4>0 |  res2>(gres+6) | res2<(gres-6)){
                  //      cout << "bad interpolation - missing orbital data at near nadir grid cell..NAN left side of line.." << i + 1 << endl;
                  //      cout<<"gd.."<<gd+1<<"XBIN.."<<k+1<<"lati_l[k]"<<lati_l[k]<<"loni_l[k]"<<loni_l[k]<<"lati_r[k]"<<lati_r[k]<<"loni_r[k]"<<loni_r[k]<<"lati2[i].."<<lati2[i]<<endl;
                        lati2[i] = NAN;
                        loni2[i] = NAN;
                        for (int k = 0;k < (nbinx-1) / 2;k++) {
                            lat_gd[gd][(nbinx-1) / 2 - 1 - k] = NAN;//as k increases we are farther away from nadpix or along-track position---
                            lon_gd[gd][(nbinx-1) / 2 - 1 - k] = NAN;                        
                         }   
                        //       cout<<"nan value.left.."<<endl;
                    }
                    else {
                        lat_gd[gd][(nbinx-1) / 2 - 1 - k] = latih2;//as k increases we are farther away from nadpix or along-track position---
                        lon_gd[gd][(nbinx-1) / 2 - 1 - k] = lonih2;                    
                       }  
                fnan1=0;
                fnan2=0;
                fnan3=0;
                fnan4=0;

                }//end pixel loop
              

                htot = 0;
                gd++;
            }//end if every 100 gd
        }//end gridlines loop
        cout << "`L1C gridding successful!......................................." << endl;

        cout<<"cutoff_r.size().."<<cutoff_r.size()<<"cutoff_l.size().."<<cutoff_l.size()<<endl;
        cout<<"first line index..."<<cutoff_l.size()-1<<"last_line.index."<<cutoff_r[0]<<endl;


//writing results ***************************************************************************
        float lontemp = 0.;
        selday = l1cinput->selday;
        selmon = l1cinput->selmon;
        selyear = l1cinput->selyear;

        std::string timestr,missionstr,fname_out, pathstr, senstr, monstr, daystr, yearstr, prodstr, gdstr, swtstr, swtnum, extstr,ofilestr;
        pathstr = "out/";
        missionstr="PACE";
        senstr = "OCIS";
        monstr = std::to_string(selmon);
        daystr = std::to_string(selday);
        yearstr = std::to_string(selyear);
        timestr="T00:00:00Z";
        prodstr = "L1Cgrid";
        //   gdstr="_gdFULL_";
        swtstr = "_swt";

        if (l1cfile->l1c_pflag >= 3) {
            swtd = l1cfile->swtnum;
            swtnum = std::to_string(swtd);
            cout << "swtnum.." << swtnum << endl;
        }
        else swtnum = std::to_string(swtd);

        extstr = ".nc";
        ofilestr=std::string(l1cinput->ofile);
            
        fname_out = pathstr + ofilestr+"FB_"+ prodstr+extstr;

        string  ATT_NAME="Units", ATT_VAL="degrees",GATT_NAME1="title",GATT_VAL1="PACE OCI Level-1C Data",GATT_NAME2="instrument",GATT_VAL2="OCI",GATT_NAME3="processing_version",GATT_VAL3="V1.0",GATT_NAME4="Conventions",GATT_VAL4="CF-1.6";
        string GATT_NAME5="institution",GATT_VAL5="NASA Goddard Space Flight Center, Ocean Biology Processing Group",GATT_NAME6="license",GATT_VAL6="http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
        string GATT_NAME7="naming_authority",GATT_VAL7="gov.nasa.gsfc.sci.oceancolor",GATT_NAME8="keywords_vocabulary",GATT_VAL8="NASA Global Change Master Directory (GCMD) Science Keywords";
        string GATT_NAME9="stdname_vocabulary",GATT_VAL9="NetCDF Climate and Forecast (CF) Metadata Convention",GATT_NAME10="creator_name",GATT_VAL10="NASA/GSFC",GATT_NAME11="creator_email",GATT_VAL11="data@oceancolor.gsfc.nasa.gov";
        string GATT_NAME12="creator_url",GATT_VAL12="http://oceancolor.gsfc.nasa.gov",GATT_NAME13="project",GATT_VAL13="PACE Project",GATT_NAME14="publisher_name",GATT_VAL14="NASA/GSFC";
        string GATT_NAME15="publisher_email",GATT_VAL15="data@oceancolor.gsfc.nasa.gov",GATT_NAME16="publisher_url",GATT_VAL16="http://oceancolor.gsfc.nasa.gov",GATT_NAME17="processing_level",GATT_VAL17="L1C";
        string GATT_NAME18="cdm_data_type",GATT_VAL18="swath",GATT_NAME19="orbit_number",GATT_VAL19="12345",GATT_NAME20="history",GATT_VAL20="",GATT_NAME21="CDL_version_date",GATT_VAL21="2021-09-10",GATT_NAME22="product_name",GATT_VAL22=fname_out;
        string GATT_NAME23="startDirection",GATT_VAL23="Ascending",GATT_NAME24="endDirection",GATT_VAL24="Ascending",GATT_NAME25="time_coverage_start",GATT_VAL25=yearstr+"-"+monstr+"-"+daystr+"-"+timestr,GATT_NAME26="time_coverage_end",GATT_VAL26=yearstr+"-"+monstr+"-"+daystr+"-"+timestr,GATT_NAME27="date{_created",GATT_VAL27="2021-09-10T15:12:41Z",GATT_NAME28="sun_earth_distance",GATT_VAL28="0.990849042172323",GATT_NAME29="terrain_data_source",GATT_VAL29="",GATT_NAME30="spectral_response_function",GATT_VAL30="",GATT_NAME31="systematic_uncertainty_model",GATT_VAL31="",GATT_NAME32="nadir_bin",GATT_VAL32="12345",GATT_NAME33="bin_size_at_nadir",GATT_VAL33="2.6km2"; 


        l1cfile->gridname = fname_out.c_str();
        filename_lt = fname_out.c_str();
        cout<<"creating file for L1C coor..."<<filename_lt<<endl;
//        if ((status = nc_create(filename_lt, NC_CLOBBER, &ncid_out))) 
        if ((status = nc_create(filename_lt, NC_CLOBBER | NC_NETCDF4, &ncid_out)))
            check_err(status, __LINE__, __FILE__);
        //define dims
     // Define the dimensions,vars and attributes at the root level
        NY = num_gridlines;
        NX = nbinx-1;

//        NY1=cutoff_l.size()-1;
//        NY2=cutoff_r[0];

        cutoff_r.clear();
        cutoff_l.clear();


//        NY1=58-1;
//        NY2=3871-1;
 

        for (int i = 0; i < NY; i++) {
                if (lat_gd[i][nadpix] >=-82){                
                    NY1 = i;}                            
            if (NY1 >= 0) break;           
        }

        for (int i = 0;i < NY; i++) {
                if (lat_gd[i][nadpix] >=82) 
                    NY2 = i - 1;                            
            if (NY2 >= 0) break;
        }

        l1cfile->NY1 = NY1;
        l1cfile->NY2 = NY2;
        NY = NY2 - NY1 + 1;

        cout<<"NY1.."<<NY1<<"NY2.."<<NY2<<"NX.."<<NX<<endl;
   //     exit(1);
        

        lat_out = allocate2d_float(NY, NX);
        lon_out = allocate2d_float(NY, NX);

        int c=0;
        for (int i = NY1; i < NY2+1; i++) {
            for (int j = 0; j < NX; j++) {
                lon_out[c][j]=lon_gd[i][j];
                lat_out[c][j] = lat_gd[i][j];
            }
            cout<<"gd#..."<<c+1<<"latgd..at nadir"<<lat_gd[i][NX-1]<<"longd.."<<lon_gd[i][NX-1]<<endl;
            c++;
          }

        NVIEWS= 2;//for OCI
        NBANDS=249;
//        lt_out = allocate2d_float(NY, NX);
  //DEF DIMENSIONS
        if ((status = nc_def_dim(ncid_out, "bins_across_track", NX, &x_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(ncid_out, "bins_along_track", NY, &y_dimid)))
            check_err(status, __LINE__, __FILE__);
        //dims for output var
        dimids[0] = y_dimid;
        dimids[1] = x_dimid;
        NDIMS=2;
        if ((status = nc_def_dim(ncid_out, "number_of_views", NVIEWS, &v_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(ncid_out, "intensity_bands_per_view", NBANDS, &b_dimid)))
            check_err(status, __LINE__, __FILE__);

//define attributes
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME1.c_str(), strlen(GATT_VAL1.c_str()),
          GATT_VAL1.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME2.c_str(), strlen(GATT_VAL2.c_str()),
          GATT_VAL2.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME3.c_str(), strlen(GATT_VAL3.c_str()),
          GATT_VAL3.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME4.c_str(), strlen(GATT_VAL4.c_str()),
          GATT_VAL4.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME5.c_str(), strlen(GATT_VAL5.c_str()),
          GATT_VAL5.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME6.c_str(), strlen(GATT_VAL6.c_str()),
          GATT_VAL6.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME7.c_str(), strlen(GATT_VAL7.c_str()),
          GATT_VAL7.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME8.c_str(), strlen(GATT_VAL8.c_str()),
          GATT_VAL8.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME9.c_str(), strlen(GATT_VAL9.c_str()),
          GATT_VAL9.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME10.c_str(), strlen(GATT_VAL10.c_str()),
          GATT_VAL10.c_str()))
           check_err(status, __LINE__, __FILE__);
        if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME11.c_str(), strlen(GATT_VAL11.c_str()),
          GATT_VAL11.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME12.c_str(), strlen(GATT_VAL12.c_str()),
          GATT_VAL12.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME13.c_str(), strlen(GATT_VAL13.c_str()),
          GATT_VAL13.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME14.c_str(), strlen(GATT_VAL14.c_str()),
          GATT_VAL14.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME15.c_str(), strlen(GATT_VAL15.c_str()),
          GATT_VAL15.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME16.c_str(), strlen(GATT_VAL16.c_str()),
          GATT_VAL16.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME17.c_str(), strlen(GATT_VAL17.c_str()),
          GATT_VAL17.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME18.c_str(), strlen(GATT_VAL18.c_str()),
          GATT_VAL18.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME19.c_str(), strlen(GATT_VAL19.c_str()),
          GATT_VAL19.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME20.c_str(), strlen(GATT_VAL20.c_str()),
          GATT_VAL20.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME21.c_str(), strlen(GATT_VAL21.c_str()),
          GATT_VAL21.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME22.c_str(), strlen(GATT_VAL22.c_str()),
          GATT_VAL22.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME23.c_str(), strlen(GATT_VAL23.c_str()),
          GATT_VAL23.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME24.c_str(), strlen(GATT_VAL24.c_str()),
          GATT_VAL24.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME25.c_str(), strlen(GATT_VAL25.c_str()),
          GATT_VAL25.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME26.c_str(), strlen(GATT_VAL26.c_str()),
          GATT_VAL26.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME27.c_str(), strlen(GATT_VAL27.c_str()),
          GATT_VAL27.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME28.c_str(), strlen(GATT_VAL28.c_str()),
          GATT_VAL28.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME29.c_str(), strlen(GATT_VAL29.c_str()),
          GATT_VAL29.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME30.c_str(), strlen(GATT_VAL30.c_str()),
          GATT_VAL30.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME31.c_str(), strlen(GATT_VAL31.c_str()),
          GATT_VAL31.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME32.c_str(), strlen(GATT_VAL32.c_str()),
          GATT_VAL32.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME33.c_str(), strlen(GATT_VAL33.c_str()),
          GATT_VAL33.c_str()))
          check_err(status, __LINE__, __FILE__);

      //define groups                  check_err(status, __LINE__, __FILE__);
       if ((status = nc_def_grp(ncid_out, "geolocation_data", &grp_coor)))  //netcdf-4
          check_err(status, __LINE__, __FILE__);
//        if ((status = nc_def_grp(ncid_out, "observation_data", &grp_obs)))  //netcdf-4
//        check_err(status, __LINE__, __FILE__);

//DEF DIMENSIONS
  /*      if ((status = nc_def_dim(grp_coor, "x", NX, &x_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(grp_coor, "y", NY, &y_dimid)))
            check_err(status, __LINE__, __FILE__);
        dimids[0] = y_dimid;
        dimids[1] = x_dimid;

        NDIMS=2;
*/
        //def var      
        if ((status = nc_def_var(grp_coor, "latitude", NC_FLOAT, NDIMS,
                dimids, &varid1)))
                check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_var(grp_coor, "longitude", NC_FLOAT, NDIMS,
                dimids, &varid2)))
                check_err(status, __LINE__, __FILE__);
    //leave define mode-----------------------
        if ((status = nc_enddef(grp_coor))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);
//        if ((status = nc_def_var(grp_obs, "lt", NC_FLOAT, NDIMS,
//                dimids, &varid3)))
//                check_err(status, __LINE__, __FILE__);
                    //leave define mode-----------------------
//        if ((status = nc_enddef(grp_obs))) //done def vars etc
//                        check_err(status, __LINE__, __FILE__);

        //writing the whole thing
        if ((status = nc_put_var_float(grp_coor, varid1, &lat_out[0][0])))
            check_err(status, __LINE__, __FILE__);

        if ((status = nc_put_var_float(grp_coor, varid2, &lon_out[0][0])))
            check_err(status, __LINE__, __FILE__);
//        if ((status = nc_put_var_float(grp_obs, varid3, &lt_out[0][0])))
//            check_err(status, __LINE__, __FILE__);
//define attributes of vars in groups
//UNITS in degrees ----
        if (nc_put_att_text(grp_coor, varid1, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
             check_err(status, __LINE__, __FILE__);
        if (nc_put_att_text(grp_coor, varid2, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str())) 
             check_err(status, __LINE__, __FILE__);


//close file
        if ((status = nc_close(ncid_out)))
           check_err(status, __LINE__, __FILE__);

         delete [] (lat_out);
         delete [] (lon_out);
  
 

//CASE ONLY LON POS VALUES***************************************************************************

        //Only positive longitude 0-360 degrees for visualization ****************************************************************
        //lat lon grid for visualization--
        prodstr = "L1Cgrid_lonpos";
        //   gdstr="_gdFULL_";
        //recompute number of gridlines based on lat limits -80 to 80 degrees  
        //asc orbit goes from negative to positive latitude....

        NY = num_gridlines-2;
        NX = nbinx-1;

        for (int i = 0; i < NY; i++) {
            for (int j = 0; j < NX; j++) {
 //               if (lat_gd[i][j] >= -80.) c1++;
                if (lat_gd[i][j] >=-90) c1++;
                if (c1 == nbinx-1) {
                    NY1 = i;
                }
            }
            if (NY1 >= 0) break;
            c1 = 1;
        }

        for (int i = 0;i < NY; i++) {
            for (int j = 0; j < NX; j++) {
                if (lat_gd[i][j] >90.) {
                    NY2 = i - 1;
                }
            }
            if (NY2 >= 0) break;
        }

       lat_out = allocate2d_float(NY, NX);
       lon_out = allocate2d_float(NY, NX);

        c=0;
        for (int i = NY1; i < NY2 + 1; i++) {
            for (int j = 0; j < NX; j++) {
                if(lon_gd[i][j]<0) lon_out[c][j]=lon_gd[i][j]+360;
                else lon_out[c][j]=lon_gd[i][j];
                lat_out[c][j] = lat_gd[i][j];
            }
     //       cout<<"gd#..."<<c+1<<"latgd.."<<lat_gd[i][257]<<"longd.."<<lon_gd[i][257]<<endl;
            c++;
        }


        extstr = ".nc";
//        fname_out = pathstr + missionstr + "_" + senstr + "." + yearstr + monstr + daystr + timestr + prodstr+extstr;
        fname_out = pathstr + ofilestr+"_"+ prodstr+extstr;

        l1cfile->gridname = fname_out.c_str();
        filename_lt = fname_out.c_str();
        cout<<"creating file for L1C coor.with only POS LONGITUDE.."<<filename_lt<<endl;
//        if ((status = nc_create(filename_lt, NC_CLOBBER, &ncid_out))) 
        if ((status = nc_create(filename_lt, NC_CLOBBER | NC_NETCDF4, &ncid_out)))
            check_err(status, __LINE__, __FILE__);
        //define dims
     // Define the dimensions,vars and attributes at the root level


  //DEF DIMENSIONS
        if ((status = nc_def_dim(ncid_out, "bins_across_track", NX, &x_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(ncid_out, "bins_along_track", NY, &y_dimid)))
            check_err(status, __LINE__, __FILE__);
        //dims for output var
        dimids[0] = y_dimid;
        dimids[1] = x_dimid;
        NDIMS=2;
        if ((status = nc_def_dim(ncid_out, "number_of_views", NVIEWS, &v_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(ncid_out, "intensity_bands_per_view", NBANDS, &b_dimid)))
            check_err(status, __LINE__, __FILE__);
/*
//def variables at the root level!
        if ((status = nc_def_var(ncid_out, "latitude", NC_FLOAT, NDIMS,
                dimids, &varid1)))
                check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_var(ncid_out, "longitude", NC_FLOAT, NDIMS,
                dimids, &varid2)))
                check_err(status, __LINE__, __FILE__);

        if ((status = nc_enddef(ncid_out))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);
        if ((status = nc_put_var_float(ncid_out, varid1, &lat_out[0][0])))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_put_var_float(ncid_out, varid2, &lon_out[0][0])))
            check_err(status, __LINE__, __FILE__);

//define attributes
       if (nc_put_att_text(ncid_out, varid1, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, varid2, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str())) 
           check_err(status, __LINE__, __FILE__);
*/
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME1.c_str(), strlen(GATT_VAL1.c_str()),
          GATT_VAL1.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME2.c_str(), strlen(GATT_VAL2.c_str()),
          GATT_VAL2.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME3.c_str(), strlen(GATT_VAL3.c_str()),
          GATT_VAL3.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME4.c_str(), strlen(GATT_VAL4.c_str()),
          GATT_VAL4.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME5.c_str(), strlen(GATT_VAL5.c_str()),
          GATT_VAL5.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME6.c_str(), strlen(GATT_VAL6.c_str()),
          GATT_VAL6.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME7.c_str(), strlen(GATT_VAL7.c_str()),
          GATT_VAL7.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME8.c_str(), strlen(GATT_VAL8.c_str()),
          GATT_VAL8.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME9.c_str(), strlen(GATT_VAL9.c_str()),
          GATT_VAL9.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME10.c_str(), strlen(GATT_VAL10.c_str()),
          GATT_VAL10.c_str()))
           check_err(status, __LINE__, __FILE__);
        if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME11.c_str(), strlen(GATT_VAL11.c_str()),
          GATT_VAL11.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME12.c_str(), strlen(GATT_VAL12.c_str()),
          GATT_VAL12.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME13.c_str(), strlen(GATT_VAL13.c_str()),
          GATT_VAL13.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME14.c_str(), strlen(GATT_VAL14.c_str()),
          GATT_VAL14.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME15.c_str(), strlen(GATT_VAL15.c_str()),
          GATT_VAL15.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME16.c_str(), strlen(GATT_VAL16.c_str()),
          GATT_VAL16.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME17.c_str(), strlen(GATT_VAL17.c_str()),
          GATT_VAL17.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME18.c_str(), strlen(GATT_VAL18.c_str()),
          GATT_VAL18.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME19.c_str(), strlen(GATT_VAL19.c_str()),
          GATT_VAL19.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME20.c_str(), strlen(GATT_VAL20.c_str()),
          GATT_VAL20.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME21.c_str(), strlen(GATT_VAL21.c_str()),
          GATT_VAL21.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME22.c_str(), strlen(GATT_VAL22.c_str()),
          GATT_VAL22.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME23.c_str(), strlen(GATT_VAL23.c_str()),
          GATT_VAL23.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME24.c_str(), strlen(GATT_VAL24.c_str()),
          GATT_VAL24.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME25.c_str(), strlen(GATT_VAL25.c_str()),
          GATT_VAL25.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME26.c_str(), strlen(GATT_VAL26.c_str()),
          GATT_VAL26.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME27.c_str(), strlen(GATT_VAL27.c_str()),
          GATT_VAL27.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME28.c_str(), strlen(GATT_VAL28.c_str()),
          GATT_VAL28.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME29.c_str(), strlen(GATT_VAL29.c_str()),
          GATT_VAL29.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME30.c_str(), strlen(GATT_VAL30.c_str()),
          GATT_VAL30.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME31.c_str(), strlen(GATT_VAL31.c_str()),
          GATT_VAL31.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME32.c_str(), strlen(GATT_VAL32.c_str()),
          GATT_VAL32.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME33.c_str(), strlen(GATT_VAL33.c_str()),
          GATT_VAL33.c_str()))
          check_err(status, __LINE__, __FILE__);




      //define groups                  check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_grp(ncid_out, "geolocation_data", &grp_coor)))  //netcdf-4
             check_err(status, __LINE__, __FILE__);
/*
//DEF DIMENSIONS
        if ((status = nc_def_dim(grp_coor, "x", NX, &x_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(grp_coor, "y", NY, &y_dimid)))
            check_err(status, __LINE__, __FILE__);
        dimids[0] = y_dimid;
        dimids[1] = x_dimid;
        NDIMS=2;
*/
        //def var      
        if ((status = nc_def_var(grp_coor, "latitude", NC_FLOAT, NDIMS,
                dimids, &varid1)))
                check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_var(grp_coor, "longitude", NC_FLOAT, NDIMS,
                dimids, &varid2)))
                check_err(status, __LINE__, __FILE__);

                    //leave define mode-----------------------
        if ((status = nc_enddef(grp_coor))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);

        //writing the whole thing
        if ((status = nc_put_var_float(grp_coor, varid1, &lat_out[0][0])))
            check_err(status, __LINE__, __FILE__);

        if ((status = nc_put_var_float(grp_coor, varid2, &lon_out[0][0])))
            check_err(status, __LINE__, __FILE__);

//define attributes of vars in groups
        if (nc_put_att_text(grp_coor, varid1, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
             check_err(status, __LINE__, __FILE__);
        if (nc_put_att_text(grp_coor, varid2, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
             check_err(status, __LINE__, __FILE__);


        if ((status = nc_close(ncid_out)))
           check_err(status, __LINE__, __FILE__);


        delete [] (lat_out);
        delete [] (lon_out);
        lat_out=nullptr;
        lon_out=nullptr;


//CASE INVERTED ALONG/ACROSS INDEXES ***************************************************************************
        prodstr = "L1Cgrid_inverted";


        lat_out = allocate2d_float(NY, NX);
        lon_out = allocate2d_float(NY, NX);

        c = 0;
        for (int i = NY1; i < NY2 + 1; i++) {
            for (int j = 0; j < NX; j++) {
                lat_out[NY - 1 - c][NX - 1 - j] = lat_gd[i][j];
                lon_out[NY - 1 - c][j] = lon_gd[i][j];
            }
            c++;
        }
       
        extstr = ".nc";
//        fname_out = pathstr + missionstr + "_" + senstr + "." + yearstr + monstr + daystr + timestr + prodstr+extstr;
        fname_out = pathstr + ofilestr+"_"+ prodstr+extstr;
        l1cfile->gridname = fname_out.c_str();
        filename_lt = fname_out.c_str();
        cout<<"creating file for L1C coor.with x/y inverted indexes.."<<filename_lt<<endl;
//        if ((status = nc_create(filename_lt, NC_CLOBBER, &ncid_out))) 
        if ((status = nc_create(filename_lt, NC_CLOBBER | NC_NETCDF4, &ncid_out)))
            check_err(status, __LINE__, __FILE__);
        //define dims
     // Define the dimensions,vars and attributes at the root level


  //DEF DIMENSIONS
        if ((status = nc_def_dim(ncid_out, "bins_across_track", NX, &x_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(ncid_out, "bins_along_track", NY, &y_dimid)))
            check_err(status, __LINE__, __FILE__);
        //dims for output var
        dimids[0] = y_dimid;
        dimids[1] = x_dimid;
        NDIMS=2;
        if ((status = nc_def_dim(ncid_out, "number_of_views", NVIEWS, &v_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(ncid_out, "intensity_bands_per_view", NBANDS, &b_dimid)))
            check_err(status, __LINE__, __FILE__);
/*
//def variables at the root level!
        if ((status = nc_def_var(ncid_out, "latitude", NC_FLOAT, NDIMS,
                dimids, &varid1)))
                check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_var(ncid_out, "longitude", NC_FLOAT, NDIMS,
                dimids, &varid2)))
                check_err(status, __LINE__, __FILE__);

        if ((status = nc_enddef(ncid_out))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);
        if ((status = nc_put_var_float(ncid_out, varid1, &lat_out[0][0])))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_put_var_float(ncid_out, varid2, &lon_out[0][0])))
            check_err(status, __LINE__, __FILE__);

//define attributes
       if (nc_put_att_text(ncid_out, varid1, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, varid2, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str())) 
           check_err(status, __LINE__, __FILE__);
*/

       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME1.c_str(), strlen(GATT_VAL1.c_str()),
          GATT_VAL1.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME2.c_str(), strlen(GATT_VAL2.c_str()),
          GATT_VAL2.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME3.c_str(), strlen(GATT_VAL3.c_str()),
          GATT_VAL3.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME4.c_str(), strlen(GATT_VAL4.c_str()),
          GATT_VAL4.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME5.c_str(), strlen(GATT_VAL5.c_str()),
          GATT_VAL5.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME6.c_str(), strlen(GATT_VAL6.c_str()),
          GATT_VAL6.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME7.c_str(), strlen(GATT_VAL7.c_str()),
          GATT_VAL7.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME8.c_str(), strlen(GATT_VAL8.c_str()),
          GATT_VAL8.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME9.c_str(), strlen(GATT_VAL9.c_str()),
          GATT_VAL9.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME10.c_str(), strlen(GATT_VAL10.c_str()),
          GATT_VAL10.c_str()))
           check_err(status, __LINE__, __FILE__);
        if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME11.c_str(), strlen(GATT_VAL11.c_str()),
          GATT_VAL11.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME12.c_str(), strlen(GATT_VAL12.c_str()),
          GATT_VAL12.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME13.c_str(), strlen(GATT_VAL13.c_str()),
          GATT_VAL13.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME14.c_str(), strlen(GATT_VAL14.c_str()),
          GATT_VAL14.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME15.c_str(), strlen(GATT_VAL15.c_str()),
          GATT_VAL15.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME16.c_str(), strlen(GATT_VAL16.c_str()),
          GATT_VAL16.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME17.c_str(), strlen(GATT_VAL17.c_str()),
          GATT_VAL17.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME18.c_str(), strlen(GATT_VAL18.c_str()),
          GATT_VAL18.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME19.c_str(), strlen(GATT_VAL19.c_str()),
          GATT_VAL19.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME20.c_str(), strlen(GATT_VAL20.c_str()),
          GATT_VAL20.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME21.c_str(), strlen(GATT_VAL21.c_str()),
          GATT_VAL21.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME22.c_str(), strlen(GATT_VAL22.c_str()),
          GATT_VAL22.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME23.c_str(), strlen(GATT_VAL23.c_str()),
          GATT_VAL23.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME24.c_str(), strlen(GATT_VAL24.c_str()),
          GATT_VAL24.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME25.c_str(), strlen(GATT_VAL25.c_str()),
          GATT_VAL25.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME26.c_str(), strlen(GATT_VAL26.c_str()),
          GATT_VAL26.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME27.c_str(), strlen(GATT_VAL27.c_str()),
          GATT_VAL27.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME28.c_str(), strlen(GATT_VAL28.c_str()),
          GATT_VAL28.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME29.c_str(), strlen(GATT_VAL29.c_str()),
          GATT_VAL29.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME30.c_str(), strlen(GATT_VAL30.c_str()),
          GATT_VAL30.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME31.c_str(), strlen(GATT_VAL31.c_str()),
          GATT_VAL31.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME32.c_str(), strlen(GATT_VAL32.c_str()),
          GATT_VAL32.c_str()))
           check_err(status, __LINE__, __FILE__);
       if (nc_put_att_text(ncid_out, NC_GLOBAL, GATT_NAME33.c_str(), strlen(GATT_VAL33.c_str()),
          GATT_VAL33.c_str()))
          check_err(status, __LINE__, __FILE__);

      //define groups                  check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_grp(ncid_out, "geolocation_data", &grp_coor)))  //netcdf-4
           check_err(status, __LINE__, __FILE__);
/*
//DEF DIMENSIONS
        if ((status = nc_def_dim(grp_coor, "x", NX, &x_dimid)))
            check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_dim(grp_coor, "y", NY, &y_dimid)))
            check_err(status, __LINE__, __FILE__);
        dimids[0] = y_dimid;
        dimids[1] = x_dimid;
        NDIMS=2;
*/
        //def var      
        if ((status = nc_def_var(grp_coor, "latitude", NC_FLOAT, NDIMS,
                dimids, &varid1)))
                check_err(status, __LINE__, __FILE__);
        if ((status = nc_def_var(grp_coor, "longitude", NC_FLOAT, NDIMS,
                dimids, &varid2)))
                check_err(status, __LINE__, __FILE__);

                    //leave define mode-----------------------
        if ((status = nc_enddef(grp_coor))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);

        //writing the whole thing
        if ((status = nc_put_var_float(grp_coor, varid1, &lat_out[0][0])))
            check_err(status, __LINE__, __FILE__);

        if ((status = nc_put_var_float(grp_coor, varid2, &lon_out[0][0])))
            check_err(status, __LINE__, __FILE__);


//define attributes of vars in groups
        if (nc_put_att_text(grp_coor, varid1, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
             check_err(status, __LINE__, __FILE__);
        if (nc_put_att_text(grp_coor, varid2, ATT_NAME.c_str(), strlen(ATT_VAL.c_str()), ATT_VAL.c_str()))
             check_err(status, __LINE__, __FILE__);



        if ((status = nc_close(ncid_out)))
           check_err(status, __LINE__, __FILE__);


        delete [] (lat_out);
        delete [] (lon_out);
        lat_out=nullptr;
        lon_out=nullptr;



//******** save az_east param as nc ***************************************************************
//az_east[count]
        prodstr="az_east";
        extstr = ".nc";
        fname_out = pathstr + ofilestr+"_"+ prodstr+extstr;
        filename_lt = fname_out.c_str();
        cout<<"creating file for az_east values in degrees as nc format.."<<filename_lt<<endl;
//        if ((status = nc_create(filename_lt, NC_CLOBBER, &ncid_out))) 
        if ((status = nc_create(filename_lt, NC_CLOBBER | NC_NETCDF4, &ncid_out)))
            check_err(status, __LINE__, __FILE__);
        //define dims
     // Define the dimensions,vars and attributes at the root level


  //DEF DIMENSIONS
        if ((status = nc_def_dim(ncid_out, "bins_along_track", NY, &y_dimid)))
            check_err(status, __LINE__, __FILE__);
//DEF VARIABLE      
        if ((status = nc_def_var(ncid_out, "az_east", NC_FLOAT, 1,
                &y_dimid, &varid1)))
                check_err(status, __LINE__, __FILE__);

        if ((status = nc_enddef(ncid_out))) //done def vars etc
                        check_err(status, __LINE__, __FILE__);

          for (int count = 0; count < num_gridlines - 2; count++) {
                if (az_east[count] < 0.001) az_east[count] = 0.0;}
        //writing the whole thing
        if ((status = nc_put_var_float(ncid_out, varid1, &az_east[0])))
            check_err(status, __LINE__, __FILE__);
//close file
        if ((status = nc_close(ncid_out)))
           check_err(status, __LINE__, __FILE__);


//***********************************************************************************************
        prodstr = "L1CgridV";
        extstr=".csv";

        fname_out = pathstr + ofilestr+"_"+ prodstr+extstr;
        //   std::ofstream myl1c("/accounts/mamontes/images/OCIS/swt1_ALLgd_L1C500_swtd13again.csv");
        std::ofstream myl1cV(fname_out);
        // for(int count = 0; count <num_gridlines; count ++){
             /*    for(int gp = 0; gp <nbinx/2; gp ++){
                    if (loni_l[gp]<0.) loni_l[gp]+=360;
                    if (loni_r[gp]<0.) loni_r[gp]+=360;
                        myl1c << loni_l[gp] << ","<<lati_l[gp]<<endl;
                        myl1c << loni_r[gp] << ","<<lati_r[gp]<<endl;}*/

        for (int count = 0; count < num_gridlines - 2;count++) {
            //         for(int count = 0; count <num_gridlines/100;count ++){
            for (int gp = 0; gp < nbinx; gp++) {
                lontemp = lon_gd[count][gp];
                if (lon_gd[count][gp] < -180.) lontemp = lon_gd[count][gp] + 360;
                if (lon_gd[count][gp] > +180.) lontemp = lon_gd[count][gp] - 360;
                if (lontemp < 0.) lontemp += 360.;//make better visualization, all lon values positive

                myl1cV << lontemp << "," << lat_gd[count][gp] << endl;            
           }        
         }


        myl1cV.close();


        prodstr = "lat_nadpix";
        fname_out = pathstr + ofilestr+"_"+ prodstr+extstr;
        //   fname_out= "/accounts/mamontes/images/OCIS/out/"+"L1Cgrid_"+"lat_nadpix"+"_"+swtstr+".txt";
        //   ofstream myfile1 ("/accounts/mamontes/images/OCIS/lat_nadpix_swtd13again.txt");
        ofstream myfile1(fname_out);
        if (myfile1.is_open())
        {
            myfile1 << "lat_nadpix.\n";
            for (int count = 0; count < num_gridlines - 2; count++) {
                //   for(int gp = 0; gp <nbinx; gp ++){
                myfile1 << lati2[count] << "," << endl;
            }
            myfile1.close();
        }
        else cout << "Unable to open file";
        prodstr = "lon_nadpix";
        fname_out = pathstr + ofilestr+"_"+ prodstr+extstr;
        //  fname_out= "/accounts/mamontes/images/OCIS/out/"+"L1Cgrid_"+"lon_nadpix"+"_"+swtstr+".txt";
        //  ofstream myfile2 ("/accounts/mamontes/images/OCIS/lon_nadpix_swtd13again.txt");
        ofstream myfile2(fname_out);
        if (myfile2.is_open())
        {
            myfile2 << "lon_nadpix.\n";
            for (int count = 0; count < num_gridlines - 2; count++) {
                //   for(int gp = 0; gp <nbinx; gp ++){
                myfile2 << loni2[count] << "," << endl;
            }
            myfile2.close();
        }
        else cout << "Unable to open file";

        prodstr = "az_east";
        fname_out = pathstr + ofilestr+"_"+ prodstr+extstr;

        l1cfile->azeast_name = fname_out.c_str();



        //  fname_out= "/accounts/mamontes/images/OCIS/out/"+"L1Cgrid_"+"lon_nadpix"+"_"+swtstr+".txt";
        //  ofstream myfile2 ("/accounts/mamontes/images/OCIS/lon_nadpix_swtd13again.txt");
        ofstream myfile3(fname_out);
        if (myfile3.is_open())
        {
            myfile3 << "az_east.\n";
            for (int count = 0; count < num_gridlines - 2; count++) {
                if (az_east[count] < 0.001) az_east[count] = 0.0;
                myfile3 << az_east[count] << "," << endl;
            }
            myfile3.close();
        }
        else cout << "Unable to open file";
        if(lati_l!=nullptr)
        delete[](lati_l);
        if(lati_r!=nullptr)
        delete[](lati_r);
        if(loni_l!=nullptr)
        delete[](loni_l);
        if(loni_r!=nullptr)
        delete[](loni_r);


        cout << " SUCCESS writing results full L1C grid center positions for--swath #..." << swtd << "\n";

        return 0;
    }




int32_t L1C::mov_SOCEA(l1c_filehandle* l1cfile, L1C_input* l1cinput, double* tcross, int16_t* file_id, int16_t* swtd_id, int16_t* nfiles_swath, double* ect_swtd, int16_t* tod, int16_t* orbdir, float* mgv_swath,size_t *num_scans_swt,double *time_swt,float **pos_swt,float **vel_swt) {
        int dimid, status;
        int32_t n_files, nadpix;
        std::string str;
        const char* ptstr;
        float lat1, lat2,***lat_3d=nullptr;
        short ***solz_3d=nullptr;
        int16_t syear, smon, sday,shour,smin, selyear = -1, selmon = -1, selday = -1,syear2, smon2, sday2;
        double seconds;
        size_t att_len, swt = 1, nite = 0;
        float* lat_id=nullptr, ** velr=nullptr,**posr=nullptr,***posr2=nullptr,***velr2=nullptr;
        int16_t* swath_id, * tod_id, * orbdir_id;
        int asc_mode, day_mode;//per file--asceding mode: 1 ascending 0 descending, day_mode: 1 day, 0 ys night
        size_t n_orb_rec;
        int ixfile = 0;
        double sum3=0.,rl2,pos_norm,clatg2,fe=1/298.257;
        double *stime=nullptr;
        float **stime2=nullptr;
        size_t nscans_swt,first_line=5000,last_line=5000; 
        double start_time=-1, end_time=-1;
        string y_swt,mo_swt,d_swt,h_swt,mi_swt,s_swt,tswt_ini,tswt_end,s_swt2,tswt_ini_file,fdatcreate;
        int logoff=-1;
        std::string  ifile_str;
        char* ifile_char;
        int blueGrp,nav1Grp;
        float tmp = 0.0;

        ifile_str = l1cinput->files[0];
        ifile_char = &ifile_str[0];
        file_format format = getFormat(ifile_char);
        l1cfile->format=format.type;

        nframes=l1cfile->nframes;
        nviews = l1cfile->n_views;
        if(format.type==FT_OCIS) nviews=2;

        num_pixels = l1cfile->npix;
        num_scans = l1cfile->nscan;        
        nadpix = l1cfile->nadpix;
        n_files = l1cfile->ifiles.size();

        int16_t nswath = l1cfile->nswath;
        nscans_swt=num_scans*n_files;

        cout<<"format.type..."<<format.type<<"nadpix..."<<nadpix<<"num_pixels.."<<num_pixels<<"nviews.."<<nviews<<"num_scans.."<<num_scans<<"nfiles.."<<n_files<<"scans x files.."<<nscans_swt<<endl;
    
        swath_id = (int16_t*)calloc(n_files,sizeof(int16_t));
        lat_id = (float*)calloc(n_files,sizeof(float));
        tod_id = (int16_t*)calloc(n_files,sizeof(int16_t));
        orbdir_id = (int16_t*)calloc(n_files,sizeof(int16_t));

        //big loop
        for (int i = 0;i < n_files;i++) {
            str = l1cfile->ifiles[i];
            ptstr = str.c_str();
            swath_id[i] = 0;

            // Open the netcdf4 input file
            status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- Error failed nc_open.\n");
                        exit(EXIT_FAILURE);
            }

            
            status = nc_inq_attlen(ncid_L1B, NC_GLOBAL, "date_created", &att_len);
            check_err(status, __LINE__, __FILE__);
            // allocate required space before retrieving values
            char* fdate_create = (char*)malloc(att_len + 1); // + 1 for trailing null
         // get global attribute values
            status = nc_get_att_text(ncid_L1B, NC_GLOBAL, "date_created",fdate_create);
            check_err(status, __LINE__, __FILE__);
            fdate_create[att_len] = '\0';

            l1cfile->date_created=fdate_create;

             //make sure the date (year/month.day) is the same for the bunch of files--
            status = nc_inq_attlen(ncid_L1B, NC_GLOBAL, "time_coverage_start", &att_len);
            check_err(status, __LINE__, __FILE__);
            // allocate required space before retrieving values
            char* time_str = (char*)malloc(att_len + 1); // + 1 for trailing null
         // get global attribute values
            status = nc_get_att_text(ncid_L1B, NC_GLOBAL, "time_coverage_start", time_str);
            check_err(status, __LINE__, __FILE__);
            time_str[att_len] = '\0';
            double start_time = isodate2unix(time_str);
   //         unix2ymds(start_time, &syear, &smon, &sday, &secs);
            unix2ymdhms(start_time, &syear,&smon,&sday,&shour,&smin,&seconds);
        
            delete [](time_str);


//*********************************************************************************************************
//            lat = allocate2d_float(num_scans, num_pixels);

           if(format.type==FT_SPEXONE){
              status = nc_inq_grp_ncid(ncid_L1B, "GEOLOCATION_DATA", &geoGrp);
              check_err(status, __LINE__, __FILE__);
              lat = allocate2d_float(num_scans, num_pixels);
              status = nc_inq_varid(geoGrp, "latitude", &latId);//scans x velements
              check_err(status, __LINE__, __FILE__);
              status = nc_get_var_float(geoGrp, latId, &lat[0][0]);
              check_err(status, __LINE__, __FILE__);
              lat_id[i] = lat[0][nadpix];
              cout<<"granule#..."<<i+1<<"num_pixels.."<<num_pixels<<"lat first pix.."<<lat[0][0]<<"lat last pix.."<<lat[num_scans-1][num_pixels-1]<<endl;
             }
           else if(format.type==FT_HARP2){
              status = nc_inq_grp_ncid(ncid_L1B, "blue", &blueGrp);
              check_err(status, __LINE__, __FILE__);
              lat_3d = allocate3d_float(nviews,num_scans, num_pixels);
              status = nc_inq_varid(blueGrp, "Latitude", &latId);//scans x velements
              check_err(status, __LINE__, __FILE__);
              status = nc_get_var_float(blueGrp, latId, &lat_3d[0][0][0]);
              check_err(status, __LINE__, __FILE__);
              lat_id[i] = lat_3d[4][0][nadpix];
              cout<<"granule#..."<<i+1<<"num_pixels.."<<num_pixels<<"lat first pix.."<<lat_3d[4][0][0]<<"lat last pix.."<<lat_3d[4][num_scans-1][num_pixels-1]<<endl;
             }
            else{//OCIS 
              status = nc_inq_grp_ncid(ncid_L1B, "geolocation_data", &geoGrp);
              check_err(status, __LINE__, __FILE__);
              lat = allocate2d_float(num_scans, num_pixels);
              status = nc_inq_varid(geoGrp, "latitude", &latId);//scans x velements
              check_err(status, __LINE__, __FILE__);
              status = nc_get_var_float(geoGrp, latId, &lat[0][0]);
              check_err(status, __LINE__, __FILE__);
              lat_id[i] = lat[0][nadpix];
              cout<<"granule#..."<<i+1<<"num_pixels.."<<num_pixels<<"lat first pix.."<<lat[0][0]<<"lat last pix.."<<lat[num_scans-1][num_pixels-1]<<endl;
            }


            //opening solz, solza must be <= 90 for daylight
            if(format.type == FT_HARP2){
               solz_3d = allocate3d_short(nviews,num_scans, num_pixels);
               status = nc_inq_varid(blueGrp, "Solar_Zenith", &solzId);//scans x velements
               check_err(status, __LINE__, __FILE__);
               status = nc_get_var_short(blueGrp, solzId, &solz_3d[0][0][0]);

               //check daylight conditions--only checking first scan!! and central pixel!         
               tmp = solz_3d[4][0][nadpix] * 0.005555556 + 0.0;
            }
            else{//OCI/SPEX
               solz = allocate2d_short(num_scans, num_pixels);
               status = nc_inq_varid(geoGrp, "solar_zenith", &solzId);//scans x velements
               check_err(status, __LINE__, __FILE__);
               status = nc_get_var_short(geoGrp, solzId, &solz[0][0]);

               //check daylight conditions--only checking first scan!! and central pixel!
           
               tmp = solz[0][nadpix] * 0.01 + 0.0;
            }


            if (tmp > 93.0) {//nightime, due to non-nadir viewing angle 
                day_mode = 0;
            }
            else {

                day_mode = 1;            
              }//daytime

            day_mode=1;//we only want daytime ---
            tod_id[i] = day_mode;

            for (unsigned int k = 0;k < num_scans;k++) {
                for (unsigned int j = 0;j < num_pixels;j++) {
                    if (tmp < 0.0) {
                        cout << "Error--negative solz[k][j]..." << tmp << endl;
                        exit(1);                    
                      }
//                    if(tmp>90.){ cout<<"nightime.....................at .solar zenith>90..."<<"granule#..."<<i+1<<"scan#.."<<k+1<<"pix#...."<<j+1<<endl;}   
                }            
            }   


            asc_mode = l1cfile->orb_dir;
            orbdir_id[i] = asc_mode;



           if(format.type == FT_HARP2){
//ascending
            for (size_t k =num_scans-1;k>0 ;k--) {
                    //checking partial granules--          
                   if(asc_mode==1 && i==0 && lat_3d[4][k-1][nadpix]>lat_3d[4][k][nadpix] && first_line==5000){
                         first_line=k;
                         cout<<"first line.."<<first_line<<endl;                    
                         }
                    }
             for (size_t k = 0;k<num_scans-1 ;k++) {
                     if(asc_mode==1 && lat_3d[4][k+1][nadpix]<lat_3d[4][k][nadpix] && i==n_files-1 && last_line==5000){//assuming change of direction NH, last granule
                         last_line=k;
                         cout<<"last line.."<<last_line<<endl;
                         }
                     }
//descending
              for (size_t k =num_scans-1;k>0 ;k--) {                     //descending
                      if(asc_mode==0 && lat_3d[4][k-1][nadpix]<lat_3d[4][k][nadpix] && i==0 && first_line==5000){//assuming change of direction SH, first granule
                         first_line=k;
                         cout<<"first line.."<<first_line<<endl;
                         }
              }
              
              for (size_t k = 0;k<num_scans-1 ;k++) {
                     if(asc_mode==0 && lat_3d[4][k+1][nadpix]>lat_3d[4][k][nadpix] && i==n_files-1 && last_line==5000){//assuming change of direction NH, last granule
                         last_line=k;
                         cout<<"last line.."<<last_line<<endl;
                         }
                   }


            //checking min max lat limits------------
              for (size_t k =0;k < num_scans;k++) {
                for (size_t j = 0;j < num_pixels;j++) {
                    //check north/south boundings
                    if (lat_3d[4][k][j] > l1cfile->latnorth) {
                        lat_3d[4][k][j] = l1cfile->latnorth;
                        cout << "lat values above 90 degrees..set value to 90" << endl;                    
                      }
                    else if (lat_3d[4][k][j] < l1cfile->latsouth) {
                        lat_3d[4][k][j] = l1cfile->latsouth;
                        cout << "lat values below 90 degrees.set value to -90." << endl;                    
                        }                
                  }

              } 
              
                 
            //loading selected day,month and year-----    
            if (i == 0) {
                selday = l1cinput->selday;
                selmon = l1cinput->selmon;
                selyear = l1cinput->selyear;
                lat1 = lat_3d[4][0][nadpix];
            }

   
  
            if (syear == selyear && smon == selmon && sday == selday) {
                //segment swath based on orbit direction and lat 
                //check if lat between min and max
                //check solar zenith for only daylight granules >90 degrees is nighttime
                //assuming files are ordered by time when creating the initial list
                //we assume only daytime files
                //we compare the nadir pixel

                 //ascending or descending  orbit---------
                if (i >= 0 && day_mode == 0) { //nightime--
                    swt++;
                    nite++;
                    swath_id[i] = swt;
                    cout << "nightime file....swath #.." << swt << endl;
                    if (tcross[i] > 0.0) {
                        cout << "nightime eq crossing..for file [" << i << "]" << endl;

                    }
                }

                if (i == 0 && day_mode == 1) { //first file---
                    lat1 = lat_3d[4][0][nadpix];
                    lat2 = lat_3d[4][num_scans - 1][nadpix];
                    cout << "FIRST filename.." << ptstr << "file number..." << i + 1 << "asc_mode.." << asc_mode << "daymode.." << day_mode << "swt..." << swt << "lat1..(NADPIX)." << lat1 << "lat2...(NADPIX)" << lat2 << endl;
                    swath_id[i] = swt;
                    lat1 = lat2;
                }

                //ascending orbit
                if (i >= 1 && asc_mode == 1 && day_mode == 1)//at least two files to compare
                {   lat1=lat2;
                    lat2 = lat_3d[4][0][nadpix];
                    if (tod_id[i - 1] == 0) swt++;
                    cout << "filename.." << ptstr << "file number..." << i + 1 << "asc_mode.." << asc_mode << "daymode.." << day_mode << "swt..." << swt << "lat1..(NADPIX)." << lat1 << "lat2...(NADPIX)" << lat2 << endl;
                    swath_id[i] = swt;
                }

            }//if same day


           }
//--------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------           
           else if(format.type == FT_SPEXONE){
                swt=0;
                //I need to count the number of swaths                
                 //ascending or descending  orbit---------
                if (i >= 0 && day_mode == 0) { //nightime--
                    swt++;
                    nite++;
                    swath_id[i] = swt;
                    cout << "nightime file....swath #.." << swt << endl;
                    if (tcross[i] > 0.0) {
                        cout << "nightime eq crossing..for file [" << i << "]" << endl;

                    }
                }
                //only 1 orbit file!!!!
                if (i == 0 && day_mode == 1) { //first file---
                    lat1 = lat[0][nadpix];
                    lat2 = lat[num_scans - 1][nadpix];
                    cout << "FIRST filename.." << ptstr << "file number..." << i + 1 << "asc_mode.." << asc_mode << "daymode.." << day_mode << "swt..." << swt << "lat1..(NADPIX)." << lat1 << "lat2...(NADPIX)" << lat2 << endl;
                    swath_id[i] = swt;
                    lat1 = lat2;
                }

                //ascending orbit
                if (i >= 1 && asc_mode == 1 && day_mode == 1)//at least two files to compare
                {   lat1=lat2;
                    lat2 = lat[0][nadpix];
                    if (tod_id[i - 1] == 0) swt++;
                    cout << "filename.." << ptstr << "file number..." << i + 1 << "asc_mode.." << asc_mode << "daymode.." << day_mode << "swt..." << swt << "lat1..(NADPIX)." << lat1 << "lat2...(NADPIX)" << lat2 << endl;
                    swath_id[i] = swt;
                }
            

           } 
           else{   //OCI
//ascending
              for (size_t k =num_scans-1;k>0 ;k--) {
                    //checking partial granules--          
                   if(asc_mode==1 && i==0 && lat[k-1][nadpix]>lat[k][nadpix] && first_line==5000){
                         first_line=k;
                         cout<<"first line.."<<first_line<<endl;                    
                         }
                    }
              for (size_t k = 0;k<num_scans-1 ;k++) {
                     if(asc_mode==1 && lat[k+1][nadpix]<lat[k][nadpix] && i==n_files-1 && last_line==5000){//assuming change of direction NH, last granule
                         last_line=k;
                         cout<<"last line.."<<last_line<<endl;
                         }
                     }
//descending
              for (size_t k =num_scans-1;k>0 ;k--) {                     //descending
                      if(asc_mode==0 && lat[k-1][nadpix]<lat[k][nadpix] && i==0 && first_line==5000){//assuming change of direction SH, first granule
                         first_line=k;
                         cout<<"first line.."<<first_line<<endl;
                         }
                   }
              
              for (size_t k = 0;k<num_scans-1 ;k++) {
                     if(asc_mode==0 && lat[k+1][nadpix]>lat[k][nadpix] && i==n_files-1 && last_line==5000){//assuming change of direction NH, last granule
                         last_line=k;
                         cout<<"last line.."<<last_line<<endl;
                         }
                   }


            //checking min max lat limits------------
              for (size_t k =0;k < num_scans;k++) {
                for (size_t j = 0;j < num_pixels;j++) {
                    //check north/south boundings
                    if (lat[k][j] > l1cfile->latnorth) {
                        lat[k][j] = l1cfile->latnorth;
                        cout << "lat values above 90 degrees..set value to 90" << endl;                    
                      }
                    else if (lat[k][j] < l1cfile->latsouth) {
                        lat[k][j] = l1cfile->latsouth;
                        cout << "lat values below 90 degrees.set value to -90." << endl;                    
                        }                
                  }

              } 
              
                 
            //loading selected day,month and year-----    
            if (i == 0) {
                selday = l1cinput->selday;
                selmon = l1cinput->selmon;
                selyear = l1cinput->selyear;
                lat1 = lat[0][nadpix];
            }

    
               
            if (syear == selyear && smon == selmon && sday == selday) {
                //segment swath based on orbit direction and lat 
                //check if lat between min and max
                //check solar zenith for only daylight granules >93 degrees is nighttime
                //assuming files are ordered by time when creating the initial list
                //we assume only daytime files
                //we compare the nadir pixel
                
     
                 //ascending or descending  orbit---------
                if (i >= 0 && day_mode == 0) { //nightime--
                    swt++;
                    nite++;
                    swath_id[i] = swt;
                    cout << "nightime file....swath #.." << swt << endl;
                    if (tcross[i] > 0.0) {
                        cout << "nightime eq crossing..for file [" << i << "]" << endl;

                    }
                }

                if (i == 0 && day_mode == 1) { //first file---
                    lat1 = lat[0][nadpix];
                    lat2 = lat[num_scans - 1][nadpix];
                    cout << "FIRST filename.." << ptstr << "file number..." << i + 1 << "asc_mode.." << asc_mode << "daymode.." << day_mode << "swt..." << swt << "lat1..(NADPIX)." << lat1 << "lat2...(NADPIX)" << lat2 << endl;
                    swath_id[i] = swt;
                    lat1 = lat2;
                }

                //ascending orbit
                if (i >= 1 && asc_mode == 1 && day_mode == 1)//at least two files to compare
                {   lat1=lat2;
                    lat2 = lat[0][nadpix];
                    if (tod_id[i - 1] == 0) swt++;
                    cout << "filename.." << ptstr << "file number..." << i + 1 << "asc_mode.." << asc_mode << "daymode.." << day_mode << "swt..." << swt << "lat1..(NADPIX)." << lat1 << "lat2...(NADPIX)" << lat2 << endl;
                    swath_id[i] = swt;
                }

              }//if same day

          }


           //close file   
            status = nc_close(ncid_L1B);
            check_err(status, __LINE__, __FILE__);
            if (lat != nullptr)
               delete[](lat);
            if (solz != nullptr)
               delete[](solz);
            if (lat_3d != nullptr)
               delete [] (lat_3d);
            if (solz_3d != nullptr)
               delete [] (solz_3d);

        }//end big loop of granules



//-------------------------------------------------------------------------------------------        
        cout << "#nightime swaths..." << nite << "daytime swaths..." << swt - nite << endl;


        //obtain swath files qith equat crossing-- one swath at the time
        int16_t swtc = 1;
        float vxyz = 0.0, sum = 0.0, mov = 0.0, mgv = 0.0;//Re earth radius in km, hsat sat height above earth surface in km
        size_t num_scans_tot = 0;
        int16_t swtc_id[nswath];
        int16_t swt_id = 0, cross_id = 1;

        cout << "nswath..or total (day/night)  crossing swaths..." << nswath << endl;

        //find crossing swath ids
        int k = 0;
        for (int i = 0;i < n_files;i++) {
            if (tcross[i] > 0.0 && tod_id[i] == 1) {
                swtc_id[k] = swath_id[i];
                k++;            
           }        
         }  
        ndswaths = k;


        //dynamic vector to store file_id
        std::vector<int16_t> swt_files, dswaths;

        //compute mgv for each swath--(day/night or asc/desc) but always crossing---
        int j = 0, tc = 0;
        size_t r=0;
        int16_t nfiles_swt = 0;
        while (j < nswath) {
            cout<<"computing mean velocity for the swath #..........................................................................................................."<<swtc_id[j]<<endl;
            swt_id = swtc_id[j];
            for (int i = 0;i < n_files;i++) {
                if (swt_id == swath_id[i] && tod_id[i] == 1) {
                    swt_files.push_back(i + 1);
                    ect_swtd[tc] = tcross[i];
                    // dswaths.push_back(swt_id);
                    dswaths.push_back(cross_id);
                    ixfile = i;
                    cout << "processing daytime swath with eq cross #...." << swtc << "file #..." << i + 1 << "swath_id.." << swath_id[i] << endl;
                    nfiles_swt++;

                    // Open the netcdf4 input file
                    str = l1cfile->ifiles[i];
                    ptstr = str.c_str();

                    status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
                    if (status != NC_NOERR) {
                         fprintf(stderr, "-E- failed nc_open.\n");
                        exit(EXIT_FAILURE);
                    }

         //the number of scans can change from 1 granule to another.....
                   if(format.type == FT_SPEXONE){
                       status = nc_inq_dimid(ncid_L1B, "bins_along_track", &dimid);
                       if (status != NC_NOERR) {
                           fprintf(stderr, "-E- Error reading number_of_scans.\n");
                           exit(EXIT_FAILURE);
                         }
                        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

                       status = nc_inq_dimid(ncid_L1B, "SC_records", &dimid);
                       if (status != NC_NOERR) {
                           fprintf(stderr, "-E- Error reading number_of_orbital records.\n");
                           exit(EXIT_FAILURE);
                         }
                        nc_inq_dimlen(ncid_L1B, dimid, &n_orb_rec);

                        stime = new double[n_orb_rec]();
                        status = nc_inq_grp_ncid(ncid_L1B, "NAVIGATION_DATA", &scGrp);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_inq_varid(scGrp, "orb_time", &otId);
                        status = nc_get_var_double(scGrp, otId, stime);
                        check_err(status, __LINE__, __FILE__);
                       
                        posr = allocate2d_float(n_orb_rec, 3);
                        velr = allocate2d_float(n_orb_rec, 3);

                       ////open velocity vectors
                        status = nc_inq_varid(scGrp, "orb_pos", &opId);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_get_var_float(scGrp, opId, &posr[0][0]);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_inq_varid(scGrp, "orb_vel", &ovId);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_get_var_float(scGrp, ovId, &velr[0][0]);
                        check_err(status, __LINE__, __FILE__);
                      }
                   else if(format.type == FT_HARP2){
                       status = nc_inq_dimid(ncid_L1B, "Swath_Lines", &dimid);
                       if (status != NC_NOERR) {
                           fprintf(stderr, "-E- Error reading number_of_lines.\n");
                           exit(EXIT_FAILURE);
                         }
                        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

                        //reading subgroup navigation
                        status = nc_inq_grp_ncid(blueGrp, "navigation", &nav1Grp);
                        if (status != NC_NOERR) {
                             fprintf(stderr, "-E- Error reading subgroup navigation --blue.\n");
                            exit(EXIT_FAILURE);
                         }
                                                                                                         
                        stime2 =allocate2d_float(nviews, nframes);
                        status = nc_inq_varid(nav1Grp, "Sec_of_Day", &otId);
                        status = nc_get_var_float(nav1Grp, otId, &stime2[0][0]);
                        check_err(status, __LINE__, __FILE__);

                        posr2 = allocate3d_float(nviews,nframes,3);
                        velr2 = allocate3d_float(nviews,nframes,3);

                       ////open velocity vectors
                        status = nc_inq_varid(nav1Grp, "Orb_Position", &opId);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_get_var_float(nav1Grp, opId, &posr2[0][0][0]);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_inq_varid(nav1Grp, "Orb_Velocity", &ovId);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_get_var_float(nav1Grp, ovId, &velr2[0][0][0]);
                        check_err(status, __LINE__, __FILE__);
                      }
                   else{    //OCI                   
                       status = nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid);                      
                       if (status != NC_NOERR) {
                           fprintf(stderr, "-E- Error reading number_of_scans.\n");
                           exit(EXIT_FAILURE);
                         }
                        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

                        stime = new double[num_scans]();
                        status = nc_inq_grp_ncid(ncid_L1B, "scan_line_attributes", &scGrp);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_inq_varid(scGrp, "time", &otId);
                        status = nc_get_var_double(scGrp, otId, stime);
                        check_err(status, __LINE__, __FILE__);

                        n_orb_rec = num_scans;
                        posr = allocate2d_float(n_orb_rec, 3);
                        velr = allocate2d_float(n_orb_rec, 3);

                     ////open velocity vectors
                        status = nc_inq_grp_ncid(ncid_L1B, "navigation_data", &navGrp);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_inq_varid(navGrp, "orb_pos", &opId);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_get_var_float(navGrp, opId, &posr[0][0]);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_inq_varid(navGrp, "orb_vel", &ovId);
                        check_err(status, __LINE__, __FILE__);
                        status = nc_get_var_float(navGrp, ovId, &velr[0][0]);
                        check_err(status, __LINE__, __FILE__);
                      }
                         
                    //compute vx velocity for each scanline or orb_vector---
                    if(i==0 && first_line==5000) first_line=0;
                    if(i==n_files-1 && last_line==5000) last_line=num_scans-1;

                    cout<<"granule#..."<<i+1<<"first line..."<<first_line<<"last line.."<<last_line<<endl;


//********* time metadata **************************************************************
                    if(format.type == FT_HARP2){
                        if(i==0) start_time=stime2[4][first_line];//view 5 or smalles nadir angle
                        if(i==n_files-1) end_time=stime2[4][last_line];
                    }
                    else{//OCI/SPEX
                          if(i==0) start_time=stime[first_line];//view 5 or smalles nadir angle
                          if(i==n_files-1) end_time=stime[last_line];
                    }




//ini time swath
               //     unix2ymds(start_time, &syear, &smon, &sday,&secs);
  //                  double tfile_ini_sec=ymds2unix(syear,smon,sday,secs);  
//                    unix2ymdhms(tfile_ini_sec,&syear,&smon,&sday, &shour, &smin, &seconds); 
//                    
                    unix2ymdhms(start_time, &syear2,&smon2,&sday2,&shour,&smin,&seconds);

                    y_swt=std::to_string(syear);
                    mo_swt=std::to_string(smon);
                    d_swt=std::to_string(sday);
                    h_swt=std::to_string(shour);
                    mi_swt=std::to_string(smin);
                    s_swt2=std::to_string(round(seconds));

                    length = (int) floor( log10 (smon) ) + 1;
                    if(length==1) mo_swt="0"+mo_swt;
                    length = (int) floor( log10 (sday) ) + 1;
                    if(length==1) d_swt="0"+d_swt;

                    if(shour==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (shour+logoff)) + 1;
                    if(length==1) h_swt="0"+h_swt;

                    if(smin==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (smin+logoff)) + 1;
                    if(length==1) mi_swt="0"+mi_swt;

                    if(seconds==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (round(seconds)+logoff)) + 1;
                    if(length==1) s_swt2="0"+s_swt2;

                    tswt_ini=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);
                    tswt_ini_file=y_swt+mo_swt+d_swt+"T"+h_swt+mi_swt+s_swt2.substr(0,2);

                    l1cfile->tswt_ini=tswt_ini;
                    l1cfile->tswt_ini_file=tswt_ini_file;

//                    cout<<"tswt_ini_file.."<<tswt_ini_file<<"s_swt2.substr(0,2)..."<<s_swt2.substr(0,2)<<endl;
//                    exit(1);


//end time swath
              
                    unix2ymdhms(end_time, &syear2,&smon2,&sday2,&shour,&smin,&seconds);

                    y_swt=std::to_string(syear);
                    mo_swt=std::to_string(smon);
                    d_swt=std::to_string(sday);
                    h_swt=std::to_string(shour);
                    mi_swt=std::to_string(smin);
                    s_swt2=std::to_string(round(seconds));

                    length = (int) floor( log10 (smon) ) + 1;
                    if(length==1) mo_swt="0"+mo_swt;
                    length = (int) floor( log10 (sday) ) + 1;
                    if(length==1) d_swt="0"+d_swt;

                    if(shour==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (shour+logoff)) + 1;
                    if(length==1) h_swt="0"+h_swt;

                    if(smin==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (smin+logoff)) + 1;
                    if(length==1) mi_swt="0"+mi_swt;

                    if(seconds==0) logoff=1; else logoff=0;
                    length = (int) floor( log10 (round(seconds)+logoff)) + 1;
                    if(length==1) s_swt2="0"+s_swt2;

                    tswt_end=y_swt+"-"+mo_swt+"-"+d_swt+"T"+h_swt+":"+mi_swt+":"+s_swt2.substr(0,2);

                    l1cfile->tswt_end=tswt_end;
                                            
//********************i*************************************************************************************
                    if(format.type == FT_SPEXONE){
                      num_records=n_orb_rec;
                       }
                    else if(format.type== FT_HARP2){
                      num_records=nframes;
                      }  
                    else {num_records=num_scans;}//OCI

//compute sum3 and swath vel/pos
//
                   if(format.type== FT_HARP2){
                    for (size_t j = 0;j < num_records;j++) {
                      if((i==0 && j>=first_line) || (i==n_files-1 && j<=last_line) || (i>0 && i<n_files-1 && j>=0 && j<=num_records-1)){
                        vxyz = sqrt(velr2[4][j][0] * velr2[4][j][0] + velr2[4][j][1] * velr2[4][j][1] + velr2[4][j][2] * velr2[4][j][2]);//units m/s 
                   
                        sum += + vxyz;//ORB VELOCITY

                        pos_norm=sqrt(posr2[4][j][0]*posr2[4][j][0]+posr2[4][j][1]*posr2[4][j][1]+posr2[4][j][2]*posr2[4][j][2]);
                        clatg2=sqrt(posr2[4][j][0]*posr2[4][j][0]+posr2[4][j][1]*posr2[4][j][1])/pos_norm;
                        rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                        sum3+=vxyz*rl2/pos_norm;//ground velocity

                        time_swt[r]=stime2[4][j];//orb info                  

                        pos_swt[r][0]=posr2[4][j][0];
                        pos_swt[r][1]=posr2[4][j][1];    
                        pos_swt[r][2]=posr2[4][j][2];

                        vel_swt[r][0]=velr2[4][j][0];
                        vel_swt[r][1]=velr2[4][j][1];
                        vel_swt[r][2]=velr2[4][j][2];

                        r++;
                      }
                    }  
                   }
                   else{ //OCI/SPEX
                     for (size_t j = 0;j < num_records;j++) {
                      if((i==0 && j>=first_line) || (i==n_files-1 && j<=last_line) || (i>0 && i<n_files-1 && j>=0 && j<=num_records-1)){
                        vxyz = sqrt(velr[j][0] * velr[j][0] + velr[j][1] * velr[j][1] + velr[j][2] * velr[j][2]);//units m/s

                        sum += + vxyz;//ORB VELOCITY

                        pos_norm=sqrt(posr[j][0]*posr[j][0]+posr[j][1]*posr[j][1]+posr[j][2]*posr[j][2]);
                        clatg2=sqrt(posr[j][0]*posr[j][0]+posr[j][1]*posr[j][1])/pos_norm;
                        rl2=Re*(1-fe)/(sqrt(1-(2-fe)*fe*clatg2*clatg2));

                        sum3+=vxyz*rl2/pos_norm;//ground velocity

                        time_swt[r]=stime[j];//orb info

                        pos_swt[r][0]=posr[j][0];
                        pos_swt[r][1]=posr[j][1];
                        pos_swt[r][2]=posr[j][2];

                        vel_swt[r][0]=velr[j][0];
                        vel_swt[r][1]=velr[j][1];
                        vel_swt[r][2]=velr[j][2];

                        r++;
                      }
                    }
                   }


                    if(i==0) num_scans_tot=num_scans_tot+num_records-first_line;
                    else if(i==n_files-1) num_scans_tot=num_scans_tot+last_line+1;                            
                    else num_scans_tot += num_records;

                    tc++;

                    status = nc_close(ncid_L1B);
                    check_err(status, __LINE__, __FILE__);

                    delete[](posr);
                    delete[](velr);
                    delete[](posr2);
                    delete[](velr2);
                    delete[](stime);
                    delete[](stime2);
                }//end if
            }//end for


            tod[swtc - 1] = tod_id[ixfile];//using the last file of the swath
            orbdir[swtc - 1] = orbdir_id[ixfile];
            if(format.type==FT_SPEXONE || format.type==FT_HARP2){
               mov = (sum / num_scans_tot);//mean velocity in m/s x 10^3
               mgv = (sum3 / num_scans_tot)*1000;
              }
            else{//OCI
               mov = (sum / num_scans_tot) * 1000;//mean velocity in m/s x 10^3
               mgv = (sum3 / num_scans_tot) * 1000;
            }


            cout<<"# scans total..swath."<<num_scans_tot<<"first line.."<<first_line<<"last line.."<<last_line<<endl;

            *num_scans_swt=num_scans_tot;

            mgv_swath[swtc - 1] = mgv;
            nfiles_swath[swtc - 1] = nfiles_swt;


            cout<<"# of orbital elements for the swath........"<<r<<"mov..in m/s : "<<mov<<"mgv..."<<mgv<<endl;

            num_scans_tot = 0;
            vxyz = 0.0;
            sum = 0.0;
            sum3=0.;
            j++;
            swtc++;
            nfiles_swt = 0;
            cross_id++;
            r=0;
        }//while end for vel processing



        nswtfiles = swt_files.size();

        l1cfile->nswt_files = nswtfiles;
        l1cfile->ndswaths = ndswaths;

        k = 0; //day swaths only!! 
        for (std::vector<int16_t>::const_iterator c = swt_files.begin(); c != swt_files.end(); ++c) {
            file_id[k] = *c;
            k++;        
           }

        k = 0;
        for (std::vector<int16_t>::const_iterator c = dswaths.begin(); c != dswaths.end(); ++c) {
            swtd_id[k] = *c;
            k++;        
           }  

        swt_files.clear();
        dswaths.clear();

        delete[](tod_id);
        delete[](orbdir_id);
        delete[](swath_id);
        delete[](lat_id);

        return 0;
    }



//same as ect_sf but including HARP AND SPEX sensors---
int32_t L1C::ect_sf2(const char* filename, L1C_input* l1cinput,l1c_filehandle* l1cfile) { //ect version for single file
        int dimid,dimidz,dimidx,dimidy,status,blueGrp,v1Grp,nav1Grp,dimidv;
        int asc_mode = 0, scan_brack[2] = { -1,-1 };
        float lat1, lat2, lon1, lon2;//lat,lon defined globally
        float ***latframe=nullptr,***lonframe=nullptr,**lat=nullptr,**lon=nullptr;
        double t1 = -1, t2 = -1;
        double tcross=-1,secs;
        float loncross=-999.;
        std::string  ifile_str;
        char* ifile_char;
        const char *nscan_str,*npix_str;
        double *stime=nullptr;
        float **stime2=nullptr;
        int16_t syear, smon, sday, selyear = -1, selmon = -1, selday = -1;
        size_t att_len;

        l1cfile->l1b_name = filename;
        selday = l1cinput->selday;
        selmon = l1cinput->selmon;
        selyear = l1cinput->selyear;


        l1cfile->l1b_name = filename;

        // Open the netcdf4 input file
        status = nc_open(filename, NC_NOWRITE, &ncid_L1B);

        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error failed nc_open...\n");
            exit(EXIT_FAILURE);
        }

        //Open dimensions
        // num_scans
       ifile_str = l1cinput->files[0];
       ifile_char = &ifile_str[0];
       file_format format = getFormat(ifile_char);

       format.type=l1cfile->format;
       cout<<"format type..in ect_sf2"<<format.type<<endl;

        //checking input date
       status = nc_inq_attlen(ncid_L1B, NC_GLOBAL, "time_coverage_start", &att_len);
            check_err(status, __LINE__, __FILE__);
            // allocate required space before retrieving values
            char* time_str = (char*)malloc(att_len + 1); // + 1 for trailing null
         // get global attribute values
            status = nc_get_att_text(ncid_L1B, NC_GLOBAL, "time_coverage_start", time_str);
            check_err(status, __LINE__, __FILE__);
            time_str[att_len] = '\0';
            double start_time = isodate2unix(time_str);
            unix2ymds(start_time, &syear, &smon, &sday, &secs);

       if (syear == selyear && smon == selmon && sday == selday) {
          cout<<"computing equatorial crossing time...............................................................for file."<<filename<<endl;
       }
       else{
         cout<<"input date does not match year.."<<syear<<"month.."<<smon<<"day.."<<sday<<" extracted from the l1b file....."<<endl;
      //   exit(1);
       }

//*******legacy sensors including OCI **************
      if (format.type != FT_SPEXONE && format.type != FT_HARP2 && format.type != FT_HARP){
        cout<<"reading num scan for legacy sensors"<<endl;  
        status = nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_scans.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

        //across-track bins
        status = nc_inq_dimid(ncid_L1B, "ccd_pixels", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_pixels per line...\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);

        }

// *********** SPEXONE **************************
       if (format.type == FT_SPEXONE){
          cout<<"Processing SPEXone equat crossing time......"<<endl;

          nscan_str="bins_along_track";
          npix_str="spatial_samples_per_image";
          status = nc_inq_dimid(ncid_L1B,nscan_str, &dimidy);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_scans.\n");
            exit(EXIT_FAILURE);           
            }
           status = nc_inq_dimid(ncid_L1B,npix_str, &dimidx);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_pixels per line...\n");
            exit(EXIT_FAILURE);
           }
           status = nc_inq_dimid(ncid_L1B,"number_of_views", &dimidv);
          if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading dimensions number_of views.\n");
            exit(EXIT_FAILURE);
          }
          nc_inq_dimlen(ncid_L1B, dimidy, &num_scans);
          nc_inq_dimlen(ncid_L1B, dimidx, &num_pixels);
          nc_inq_dimlen(ncid_L1B, dimidv, &nviews);
         } 

 //********** HARP2 ***************************      
       if (format.type == FT_HARP2){
          status = nc_inq_grp_ncid(ncid_L1B, "blue", &blueGrp);
          if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading group blue.\n");
            exit(EXIT_FAILURE);
          }
          status = nc_inq_dimid(blueGrp,"Views", &dimidv);
          if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading dimensions number_of views.\n");
            exit(EXIT_FAILURE);
          }
       //reading subgroup navigation
          status = nc_inq_grp_ncid(blueGrp, "navigation", &nav1Grp);
          if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading subgroup navigation --blue.\n");
            exit(EXIT_FAILURE);
          }
          status = nc_inq_dimid(nav1Grp,"Frames", &dimidz);
          if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of frames.\n");
            exit(EXIT_FAILURE);
          }  
          
          status = nc_inq_dimid(ncid_L1B,"Swath_Lines", &dimidy);
          if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_scans.\n");
            exit(EXIT_FAILURE);
          }
          status = nc_inq_dimid(ncid_L1B,"Swath_Pixels", &dimidx);
          if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_pixels per line...\n");
            exit(EXIT_FAILURE);
            }
         nc_inq_dimlen(ncid_L1B, dimidz, &num_frames);
         nc_inq_dimlen(ncid_L1B, dimidy, &num_scans);
         nc_inq_dimlen(ncid_L1B, dimidx, &num_pixels);
         nc_inq_dimlen(ncid_L1B, dimidv, &nviews);

         l1cfile->nframes=num_frames;

         cout<<"#views.."<<nviews<<"#Frames.."<<num_frames<<"#lines.."<<num_scans<<"#pixels.."<<num_pixels<<endl;
        }

 //********** HARP ***************************      
       if (format.type == FT_HARP){
         status = nc_inq_grp_ncid(ncid_L1B, "blue", &blueGrp);
       if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading group blue.\n");
            exit(EXIT_FAILURE);
          }
         status = nc_inq_dimid(blueGrp,"Views", &dimidv);
       if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading dimensions number_of views.\n");
            exit(EXIT_FAILURE);
          }
       //reading subgroup view1
          status = nc_inq_grp_ncid(blueGrp, "View01", &v1Grp);
       if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading subgroup view 1 --blue.\n");
            exit(EXIT_FAILURE);
          }
   //lines and pixels    
          status = nc_inq_dimid(v1Grp,"Swath_Pixels", &dimidx);
       if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading # pixel dimensions view01.\n");
            exit(EXIT_FAILURE);
          }
         status = nc_inq_dimid(v1Grp,"Swath_Lines", &dimidy);
       if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading # lines dimensions view01.\n");
            exit(EXIT_FAILURE);
          }
         nc_inq_dimlen(ncid_L1B, dimidv, &nviews);
         nc_inq_dimlen(ncid_L1B, dimidy, &num_scans);
         nc_inq_dimlen(ncid_L1B, dimidx, &num_pixels);
        }

        l1cfile->nscan = num_scans;
        l1cfile->npix = num_pixels;
        l1cfile->n_views = nviews;

  //time line and geolocation attributes------     
  //************************************************************
  //
       if(format.type == FT_SPEXONE){
        stime = new double[num_scans]();    
        lat = allocate2d_float(num_scans, num_pixels);
        lon = allocate2d_float(num_scans, num_pixels);              
        status = nc_inq_grp_ncid(ncid_L1B, "BIN_ATTRIBUTES", &scGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(scGrp, "image_time", &otId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(scGrp, otId, stime);
        check_err(status, __LINE__, __FILE__);
        lat = allocate2d_float(num_scans, num_pixels);
        lon = allocate2d_float(num_scans, num_pixels);
        status = nc_inq_grp_ncid(ncid_L1B, "GEOLOCATION_DATA", &navGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(navGrp, "latitude", &latId);//scans x velements
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(navGrp, latId, &lat[0][0]);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(navGrp, "longitude", &lonId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(navGrp, lonId, &lon[0][0]);
        check_err(status, __LINE__, __FILE__);
        }

       else if(format.type == FT_HARP2){
        stime2 = allocate2d_float(nviews, num_frames);
        status = nc_inq_varid(nav1Grp, "Sec_of_Day", &otId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(nav1Grp, otId, &stime2[0][0]);
        check_err(status, __LINE__, __FILE__);

        latframe = allocate3d_float(nviews,num_scans, num_pixels);
        lonframe = allocate3d_float(nviews,num_scans, num_pixels);   
   
        status = nc_inq_varid(blueGrp, "Latitude", &latId);//scans x velements
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(blueGrp, latId, &latframe[0][0][0]);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(blueGrp, "Longitude", &lonId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(blueGrp, lonId, &lonframe[0][0][0]);
        check_err(status, __LINE__, __FILE__);
       }
        else if(format.type == FT_HARP){
         cout<<"checking lat/lon arrays...HARP.."<<endl;   
         stime = new double[num_scans]();
         lat = allocate2d_float(num_scans, num_pixels);
         lon = allocate2d_float(num_scans, num_pixels);
   /*     stime = new double[num_frames]();
        status = nc_inq_grp_ncid(ncid_L1B, "NAVIGATION", &scGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(scGrp, "msec", &otId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(scGrp, otId, stime);
        check_err(status, __LINE__, __FILE__);
        latframe = allocate3d_float(num_frames,num_scans, num_pixels);
        lonframe = allocate3d_float(num_frames,num_scans, num_pixels);
        status = nc_inq_grp_ncid(ncid_L1B, "Sensor1", &sen1Grp);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_varid(sen1Grp, "Latitude", &latId);//scans x velements
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(sen1Grp, latId, &latframe[0][0][0]);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(sen1Grp, "Longitude", &lonId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(sen1Grp, lonId, &lonframe[0][0][0]);
        check_err(status, __LINE__, __FILE__);
        */
       }
       else{ //legacy sensors including OCI
        stime = new double[num_scans]();   
        status = nc_inq_grp_ncid(ncid_L1B, "scan_line_attributes", &scGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(scGrp, "time", &otId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(scGrp, otId, stime);
        check_err(status, __LINE__, __FILE__); 
        lat = allocate2d_float(num_scans, num_pixels);
        lon = allocate2d_float(num_scans, num_pixels);   
        status = nc_inq_grp_ncid(ncid_L1B, "geolocation_data", &navGrp);
        check_err(status, __LINE__, __FILE__);         
        status = nc_inq_varid(navGrp, "latitude", &latId);//scans x velements
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(navGrp, latId, &lat[0][0]);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(navGrp, "longitude", &lonId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(navGrp, lonId, &lon[0][0]);
        check_err(status, __LINE__, __FILE__);
        }

   
        //find scans bracketing equator
        int nad_pix = floor(num_pixels / 2) + 1;
        int timeix=0;
        nad_pix = nad_pix - 1;
        l1cfile->nadpix = nad_pix;

       if(format.type == FT_HARP2 | format.type == FT_HARP){
           //blue channel
        if (latframe[4][0][nad_pix] > latframe[4][1][nad_pix])
            asc_mode = 0;//descending
        else
            asc_mode = 1;//ascending
           }
        else{//OCIS/SPEXONE
         if (lat[0][nad_pix] > lat[1][nad_pix])
            asc_mode = 0;//descending
         else
            asc_mode = 1;//ascending
           }    


        l1cfile->orb_dir = asc_mode;


 //       for (unsigned i = 0;i < num_scans - 1;i++) { 
//             cout<<"lat..."<<latframe[4][i][nad_pix]<<endl;}

//HARP2-----------
//use view nadir to compute crossing time!!!!
//
//--------------------------------------------------
//blue
//Views = -49.81, -40.21, -28.92, -15.83, -1.6, 5.55, 19.6, 32.13, 42.78,
//      51.82
//
//
//green
//Views = -52.68, -43.63, -32.92, -20.33, -6.34, 0.81, 15.04, 28.2, 39.47, 
//      49.03 ;
//
////red
//      Views = -55.85, -54.54, -53.24, -51.87, -50.46, -48.99, -47.46, -45.94,
//      -44.36, -42.72, -41.02, -39.27, -37.52, -35.7, -33.83, -31.85, -29.86,
//      -27.86, -25.75, -23.64, -21.4, -19.16, -16.92, -14.61, -12.25, -9.88,
//      -7.51, -5.14, -2.77, -0.43, 1.99, 4.35, 6.72, 9.09, 11.46, 13.83,
//      16.19, 18.43, 20.67, 22.85, 24.97, 27.09, 29.14, 31.13, 33.01, 34.88,
//      36.75, 38.5, 40.26, 41.9, 43.54, 45.17, 46.7, 48.23, 49.64, 51.06,
//      52.47, 53.78, 55.09, 56.4 
//
//nir-----------------------------------------------
//      -50.98, -41.74, -30.73, -17.97, -3.94, 3.16, 17.24, 30.01, 40.92,
//      50.16
//---------------------------------------------------------------------
     if(format.type == FT_HARP2 || format.type == FT_HARP){
  //      for (unsigned k = 0;k < num_frames;k++) {

         //BLUE GROUP
         for (unsigned i = 0;i < num_scans - 1;i++) { //improve with binary search
            if (asc_mode == 1 && latframe[4][i][nad_pix] < 0. && latframe[4][i + 1][nad_pix]>0.) {
                scan_brack[0] = i + 1;
                scan_brack[1] = i + 2;
                lat1 = latframe[4][i][nad_pix];
                lat2 = latframe[4][i + 1][nad_pix];
                lon1 = lonframe[4][i][nad_pix];
                lon2 = lonframe[4][i + 1][nad_pix];
//                timeix=k;
                timeix=i;//frames and scans are the same in this case
                i=num_scans;
//                k=num_frames;
//
               cout<<"lat1..."<<lat1<<"lat2..."<<lat2<<"asc_mode..."<<asc_mode<<endl;                
                break;
            }

            else if (asc_mode == 0 && latframe[4][i][nad_pix] > 0. && latframe[4][i + 1][nad_pix] < 0) { //descending orbit
                scan_brack[0] = i + 2;//negative lat first convention
                scan_brack[1] = i + 1;
                lat1 = latframe[4][i + 1][nad_pix];
                lat2 = latframe[4][i][nad_pix];
                lon1 = lonframe[4][i + 1][nad_pix];
                lon2 = lonframe[4][i][nad_pix];
                timeix=i;
                i=num_scans;
                break;
            }
           }
        // }//end for
         }
     else{
//OCIS/SPEXONE--------    
        for (unsigned i = 0;i < num_scans - 1;i++) { //improve with binary search
            if (asc_mode == 1 && lat[i][nad_pix] < 0. && lat[i + 1][nad_pix]>0.) {                
                scan_brack[0] = i + 1;
                scan_brack[1] = i + 2;
                lat1 = lat[i][nad_pix];
                lat2 = lat[i + 1][nad_pix];
                lon1 = lon[i][nad_pix];
                lon2 = lon[i + 1][nad_pix];
                timeix=i;
                i=num_scans;       
                break;            
            }

            else if (asc_mode == 0 && lat[i][nad_pix] > 0. && lat[i + 1][nad_pix] < 0) { //descending orbit
                scan_brack[0] = i + 2;//negative lat first convention
                scan_brack[1] = i + 1;
                lat1 = lat[i + 1][nad_pix];
                lat2 = lat[i][nad_pix];
                lon1 = lon[i + 1][nad_pix];
                lon2 = lon[i][nad_pix];
                timeix=i;
                i=num_scans;
                break;            
            }     
           } //end for
          }


        if (scan_brack[0] > -1) {//if file with equat crossing

          if(format.type == FT_HARP2 || format.type == FT_HARP){
              //blue group, view #5 closest to nadir!!
               t1 = stime2[4][timeix];
               t2 = stime2[4][timeix + 1];
               //interpolate time--linear--this can be improved using Harversine!!
               tcross = t1 - lat1 * (t2 - t1) / (lat2 - lat1);//equat crossing time
               float dtcross = (tcross - t1) / (t2 - t1);
               loncross = lon1 + (lon2 - lon1) * dtcross;//latcross is zero
               l1cfile->eqt = tcross;
               l1cfile->orbit_node_lon = loncross;
               cout<<"tcross..in ect_sf2.."<<tcross<<"loncross.."<<loncross<<"asc_mode..."<<asc_mode<<endl;
            }
            else{ //OCI/SPEXONE          
               t1 = stime[timeix];
               t2 = stime[timeix + 1];
               //interpolate time--linear--this can be improved using Harversine!!
               tcross = t1 - lat1 * (t2 - t1) / (lat2 - lat1);//equat crossing time
               float dtcross = (tcross - t1) / (t2 - t1);
               loncross = lon1 + (lon2 - lon1) * dtcross;//latcross is zero
               l1cfile->eqt = tcross;
               l1cfile->orbit_node_lon = loncross;           
               cout<<"tcross..in ect_sf2.."<<tcross<<"loncross.."<<loncross<<"asc_mode..."<<asc_mode<<endl;
            }}
        else {
            l1cfile->eqt = -999.0;
            l1cfile->orbit_node_lon = -999.0;    
            cout<<"tcross..in ect_sf2.."<<tcross<<"loncross.."<<loncross<<"asc_mode..."<<asc_mode<<endl;      
        }
        
        status = nc_close(ncid_L1B);
        check_err(status, __LINE__, __FILE__);

        if (time_str != nullptr)
           delete [](time_str);
        if (stime != nullptr)
           delete [](stime);
        if (stime2 != nullptr)
           delete [](stime2);
        if (lat != nullptr)
           delete [](lat);
        if (lon != nullptr)
           delete [](lon);
        if (latframe != nullptr)
           delete [](latframe);
        if (lonframe != nullptr)
           delete [](lonframe);


        return status;
    }



    //this version includes info coming from Don open file routines--

    int32_t L1C::load_l1c_filehandle4(l1c_filehandle* l1cfile, L1C_input* l1cinput) {
        int nfiles;
        string filename;
        nfiles = l1cinput->files.size();

        for (int j = 0;j < nfiles;j++) {       
            filename = l1cinput->files[j];
            l1cfile->ifiles.push_back(filename);        
        }
  
        l1cfile->l1c_pflag = l1cinput->l1c_pflag;
        l1cfile->swtnum = l1cinput->swtnum;
  

        for (int j = 0;j < 10;j++) {  //up to 10 granules
            l1cfile->selgran[j] = l1cinput->selgran[j];
        }
        //projection params
        l1cfile->gres = l1cinput->grid_resolution;
        l1cfile->proj_type = l1cinput->projection;
        l1cfile->cloud_corrected = l1cinput->cloud_correct;
        //multi  attributes (view, pol, bands) 
        l1cfile->overlap_vflag = l1cinput->overlap_vflag;
        l1cfile->overlap_pflag = l1cinput->overlap_pflag;
        l1cfile->overlap_bflag = l1cinput->overlap_bflag;
        //uncertainty params l1c merged products
        l1cfile->unc_meth = l1cinput->unc_meth;
        l1cfile->unc_thres_v = l1cinput->unc_thres_v;
        l1cfile->unc_thres_p = l1cinput->unc_thres_p;
        l1cfile->unc_thres_b = l1cinput->unc_thres_b;

        //calibration attributes
        //l1cfile->Fobar=file->Fobar;//nc

        cout << "ok transfering filehandle to l1c_filehandle info.." << endl;
        return 0;
    }


}//close l1c namespace 
