#include <fstream>
#include <iostream>
#include <sstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <string.h>
#include <string>
#include "hawkeye_methods.h"
//#include "hawkeyeUtil.h"
#include "netcdf.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"
#include <math.h>
#include <timeutils.h>

#define RADEG 57.29577951

using namespace std;

//this function ONLY WORKS WITH MONOTONIC INCREASING X VARIABLE!!!
int latlon_interp_dist(size_t num_gridlines, double *dist_ai, float *lati, float *loni,double *dist, float *lati2, float *loni2){
    std::vector<double> dscan;
    std::vector<double> latpix,lonpix;
  
      for(size_t i=0;i<num_gridlines-1;i++){
              dscan.push_back(dist_ai[i]);
              lati[i]+=90;
              latpix.push_back(lati[i]);
              if (loni[i]<0.0) loni[i]+=360;
              lonpix.push_back(loni[i]);
           }
      //interpolate pix time series--
      for(size_t k=0;k<num_gridlines-1;k++){
           tk::spline s1,s2;
           s1.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);
           s2.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);
         
           s1.set_points(dscan,latpix);
           s2.set_points(dscan,lonpix);
           lati2[k]=s1(dist[k]);
           loni2[k]=s2(dist[k]);

           if (loni2[k]>180.) loni2[k]-=360;
           lati2[k]-=90;
           }

     dscan.clear();
     latpix.clear();
     lonpix.clear();
     return(0);
}


int latlon_interp2(size_t num_gridlines,size_t num_pixels, double *gtime, float *lati, float *loni,double *gtime2, float *lati2, float *loni2){
    std::vector<double> tscan;
    std::vector<double> latpix,lonpix;

      for(size_t i=0;i<num_gridlines;i++){
              tscan.push_back(gtime[i]);
              lati[i]+=90;
              latpix.push_back(lati[i]);
              if (loni[i]<0.0) loni[i]+=360;
              lonpix.push_back(loni[i]);
           }
      //interpolate pix time series--
      for(size_t k=0;k<num_gridlines;k++){
           tk::spline s1,s2;
           s1.set_points(tscan,latpix);
           s2.set_points(tscan,lonpix);
           lati2[k]=s1(gtime2[k]);
           loni2[k]=s2(gtime2[k]);
           if (loni2[k]>180.) loni2[k]-=360;
           lati2[k]-=90;
           }

     tscan.clear();
     latpix.clear();
     lonpix.clear();
     return(0);
}



int latlon_interp1pix(size_t num_gridlines, size_t gd_index,double *tgrid, float *lat_nad, float *lon_nad, float *latipix, float *lonipix){
    std::vector<double> tscan;
    std::vector<double> latpix,lonpix;

      for(size_t i=0;i<num_gridlines;i++){
              tscan.push_back(tgrid[i]);
              lat_nad[i]+=90;
              latpix.push_back(lat_nad[i]);
              if (lon_nad[i]<0.0) lon_nad[i]+=360;
              lonpix.push_back(lon_nad[i]);
           }
    
           tk::spline s1,s2;
           s1.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);
           s2.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);

           s1.set_points(tscan,latpix);
           s2.set_points(tscan,lonpix);
                      
           *latipix=s1(tgrid[gd_index]);
           *lonipix=s2(tgrid[gd_index]);
           if (*lonipix>180.) *lonipix-=360;
           *latipix-=90;

     tscan.clear();
     latpix.clear();
     lonpix.clear();
  return 0;
  }



int latlon_interp_vec(size_t n_orb_rec, size_t num_gridlines,double *torb, double *latorb, double *lonorb,double *tmgv, float *lati, float *loni){
    double elon1,elon2,elat1,elat2,thres=5;
    std::vector<double> tscan,tlat,tlon;
    std::vector<double> latpix,lonpix,latclean,lonclean;
    std::vector<int> latidx,lonidx;
    double latemp=0.,lontemp=0.;

      for(size_t i=0;i<n_orb_rec;i++){
              tscan.push_back(torb[i]);
              latemp=latorb[i]+90;
              latpix.push_back(latemp);

              if (lonorb[i]<0.0){ 
                  lontemp=lonorb[i]+360;}
              else lontemp=lonorb[i];

              lonpix.push_back(lontemp);
           }

      for(size_t k=0;k<num_gridlines;k++){
           tk::spline s1,s2;
           s1.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);
           s2.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);

           s1.set_points(tscan,latpix);
           s2.set_points(tscan,lonpix);
           lati[k]=s1(tmgv[k]);
           loni[k]=s2(tmgv[k]);
           if (loni[k]>180.) loni[k]-=360;
           lati[k]-=90;


          if(k>=1 && k+1<num_gridlines){
            elon1=abs(loni[k]-loni[k-1]);
            elon2=abs(loni[k+1]-loni[k]);
            elat1=abs(lati[k]-lati[k-1]);
            elat2=abs(lati[k+1]-lati[k]);

              if(abs(100*elon2/elon1)<=100+thres && abs(100*elon2/elon1)>=100-thres){
                 tlon.push_back(tmgv[k]);
                 lonidx.push_back(k);
                 lonclean.push_back(loni[k]);}
              if(abs(100*elat2/elat1)<=100+thres && abs(100*elat2/elat1)>=100-thres){
                 tlat.push_back(tmgv[k]);
                 latidx.push_back(k);
                 latclean.push_back(lati[k]);}
             }

          }
 
     tscan.clear();
     tlat.clear();
     tlon.clear();
     latpix.clear();
     lonpix.clear();
     latidx.clear();
     lonidx.clear();
     latclean.clear();
     lonclean.clear();
      return 0;
  }





     


int latlon_interp(size_t n_orb_rec, size_t num_gridlines,size_t num_pixels, double *torb, float **lat, float **lon,double *time, float *lati, float *loni){
    double elon1,elon2,elat1,elat2,thres=5;

    std::vector<double> tscan,tlat,tlon;
    std::vector<double> latpix,lonpix,latclean,lonclean;
    std::vector<int> latidx,lonidx;

    int j=(num_pixels-1)/2;//nad pixel
    cout<<"interpolating based on nad pixels....index #.."<<j<<endl; 
      for(size_t i=0;i<n_orb_rec;i++){
              tscan.push_back(torb[i]);
              lat[i][j]+=90;
              latpix.push_back(lat[i][j]);
       
              if (lon[i][j]<0.0) lon[i][j]+=360;
              lonpix.push_back(lon[i][j]);
           }
        

      for(size_t k=0;k<num_gridlines;k++){
           tk::spline s1,s2;
           s1.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);
           s2.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);

           s1.set_points(tscan,latpix);
           s2.set_points(tscan,lonpix);
           lati[k]=s1(time[k]);
           loni[k]=s2(time[k]);
           if (loni[k]>180.) loni[k]-=360;
           lati[k]-=90;
          
          if(k>=1 && k+1<num_gridlines){
            elon1=abs(loni[k]-loni[k-1]);
            elon2=abs(loni[k+1]-loni[k]);
            elat1=abs(lati[k]-lati[k-1]);
            elat2=abs(lati[k+1]-lati[k]);

              if(abs(100*elon2/elon1)<=100+thres && abs(100*elon2/elon1)>=100-thres){
                 tlon.push_back(time[k]);
                 lonidx.push_back(k);
                 lonclean.push_back(loni[k]);}
              if(abs(100*elat2/elat1)<=100+thres && abs(100*elat2/elat1)>=100-thres){
                 tlat.push_back(time[k]);
                 latidx.push_back(k);
                 latclean.push_back(lati[k]);}

             }
            
          }
   
     tscan.clear();
     tlat.clear();
     tlon.clear();
     latpix.clear(); 
     lonpix.clear();
     latidx.clear();
     lonidx.clear();
     latclean.clear();
     lonclean.clear();

 return 0;
  }



int orb_interp2(size_t n_orb_rec, size_t sdim, double *torb, orb_array2 *p,
               orb_array2 *v, double *time, orb_array2 *posi, orb_array2 *veli) {
  double a0[3], a1[3], a2[3], a3[3];

  size_t ind = 0;
  for (size_t i = 0; i < sdim; i++) {
    //  Find input orbit vectors bracketing scan
    for (size_t j = ind; j < n_orb_rec; j++) {
      if (time[i] > torb[j] && time[i] <= torb[j + 1]) {
        ind = j;
        break;
      }
    }

    //  Set up cubic interpolation
    double dt = torb[ind + 1] - torb[ind];
    for (size_t j = 0; j < 3; j++) {
      a0[j] = p[ind][j];
      a1[j] = v[ind][j] * dt;
      a2[j] = 3 * p[ind + 1][j] - 3 * p[ind][j] - 2 * v[ind][j] * dt -
              v[ind + 1][j] * dt;
      a3[j] = 2 * p[ind][j] - 2 * p[ind + 1][j] + v[ind][j] * dt +
              v[ind + 1][j] * dt;
    }

    //  Interpolate orbit position and velocity components to scan line time
    double x = (time[i] - torb[ind]) / dt;
    double x2 = x * x;
    double x3 = x2 * x;
    for (size_t j = 0; j < 3; j++) {
      posi[i][j] = a0[j] + a1[j] * x + a2[j] * x2 + a3[j] * x3;
      veli[i][j] = (a1[j] + 2 * a2[j] * x + 3 * a3[j] * x2) / dt;
    }
  }  // i-loop

  return 0;
}



int orb_interp(size_t n_orb_rec, size_t sdim, double *torb, orb_array *p,
               orb_array *v, double *time, orb_array *posi, orb_array *veli) {
  double a0[3], a1[3], a2[3], a3[3];

  size_t ind = 0;
  for (size_t i = 0; i < sdim; i++) {
    //  Find input orbit vectors bracketing scan
    for (size_t j = ind; j < n_orb_rec; j++) {
      if (time[i] > torb[j] && time[i] <= torb[j + 1]) {
        ind = j;
        break;
      }
    }

    //  Set up cubic interpolation
    double dt = torb[ind + 1] - torb[ind];
    for (size_t j = 0; j < 3; j++) {
      a0[j] = p[ind][j];
      a1[j] = v[ind][j] * dt;
      a2[j] = 3 * p[ind + 1][j] - 3 * p[ind][j] - 2 * v[ind][j] * dt -
              v[ind + 1][j] * dt;
      a3[j] = 2 * p[ind][j] - 2 * p[ind + 1][j] + v[ind][j] * dt +
              v[ind + 1][j] * dt;
    }

    //  Interpolate orbit position and velocity components to scan line time
    double x = (time[i] - torb[ind]) / dt;
    double x2 = x * x;
    double x3 = x2 * x;
    for (size_t j = 0; j < 3; j++) {
      posi[i][j] = a0[j] + a1[j] * x + a2[j] * x2 + a3[j] * x3;
      veli[i][j] = (a1[j] + 2 * a2[j] * x + 3 * a3[j] * x2) / dt;
    }
  }  // i-loop

  return 0;
}


int j2000_to_ecr(int32_t iyr, int32_t idy, double sec, double ecmat[3][3]) {
  // Get J2000 to ECEF transformation matrix

  // Arguments:

  // Name               Type    I/O     Description
  // --------------------------------------------------------
  // iyr         I       I      Year, four digits
  // idy         I       I      Day of year
  // sec        R*8      I      Seconds of day
  // ecmat(3,3)  R       O      J2000 to ECEF matrix

  // Get transformation from J2000 to mean-of-date inertial
  double j2mod[3][3];
  j2000_to_mod(iyr, idy, sec, j2mod);

  // Get nutation and UT1-UTC (once per run)
  double xnut[3][3], ut1utc;
  get_nut(iyr, idy, xnut);
  get_ut1(iyr, idy, ut1utc);

  // Compute Greenwich hour angle for time of day
  double day = idy + (sec + ut1utc) / 86400;
  double gha, gham[3][3];
  gha2000(iyr, day, gha);

  gham[0][0] = cos(gha / RADEG);
  gham[1][1] = cos(gha / RADEG);
  gham[2][2] = 1;
  gham[0][1] = sin(gha / RADEG);
  gham[1][0] = -sin(gha / RADEG);

  gham[0][2] = 0;
  gham[2][0] = 0;
  gham[1][2] = 0;
  gham[2][1] = 0;

  // Combine all transformations
  gsl_matrix_view A = gsl_matrix_view_array(&gham[0][0], 3, 3);  // gham
  gsl_matrix_view B = gsl_matrix_view_array(&xnut[0][0], 3, 3);  // xnut
  gsl_matrix *C = gsl_matrix_alloc(3, 3);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &A.matrix, &B.matrix, 0.0, C);

  gsl_matrix_view D = gsl_matrix_view_array(&j2mod[0][0], 3, 3);  // j2mod
  gsl_matrix *E = gsl_matrix_alloc(3, 3);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C, &D.matrix, 0.0, E);
  double *ptr_E = gsl_matrix_ptr(E, 0, 0);

  memcpy(ecmat, ptr_E, 9 * sizeof(double));

  gsl_matrix_free(C);
    gsl_matrix_free(E);

  return 0;
}

int j2000_to_mod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]) {
  // Get J2000 to MOD (precession) transformation

  // Arguments:

  // Name               Type    I/O     Description
  // --------------------------------------------------------
  // iyr         I       I      Year, four digits
  // idy         I       I      Day of year
  // sec        R*8      I      Seconds of day
  // j2mod(3,3)  R       O      J2000 to MOD matrix

  int16_t iyear = iyr;
  int16_t iday = idy;

  double t = jday(iyear, 1, iday) - (double)2451545.5 + sec / 86400;
  t /= 36525;

  double zeta0 = t * (2306.2181 + 0.302 * t + 0.018 * t * t) / RADEG / 3600;
  double thetap = t * (2004.3109 - 0.4266 * t - 0.04160 * t * t) / RADEG / 3600;
  double xip = t * (2306.2181 + 1.095 * t + 0.018 * t * t) / RADEG / 3600;

  j2mod[0][0] = -sin(zeta0) * sin(xip) + cos(zeta0) * cos(xip) * cos(thetap);
  j2mod[0][1] = -cos(zeta0) * sin(xip) - sin(zeta0) * cos(xip) * cos(thetap);
  j2mod[0][2] = -cos(xip) * sin(thetap);
  j2mod[1][0] = sin(zeta0) * cos(xip) + cos(zeta0) * sin(xip) * cos(thetap);
  j2mod[1][1] = cos(zeta0) * cos(xip) - sin(zeta0) * sin(xip) * cos(thetap);
  j2mod[1][2] = -sin(xip) * sin(thetap);
  j2mod[2][0] = cos(zeta0) * sin(thetap);
  j2mod[2][1] = -sin(zeta0) * sin(thetap);
  j2mod[2][2] = cos(thetap);

  return 0;
}

int get_nut(int32_t iyr, int32_t idy, double xnut[3][3]) {
  int16_t iyear = iyr;
  int16_t iday = idy;

  double t = jday(iyear, 1, iday) - (double)2451545.5;

  double xls, gs, xlm, omega;
  double dpsi, eps, epsm;
  ephparms(t, xls, gs, xlm, omega);
  nutate(t, xls, gs, xlm, omega, dpsi, eps, epsm);

  xnut[0][0] = cos(dpsi / RADEG);
  xnut[1][0] = -sin(dpsi / RADEG) * cos(epsm / RADEG);
  xnut[2][0] = -sin(dpsi / RADEG) * sin(epsm / RADEG);
  xnut[0][1] = sin(dpsi / RADEG) * cos(eps / RADEG);
  xnut[1][1] = cos(dpsi / RADEG) * cos(eps / RADEG) * cos(epsm / RADEG) +
               sin(eps / RADEG) * sin(epsm / RADEG);
  xnut[2][1] = cos(dpsi / RADEG) * cos(eps / RADEG) * sin(epsm / RADEG) -
               sin(eps / RADEG) * cos(epsm / RADEG);
  xnut[0][2] = sin(dpsi / RADEG) * sin(eps / RADEG);
  xnut[1][2] = cos(dpsi / RADEG) * sin(eps / RADEG) * cos(epsm / RADEG) -
               cos(eps / RADEG) * sin(epsm / RADEG);
  xnut[2][2] = cos(dpsi / RADEG) * sin(eps / RADEG) * sin(epsm / RADEG) +
               cos(eps / RADEG) * cos(epsm / RADEG);

  return 0;
}

 int ephparms(double t, double &xls, double &gs, double &xlm, double &omega) {
  //  This subroutine computes ephemeris parameters used by other Mission
  //  Operations routines:  the solar mean longitude and mean anomaly, and
  //  the lunar mean longitude and mean ascending node.  It uses the model
  //  referenced in The Astronomical Almanac for 1984, Section S
  //  (Supplement) and documented for the SeaWiFS Project in "Constants
  //  and Parameters for SeaWiFS Mission Operations", in TBD.  These
  //  parameters are used to compute the solar longitude and the nutation
  //  in longitude and obliquity.

  //  Sun Mean Longitude
  xls = (double)280.46592 + ((double)0.9856473516) * t;
  xls = fmod(xls, (double)360);

  //  Sun Mean Anomaly
  gs = (double)357.52772 + ((double)0.9856002831) * t;
  gs = fmod(gs, (double)360);

  //  Moon Mean Longitude
  xlm = (double)218.31643 + ((double)13.17639648) * t;
  xlm = fmod(xlm, (double)360);

  //  Ascending Node of Moon's Mean Orbit
  omega = (double)125.04452 - ((double)0.0529537648) * t;
  omega = fmod(omega, (double)360);

  return 0;
}

int nutate(double t, double xls, double gs, double xlm, double omega,
           double &dpsi, double &eps, double &epsm) {
  //  This subroutine computes the nutation in longitude and the obliquity
  //  of the ecliptic corrected for nutation.  It uses the model referenced
  //  in The Astronomical Almanac for 1984, Section S (Supplement) and
  //  documented for the SeaWiFS Project in "Constants and Parameters for
  //  SeaWiFS Mission Operations", in TBD.  These parameters are used to
  //  compute the apparent time correction to the Greenwich Hour Angle and
  //  for the calculation of the geocentric Sun vector.  The input
  //  ephemeris parameters are computed using subroutine ephparms.  Terms
  //  are included to 0.1 arcsecond.

  //  Nutation in Longitude
  dpsi = -17.1996 * sin(omega / RADEG) + 0.2062 * sin(2. * omega / RADEG) -
         1.3187 * sin(2. * xls / RADEG) + 0.1426 * sin(gs / RADEG) -
         0.2274 * sin(2. * xlm / RADEG);

  //  Mean Obliquity of the Ecliptic
  epsm = (double)23.439291 - ((double)3.560e-7) * t;

  //  Nutation in Obliquity
  double deps = 9.2025 * cos(omega / RADEG) + 0.5736 * cos(2. * xls / RADEG);

  //  True Obliquity of the Ecliptic
  eps = epsm + deps / 3600;

  dpsi = dpsi / 3600;

  return 0;
}

int get_ut1(int32_t iyr, int32_t idy, double &ut1utc) {
  int16_t iyear = iyr;
  int16_t iday = idy;

  static int32_t ijd[25000];
  static double ut1[25000];
  string utcpole = "$OCVARROOT/var/modis/utcpole.dat";
  static bool first = true;
  int k = 0;
  if (first) {
    string line;
    expandEnvVar(&utcpole);
    istringstream istr;

    ifstream utcpfile(utcpole.c_str());
    if (utcpfile.is_open()) {
      getline(utcpfile, line);
      getline(utcpfile, line);
      while (getline(utcpfile, line)) {

        istr.clear();
        istr.str(line.substr(0, 5));
        istr >> ijd[k];
        istr.clear();
        istr.str(line.substr(44, 9));
        istr >> ut1[k];
        k++;
      }
      ijd[k] = 0;
      utcpfile.close();
      first = false;
    } else {
      cout << utcpole.c_str() << " not found" << endl;
      exit(1);
    }
  }  // if (first)

  k = 0;
  int mjd = jday(iyear, 1, iday) - 2400000;
  while (ijd[k] > 0) {
    if (mjd == ijd[k]) {
      ut1utc = ut1[k];
      break;
    }
    k++;
  }

  return 0;
}

int gha2000(int32_t iyr, double day, double &gha) {
  //  This subroutine computes the Greenwich hour angle in degrees for the
  //  input time.  It uses the model referenced in The Astronomical Almanac
  //  for 1984, Section S (Supplement) and documented for the SeaWiFS
  //  Project in "Constants and Parameters for SeaWiFS Mission Operations",
  //  in TBD.  It includes the correction to mean sidereal time for nutation
  //  as well as precession.
  //

  //  Compute days since J2000
  int16_t iday = day;
  double fday = day - iday;
  int jd = jday(iyr, 1, iday);
  double t = jd - (double)2451545.5 + fday;

  //  Compute Greenwich Mean Sidereal Time      (degrees)
  double gmst = (double)100.4606184 + (double)0.9856473663 * t +
                (double)2.908e-13 * t * t;

  //  Check if need to compute nutation correction for this day
  double xls, gs, xlm, omega;
  double dpsi, eps, epsm;
  ephparms(t, xls, gs, xlm, omega);
  nutate(t, xls, gs, xlm, omega, dpsi, eps, epsm);

  //  Include apparent time correction and time-of-day
  gha = gmst + dpsi * cos(eps / RADEG) + fday * 360;
  gha = fmod(gha, (double)360);

  return 0;
}


int expandEnvVar(std::string *sValue) {
  if ( (*sValue).find_first_of( "$" ) == string::npos) return 0;
  string::size_type posEndIdx = (*sValue).find_first_of( "/" );
  if ( posEndIdx == string::npos) return 0;
  char *envVar_str = getenv((*sValue).substr( 1, posEndIdx-1 ).c_str());
  if (envVar_str == 0x0) {
    printf("Environment variable: %s not defined.\n", envVar_str);
    exit(1);
  }
  *sValue = envVar_str + (*sValue).substr( posEndIdx);

  return 0;
}


























































