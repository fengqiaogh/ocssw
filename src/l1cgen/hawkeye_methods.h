#ifndef _HAWKEYE_METHODS_H_
#define _HAWKEYE_METHODS_H_

#include <string.h>
#include <string>

typedef float quat_array[4];
typedef float orb_array[3];
typedef double quat_array2[4];
typedef double orb_array2[3];

int orb_interp(size_t n_SC_rec, size_t sdim,
               double *torb, orb_array *p, orb_array *v,
               double *time, orb_array *posi, orb_array *veli);
int orb_interp2(size_t n_SC_rec, size_t sdim,
               double *torb, orb_array2 *p, orb_array2 *v,
               double *time, orb_array2 *posi, orb_array2 *veli);


int latlon_interp(size_t n_SC_rec,size_t num_gridlines,size_t num_pixels,
               double *torb, float **lat, float **lon,
               double *gtime, float *lati, float *loni);
int latlon_interp2(size_t num_gridlines,size_t num_pixels,
               double *gtime, float *lati, float *loni,
               double *gtime2, float *lati2, float *loni2);
int latlon_interp_vec(size_t n_orb_rec, size_t num_gridlines,double *torb, double *latorb, double *lonorb,double *tmgv, float *lati, float *loni);

int latlon_interp_dist(size_t num_gridlines,double *dist_ai,float *lati,float *loni,double *dist_a,float *lati2,float *loni2);
int latlon_interp1pix(size_t num_gridlines, size_t gd_index,double *tgrid, float *lat_nad, float *lon_nad, float *latpix, float *lonpix);

int j2000_to_ecr(int32_t iyr, int32_t idy, double sec, double ecmat[3][3]);
int j2000_to_mod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]);
int get_nut(int32_t iyr, int32_t idy, double xnut[3][3]);
int get_ut1(int32_t iyr, int32_t idy, double &ut1utc);
int ephparms(double t, double &xls, double &gs, double &xlm, double &omega);
int nutate(double t, double xls, double gs, double xlm, double omega,
           double &dpsi, double &eps, double &epsm);
int gha2000(int32_t iyr, double day, double &gha);


int expandEnvVar(std::string *sValue);


#endif  // _GEOLOCATE_HAWKEYE_H_
