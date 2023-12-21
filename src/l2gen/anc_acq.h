#ifndef _ANC_ACQ_PROTO
#define _ANC_ACQ_PROTO

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/*  define calls for the anc_acq here  */
/* interpolation structure for one level, param, time */
typedef struct gen_int_str_def {
    gsl_spline2d *int_id;
    gsl_interp_accel *accel_lon, *accel_lat;
    double anc_time;
    double *lat_coord, *lon_coord;
    unsigned char *qual;
    int32_t nlon, nlat; /* sizes involved */
} gen_int_str;

int anc_acq_init(instr *, l1str *, int32_t *);
int32_t anc_acq_ck(char *, char *);
int32_t anc_acq_ecmwf_init(char **, char **, int, int32_t);
int anc_acq_lin_olci(int, char *, l1str *);
int anc_acq_lin(int32_t, l1str *);
int32_t anc_acq_f_stat( char **, char, int32_t );
float anc_miss_fill(int32_t);
float bilin_interp(float *, int, int, int, float, float);
int64_t jd4713bc_get_jd(int32_t, int32_t, int32_t);
int jd4713bc_get_date(int64_t, int32_t *, int32_t *, int32_t *);
int32_t anc_acq_read_gmao( char *, char *, char *, float **, unsigned char **,
  double *, int32_t *, int32_t *, int32_t *, double **, double ** );
int32_t anc_acq_lin_met( l1str * );
int32_t anc_acq_lin_prof( l1str * );
int32_t anc_acq_lin_aerosol( l1str * );
int32_t anc_acq_lin_oz( l1str * );
int32_t anc_acq_gmao_met_prep( char *, gen_int_str * );
int32_t anc_acq_gmao_prof_prep( char *, gen_int_str *, int32_t );
int32_t anc_acq_gmao_aer_prep(char *file, gen_int_str *aer_int);
int32_t anc_acq_gmao_oz_prep( char *, gen_int_str * );
int32_t anc_acq_fnd_t_interp( double, double *, int32_t, int32_t *, 
  int32_t *, float * );
int32_t anc_acq_eval_pt( gen_int_str *, int32_t, int32_t, float, float, 
  int32_t, int32_t *, float, int32_t, int32_t, int32_t, float *, float * );
int32_t init_anc_add( l1str * );
int32_t init_anc_aerosol(l1str *l1rec);
#endif
