#include <stdint.h>
#include <fstream>
#include <timeutils.h>
#include <netcdf>
#include <vector>

#include "common.h"

#define NBWAVE 512
#define NRWAVE 512
#define NIWAVE 9
#define NTIMES 2
#define NTEMPS 30
/*
typedef struct {
  double master_clock;
  double MCE_clock;

  double sc_to_tilt[3][3];
  double tilt_axis[3];
  double tilt_angles[2];
  double tilt_to_oci_mech[3][3];
  double oci_mech_to_oci_opt[3][3];
  double rta_axis[3];
  double ham_axis[3];
  double ham_at_angles[2];
  double ham_ct_angles[2];
  double rta_enc_scale;
  double ham_enc_scale;

  int32_t rta_nadir;
} geo_struct;
*/
typedef struct {
  uint16_t ldims[6];
  float **K1;
  float ***K2;
  float ***K3_coef;
  float ***K4_coef;
  double **K5;
} cal_lut_struct;

typedef struct {
  uint16_t ldims[6];
  float **K1K2;
  float ***K3_coef;
  float ***K4_coef;
  double **K5;
} gains_struct;


#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )

#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )


class l1bFile {
  std::string fileName;
  
  int ngrps;
  int ndims;
  
  netCDF::NcDim ncDims[1000];
  
 public:
  l1bFile();
  ~l1bFile();

  netCDF::NcFile *l1bfile;

  netCDF::NcGroup ncGrps[10];

//  std::string platform;
// int apktsize;
// int bpktsize;
// int EV_APIDs;
/*
  int createl1b( char* l1b_filename, uint16_t nscan_good,
                   uint16_t pcdim, uint16_t bbb, uint16_t rbb,
                   uint16_t psdim, uint16_t swb);
*/

  int parseDims( std::string dimString, std::vector<netCDF::NcDim>& varDims);

  int write_oci_science_data( uint32_t isc,
                              uint16_t nbbs, uint16_t nrbs, uint16_t nswb,
                              uint16_t ncps, uint16_t nsps,
                              uint16_t **bsci, uint16_t **rsci,
                              uint32_t **ssci, int8_t *sfrms);


  int write_granule_metadata( std::string tstart, std::string tend,
                              std::string l1b_name);
  
  int close();
};

int get_agg_mat( size_t *ia, int16_t iagg, int16_t jagg[16], uint16_t nib,
                 uint32_t& nbb, uint32_t ntb[16], float **amat, float **gmat);

int read_oci_cal_lut( netCDF::NcFile *calLUTfile, std::string tag,
                      netCDF::NcGroup gidLUT, uint32_t& banddim,
                      cal_lut_struct& cal_lut);

int make_oci_gains( uint32_t nib, uint32_t banddim, uint16_t iyr, uint16_t idom,
                    double stime, double K2t[NTIMES], cal_lut_struct& cal_lut,
                    float **gmat, gains_struct& gains);


int createField( netCDF::NcGroup &ncGrp, const char *sname, const char *lname, 
                   const char *standard_name, const char *units,
                   void *fill_value, const char *flag_values,
                   const char *flag_meanings,
                   double low, double high, int nt,
                   std::vector<netCDF::NcDim>& varVec);

template <typename T>
int get_oci_dark( size_t iscn, uint32_t nscan, uint8_t *hside, uint16_t ndsc,
                  uint16_t nskp, int16_t iags, int16_t iagd, uint32_t ntaps,
                  int16_t *jagg, uint32_t dfill, int16_t ndc, T ***dark,
                  uint32_t nib, float *dc, int16_t& iret);

int get_oci_temp_corr( uint32_t nib, gains_struct gains, float K3T[NTEMPS],
                       float *caltemps, uint32_t nscan, float *k3);

int get_oci_rvs_corr( uint32_t nib, uint16_t pdim, uint8_t hside,
                      gains_struct gains, double *theta, float **k4);

int get_oci_lin_corr( uint32_t nib, uint16_t pdim, gains_struct gains,
                      float K3T[NTEMPS], float *caltemps, float **dn,
                      float **k5);

template <typename T> T*** make3dT( size_t dims[3]);


inline
int expandEnvVar( std::string *sValue) {
  if ( (*sValue).find_first_of( "$" ) == std::string::npos) return 0;
  std::string::size_type posEndIdx = (*sValue).find_first_of( "/" );
  if ( posEndIdx == std::string::npos) return 0;
  const std::string envVar = sValue->substr (1, posEndIdx - 1);
  char *envVar_str = getenv(envVar.c_str());
  if (envVar_str == 0x0) {
    printf("Environment variable: %s not defined.\n", sValue->c_str());
    exit(1);
  }
  *sValue = envVar_str + (*sValue).substr( posEndIdx);

  return 0;
}

l1bFile::l1bFile() {

}


l1bFile::~l1bFile() {
  
}

