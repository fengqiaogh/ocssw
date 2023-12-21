#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <netcdf>
#include <sstream>

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

  int32_t rta_nadir[2];
} geo_struct;

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )

#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )

int read_mce_tlm( netCDF::NcFile *l1afile, netCDF::NcGroup egid,
                  uint32_t nmcescan, uint32_t nenc, int32_t& ppr_off,
                  double& revpsec, double&secpline, int16_t& board_id,
                  int32_t *mspin, int32_t *ot_10us,
                  uint8_t *enc_count, float **hamenc, float **rtaenc);

int get_ev( double secpline, int16_t *dtype, int16_t *lines, int16_t *iagg,
            uint16_t& pcdim, uint16_t& psdim, double& ev_toff,
            float *clines, float *slines, double *deltc, double *delts,
            int16_t& iret);

int get_oci_vecs( uint32_t nscan, uint16_t pdim, int32_t *rta_nadir,
                  double ev_toff, int32_t spin, float *clines, double *delt,
                  double revpsec, int32_t ppr_off, int16_t board_id,
                  int32_t *mspin, int32_t *ot_10us, uint8_t *enc_count,
                  float **hamenc, float **rtaenc, float **pview,
                  double *theta, int16_t& iret);

int createField( netCDF::NcGroup &ncGrp, const char *sname, const char *lname, 
                   const char *standard_name, const char *units,
                   void *fill_value, const char *flag_values,
                   const char *flag_meanings,
                   double low, double high, int nt,
                   std::vector<netCDF::NcDim>& varVec);
