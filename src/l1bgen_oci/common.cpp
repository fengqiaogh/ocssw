#include "common.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


int read_mce_tlm( NcFile *l1afile, NcGroup egid, uint32_t nmcescan,
                  uint32_t nenc, int32_t& ppr_off,
                  double& revpsec, double&secpline,
                  int16_t& board_id, int32_t *mspin, int32_t *ot_10us,
                  uint8_t *enc_count, float **hamenc, float **rtaenc) {

  NcDim mce_dim = l1afile->getDim("MCE_block");
  uint32_t mce_blk = mce_dim.getSize();
  NcDim ddc_dim = l1afile->getDim("DDC_tlm");
  uint32_t nddc = ddc_dim.getSize();
  NcDim tlm_dim = l1afile->getDim("tlm_packets");
  uint32_t ntlmpack = tlm_dim.getSize();
  
  uint8_t **mtlm = new uint8_t *[nmcescan];
  mtlm[0] = new uint8_t[mce_blk*nmcescan];
  for (size_t i=1; i<nmcescan; i++) mtlm[i] = mtlm[i-1] + mce_blk;

  int16_t **enc = new int16_t *[nmcescan];
  enc[0] = new int16_t[4*nenc*nmcescan];
  for (size_t i=1; i<nmcescan; i++) enc[i] = enc[i-1] + 4*nenc;
  
  uint8_t **ddctlm = new uint8_t *[ntlmpack];
  ddctlm[0] = new uint8_t[nddc*ntlmpack];
  for (size_t i=1; i<ntlmpack; i++) ddctlm[i] = ddctlm[i-1] + nddc;

  NcVar var;
  vector<size_t> start, count;
  start.push_back(0);
  start.push_back(0);
  count.push_back(1);

  var = egid.getVar( "MCE_telemetry");
  var.getVar( &mtlm[0][0]);
  
  var = egid.getVar( "MCE_encoder_data");
  count.pop_back();
  count.push_back(4*nenc);
  var.getVar( &enc[0][0]);

  var = egid.getVar( "MCE_spin_ID");
  var.getVar( mspin);

  var = egid.getVar( "DDC_telemetry");
  var.getVar( &ddctlm[0][0]);

  int32_t max_enc_cts = 131072; // 2^17
  double clock[2] = {1.36e8,1.0e8};

  // Get ref_pulse_divider and compute commanded rotation rate
  uint32_t ui32;
  uint32_t ref_pulse_div[2];
  memcpy( &ui32, &mtlm[0][0], 4);
  ref_pulse_div[0] = SWAP_4( ui32) % 16777216; // 2^24
  memcpy( &ui32, &mtlm[0][4], 4);
  ref_pulse_div[1] = SWAP_4( ui32) % 16777216; // 2^24

  int32_t ref_pulse_sel = mtlm[0][9] / 128;
  
  revpsec = clock[ref_pulse_sel] / 2 / max_enc_cts /
    (ref_pulse_div[ref_pulse_sel]/256.0 + 1);

  // Get PPR offset and on-time_10us
  memcpy( &ui32, &mtlm[0][8], 4);
  ppr_off = SWAP_4( ui32) % max_enc_cts;
  for (size_t i=0; i<nmcescan; i++) {
    memcpy( &ui32, &mtlm[i][368], 4);
    ot_10us[i] = SWAP_4( ui32) % 4;
  }

  // Get MCE board ID  
  board_id = (int16_t) mtlm[0][322] / 16;

  // Get TDI time and compute time increment per line
  uint16_t ui16;
  memcpy( &ui16, &ddctlm[0][346], 2);
  int32_t tditime = SWAP_2( ui16);
  secpline = (tditime+1) / clock[0]; 

  // Get valid encoder count, HAM and RTA encoder data
  for (size_t i=0; i<nmcescan; i++) enc_count[i] = mtlm[i][475];
  float enc_s = 81.0 / 2560;
  for (size_t i=0; i<nmcescan; i++) {
    for (size_t j=0; j<nenc; j++) {
      hamenc[i][j] = enc[i][4*j+0] * enc_s;
      rtaenc[i][j] = enc[i][4*j+1] * enc_s;
    }
  }

  delete[] mtlm[0];
  delete[] mtlm;

  delete [] enc[0];
  delete [] enc;

  delete [] ddctlm[0];
  delete [] ddctlm;

  return 0;
}


int get_ev( double secpline, int16_t *dtype, int16_t *lines, int16_t *iagg,
            uint16_t& pcdim, uint16_t& psdim, double& ev_toff,
            float *clines, float *slines, double *deltc, double *delts,
            int16_t &iret) {

  // Find end of no-data zone
  int16_t iz=0, line0=0;
  iret = -1;
  while ( dtype[iz] == 0) {
    line0 += lines[iz];
    iz++;
  }
  if (iz == 10) return 0;
  
  // Find number of pixels in Earth views
  pcdim = 0;
  psdim = 0;
  int16_t linen = line0;
  for (size_t i=iz; i<9; i++) {
    // Check for not dark or no-data
    if ( dtype[i] != 0 && dtype[i] != 2 && dtype[i] < 10) {

      //if ( dtype[i] == 1 || dtype[i] == 3 || dtype[i] == 4 || dtype[i] == 6 ||
      // dtype[i] == 7 || dtype[i] == 9) {

      uint16_t np = lines[i] / iagg[i];
      for (size_t j=0; j<np; j++) {
        clines[pcdim+j] = linen + j*iagg[i] + 0.5*iagg[i] - 64;
      }
      pcdim += np;
      uint16_t ns = lines[i] / 8;
      for (size_t j=0; j<ns; j++) {
        slines[psdim+j] = linen + j*8 + 3.5;
      }
      psdim += ns;
      iret = 0;
    }
    linen += lines[i];
  }

  // Calculate times
  for (size_t i=0; i<(size_t) pcdim; i++) deltc[i] = secpline * clines[i];
  ev_toff = 0.5 * (deltc[0] + deltc[pcdim-1]);
  for (size_t i=0; i<(size_t) pcdim; i++) deltc[i] -= ev_toff;
    
  for (size_t i=0; i<(size_t) psdim; i++)
    delts[i] = secpline * slines[i] - ev_toff;

  return 0;
}

int get_oci_vecs( uint32_t nscan, uint16_t pdim, int32_t *rta_nadir,
                  double ev_toff, int32_t spin, float *clines, double *delt,
                  double revpsec, int32_t ppr_off, int16_t board_id,
                  int32_t *mspin, int32_t *ot_10us, uint8_t *enc_count,
                  float **hamenc, float **rtaenc, float **pview,
                  double *theta, int16_t& iret) {

  // This program generates the OCI Earth view vectors for one spin.  
  // It uses MCE telemetry and encoder data.  Further refinements will be made
  // as the instrument optics model and test results become available.  
  // Reference: "Use of OCI Telemetry to Determine Pixel Line-of-Sight",
  // F. Patt, 2020-05-18

  int32_t max_enc_cts = 131072; // 2^17
  double dtenc = 0.001;
  constexpr double pi = 3.14159265358979323846;
  int16_t bd_id = board_id % 2;
  
  double rad2asec = (180/pi) * 3600;

  // Compute scan angle corresponding to PPR
  //  float pprang = 2 * pi * (ppr_off - geoLUT.rta_nadir[bd_id]) / max_enc_cts;
  float pprang = 2 * pi * (ppr_off - rta_nadir[bd_id]) / max_enc_cts;
  if (pprang > pi) pprang -= 2*pi;
  
  // Compute ideal scan angles for science pixels
  double *toff = new double[pdim];
  for (size_t i=0; i<pdim; i++) toff[i] = delt[i] + ev_toff;

  for (size_t i=0; i<pdim; i++)
    theta[i] = pprang + 2 * pi * revpsec * toff[i];
  // Interpolate encoder data to pixel times and add to scan angles
  // RTA only for now, include HAM when we know how.

  double *thetacor = new double[pdim];

  int isp = -1;
  for (size_t i=0; i<nscan; i++) {
    if ( mspin[i] == spin) {
      isp = (int) i;
      break;
    }
  }

  // Interpolate encoder data to pixel times and add to scan angles
  // RTA only for now, include HAM when we know how.
  if ( isp == -1) {
    cout << "No MCE encoder data for spin: " << spin << endl;
    iret = 1;
  } else {
    size_t ip = 0, ke;
    double tenc_ke;
    while ( ip < pdim) {
      // Encoder sample times at 1 KHz
      for (size_t j=0; j<enc_count[isp]; j++) {
        double tenc = j*dtenc;
        if ( tenc < toff[ip] && (tenc+dtenc) > toff[ip]) {
          ke = j;
          tenc_ke = tenc;
          break;
        }
      }
      size_t njp=0;
      for (size_t i=0; i<pdim; i++) {
        if ( toff[i] >= tenc_ke && toff[i] < tenc_ke+dtenc) {
          double ft = (toff[i] - tenc_ke) / dtenc;
          thetacor[i] =
            (1-ft)*rtaenc[isp][ke] + ft*rtaenc[isp][ke+1] -
            ((1-ft)*hamenc[isp][ke] + ft*hamenc[isp][ke+1])*0.25;
          njp++;
        }
      }
      ip += njp;
    }
  }
  
  // Simple, planar view vector model to start
  // Update when optical model is available.
  for (size_t i=0; i<pdim; i++) {
    theta[i] = theta[i] - thetacor[i] / rad2asec;
    pview[i][0] = 0.0;
    pview[i][1] = sin(theta[i]);
    pview[i][2] = cos(theta[i]);
  }

  delete [] thetacor;
  delete [] toff;

  return 0;
}

int createField( NcGroup &ncGrp, const char *sname, const char *lname, 
                 const char *standard_name, const char *units,
                 void *fill_value, const char *flag_values,
                 const char *flag_meanings,
                 double low, double high, int nt, vector<NcDim>& varVec) {

  size_t dimlength;


  /* Create the NCDF dataset */
  NcVar ncVar;
  try {
    ncVar = ncGrp.addVar(sname, nt, varVec);   
  }
  catch ( NcException& e) {
    cout << e.what() << endl;
    cerr << "Failure creating variable: " << sname << endl;
    exit(1);
  }

  // Set fill value
  double fill_value_dbl;
  memcpy( &fill_value_dbl, fill_value, sizeof(double));

  int8_t i8;
  uint8_t ui8;
  int16_t i16;
  uint16_t ui16;
  int32_t i32;
  uint32_t ui32;
  float f32;

  if ( low != fill_value_dbl) {
    if ( nt == NC_BYTE) {
      i8 = fill_value_dbl;
      ncVar.setFill(true, (void *) &i8);
    } else if ( nt == NC_UBYTE) {
      ui8 = fill_value_dbl;
      ncVar.setFill(true, (void *) &ui8);
    } else if ( nt == NC_SHORT) {
      i16 = fill_value_dbl;
      ncVar.setFill(true, (void *) &i16);
    } else if ( nt == NC_USHORT) {
      ui16 = fill_value_dbl;
      ncVar.setFill(true, (void *) &ui16);
    } else if ( nt == NC_INT) {
      i32 = fill_value_dbl;
      ncVar.setFill(true, (void *) &i32);
    } else if ( nt == NC_UINT) {
      ui32 = fill_value_dbl;
      ncVar.setFill(true, (void *) &ui32);
    } else if ( nt == NC_FLOAT) {
      f32 = fill_value_dbl;
      ncVar.setFill(true, (void *) &f32);
    } else {
      ncVar.setFill(true, (void *) &fill_value_dbl);
    }
  }

  /* vary chunck size based on dimensions */ 
  int do_deflate = 0;
  vector<size_t> chunkVec;
  if ( varVec.size() == 3 && (strncmp(sname, "EV_", 3) == 0)) {
    dimlength = varVec[2].getSize();

    chunkVec.push_back(1);
    chunkVec.push_back(16);
    chunkVec.push_back(dimlength/10);

    do_deflate = 1;
  }

  /* Set compression */
  if ( do_deflate) {
    /* First set chunking */
    try {
      ncVar.setChunking(ncVar.nc_CHUNKED, chunkVec);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure setting chunking: " << sname << endl;
      exit(1);
    }

    try {
      ncVar.setCompression(true, true, 5);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure setting compression: " << sname << endl;
      exit(1);
    }
  }


  /* Add a "long_name" attribute */
  try {
    ncVar.putAtt("long_name", lname);
  }
  catch ( NcException& e) {
    e.what();
    cerr << "Failure creating 'long_name' attribute: " << lname << endl;
    exit(1);
  }

  if ( strcmp( flag_values, "") != 0) {

    size_t curPos=0;

    string fv;
    fv.assign( flag_values);
    size_t pos = fv.find("=", curPos);
    fv = fv.substr(pos+1);

    size_t semicln = fv.find(";");
    pos = 0;

    int8_t vec[1024];
    int n = 0;
    while(pos != semicln) {
      pos = fv.find(",", curPos);
      if ( pos == string::npos) 
        pos = semicln;

      string flag_value;
      istringstream iss(fv.substr(curPos, pos-curPos));
      iss >> skipws >> flag_value;
      vec[n++] = atoi( flag_value.c_str());
      curPos = pos + 1;
    }

    try {
      ncVar.putAtt("flag_values", nt, n, vec);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'flag_values' attribute: " << lname << endl;
      exit(1);
    }
  }

  /* Add a "flag_meanings" attribute if specified*/
  if ( strcmp( flag_meanings, "") != 0) {

    try {
      ncVar.putAtt("flag_meanings", flag_meanings);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'flag_meanings' attribute: "
           << flag_meanings << endl;
      exit(1);
    }
  }

  /* Add "valid_min/max" attributes if specified */
  if (low < high) {
    switch(nt) {              /* Use the appropriate number type */
    case NC_BYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;

        try {
          ncVar.putAtt("valid_min", NC_BYTE, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_BYTE, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    case NC_UBYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;

        try {
          ncVar.putAtt("valid_min", NC_UBYTE, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_UBYTE, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    case NC_SHORT:
      {
	int16_t vr[2];
	vr[0] = (int16_t)low;
	vr[1] = (int16_t)high;

        try {
          ncVar.putAtt("valid_min", NC_SHORT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_SHORT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    case NC_USHORT:
      {
	uint16_t vr[2];
	vr[0] = (uint16_t)low;
	vr[1] = (uint16_t)high;

        try {
          ncVar.putAtt("valid_min", NC_USHORT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_USHORT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
              }
      break;
    case NC_INT:
      {
	int32_t vr[2];
	vr[0] = (int32_t)low;
	vr[1] = (int32_t)high;

        try {
          ncVar.putAtt("valid_min", NC_INT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_INT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }

      }
      break;
    case NC_UINT:
      {
	uint32_t vr[2];
	vr[0] = (uint32_t)low;
	vr[1] = (uint32_t)high;

        try {
          ncVar.putAtt("valid_min", NC_UINT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_UINT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
        
      }
      break;
    case NC_FLOAT:
      {
	float vr[2];
	vr[0] = (float)low;
	vr[1] = (float)high;

        try {
          ncVar.putAtt("valid_min", NC_FLOAT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_FLOAT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }        
      }
      break;
    case NC_DOUBLE:
      {
	double vr[2];
	vr[0] = low;
	vr[1] = high;

        try {
          ncVar.putAtt("valid_min", NC_DOUBLE, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_DOUBLE, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    default:
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"Got unsupported number type (%d) ",nt);
      fprintf(stderr,"while trying to create NCDF variable, \"%s\", ",sname);
      return(1);
    }
  }           
    
  /* Add a "units" attribute if one is specified */
  if(units != NULL && *units != 0) {

    try {
      ncVar.putAtt("units", units);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'units' attribute: " << units << endl;
      exit(1);
    }
  }

  /* Add a "standard_name" attribute if one is specified */
  if(standard_name != NULL && *standard_name != 0) {
    try {
      ncVar.putAtt("standard_name", units);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'standard_name' attribute: "
           << standard_name << endl;
      exit(1);
    }
  }
  
  return 0;
}

