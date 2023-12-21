#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <libgen.h>

#include "nc4utils.h"
#include "global_attrs.h"
#include "l1bgen_oci.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#define VERSION "0.045"

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           04/29/20 0.01  Original development
//  Joel Gales     SAIC           01/13/21 0.02  Add support for SWIR
//  Joel Gales     SAIC           08/11/21 0.03  Add support for handling
//                                               fill values in science data
//  Joel Gales     SAIC           08/12/21 0.04  Add support for hysteresis
//                                               correction
//  Joel Gales     SAIC           09/23/21 0.045 Initialize uninitialized
//                                               variables

ofstream tempOut;

int main (int argc, char* argv[])
{

  cout << "l1bgen_oci " << VERSION << " (" 
       <<  __DATE__ << " " << __TIME__ << ")" << endl;

  if ( argc == 1) {
    cout << endl << 
      "l1bgen_oci OCI_L1A_file cal_LUT_file OCI_L1B_file"
	 << endl;
    return 0;
  }
  
  // ********************************* //
  // *** Read calibration LUT file *** //
  // ********************************* //
  NcFile *calLUTfile = new NcFile( argv[2], NcFile::read);

  NcGroup gidCommon, gidBlue, gidRed, gidSWIR;
  gidCommon = calLUTfile->getGroup( "common");
  gidBlue = calLUTfile->getGroup( "blue");
  gidRed = calLUTfile->getGroup( "red");
  gidSWIR = calLUTfile->getGroup( "SWIR");
    
  float bwave[NBWAVE];
  float rwave[NRWAVE];
  float swave[NIWAVE];
  float spass[NIWAVE] = {45, 80, 30, 30, 15, 75, 75, 50, 75};
  double K2t[NTIMES];
  float K3T[NTEMPS];
  
  NcVar var;

  var = gidCommon.getVar( "blue_wavelength");
  var.getVar( bwave);

  var = gidCommon.getVar( "red_wavelength");
  var.getVar( rwave);

  var = gidCommon.getVar( "SWIR_wavelength");
  var.getVar( swave);

  var = gidCommon.getVar( "K2t");
  var.getVar( K2t);

  var = gidCommon.getVar( "K3T");
  var.getVar( K3T);

  string tag;

  cal_lut_struct blue_lut;
  cal_lut_struct red_lut;
  cal_lut_struct swir_lut;

  uint32_t bbanddim, rbanddim, sbanddim;
  
  tag.assign("blue");
  read_oci_cal_lut( calLUTfile, tag, gidBlue, bbanddim, blue_lut);

  tag.assign("red");
  read_oci_cal_lut( calLUTfile, tag, gidRed, rbanddim, red_lut);

  tag.assign("SWIR");
  read_oci_cal_lut( calLUTfile, tag, gidSWIR, sbanddim, swir_lut);

  // Read hysterisis parameters
  float hysttime[9][3];
  float hystamp[9][3];
  
  var = gidSWIR.getVar( "hyst_time_const");
  var.getVar( &hysttime[0][0]);
  var = gidSWIR.getVar( "hyst_amplitude");
  var.getVar( &hystamp[0][0]);
  
  calLUTfile->close();
  
  //  for (size_t i=0; i<NBWAVE; i++) cout << bwave[i] << endl;

  NcGroupAtt att;

  static l1bFile outfile;
  outfile.l1bfile = new NcFile( argv[3], NcFile::write);

  outfile.ncGrps[0] = outfile.l1bfile->getGroup( "sensor_band_parameters");
  outfile.ncGrps[1] = outfile.l1bfile->getGroup( "scan_line_attributes");
  outfile.ncGrps[2] = outfile.l1bfile->getGroup( "geolocation_data");
  outfile.ncGrps[3] = outfile.l1bfile->getGroup( "navigation_data");
  outfile.ncGrps[4] = outfile.l1bfile->getGroup( "observation_data");
  /*
  group: sensor_band_parameters {
group: scan_line_attributes {
group: geolocation_data {
group: navigation_data {
group: observation_data {
  */
  
  // Append call sequence to existing history
  string history = get_history(outfile.l1bfile);
  history.append(call_sequence(argc, argv));

  int32_t rta_nadir[2];
  att = outfile.l1bfile->getAtt("rta_nadir");
  att.getValues(rta_nadir);

  NcDim nscanl1b_dim = outfile.l1bfile->getDim("number_of_scans");
  uint32_t nscanl1b = nscanl1b_dim.getSize();
  
  double *evtime = new double[nscanl1b];

  var = outfile.ncGrps[1].getVar( "time");
  var.getVar( evtime);

  
  // Open a read data from L1Afile
  NcFile *l1afile = new NcFile( argv[1], NcFile::read);
   
  NcGroup ncGrps[4];

  ncGrps[0] = l1afile->getGroup( "scan_line_attributes");
  ncGrps[1] = l1afile->getGroup( "spatial_spectral_modes");
  ncGrps[2] = l1afile->getGroup( "engineering_data");
  ncGrps[3] = l1afile->getGroup( "science_data");

  // Get date (change this when year and day are added to time field)
  string tstart, tend;
  att = l1afile->getAtt("time_coverage_start");
  att.getValues(tstart);
  cout << tstart << endl;

  att = l1afile->getAtt("time_coverage_end");
  att.getValues(tend);
  cout << tend << endl;

  uint16_t iyr, imn, idom;
  istringstream iss;

  iss.str(tstart.substr(0,4));
  iss >> iyr; iss.clear();
  iss.str(tstart.substr(5,2));
  iss >> imn; iss.clear();
  iss.str(tstart.substr(8,2));
  iss >> idom;
  int32_t jd = jday(iyr, imn, idom);

  int32_t iyr32, idy32;
  jdate( jd, &iyr32, &idy32);

  // Get numbers of blue and red bands
  NcDim blue_dim = l1afile->getDim("blue_bands");
  uint32_t bbands = blue_dim.getSize();
  NcDim red_dim = l1afile->getDim("red_bands");
  uint32_t rbands = red_dim.getSize();

  // Scan time, spin ID and HAM side
  NcDim nscan_dim = l1afile->getDim("number_of_scans");
  uint32_t nscan = nscan_dim.getSize();

  NcDim mcescan_dim = l1afile->getDim("number_of_mce_scans");
  uint32_t nmcescan = mcescan_dim.getSize();
  
  double *sstime = new double[nscan];
  var = ncGrps[0].getVar( "scan_start_time");
  var.getVar( sstime);

  int32_t *spin = new int32_t[nscan];
  var = ncGrps[0].getVar( "spin_ID");
  var.getVar( spin);

  uint8_t *hside = new uint8_t[nscan];
  var = ncGrps[0].getVar( "HAM_side");
  var.getVar( hside);
  
  uint32_t nscan_good=0;
  for (size_t i=0; i<nscan; i++) {
    if ( spin[i] > 0) {
      sstime[nscan_good] = sstime[i];
      hside[nscan_good] = hside[i];
      nscan_good++;
    }
  }

  // Check for and fill in missing scan times
  ////////////////// check_scan_times, sstime, sfl JMG


  // ******************************************** //
  // *** Get spatial and spectral aggregation *** //
  // ******************************************** //
  NcDim ntaps_dim = l1afile->getDim("number_of_taps");
  uint32_t ntaps = ntaps_dim.getSize();
  NcDim spatzn_dim = l1afile->getDim("spatial_zones");
  uint32_t spatzn = spatzn_dim.getSize();

  int16_t *dtype = new int16_t[spatzn];
  var = ncGrps[1].getVar( "spatial_zone_data_type");
  var.getVar( dtype);

  int16_t *lines = new int16_t[spatzn];
  var = ncGrps[1].getVar( "spatial_zone_lines");
  var.getVar( lines);

  int16_t *iagg = new int16_t[spatzn];
  var = ncGrps[1].getVar( "spatial_aggregation");
  var.getVar( iagg);

  int16_t *bagg = new int16_t[ntaps];
  var = ncGrps[1].getVar( "blue_spectral_mode");
  var.getVar( bagg);

  int16_t *ragg = new int16_t[ntaps];
  var = ncGrps[1].getVar( "red_spectral_mode");
  var.getVar( ragg);

  
  // ********************************************************************* //
  // *** Get # of EV lines/offset from scan start time to EV mid-time  *** //
  // ********************************************************************* //
  // This will be done by geolocation when integrated into L1B
  int32_t ppr_off;
  double revpsec, secpline;

  int32_t *mspin = new int32_t[nmcescan];
  int32_t *ot_10us = new int32_t[nmcescan];
  uint8_t *enc_count = new uint8_t[nmcescan];
  int16_t board_id;
  
  NcDim nenc_dim = l1afile->getDim("encoder_samples");
  uint32_t nenc = nenc_dim.getSize();
  
  float **hamenc = new float *[nmcescan];
  hamenc[0] = new float[nenc*nmcescan];
  for (size_t i=1; i<nmcescan; i++) hamenc[i] = hamenc[i-1] + nenc;

  float **rtaenc = new float *[nmcescan];
  rtaenc[0] = new float[nenc*nmcescan];
  for (size_t i=1; i<nmcescan; i++) rtaenc[i] = rtaenc[i-1] + nenc;

  read_mce_tlm( l1afile, ncGrps[2], nmcescan, nenc, ppr_off, revpsec, secpline,
                board_id, mspin, ot_10us, enc_count, &hamenc[0], rtaenc);

  float clines[32400], slines[4050];
  uint16_t pcdim, psdim;
  int16_t iret;
  double ev_toff, deltc[32400], delts[4050];
  get_ev( secpline, dtype, lines, iagg, pcdim, psdim, ev_toff, clines, slines,
          deltc, deltc, iret);
  if ( iret < 0) {
    cout << "No science collect in file: " << argv[1] << endl;
    l1afile->close();
    return 1;
  }

  size_t ka;
  for (size_t i=0; i<spatzn; i++) {
    if ( dtype[i] != 0 && dtype[i] != 2 && dtype[i] != 10) {
      ka = i;
      break;
    }
  }


  // *********************************************************** //
  // *** Generate matrices for spectral and gain aggregation *** //
  // *********************************************************** //
  size_t *ia;
  uint32_t iia;
  uint32_t ntb[16];
  ia = new size_t[ntaps];
  
  // Blue bands
  iia=0;
  for (size_t i=0; i<ntaps; i++) {
    if ( bagg[i] > 0) {
      ia[iia] = i;
      iia++;
    }
  }

  uint32_t bib=1, bbb=1;
  if ( iia == 0) {
    cout << "All blue taps disabled" << endl;
  } else {
    for (size_t i=0; i<16; i++) ntb[i] = 0;
    for (size_t i=0; i<iia; i++)
      ntb[ia[i]] = 32 / bagg[ia[i]];
    bib = 0;
    for (size_t i=0; i<16; i++) bib += ntb[i];
  }
  
  // Note: bgmat & rgmat not necessarily contiguous
  float **bamat = new float*[bib];  
  float **bgmat = new float*[512];

  if (bib != 1)
    get_agg_mat( ia, iagg[ka], bagg, bib, bbb, ntb, bamat, bgmat);

  if (bib != bbands) {
    cout << "Number of blue bands in file: " << argv[1] <<
      " not consistent with spectral aggregation" << endl;
    l1afile->close();
    return 1;
  } else if ( bib < 4) cout << "No blue bands in file: " << argv[1] << endl;
    
  // Red bands
  iia=0;
  for (size_t i=0; i<ntaps; i++) {
    if ( ragg[i] > 0) {
      ia[iia] = i;
      iia++;
    }
  }

  uint32_t rib=1, rbb=1;
  if ( iia == 0) {
    cout << "All red taps disabled" << endl;
  } else {
    for (size_t i=0; i<16; i++) ntb[i] = 0;
    for (size_t i=0; i<iia; i++) ntb[ia[i]] = 32 / ragg[ia[i]];
    rib = 0;
    for (size_t i=0; i<16; i++) rib += ntb[i];
  }
  
  float **ramat = new float*[rib];  
  float **rgmat = new float*[512];

  if (rib != 1)
    get_agg_mat( ia, iagg[ka], ragg, rib, rbb, ntb, ramat, rgmat);

  if (rib != rbands) {
    cout << "Number of red bands in file: " << argv[1] <<
      " not consistent with spectral aggregation" << endl;
    l1afile->close();
    return 1;
  } else if ( rib < 4) cout << "No red bands in file: " << argv[1] << endl;

  uint16_t swb = 9;


  // ********************************* //
  // *** Get dark collect location *** //
  // ********************************* //
  int16_t kd=-1;
  for (size_t i=0; i<spatzn; i++) {
    if ( dtype[i] == 2) {
      kd = (int16_t) i;
    }
  }
  if ( kd == -1) {
    cout << "No dark collect in file: " << argv[1] << endl;
    l1afile->close();
    return 1;
  }

  int16_t ldc=0, lds=0;
  for (size_t i=0; i<(size_t) kd; i++) {
    if ( dtype[i] != 0 && dtype [1] != 10) {
      ldc += lines[i] / iagg[i];
      lds += lines[i] / 8;
    }
  }
  int16_t ndc = lines[kd] / iagg[kd];
  int16_t nds = lines[kd] / 8;


  // *********************************************************************** //
  // *** Generate band gain structs from LUTs, date/time & gain matrices *** //
  // *********************************************************************** //
  gains_struct bgains;
  gains_struct rgains;
  gains_struct sgains;
  /*
  for (size_t i=0; i<64; i++)
    for (size_t j=0; j<512; j++)
      cout << "0: " << i << " " << j << " " << bgmat[j][i] << endl;
  */
  // Note: bgmat & rgmat not necessarily contiguous
  // That's why we pass the location of the array of pointers
  if ( bib >= 4)
    make_oci_gains( bib, bbanddim, iyr, idom, evtime[0], K2t, blue_lut,
                    &bgmat[0], bgains);

  if ( rib >= 4)
    make_oci_gains( rib, rbanddim, iyr, idom, evtime[0], K2t, red_lut,
                    &rgmat[0], rgains);

  float **sgmat = new float*[swb];
  for (size_t i=0; i<swb; i++) {
    sgmat[i] = new float[swb];
    for (size_t j=0; j<swb; j++) {
      if (i == j) sgmat[i][j] = 1.0; else sgmat[i][j] = 0.0;
    }
  }
  make_oci_gains( swb, swb, iyr, idom, evtime[0], K2t, swir_lut,
                  &sgmat[0], sgains);

  // Read selected temperature fields and interpolate to scan times
  float **caltemps = new float *[nscan_good];
  caltemps[0] = new float[30*nscan_good];
  for (size_t i=1; i<nscan_good; i++) caltemps[i] = caltemps[i-1] + 30;

  // get_oci_cal_temps JMG not yet defined
  for (size_t i=0; i<nscan_good; i++)
    for (size_t j=1; j<30; j++)
      caltemps[i][j] = 0.0;
  
  // Read dark collects from science data arrays
  vector<size_t> start, count;
  start.push_back(0);
  start.push_back(0);
  start.push_back(ldc);

  size_t dims[3];
  dims[0] = nscan_good; dims[1] = bib; dims[2] = ndc;
  uint16_t ***bdark = make3dT<uint16_t>(dims);
  dims[0] = nscan_good; dims[1] = rib; dims[2] = ndc;
  uint16_t ***rdark = make3dT<uint16_t>(dims);

  uint16_t bfill;
  uint16_t rfill;
  
  if ( bib > 4) {
    count.push_back(nscan_good);
    count.push_back(bib);
    count.push_back(ndc);
    
    var = ncGrps[3].getVar( "sci_blue");
    var.getVar( start, count, &bdark[0][0][0]);

    var.getAtt("_FillValue").getValues(&bfill);
  }

  if ( rib > 4) {
    count.clear();
    count.push_back(nscan_good);
    count.push_back(rib);
    count.push_back(ndc);
    
    var = ncGrps[3].getVar( "sci_red");
    var.getVar( start, count, &rdark[0][0][0]);

    var.getAtt("_FillValue").getValues(&rfill);
  }

  dims[0] = nscan_good; dims[1] = swb; dims[2] = nds;
  uint32_t ***sdark = make3dT<uint32_t>(dims);

  uint32_t sfill;
  
  start.clear();
  start.push_back(0);
  start.push_back(0);
  start.push_back(lds);

  count.clear();
  count.push_back(nscan_good);
  count.push_back(swb);
  count.push_back(nds);

  var = ncGrps[3].getVar( "sci_SWIR");
  var.getVar( start, count, &sdark[0][0][0]);

  var.getAtt("_FillValue").getValues(&sfill);

  uint32_t fill32;
  
  // number of scans of dark data to average; will make this an input parameter
  uint16_t ndsc = 1;
  //number of dark pixels to skip; will make this an input parameter
  uint16_t nskp = 0;
  // adjust SWIR skip factor in case of CCD aggregation LT 8
  //  uint32_t nsskp = (nskp * iagg[kd]) / 8;

  uint32_t bicount[3] = {1,bib,(uint32_t) pcdim};
  uint32_t ricount[3] = {1,rib,(uint32_t) pcdim};
  uint32_t sicount[3] = {1,swb,(uint32_t) psdim};
  uint32_t bcount[3] = {bbb,1,(uint32_t) pcdim};
  uint32_t rcount[3] = {rbb,1,(uint32_t) pcdim};
  uint32_t scount[3] = {swb,1,(uint32_t) psdim};

  // Calibrated data variables
  float **bdn = new float *[bib];
  bdn[0] = new float[pcdim*bib];
  for (size_t i=1; i<bib; i++) bdn[i] = bdn[i-1] + pcdim;

  float **rdn = new float *[rib];
  rdn[0] = new float[pcdim*rib];
  for (size_t i=1; i<rib; i++) rdn[i] = rdn[i-1] + pcdim;

  float **sdn = new float *[swb];
  sdn[0] = new float[psdim*swb];
  for (size_t i=1; i<swb; i++) sdn[i] = sdn[i-1] + psdim;

  float **bcal = new float *[bib];
  bcal[0] = new float[pcdim*bib];
  for (size_t i=1; i<bib; i++) bcal[i] = bcal[i-1] + pcdim;

  float **rcal = new float *[rib];
  rcal[0] = new float[pcdim*rib];
  for (size_t i=1; i<rib; i++) rcal[i] = rcal[i-1] + pcdim;

  float **scal = new float *[swb];
  scal[0] = new float[psdim*swb];
  for (size_t i=1; i<swb; i++) scal[i] = scal[i-1] + psdim;

  double *thetap = new double[pcdim];
  double *thetas = new double[psdim];

  float **pview = new float *[pcdim];
  pview[0] = new float[3*pcdim];
  for (size_t i=1; i<pcdim; i++) pview[i] = pview[i-1] + 3;

  float **sview = new float *[psdim];
  sview[0] = new float[3*psdim];
  for (size_t i=1; i<psdim; i++) sview[i] = sview[i-1] + 3;


  // Main loop
  // Read, calibrate and write science data
  for (size_t iscn=0; iscn<nscan_good; iscn++) {

    if ((iscn % 50) == 0) cout << "Calibrating scan: " << iscn << endl;

    // Check for valid mirror side
    if ( hside[iscn] == 0 || hside[iscn] == 1) {


      // Get scan angle
      // This will be done by geolocation when integrated into L1B
      get_oci_vecs( nscan, pcdim, rta_nadir, ev_toff, spin[iscn], clines, deltc,
                    revpsec, ppr_off, board_id, mspin, ot_10us, enc_count,
                    &hamenc[0], &rtaenc[0], pview, thetap, iret);

      //      tempOut.open ("pview.bin", ios::out | ios::trunc | ios::binary);
      //tempOut.write((char *) &pview[0][0], sizeof(float)*3*pcdim);
      //tempOut.close();

      //tempOut.open ("thetap.bin", ios::out | ios::trunc | ios::binary);
      //tempOut.write((char *) &thetap[0], sizeof(double)*pcdim);
      //tempOut.close();

      get_oci_vecs( nscan, psdim, rta_nadir, ev_toff, spin[iscn], slines, delts,
                    revpsec, ppr_off, board_id, mspin, ot_10us, enc_count,
                    &hamenc[0], &rtaenc[0], sview, thetas, iret);

      //tempOut.open ("sview.bin", ios::out | ios::trunc | ios::binary);
      //tempOut.write((char *) &sview[0][0], sizeof(float)*3*psdim);
      //tempOut.close();

      //tempOut.open ("thetas.bin", ios::out | ios::trunc | ios::binary);
      //tempOut.write((char *) &thetas[0], sizeof(double)*psdim);
      //tempOut.close();
      
      //  Blue bands
      if (bib >= 4) {

        start.clear();
        start.push_back(iscn);
        start.push_back(0);
        start.push_back(0);

        count.clear();
        count.push_back(bicount[0]);
        count.push_back(bicount[1]);
        count.push_back(bicount[2]);

        uint16_t **bsci = new uint16_t *[bib];
        bsci[0] = new uint16_t[pcdim*bib];
        for (size_t i=1; i<bib; i++) bsci[i] = bsci[i-1] + pcdim;

        var = ncGrps[3].getVar( "sci_blue");
        var.getVar( start, count, bsci[0]);

        //tempOut.open ("bsci.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) bsci[0], sizeof(uint16_t)*pcdim*bib);
        //tempOut.close();
      
        // Compute dark offset, correct data, and apply absolute and
        // temporal gain and temperature correction
        float *bdc = new float[bib];

        fill32 = bfill;
        int16_t iret;
        get_oci_dark<uint16_t>( iscn, nscan_good, hside, ndsc, nskp,
                                iagg[ka], iagg[kd], ntaps, bagg, fill32,
                                ndc, bdark, bib, bdc, iret);

        float *k3 = new float[bib];
        get_oci_temp_corr( bib, bgains, K3T, caltemps[iscn], nscan_good, k3);
        for (size_t j=0; j<bib; j++) {
          for (size_t k=0; k<pcdim; k++) {

            // Handle fill value
            if (bsci[j][k] == bfill) {
              bdn[j][k] = -999;
              bcal[j][k] = -999;
              continue;
            }

            // Need to save dn for linearity correction
            bdn[j][k] = bsci[j][k] - bdc[j];
            bcal[j][k] = k3[j] * bgains.K1K2[j][hside[iscn]] * bdn[j][k];
          }
        }

        delete [] k3;
        delete [] bdc;
        
        // Compute and apply RVS and linearity
        //k4 = fltarr(pdim,nib)
        float **k4 = new float *[bib];
        k4[0] = new float[pcdim*bib];
        for (size_t i=1; i<bib; i++) k4[i] = k4[i-1] + pcdim;

        get_oci_rvs_corr( bib, pcdim, hside[iscn], bgains, thetap, k4);

        float **k5 = new float *[bib];
        k5[0] = new float[pcdim*bib];
        for (size_t i=1; i<bib; i++) k5[i] = k5[i-1] + pcdim;
        
        get_oci_lin_corr( bib, pcdim, bgains, K3T, caltemps[iscn], bdn, k5);

        for (size_t j=0; j<bib; j++) {
          for (size_t k=0; k<pcdim; k++) {
            bcal[j][k] *= k4[j][k] * k5[j][k];
          }
        }

        delete [] k4[0];
        delete [] k4;
        delete [] k5[0];
        delete [] k5;

        float **bcalb = new float *[bbb];
        bcalb[0] = new float[pcdim*bbb];
        for (size_t i=1; i<bbb; i++) bcalb[i] = bcalb[i-1] + pcdim;

        // Aggregate to L1B bands
        //bcalb = transpose(bamat#transpose(bcal))
        for (size_t j=0; j<bbb; j++) {
          for (size_t k=0; k<pcdim; k++) {
            float sum = 0.0;
            for (size_t l=0; l<bib; l++)
              if (bcal[l][k] != -999) sum += bamat[l][j]*bcal[l][k];
            bcalb[j][k] = sum;
          }
        }

        start.clear();
        start.push_back(0);
        start.push_back(iscn);
        start.push_back(0);

        count.clear();
        count.push_back(bcount[0]);
        count.push_back(bcount[1]);
        count.push_back(bcount[2]);
        
        // Output to L1B file
        var = outfile.ncGrps[4].getVar( "Lt_blue");
        var.putVar( start, count, &bcalb[0][0]);

        delete [] bcalb[0];
        delete [] bcalb;
        delete [] bsci[0];
        delete [] bsci;
      } // End blue


      //  Red bands
      if (rib >= 4) {

        start.clear();
        start.push_back(iscn);
        start.push_back(0);
        start.push_back(0);

        count.clear();
        count.push_back(ricount[0]);
        count.push_back(ricount[1]);
        count.push_back(ricount[2]);

        uint16_t **rsci = new uint16_t *[rib];
        rsci[0] = new uint16_t[pcdim*rib];
        for (size_t i=1; i<rib; i++) rsci[i] = rsci[i-1] + pcdim;

        var = ncGrps[3].getVar( "sci_red");
        var.getVar( start, count, rsci[0]);

        //tempOut.open ("rsci.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) rsci[0], sizeof(uint16_t)*pcdim*rib);
        //tempOut.close();
      
        // Compute dark offset, correct data, and apply absolute and
        // temporal gain and temperature correction
        float *rdc = new float[rib];
        fill32 = rfill;
        int16_t iret;
        get_oci_dark<uint16_t>( iscn, nscan_good, hside, ndsc, nskp,
                                iagg[ka], iagg[kd], ntaps, ragg, fill32,
                                ndc, rdark, rib, rdc, iret);

        //tempOut.open ("rdc.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) rdc, sizeof(float)*rib);
        //tempOut.close();
        float *k3 = new float[rib];
        get_oci_temp_corr( rib, rgains, K3T, caltemps[iscn], nscan_good, k3);
        //tempOut.open ("k3.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) k3, sizeof(float)*rib);
        //tempOut.close();
        for (size_t j=0; j<rib; j++) {
          for (size_t k=0; k<pcdim; k++) {

            // Handle fill value
            if (rsci[j][k] == rfill) {
              rdn[j][k] = -999;
              rcal[j][k] = -999;
              continue;
            }
            
            // Need to save dn for linearity correction
            rdn[j][k] = rsci[j][k] - rdc[j];
            rcal[j][k] = k3[j] * rgains.K1K2[j][hside[iscn]] * rdn[j][k];
          }
        }
        //tempOut.open ("rdn.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) rdn[0], sizeof(float)*pcdim*rib);
        //tempOut.close();

        delete [] k3;
        delete [] rdc;
        
        // Compute and apply RVS and linearity
        //k4 = fltarr(pdim,nib)
        float **k4 = new float *[rib];
        k4[0] = new float[pcdim*rib];
        for (size_t i=1; i<rib; i++) k4[i] = k4[i-1] + pcdim;
        get_oci_rvs_corr( rib, pcdim, hside[iscn], rgains, thetap, k4);
        //tempOut.open ("k4.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) k4[0], sizeof(float)*pcdim*rib);
        //tempOut.close();

        float **k5 = new float *[rib];
        k5[0] = new float[pcdim*rib];
        for (size_t i=1; i<rib; i++) k5[i] = k5[i-1] + pcdim;
        
        get_oci_lin_corr( rib, pcdim, rgains, K3T, caltemps[iscn], rdn, k5);
        //tempOut.open ("k5.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) k5[0], sizeof(float)*pcdim*rib);
        //tempOut.close();
        for (size_t j=0; j<rib; j++) {
          for (size_t k=0; k<pcdim; k++) {
            rcal[j][k] *= k4[j][k] * k5[j][k];
          }
        }

        //tempOut.open ("rcal.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) rcal[0], sizeof(float)*pcdim*rib);
        //tempOut.close();
        
        delete [] k4[0];
        delete [] k4;
        delete [] k5[0];
        delete [] k5;

        float **rcalb = new float *[rbb];
        rcalb[0] = new float[pcdim*rbb];
        for (size_t i=1; i<rbb; i++) rcalb[i] = rcalb[i-1] + pcdim;

        // Aggregate to L1B bands
        for (size_t j=0; j<rbb; j++) {
          for (size_t k=0; k<pcdim; k++) {
            float sum = 0.0;
            for (size_t l=0; l<rib; l++)
              if (rcal[l][k] != -999) sum += ramat[l][j]*rcal[l][k];
            rcalb[j][k] = sum;
          }
        }

        //tempOut.open ("rcalb.bin", ios::out | ios::trunc | ios::binary);
        //tempOut.write((char *) rcalb[0], sizeof(float)*pcdim*rib);
        //tempOut.close();

        start.clear();
        start.push_back(0);
        start.push_back(iscn);
        start.push_back(0);

        count.clear();
        count.push_back(rcount[0]);
        count.push_back(rcount[1]);
        count.push_back(rcount[2]);
  
        // Output to L1B file
        var = outfile.ncGrps[4].getVar( "Lt_red");
        var.putVar( start, count, &rcalb[0][0]);

        delete [] rcalb[0];
        delete [] rcalb;
        delete [] rsci[0];
        delete [] rsci;
      } // End red

      //  SWIR bands
      start.clear();
      start.push_back(iscn);
      start.push_back(0);
      start.push_back(0);

      count.clear();
      count.push_back(sicount[0]);
      count.push_back(sicount[1]);
      count.push_back(sicount[2]);

      uint32_t **ssci = new uint32_t *[swb];
      ssci[0] = new uint32_t[psdim*swb];
      for (size_t i=1; i<swb; i++) ssci[i] = ssci[i-1] + psdim;

      var = ncGrps[3].getVar( "sci_SWIR");
      var.getVar( start, count, ssci[0]);

      // Compute dark offset, correct data, and apply absolute and
      // temporal gain and temperature correction
      float *sdc = new float[swb];
      //      int16_t sagg = 8;
      int16_t sagg = 1;
      int16_t iret;
      get_oci_dark<uint32_t>( iscn, nscan_good, hside, ndsc, nskp,
                              1, 1, 1, &sagg, sfill, nds, sdark,
                              swb, sdc, iret);

      float *k3 = new float[swb];
      if ( iret != -1) {
        get_oci_temp_corr( swb, sgains, K3T, caltemps[iscn], nscan_good, k3);

        for (size_t j=0; j<swb; j++) {
          for (size_t k=0; k<psdim; k++) {

            // Handle fill value
            if (ssci[j][k] == sfill) {
              sdn[j][k] = -999;
              scal[j][k] = -999;
              continue;
            }
          
            // Need to save dn for linearity correction
            sdn[j][k] = ssci[j][k] - sdc[j];

            // Hysteresis correction
            float hc_prev[3];
            float hc[3]={0,0,0};
            float hyst = 0.0;
            for (size_t l=0; l<3; l++) {
              // Compute exponential decay constants
              float e = exp(-1.0/hysttime[j][l]);

              if ( k > 0) {
                hc[l] = hc_prev[l]*e + sdn[j][k-1]*hystamp[j][l];
                hyst += hc[l];
              }
              hc_prev[l] = hc[l];
            } // l-loop

            scal[j][k] =
              k3[j] * sgains.K1K2[j][hside[iscn]] * (sdn[j][k] - hyst);
          } // k-loop
        }
      } // iret != -1
      
      delete [] k3;
      delete [] sdc;
              
      // Compute and apply RVS and linearity
      float **k4 = new float *[swb];
      k4[0] = new float[psdim*swb];
      for (size_t i=1; i<swb; i++) k4[i] = k4[i-1] + psdim;
      if ( iret != -1)
        get_oci_rvs_corr( swb, psdim, hside[iscn], sgains, thetap, k4);

      float **k5 = new float *[swb];
      k5[0] = new float[psdim*swb];
      for (size_t i=1; i<swb; i++) k5[i] = k5[i-1] + psdim;

      if ( iret != -1)
        get_oci_lin_corr( swb, psdim, sgains, K3T, caltemps[iscn], sdn, k5);

      for (size_t j=0; j<swb; j++) {
        for (size_t k=0; k<psdim; k++) {
          if (scal[j][k] != -999)
            scal[j][k] *= k4[j][k] * k5[j][k];
        }
      }
        
      delete [] k4[0];
      delete [] k4;
      delete [] k5[0];
      delete [] k5;

      start.clear();
      start.push_back(0);
      start.push_back(iscn);
      start.push_back(0);

      count.clear();
      count.push_back(scount[0]);
      count.push_back(scount[1]);
      count.push_back(scount[2]);
  
      // Output to L1B file
      var = outfile.ncGrps[4].getVar( "Lt_SWIR");
      var.putVar( start, count, &scal[0][0]);

      delete [] ssci[0];
      delete [] ssci;

    } else {
      cout << "No mirror side index for scan: " << iscn << endl;
    } // Check for valid mirror side
  } // Scan loop

  // End Main loop
  
  // Write spectral band information
  // Calculate band centers for aggregated hyperspectral bands
  // b1bwave = bamat#bgmat#bwave

  if (bib >= 4) {
    // bgmat#bwave
    float *b1bwave = new float[bbb];
    float *sum = new float[bib];
    for (size_t i=0; i<bib; i++) {
      sum[i] = 0.0;
      for (size_t j=0; j<512; j++) {
        sum[i] += bwave[j]*bgmat[j][i];
      }
    }

    // bamat#sum
    for (size_t i=0; i<bbb; i++) {
      b1bwave[i] = 0.0;
      for (size_t j=0; j<bib; j++) {
        b1bwave[i] += bamat[j][i]*sum[j];
      }
    }
    
    start.clear();
    start.push_back(0);

    count.clear();
    count.push_back(bbb);
    var = outfile.ncGrps[0].getVar( "blue_wavelength");
    var.putVar( start, count, b1bwave);

    delete [] b1bwave;
    delete [] sum;
  }

  if (rib >= 4) {
    float *r1bwave = new float[rbb];
    float *sum = new float[rib];
    for (size_t i=0; i<rib; i++) {
      sum[i] = 0.0;
      for (size_t j=0; j<512; j++) {
        sum[i] += rwave[j]*rgmat[j][i];
      }
    }

    for (size_t i=0; i<rbb; i++) {
      r1bwave[i] = 0.0;
      for (size_t j=0; j<rib; j++) {
        r1bwave[i] += ramat[j][i]*sum[j];
      }
    }

    start.clear();
    start.push_back(0);

    count.clear();
    count.push_back(rbb);
    var = outfile.ncGrps[0].getVar( "red_wavelength");
    var.putVar( start, count, r1bwave);

    delete [] r1bwave;
    delete [] sum;
  }

  // SWIR wavelengths/bandpass
  start.clear();
  start.push_back(0);
  count.clear();
  count.push_back(NIWAVE);
  
  var = outfile.ncGrps[0].getVar( "SWIR_wavelength");
  var.putVar( start, count, swave);
  var = outfile.ncGrps[0].getVar( "SWIR_bandpass");
  var.putVar( start, count, spass);

  
  string l1b_filename;
  l1b_filename.assign(argv[3]);
  outfile.write_granule_metadata( tstart, tend, l1b_filename);

  outfile.close();
  set_global_attrs(l1b_filename, history, "");
  
  delete [] sstime;
  delete [] spin;
  delete [] hside;
  delete [] dtype;
  delete [] lines;
  delete [] iagg;
  delete [] bagg;
  delete [] ragg;
  delete [] mspin;
  delete [] ot_10us;
  delete [] enc_count;
  delete [] sgmat;
  delete [] thetap;
  delete [] thetas;
  delete [] ia;
  delete [] evtime;
  
  delete [] hamenc[0];
  delete [] hamenc;
  delete [] rtaenc[0];
  delete [] rtaenc;
  delete [] caltemps[0];
  delete [] caltemps;
  delete [] bdn[0];
  delete [] bdn;
  delete [] rdn[0];
  delete [] rdn;
  delete [] sdn[0];
  delete [] sdn;  
  delete [] bcal[0];
  delete [] bcal;
  delete [] rcal[0];
  delete [] rcal;
  delete [] scal[0];
  delete [] scal;  
  delete [] pview[0];
  delete [] pview;
  delete [] sview[0];
  delete [] sview;  

  delete [] blue_lut.K1[0];
  delete [] blue_lut.K1;
  delete [] blue_lut.K5[0];
  delete [] blue_lut.K5;

  if (bgains.K1K2 != NULL) delete [] bgains.K1K2[0];
  if (bgains.K1K2 != NULL) delete [] bgains.K1K2;
  if (bgains.K5 != NULL) delete [] bgains.K5[0];
  if (bgains.K5 != NULL) delete [] bgains.K5;

  // delete [] arrT[0][0]; delete [] arrT[0]; delete [] arrT;

  delete [] red_lut.K1[0];
  delete [] red_lut.K1;
  delete [] red_lut.K5[0];
  delete [] red_lut.K5;

  if (rgains.K1K2 != NULL) delete [] rgains.K1K2[0];
  if (rgains.K1K2 != NULL) delete [] rgains.K1K2;
  if (rgains.K5 != NULL) delete [] rgains.K5[0];
  if (rgains.K5 != NULL) delete [] rgains.K5;
  
  delete [] swir_lut.K1[0];
  delete [] swir_lut.K1;
  delete [] swir_lut.K5[0];
  delete [] swir_lut.K5;

  delete [] sgains.K1K2[0];
  delete [] sgains.K1K2;
  delete [] sgains.K5[0];
  delete [] sgains.K5;

  // Add delete for dark arrays
  
  return 0;
}

int get_agg_mat( size_t *ia, int16_t iagg, int16_t jagg[16], uint16_t nib,
                 uint32_t& nbb, uint32_t ntb[16], float **amat, float **gmat) {

  // Compute number of bands for 8x aggregation with overlapping bands
  nbb = (ntb[0]*3) / 4 + 1;
  //amat(nbb,nib)
  for (size_t i=1; i<16; i++) {
    if (jagg[i] >= jagg[i-1])
      nbb += ntb[i];
    else
      nbb += (ntb[i]*3) / 4 + ntb[i-1]/4;
  }
  
  for (size_t i=0; i<512; i++) {
    gmat[i] = new float[nib];
    for (size_t j=0; j<nib; j++) gmat[i][j] = 0.0;
  }
  
  // Populate gain aggregation matrix
  int16_t ii = 0;
  int16_t itt[16][2];
  for (size_t i=0; i<16; i++) {
    itt[i][0] = itt[i][1] = 0;
    if ( jagg[i] > 0) {
      int16_t iaf = 4;
      if ( iagg*jagg[i] < iaf) iaf = iagg*jagg[i];
      for (size_t k=0; k<ntb[i]; k++) {
        size_t ic = 32*i;
        size_t kj = k*jagg[i];
        for (size_t j=0; j<(size_t) jagg[i]; j++)
          gmat[ic+kj+j][ii+k] = (1.0/jagg[i]) * (4/iaf);
      }
      itt[i][0] = ii;
      ii += ntb[i];
      itt[i][1] = ii - 1;
    }
  }

  for (size_t i=0; i<nib; i++) {
    amat[i] = new float[nbb];
    for (size_t j=0; j<nbb; j++) amat[i][j] = 0;
  }

  // First tap
  int16_t ib;
  // for k=0,ntb(ia(0))*3/4 do amat(k,k:k+ntb(ia(0))/4-1) = jagg(ia(0))/8.0
  for (size_t k=0; k<=(ntb[ia[0]]*3)/4; k++) {
    for (size_t l=0; l<=ntb[ia[0]]/4-1; l++) amat[k+l][k] = jagg[ia[0]]/8.0;
  }
  ib = (ntb[ia[0]]*3)/4 + 1;

  // Remaining taps
  uint16_t nr;
  for (size_t i=ia[1]; i<16; i++) {
    if (ntb[i] > 0) {
      if (ntb[i] >= ntb[i-1]) {
        // Transition resolution determined by preceding tap
      	nr = ntb[i-1]/4 - 1;
        // Remaining bands using preceding tap
      	if (nr > 0) {
          for (size_t k=0; k<nr; k++) {
	    uint16_t k1 = nr - k - 1;
            uint16_t k2 = ((k+1)*ntb[i]) / ntb[i-1] - 1;
            for (size_t j=0; j<=k1; j++)
              amat[itt[i-1][1]-k1+j][ib+k] = jagg[i-1] / 8.0;
            for (size_t j=0; j<=k2; j++)
              amat[itt[i][0]+j][ib+k] = jagg[i] / 8.0;
          }
      	  ib += nr;
        }
      } else {
        // Transition resolution determined using current tap
        nr = ntb[i]/4 - 1;
        // Remaining bands using previous tap
        if (nr > 0) {
          for (size_t k=0; k<nr; k++) {
            uint16_t k1 = ((nr-k)*ntb[i-1]) / ntb[i] - 1;
	    uint16_t k2 = k;
            for (size_t j=0; j<=k1; j++)
              amat[itt[i-1][1]-k1+j][ib+k] = jagg[i-1] / 8.0;
            for (size_t j=0; j<=k2; j++)
              amat[itt[i][0]+j][ib+k] = jagg[i] / 8.0;
          }            
   	  ib += nr;
        }
      }
      // Remaining bands using this tap
      for (size_t k=0; k<=(ntb[i]*3)/4; k++) {
        for (size_t j=0; j<ntb[i]/4; j++)
        amat[itt[i][0]+k+j][ib+k] = jagg[i] / 8.0;
      }
      ib += (ntb[i]*3) / 4 + 1;
    }
  }
                            
  return 0;
}


int read_oci_cal_lut( NcFile *calLUTfile, string tag, NcGroup gidLUT,
                      uint32_t& banddim, cal_lut_struct& cal_lut) {

  size_t dims[3];
  NcDim ncDIM;
  uint32_t timedim, tempdim, tcdim, rvsdim, nldim, msdim=2;
  string bandname;
  bandname.assign( tag);
  bandname.append( "_bands");
  
  ncDIM = calLUTfile->getDim(bandname.c_str());
  banddim = ncDIM.getSize();
  
  ncDIM = calLUTfile->getDim("number_of_times");
  timedim = ncDIM.getSize();
  ncDIM = calLUTfile->getDim("number_of_temperatures");
  tempdim = ncDIM.getSize();
  ncDIM = calLUTfile->getDim("number_of_T_coefficients");
  tcdim = ncDIM.getSize();
  ncDIM = calLUTfile->getDim("number_of_RVS_coefficients");
  rvsdim = ncDIM.getSize();
  ncDIM = calLUTfile->getDim("number_of_nonlinearity_coefficients");
  nldim = ncDIM.getSize();

  float **K1 = new float *[banddim];
  K1[0] = new float[msdim*banddim];
  for (size_t i=1; i<banddim; i++) K1[i] = K1[i-1] + msdim;
  cal_lut.K1 = K1;

  dims[0] = banddim; dims[1] = msdim; dims[2] = timedim;
  float ***K2 = make3dT<float>(dims);
  cal_lut.K2 = K2;
  dims[0] = banddim; dims[1] = tempdim; dims[2] = tcdim;
  float ***K3_coef = make3dT<float>(dims);
  cal_lut.K3_coef = K3_coef;
  dims[0] = banddim; dims[1] = msdim; dims[2] = rvsdim;
  float ***K4_coef = make3dT<float>(dims);
  cal_lut.K4_coef = K4_coef;

  double **K5 = new double *[banddim];
  K5[0] = new double[nldim*banddim];
  for (size_t i=1; i<banddim; i++) K5[i] = K5[i-1] + nldim;
  cal_lut.K5 = K5;
  
  cal_lut.ldims[0] = timedim; cal_lut.ldims[1] = tempdim;
  cal_lut.ldims[2] = tcdim; cal_lut.ldims[3] = rvsdim;
  cal_lut.ldims[4] = nldim; cal_lut.ldims[5] = msdim;

  NcVar var;
  
  var = gidLUT.getVar( "K1");
  var.getVar( &cal_lut.K1[0][0]);
  var = gidLUT.getVar( "K2");
  var.getVar( &cal_lut.K2[0][0][0]);
  var = gidLUT.getVar( "K3_coef");
  var.getVar( &cal_lut.K3_coef[0][0][0]);
  var = gidLUT.getVar( "K4_coef");
  var.getVar( &cal_lut.K4_coef[0][0][0]);
  var = gidLUT.getVar( "K5");
  var.getVar( &cal_lut.K5[0][0]);
  
  return 0;
}

// float  K1K2[nib][msdim]: absolute and temporal gain
// float  K3_coef[nib][tempdim][tcdim]: temperature correction
// float  K4_coef[nib][msdim][rvsdim]: RVS correction
// double K5[nib][nldim]: linearity correction

int make_oci_gains( uint32_t nib, uint32_t banddim, uint16_t iyr, uint16_t idom,
                    double stime, double K2t[NTIMES], cal_lut_struct& cal_lut,
                    float **gmat, gains_struct& gains) {

  for (size_t i=0; i<6; i++) gains.ldims[i] = cal_lut.ldims[i];

  uint16_t timedim = gains.ldims[0];
  uint16_t tempdim = gains.ldims[1];
  uint16_t tcdim = gains.ldims[2];
  uint16_t rvsdim = gains.ldims[3];
  uint16_t nldim = gains.ldims[4];
  uint16_t msdim = gains.ldims[5];

  size_t dims[3];
  
  float **K1K2 = new float *[nib];
  K1K2[0] = new float[nib*msdim];
  for (size_t i=1; i<nib; i++) K1K2[i] = K1K2[i-1] + msdim;
  gains.K1K2 = K1K2;

  dims[0] = nib; dims[1] = tempdim; dims[2] = tcdim;
  float ***K3_coef = make3dT<float>(dims);
  gains.K3_coef = K3_coef;
  dims[0] = nib; dims[1] = msdim; dims[2] = rvsdim;
  float ***K4_coef = make3dT<float>(dims);
  gains.K4_coef = K4_coef;

  double **K5 = new double *[nib];
  K5[0] = new double[nib*nldim];
  for (size_t i=1; i<nib; i++) K5[i] = K5[i-1] + nldim;
  gains.K5 = K5;
  
  // Mirror-side dependent gains
  double *K2 = new double[banddim];
  for (size_t ms=0; ms<msdim; ms++) {

    // Get temporal gain and combine with absolute gain
    double d2 = jday(iyr, 1, idom) - 2451545 + stime/864.0;

    size_t kd=0;
    for (size_t j=NTIMES-1; j>=0; j--) {
      if ( d2 > K2t[j]) {
        kd = j;
        break;
      }
    }
    if ( kd < (size_t) (timedim-1)) {
      double ff = (d2 - K2t[kd]) / (K2t[kd+1] - K2t[kd]);
      for (size_t j=0; j<banddim; j++)
        K2[j] = cal_lut.K2[j][ms][kd]*(1.0-ff) + cal_lut.K2[j][ms][kd+1]*ff;
    } else {
      for (size_t j=0; j<banddim; j++) K2[j] = cal_lut.K2[j][ms][kd];
    }
    for (size_t j=0; j<banddim; j++) K2[j] *= cal_lut.K1[j][ms];
    for (size_t i=0; i<nib; i++) {
      gains.K1K2[i][ms] = 0;
      for (size_t j=0; j<banddim; j++)
        gains.K1K2[i][ms] += gmat[j][i] * cal_lut.K1[j][ms];
    }

    // Generate RVS coefficents
    for (size_t i=0; i<nib; i++) {
      for (size_t k=0; k<rvsdim; k++) {
        gains.K4_coef[i][ms][k] = 0;
        for (size_t j=0; j<banddim; j++)
          gains.K4_coef[i][ms][k] += gmat[j][i] * cal_lut.K4_coef[j][ms][k];
      }
    }
  }
  delete [] K2;
  
  // Generate temperature coefficients
  for (size_t i=0; i<nib; i++) {
    for (size_t k=0; k<tcdim; k++) {
      for (size_t l=0; l<tempdim; l++) {
        gains.K3_coef[i][l][k] = 0;
        for (size_t j=0; j<banddim; j++)
          gains.K3_coef[i][l][k] += gmat[j][i] * cal_lut.K3_coef[j][l][k];
      }
    }
  }

  // Generate linearity coefficients
  for (size_t i=0; i<nib; i++) {
    for (size_t k=0; k<nldim; k++) {
      gains.K5[i][k] = 0;
      for (size_t j=0; j<banddim; j++)
        gains.K5[i][k] += gmat[j][i] * cal_lut.K5[j][k];
    }
  }
  /*
  for (size_t i=0; i<nib; i++) {
    for (size_t j=0; j<banddim; j++)
      cout << "2: " << i << " " << j << " " << gmat[j][i] << endl;
  }
  */
  return 0;
}


template <typename T>
int get_oci_dark( size_t iscn, uint32_t nscan, uint8_t *hside, uint16_t ndsc,
                  uint16_t nskp, int16_t iags, int16_t iagd, uint32_t ntaps,
                  int16_t *jagg, uint32_t dfill, int16_t ndc, T ***dark,
                  uint32_t nib, float *dc, int16_t& iret) {

  // Program to generate dark corrections for OCI data by averaging the
  // dark collect data and correcting for bit shift/truncation if necessary

  // Determine number of bands per tap for hyperspectral data
  int16_t *nbndt = new int16_t[ntaps];

  if (ntaps == 16) {
    // hyperspectral bands
    for (size_t i=0; i<ntaps; i++)
      if ( jagg[i] > 0) nbndt[i] = 32 / jagg[i]; else nbndt[i] = 0;
  } else {
    for (size_t i=0; i<ntaps; i++) nbndt[i] = 9;
  }
  int16_t nbnd = 0;
  for (size_t i=0; i<ntaps; i++) nbnd += nbndt[i];

  // Select data for HAM side and determine scan indices
  int32_t *kh = new int32_t[nscan];
  int32_t nkh=0;
  for (size_t i=0; i<nscan; i++) {
    if ( hside[i] == hside[iscn]) {
      kh[i] = (int32_t) i;
      nkh++;
    } else {
      kh[i] = -1;
    }
  }

  int32_t js=0;
  for (size_t i=0; i<nscan; i++) {
    if ( kh[i] == (int32_t) iscn) {
      js = (int32_t) i;
      break;
    }
  }

  // Check for valid dark collect data within specified range
  uint16_t ndscl = ndsc;
  bool valid_dark_found = false;
  int32_t is1=js, is2=js;

  while (!valid_dark_found && ndscl <= nkh) {
    if ( ndsc > 1) {
      is1 = js - ndsc/2;
      is2 = js + ndsc/2;
      // Check for start or end of granule
      if (is1 < 0) {
        is1 = 0;
        is2 = ndsc - 1;
      }
      if (is2 >= nkh) {
        is1 = nkh - ndsc;
        is2 = nkh - 1;
      }
    }

    // If no valid dark data, expand scan range
    for (size_t i=is1; i<=(size_t) is2; i++) {
      for (size_t j=nskp; j<(size_t) ndc; j++) {
        if ( dark[kh[i]][0][j] != dfill) {
          valid_dark_found = true;
          break;
        }
      }
    }
    if ( !valid_dark_found) ndscl += 2;
  }

  if ( !valid_dark_found) {
    iret =-1;
    return 0;
  }
  
  // Loop through taps and compute dark correction
  int16_t ibnd=0;
  for (size_t i=0; i<ntaps; i++) {
    if ( jagg[i] > 0) {
      float ddiv = 1.0;
      float doff = 0.0;
      if ( iags*jagg[i] > 4) {
        ddiv = iagd * jagg[i] / 4.0;
        doff = (ddiv-1) / (2*ddiv);
      }

      for (size_t j=0; j<(size_t) nbndt[i]; j++) {
        float sum = 0.0;
        int nv = 0;
        for (size_t k=kh[is1]; k<=(size_t) kh[is2]; k++) {
          for (size_t l=nskp; l<(size_t) ndc; l++) {
            if (dark[k][ibnd+j][l] != dfill) {
              sum += dark[k][ibnd+j][l];
              nv++;
            }
          }
        }
        //        dc[ibnd+j] = ((sum/(ndc-nskp))/ddiv) - doff;
        dc[ibnd+j] = ((sum/(nv))/ddiv) - doff;
      } // j loop
    } // if ( 
    ibnd += nbndt[i];
  } // i loop

  delete [] kh;
  delete [] nbndt;

  iret = 0;
  if (ndscl > ndsc) iret = 1;
  
  return 0;
}


int get_oci_temp_corr( uint32_t nib, gains_struct gains, float K3T[NTEMPS],
                       float *caltemps, uint32_t nscan, float *k3) {

  uint16_t tempdim = gains.ldims[1];
  uint16_t tcdim = gains.ldims[2];

  for (size_t i=0; i<nib; i++) k3[i] = 1.0;

  for (size_t i=0; i<tempdim; i++) {
    float td = caltemps[i] - K3T[i];
      for (size_t j=0; j<tcdim; j++) {
        for (size_t k=0; k<nib; k++) {
          k3[k] -= gains.K3_coef[k][i][j] * powf(td, j+1);
        }
      }
  }

  return 0;
}

int get_oci_rvs_corr( uint32_t nib, uint16_t pdim, uint8_t hside,
                      gains_struct gains, double *theta, float **k4) {

  // Program to compute RVS correction from coefficients and scan angle

  uint16_t rvsdim = gains.ldims[3];

  for (size_t i=0; i<nib; i++)
    for (size_t j=0; j<pdim; j++)
      k4[i][j] = 1.0;

  for (size_t i=0; i<rvsdim; i++)
    for (size_t j=0; j<nib; j++)
      for (size_t k=0; k<pdim; k++)
        k4[j][k] += gains.K4_coef[j][hside][i] * powf(theta[k], j+1);
  
  return 0;
}

int get_oci_lin_corr( uint32_t nib, uint16_t pdim, gains_struct gains,
                      float K3T[NTEMPS], float *caltemps, float **dn,
                      float **k5) {

  // Program to compute OCI linearity correction from coefficients, dn and
  // temperatures 3rd-order polynomial of dn with 2nd-order temperature effects
  // for three temperatures.
  // This is currently a SWAG at the functional form and will be revised when
  // better known

  // Index for which temperatures to use for linearity; revisit this when known
  int16_t ibtmp[3] = {1,2,3};
  float td[3];
  for (size_t i=0; i<3; i++) td[i] = caltemps[ibtmp[i]] - K3T[ibtmp[i]];

  for (size_t i=0; i<pdim; i++) {
    // Zeroth-order correction
    for (size_t j=0; j<nib; j++) k5[j][i] = gains.K5[j][0];
    // First and second-order corrections are quadratic functions of
    // temperatures
    for (size_t k=1; k<=2; k++) {
      for (size_t l=0; l<nib; l++) {
        k5[l][i] += (gains.K5[l][k] + gains.K5[l][k+3]*powf(td[0], k) +
                     gains.K5[l][k+5]*powf(td[1], k) +
                     gains.K5[l][k+7]*powf(td[2], k))*powf(dn[l][i], k);
      }
      // Third-order correction
      for (size_t l=0; l<nib; l++)
        k5[l][i] += gains.K5[l][3]*powf(dn[l][i], 3);
    }
  }

  return 0;
}


int l1bFile::write_granule_metadata( std::string tstart, std::string tend,
                                     std::string l1b_name) {

  l1bfile->putAtt("time_coverage_start", tstart.c_str());
  l1bfile->putAtt("time_coverage_end", tend.c_str());

  // Write product file name
  l1bfile->putAtt("product_name", l1b_name.c_str());

  return 0;
}


int l1bFile::close() {

  try { 
    l1bfile->close();
  }
  catch ( NcException& e) {
    cout << e.what() << endl;
    cerr << "Failure closing: " + fileName << endl;
    exit(1);
  }
  
  return 0;
}


template <typename T>
T*** make3dT( size_t dims[3]) {

  T ***arr3d = new T **[dims[0]];

  arr3d[0] = new T*[dims[0]*dims[1]];
  arr3d[0][0] = new T[dims[0]*dims[1]*dims[2]];

  for (size_t i=1; i<dims[0]; i++) arr3d[i] = arr3d[i-1] + dims[1];

  for (size_t i=0; i<dims[0]; i++) {
    if ( i > 0) arr3d[i][0] = arr3d[i-1][0] + dims[1]*dims[2];
    for (size_t j=1; j<dims[1]; j++)
      arr3d[i][j] = arr3d[i][j-1] + dims[2];
  }

  return arr3d;
}

