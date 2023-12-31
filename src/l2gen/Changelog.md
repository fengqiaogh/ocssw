# l2gen Changelog



## 9.5.1
  - fixed sstref array size to fix memory ocerrun
  - fixed memory leaks in l1bgen_generic writer
  - forced l1 queue size to min of 5 if proc_sst is on

### Sources changed
cpl1rec.c
l1_generic_write.c
msl12_input.c
sstref.c


## 9.5.0

 - fixed a memory over run by the auto NULL termination when reading strings with the HDF5 API.  Also updated the test data for HICO.
 - added check for M13 band in VIIRS reader
 - corrected URLs to https
 - modified ice_mask Reynolds OISST read function to read CMC inputs as well
 - fixed a return value bug
 - calculate the h20 transmittance using LUTs
 - repalced the aw/bbw from sensorinfo.dat file with that from water_spectra.dat
  - added ability to read CMC SST reference data
 - mods for MISR - mostly default chl related
 - moved hdf.h out of the header
 - made sure OCTS only increments the day once
 - added a per pixel ancillary file reader
 - made landsat 7 work, somewhat
 - fixed the OCTS reader to cross a day boundary correctly
 - implemented sea surface nitrate code
 - added OCIS reader (simulated OCI data)
 - Kay Kilpatrick discovered a longstanding issue with the SST4diff usage in the SST flagging: 1. The two thresholds, SST4diff1,and SST4diff2,  should be negative, -0.8, and -1.0, respectively 2. The comparisons should be less than  eg. dSST< SST4diff1  and dSST < SSTdiff
 - fixed bug in reading tilt array for OCTS extracts
 - use scantime directly for terra electronics side effects correction if block
 - added checks for every use of the SST coefficient lut.  Made sure num lat bands is kept track of
 - Modified get_mld.cpp to handle multiple possible coordinate variable names, since they changed in latest Hycom MLD files.
 - added logic to read netCDF versions of the SST coefficient LUTs
 - added 3D products and other Paul Martinolich changes
 - added a dust correction algorihm for SST
 - minor restructuring of a CI test to allow use in get_habs_flags function
 - fixed hawkeye smile correction and radiance units conversion (yes, had to divide by 10)
 - add anc_aerosol based absorbing aerosol flag option
 - corrected MERIS cloud mask to match NOAA's implementation for CyAN
 - add support for GMAO aerosol ancillary data
 - Moved oel_util/libocssw to src/oc, updated dependent code, ignored $OCSSWROOT/include
 - set non-detect levels of CI to the minimum valid value in product.xml definition
 - gring changes for l1info
 - deleted the SOA code.
 - changes for compiling on gcc-8
 - modified to reduce calls to get_cldmask
 - cleaned up get_habs.c; made as consistent as I could to the NOAA python code
 - modified l3gen to return a 110 if that are no valid bins
 - HICO L1B reader now looks for latitude/longitude and latitudes/longitudes
 - made l3gen work with new bin class
 - fixed l1_viirs_nc leap second calculations
 - fixed compile warnings
 - fixed netcdf return value checking
 - modified reader to work with changes to the OCIA L1B file format
 - fixed undefined OCI function and cleaned up unused variables
 - Added is_l2 flag to the L1 struct, a SeaBASS reader, and various small modifications for each.
 - fixed a case where the last pixel or line is returned by lonlat2pixline.  It was off by one
 - Set NAVFAIL for any bad quality flag for Hawkeye L1A reader.
 - fixed it such that negative npp values are now returned as BAD_FLT
 - changed coefficent in get_pft_uitz.c in pico-chla zeta-max S5 0.575 -> 0.574
 - deleted unused code l2gen/fresnel.c and old hico code
 - l1a_hawkeye.c: Skip output of "dark rows".
 - l1a_hawkeye.c: Set "dark rows" to fill value.
 - added sentinel 3B and made subsensors for OLCI
 - l1a_hawkeye.c: read band co-registration, background subtraction and cropping parameters from LUT.
 - l1a_hawkeye.c: Read calibration params from NetCDF group "radiometric_calibration". Made netcdf calls less verbose.
 - added sentinel 3B and made subsensors for OLCI
 - l1a_hawkeye.c: implement band-to-band co-registration Note: this version has hard-coded offset parameters; need to implement in LUT.
 - l2gen: implement simple background subtraction for hawkeye.
 - l1a_hawkeye.c: applied formatting.
 - adding polarization correction to hawkeye L1A reader
 - hawkeye smile correction done
 - modified water_spectra code to read netCDF version of the table; added user-definable option for which water_spectra file to use
 - fixed MSI Sentinel 2A and 2B support
 - added valid input check for Lt_blue and Lt_green to complement existing Lt_red check to prevent pink clouds
 - Added OCI reader (but still a work in progress e.g. extracting only 62 bands)
 - trap bad nav so gsl interp fuction will not explode; removed some extraneous print statements
 - only process profile ancillary if provided
 - removed proc_cloud restrictions on anc_profile inputs
 - updates to support changes made to GMAO ancillary input file format
 - added VIIRS JPSS-1 support to the SDR file reader. Reader now supports the Block 2 file format. Reader wil guess the GEO filename if one is not given.
 - got rid of unused variables in /anc_acq.c
 - Modified get_l2prod_index.c and get_chl.c to support OCI (and made improvement wrt HAWKEYE)
 - Modified l1_io.c to support OCI
 - Eliminated "discarded-qualifiers" warnings introduced in commit d17cd22. ALL WARNINGS GONE!
 - Commented out all unused variables from l2gen source code; eliminated another 4778 warnings from full compile.
 - Removed unneeded variables from header files; moved a few into source code. Another 345 unused-variable warnings quashed.
 - Correctly defined constants with const in header files. Eliminated 3900 unused-variable warnings!
 - lonlat2pixline: accepts VIIRS GEO as only input file.
 - cleaned up merge from develop into hawkeye
 - Fix problem with MISR block number for multifile input
 - Add support for multi-file processing
 - src/l2gen/get_Kd.c - update Kd_lee using 2013 paper (add a bbw/bb term)
 - src/l2gen/qaa.c - correct equation in step 2 for turbid waters
 - Add support for gmpfile (MISR)
 - Fixed bug identified in olci spectrum
 - moved the replacement of the OCROOTS environment variables to a library.
 - modified attribute read to read fill, min/max val as int instead of short since the M13 band is now uint - easier to read them all as uint than test for ushort vs uint...
 - increased max NANG for polcor tables to 15 (was 7)
 - Added support for hawkeye
 - updated code from Paul Martinolich to bring QAA in line with IOCCG documentation
 - modified SeaWiFS sola/sena to be in the range of -180 to 180 for consistency with other missions and the product.xml definition
 - changed init for rho_cirrus in VIIRS readers to BAD_FLT
 - Fixed a bug on assignment of l1rec Lt
 - added sqr and Xarray to rho[tc]boxstats
 - moved meris header files inro the meris library source dir.

### Source Changed

 * aerosol.c
 * alloc_l1.c
 * alloc_l2.c
 * anc_acq.c
 * anc_acq.h
 * aph.c
 * atmocor1.c
 * atmocor2.c
 * bioOptBandShift.c
 * brdf.c
 * calcite.c
 * calc_par.c
 * calfile_utils.c
 * carder.c
 * Changelog.md
 * chl.h
 * cloud_flag.c
 * CMakeLists.txt
 * convert_band.c
 * convl12.c
 * cpl1rec.c
 * filehandle.h
 * filehandle_init.c
 * filter.c
 * filter.h
 * flags_iop.c
 * flags_iop.h
 * flags_sst.h
 * fluorescence.c
 * ftrim.c
 * fuzzy_func_v3.c
 * gas_trans.c
 * get_atrem_corl1.c
 * get_atrem_corl1v2.c
 * get_atrem_corl1v3.c
 * get_chl.c
 * get_f0.c
 * get_habs.c
 * get_ice_frac.c
 * get_Kd.c
 * get_l2prod_index.c
 * get_mld.cpp
 * get_nitrate.c
 * get_nitrate.h
 * get_npp.c
 * get_owmc.c
 * get_par.c
 * get_pft_hirata.c
 * get_pft_uitz.c
 * get_pml.c
 * get_poc.c
 * get_psd_ksm.c
 * get_qaa.c
 * get_rhos.c
 * get_rhown_nir.c
 * get_zno3.c
 * giop.c
 * giop.c.adaptive
 * glint.c
 * gringHelper.cpp
 * gringHelper.h
 * gsm.c
 * HOWTO_Add_a_sensor.txt
 * ice_mask.c
 * init_l1.c
 * init_l2.c
 * input_struc.h
 * l12_parms.h
 * l12_proto.h
 * l1_aci_hdf.c
 * l1_aci_hdf.h
 * l1a_hawkeye.c
 * l1a_hawkeye.h
 * l1a_osmi.c
 * l1a_seawifs.c
 * l1a_seawifs.h
 * l1_aviris.c
 * l1_aviris.h
 * l1_aviris_nc.c
 * l1b_misr.c
 * l1b_misr.h
 * l1b_oci.c
 * l1b_oci.h
 * l1b_ocis.c
 * l1b_ocis.h
 * l1b_viirs_nc.c
 * l1b_viirs_nc.h
 * l1c_msi.cpp
 * l1c_msi.h
 * l1c_msi_private.h
 * l1_czcs_hdf.c
 * l1_generic_write.c
 * l1_goci.c
 * l1_hdf_generic_read.c
 * l1_hdf_generic_read.h
 * l1_hdf_generic_write.c
 * l1_hico_h5.c
 * l1_hico_h5.h
 * l1_hmodis_hdf.h
 * l1_io.c
 * l1_l5tm.c
 * l1_l7etm.c
 * l1_meris_CC.c
 * l1_meris_CC.h
 * l1_meris_N1.c
 * l1_meris_N1.h
 * l1_mos_hdf.c
 * l1_mos_hdf.h
 * l1_nc_generic_read.c
 * l1_nc_generic_read.h
 * l1_ocia.c
 * l1_ocia.h
 * l1_ocm2_hdf.c
 * l1_ocm2_hdf.h
 * l1_ocmdb_hdf.c
 * l1_ocmdb_hdf.h
 * l1_ocm_hdf.c
 * l1_ocm_hdf.h
 * l1_octs_hdf.c
 * l1_octs_hdf.h
 * l1_olci.c
 * l1_olci.h
 * l1_oli.c
 * l1_oli.h
 * l1_osmi_hdf.h
 * l1_prism.c
 * l1_prism.h
 * l1_seabass.cpp
 * l1_seabass.h
 * l1_sgli.c
 * l1_sgli.h
 * l1_struc.h
 * l1subpix.c
 * l1_viirs_h5.c
 * l1_viirs_h5.h
 * l1_viirs_nc.c
 * l1_viirs_nc.h
 * l1_xcal_hdf.c
 * l1_xcal_hdf.h
 * l2binmatch_input.cpp
 * l2_flags.h
 * l2_generic.c
 * l2_hdf_generic.h
 * l2prod.h
 * l2prod_struc.h
 * l2_struc.h
 * las_iop.c
 * loadl1.c
 * lonlat2pixline.c
 * main_l1brsgen.c
 * main_l1det2det.c
 * main_l1info.cpp
 * main_l1mapgen.c
 * main_l2binmatch.cpp
 * main_l2gen.c
 * main_l3gen.cpp
 * main_lonlat2pixline.c
 * main_vcalmerge.c
 * met_cvt.c
 * met_cvt.h
 * mgiop.c
 * mph_flags.h
 * mscal_struc.h
 * msl12_input.c
 * myprod.c
 * niwa_iop.c
 * nrutil.h
 * owt.c
 * photic_depth.c
 * polcor.c
 * polcor_hawkeye.cpp
 * polcor_hawkeye.h
 * prodgen.c
 * qaa.c
 * qaa.h
 * raman.c
 * rawcal.h
 * rayleigh.c
 * rdatreminfo.c
 * read_l3bin.cpp
 * read_pixel_anc_file.cpp
 * read_pixel_anc_file.h
 * runcal.h
 * scene_meta.c
 * scene_meta.h
 * seawater.c
 * setanc.c
 * setflags.c
 * smi_climatology.c
 * smile.c
 * sssref.c
 * sst.c
 * sst_dsdi.c
 * sst_dsdi.h
 * sstref.c
 * swim.c
 * target_io.c
 * vcal.c
 * viirs_pxcvt.c
 * viirs_utls.c
 * virtual_constellation.c
 * water_spectra.c
 * xcal.c

## 9.3.0 - 2018-04-09

 - added support for OCI simulated data
 - added support for MSI (Sentinel 2A and 2B); Lansdat 5/7ETM
 - fixed VIIRS L1A file input to l1info
 - implemented new calcite CI algorithms
 - added checks for modis GEO file input into l2gen
 - modified l3gen to update the units attribute to match the output products
 - refactored opp product to npp
 - modified how the par product is read into the npp code; changed read to explicit float and added scale_factor/add_offset handling
 - added doi logic to l3gen to fix issue where products were showing incorrect doi as it was copied from the source L3b files
 - added logic to derive the default xcalfile so option need not exist in the msl12_defaults.par file
 - added NAVFAIL mask to par calculation
 - resolved uninitialized memory issue with l3gen
 - fixed OLCI input file given as one of the radiance files
 - freed memory that valgrind identified
 - implemented modifications provided by Paul Martinolich to allow l2gen to output 3D arrays for radiances and Rrs
 - changed olci reader to look in the directory of the ifile for all of the other files needed
 - renamed the modis aqua subframe correction files
 - renamed hmodis to modis and handle loading defaults files for subsensors
 - changed VIIRS to VIIRSN and implemented VIIRS J1 for SST
 - implemented new sensor info functions
 - VIIRS l1bptrs now returns reflectance instead of radiance for MOD bands.

### Source Changed
new code:
  * l1_l7etm.c
  * l1c_msi.cpp
  * getformatXML.cpp

modified code:
  * aerosol.c
  * atmocor2.c
  * bioOptBandShift.c
  * calcite.c
  * calc_par.c
  * calfile_utils.c
  * Changelog.md
  * CMakeLists.txt
  * filehandle_init.c
  * filter.c
  * fluorescence.c
  * get_atrem_corl1.c
  * get_chl.c
  * getformat.c
  * get_habs.c
  * get_l2prod_index.c
  * get_ndvi.c
  * get_npp.c
  * get_npp.h
  * get_par.c
  * init_l1.c
  * l12_proto.h
  * l1b_viirs_nc.c
  * l1_generic_write.c
  * l1_hmodis_hdf.c
  * l1_io.c
  * l1_meris_N1.c
  * l1_olci.c
  * l1_viirs_nc.c
  * l2_generic.c
  * l2prod.h
  * loadl1.c
  * main_l1brsgen.c
  * main_l1det2det.c
  * main_l1info.c
  * main_l2binmatch64.cpp
  * main_l2binmatch.cpp
  * main_l2gen.c
  * main_l3gen64.cpp
  * main_l3gen.cpp
  * main_vcalmerge.c
  * msl12_input.c
  * niwa_iop.c
  * owt.c
  * par_utils.c
  * photic_depth.c
  * pml_iop_tables.c
  * polcor.c
  * prodgen.c
  * raman.c
  * rayleigh.c
  * rdatreminfo.c
  * setanc.c
  * soa_sma.c
  * sst.c
  * sstref.c
  * virtual_constellation.c
  * water.c
  * xcal.c
removed/renamed code:
  * get_opp.c
  * opp.h

## 9.2.1 - 2017-11-21
VIIRS earth_sun_distance is calculated by esdist_() instead of being read from L1.

## 9.2.0 - 2017-11-02
Added ability to read new OLI enhanced metadata file.
normalized GCMD keywords and DOI code
added the MSI sensor
added geometry per band calculations if sensor supports it
die if a bad suite is given
revamped the generic L1B netCDF reader
added new product nKd
moved memory allocation from stack to heap
new GOCI slot assignment
updated NDVI product code
OLCI smile correction

## 9.1.0 - 2016-09-28
Add instrument SGLI

Fixed bug in sst code where nbvis was not being set
Fixed bug in l1info where second field in time output was allowed to be more than 59.999
Fixed allocation issue in calfile_utils.c and added initalization to l2_str in vcalmerge
added product chl_oci2

### Source Changed
  * sst.c
  * main_l1info.c
  * main_vcalmerge.c
  * calfile_utils.c
  * get_chl.c
  * l2prod.h
  * prodgen.c

## 9.0.1 - 2016-01-19
NBANDSIR bug fix for VIIRS L1B reader 

### Source Changed
  * l1b_viirs_nc.c

## 9.0.0 - 2016-01-12
Major reorganization of the internal structures.

### Source Changed
just about everything...

## 8.9.10 - 2015-11-12
Added support for GOCI "slot" product provided by Paul Martinolich of NRL
Modifications to OLCI reader to better mimic MERIS until OLCI specific tables are
generated - provided by Paul Martinolich
Update to rdsensorinfo.c to allow sensor_defaults.dat to optionally be indexed by 
wavelength...also from PM
Updated sst.c code for VIIRS support provided by Sue Walsh of RSMAS

### Source Changed
  *  alloc_l1.c
  *  convl12.c
  *  init_l1.c
  *  l1_goci.c
  *  l1_io.c
  *  l1_struc.h
  *  l1subpix.c
  *  l2_struc.h
  *  l2prod.h
  *  prodgen.c
  *  getformat.c
  *  l1_olci.c
  *  olci.h
  *  rdsensorinfo.c
  *  sst.c

## 8.9.9 - 2015-11-05

Added Rick Stumpf alogorithm for cloud flagging MERIS and MODIS
Use eval=2

Also added as products
use l2prod=flags_habs
    Flag = 0 Water
    Flag = 1 Cloud
    Flag = 2 Not Water
    Flag = 3 Cloud over Not Water
    
### Source Changed

  * get_habs.c
  * l12_parms.h
  * l12_proto.h
  * l2prod.h
  * mph_flags.h
  * msl12_input.c
  * prodgen.c
  * setflags.c

## 8.9.5 - 2015-11-01

Added OrbitNumber to l1info output
Added new flag, BOWTIEDEL, and function to set it

### Source Changed
   * l1_viirs_nc.c
   * main_l1info.c
   * 12_proto.h
   * l1b_viirs_nc.c
   * l2_flags.h
   * l1_viirs_h5.c

### Data Added
## 8.9.0 - 2015-10-20

Updated Atrem to do inline k-distribution method, error angle max geometry recalculation.
Added command line options, atrem_opt, atrem_full, atrem_geom, atrem_model. 
Added dynamic Atmospheric model using pixel latitude and scene day of year.
Added ability to use gas_opt (gaseous_transmission) with atrem_opt.
Updated sensor info to higher precision values of wavelengths
Normalized full calculation so that transmittances are 1.0 in blue spectrum

### Source Changed
   * get_atrem_corl1.c
   * atrem_corl1.h
   * msl12_input.c
   * input_struc.h
   * gas_trans.c
   * atmocor1.c
   * atrem_app_refl_plus_gas_removal_l2.f90

## 8.8.0 - 2015-09-13

Added HABS and KD (Harmful Algal Blooms) product
  (Kd_stumpf, CI_stumpf, MCI_stumpf, MPH_chl, MPH_flags)

### Source Changed
   * prodgen.c
   * l2_generic.c
   * l12_proto.h
   * get_Kd.c
   
 New code:
    * get_habs.c

## 8.7.0 - 2015-08-04

Added OLCI sensor for Sentinel-3

### Source Changed

  * filehandle.h (1 diff)
  * get_chl.c (1 diff)
  * get_l2prod_index.c (1 diff)
  * getformat.c (5 diffs)
  * l1_io.c (4 diffs)
  * include/sensorDefs.h (1 diff)
  * include/sensorInfo.h (4 diffs) 
 
 New code:
 
  * olci.h
  * l1_olci.c
  * l1_olci.h

### Data Added

run/data/olci

## 8.6.1 - 2015-08-03
Updated MS epsilon code and modisa eval aerosol LUTs for 4th-order poly.

### Source Changed
     modified code:
        aerosol.c

### Data Added
## 8.6.1 - 2015-07-31
Fixed bug in array allocation in gsm.c that prevented proper use with l3gen

### Source Changed
     modified code:
        gsm.c

### Data Added
## 8.6.0 - 2015-07-27
Added user contributed chlor_a l2 product (chl_abi) algorithm from P. Shanmugam

Shanmugam P. (2011), A new bio-optical algorithm for the remote sensing of algal blooms
in complex ocean waters, Journal of Geophysical Research, 116(C4), C04016.
doi:10.1029/2010JC006796

Added check in rayleigh table solz limit to prevent going outside table bounds (needed fro GOCI)

### Source Changed
     modified code:
        get_chl.c
        l12_proto.h
        l2prod.h
        prodgen.c
        rayliegh.c
 
## 8.5.0 - 2015-07-24
Added command line option to specify input PAR file for OPP calculation when
desired (e.g., when using l3gen)
    - Specify parfile=<inputfile>
    - If no file specified, attempt is made to calculate PAR  
    - File needs to be netcdf format with variables "par","lat", and "lon"
    
     modified code:
        get_opp.c
        opp.h
        msl12_input.c
        input_struc.h


## 8.4.0 - 2015-07-23
Added command line option to specfify f/Q file for brdf correction
    - Specify fqfile=<inputfile>
    - If no file specified, file is read from $OCDATAROOT/common/morel_fq.nc 
    - File needs to be netcdf format
    - morel_fq.nc created from morel_fq.dat using cre_morel_ncfile.py
    
     modified code:
        atmocor2.c
        brdf.c
        msl12_input.c
        get_rhown_nir.c
        giop.c
        input_struc.h
        l12_proto.h
       
    added code:
        cre_morel_ncfile.py


## 8.3.0 - 2015-06-30
Added PRISM sensor
    - Specify prm*.hdr file in command line ifile opt.  Will search for gain file to read.
    - Uses the first 44 bands (to 834 nm). 
    - Rayleigh files need to be updated. Using HICO files in the meantime
    - msl12_sensor_info.dat file contains pseudo values.  Needs to be updated.
made NDVI more flexible with band-centers 
modified combined 2b/3b calcite algorithm to use 3b only when inputs for 2b

### Source Changed
     modified code:
         ndvi.c
         calcite.c

## 8.2.0 - 2015-06-02

Added ORCA sensor from PACE mission
    - Used HICO parameters in msl12_sensor_info.dat file and Rayleigh files from HICO
Added AVIRIS sensor
    - Specify img.hdr file in command line ifile opt.  Will search for gain file to read.
    - Uses the first 44 bands (to 795 nm). 
    - Rayleigh files need to be updated. Using HICO files in the meantime
    - msl12_sensor_info.dat file contains pseudo values.  Needs to be updated.
    
### Source Changed
 modified code:
 
  * atmocor1.c (3 diffs)
  * filehandle.h (1 diff)
  * get_chl.c (1 diff)
  * get_l2prod_index.c (1 diff)
  * getformat.c (5 diffs)
  * l1_io.c (4 diffs)
  * msl12_input.c (1 diff)
  * raman.c (1 diff)
  * setanc.c (1 diff)
  * include/sensorDefs.h (1 diff)
  * include/sensorInfo.h (4 diffs) 
 
 New code:
 
  * aviris.h
  * l1_aviris.c
  * l1_aviris.h
  * l1_orca.c
  * l1_orca.h

## 8.1.4 - 2015-05-28

Increased solz threshold for par to 90 (needed for MODISA reprocessing 2014.0)
added support for Raman scattering correction to GIOP

### Source Changed
modified code:
  * get_par.c
  * atmocor2.c
  * giop.c
  * input_struc.h
  * l12_proto.h
  * l2_struc.h
  * msl12_input.c
new code:
  * raman.c

## 8.1.3 - 2015-05-20

Modifications to handle library reorganization;
minor code cleanup (dead wood removal);
updated calcite code to remove the 4/1.628 factor

### Source Changed
modified code:
  * atrem_corl1.h
  * brightness.c
  * calc_par.c
  * calcite.c
  * calfile_utils.c
  * filehandle.h
  * get_atrem_corl1.c
  * get_mld.c
  * get_qaa.c
  * get_zno3.c
  * getformat.c
  * input_struc.h
  * l12_parms.h
  * l12_proto.h
  * l12_seawifs.c
  * l1_aci_hdf.c
  * l1_czcs_hdf.c
  * l1_generic_write.c
  * l1_hdf_generic_read.c
  * l1_hdf_generic_write.c
  * l1_hmodis_hdf.c
  * l1_ocm2_hdf.c
  * l1_ocm_hdf.c
  * l1_ocmdb_hdf.c
  * l1_oli.c
  * l1_xcal_hdf.c
  * l1a_osmi.c
  * l1a_seawifs.c
  * l2_generic.c
  * l2_generic.h
  * l2_hdf_generic.h
  * main_l1info.c
  * msl12_input.c
  * polcor.c
  * sssref.c
  * sst.c
  * sstref.c
  * test_l1brsgen_config.cmake
  * test_l1mapgen_config.cmake
  * test_l3gen_config.cmake
  * viirs_utls.c

## 8.1.2 - 2015-05-04

Merged HICO HS and Multi-spectral data directories.  To run HICO HS copy msl12_sensor_info_HICOHS.dat
to msl12_sensor_info.dat and likewise to run MS copy msl12_sensor_info_HICOMS.dat to
msl12_sensor_info.dat in the $(OCSSW)/run/data/hico directory.

The new HICO MS uses HICO HS data files and wavelengths which are slightly different from the old HICO MS wavelengths so results will be slightly different.

Old HICO MS directory is in $(OCSSW)/run/data/hicoms.
    
## 8.1.1 - 2015-04-24
l2gen: 
Added ability to output all four aerosol models and two model ratios. 

Corrected many issues with the initial implementation of the Ahmad MSEPS
aerosol selection scheme.

Corrected issue with the calcite-specific-backscatter used in the calcite_2b
and calcite_3b products, and made explicit the inadvertently elevated thresholds
on background and switching values.

### Source Changed
	
modified code:
  * aerosol.c
  * alloc_l2.c
  * atmocor2.c
  * calcite.c
  * init_l2.c
  * l12_proto.h
  * l2_struc.h
  * prodgen.c
  * rayleigh.c

	
## 8.1.0 - 2015-04-07
l2gen: Implemented Bo-Cai Gao's Atmospheric Removal (Atrem) routines for HICO HS
(includes gases O2, O3, N2O, NO2, CO, CO2, and CH4 as well as H2O).  At the moment,
H2O is always removed.  Other gases to remove are specified in file msl12_atrem_info.dat 
in run/data/hico directory. Use command line option gas_opt=16.  Still in testing phase, 
thus it's an "undocumented program feature".

Also, made modifications so that HICO HS parameters are initialized properly.

### Source Changed
new code:
  * atrem_app_refl_plus_gas_removal_l2.f90
  * get_atrem_corl1.c
  * rdatreminfo.c
  * atrem_corl1.h
  * COMMONS_INC.f
  * tpvmr_init.f
  * cubeio.f90
  * bndprms.f
    
modified code:
  * aerosol.c
  * atmocor1.c
  * atmocor2.c
  * msl12_input.c
  * l12_parms.h
  * CMakeLists.txt
    
## 8.0.2 - 2015-03-29
l2gen: implement awhite interpolation within whitecaps.c, superceding the need
for awhite coefficients in sensorinfo.dat.

### Source Changed
modified code
  * whitecaps.c, compute awhite, don't read it sensorinfo.dat.

## 8.0.2 - 2015-03-17
l2binmatch bug fix
Generalized flag percentages for l2gen and added new flag percentages for: flags_sst,
flags_sst4, flags_giop, flags_qaa, flags_carder, flags_niwa

### Source Changed
deleted code
  * main_l2binmatch.c
  * msl2val_hdf.c
  * msl2val_struc.c
  * msl2val_struc.h
modified code
  * main_l2binmatch.cpp, corrected detnum and mside output
  * l2_generic.c, new flag cnt arrays
  * l3mapgen.c, removed an unneeded fudge factor
  * nccmp.c, added new option "NaNs are equal"

## 8.0.1 - 2015-03-10

Added new, optional correction for cirrus; bug fix for l2binmatch;
Updated OCI chl ramp range; updated test suite

### Source Changed
  * test_l2gen_config.cmake
  * test_l1info_config.cmake
  * test_l1bgen_generic_config.cmake
  * CMakeLists.txt
  * msl12_input.c
  * l1_oli.c
  * l1_oli.h
  * atmocor2.c
  * l2binmatch_input.cpp
  * getformat.c
  * bioOptBandShift.c
  * get_chl.c
  * input_struc.h
  * get_zno3.c
  * calfile_utils.c
  * main_l2binmatch.cpp


## 8.0.0 - 2015-01-09

Change to aerosol code to incorporate the multi-scattering epsilon tables for
Rayleigh, replace modis reader with new version - deprecate use of modis l1a
subset files (touches a LOT of files), add net (ocean) primary productivity
algorithms (opp), add support for new VIIRS L1A file, redesigned l2binmatch
program, implement new biooptical bandshift code.

For Ocean Primary Productivity there are three options:

    opp_befa   - Behrenfeld-Falkowski
    opp_eppley - Same as befa, but modifies pb_opt function after Eppley
    opp_cbmp2  - Primary Productivity using a chl:Carbon ratio using nine separate
                 wavelengths
  
### Source Changed
deleted code
  * get_aer_index.c
  * aerosol_ahmad.c
  * l1_modis_hdf.c
  * l1_modis_hdf.h
  * dem_s.fin

new code
  * test_l1mapgen_config.cmake
  * test_l2gen_config.cmake
  * test_l3gen_config.cmake
  * test_l1bgen_generic_config.cmake
  * test_l1brsgen_config.cmake
  * l1_viirs_nc.c
  * test_l1info_config.cmake
  * l1_viirs_nc.h
  * l2binmatch_input.h
  * freearray.h
  * bioOptBandShift.c
  * get_zno3.c
  * main_l2binmatch.cpp
  * get_opp.c
  * opp.h
  * l2binmatch_input.cpp

modified code
  * input_struc.h
  * calfile_utils.h
  * atmocor2.c
  * filehandle.h
  * main_l1brsgen.c
  * get_tricho.c
  * ipar.c
  * giop.c
  * polcor.c
  * prodgen.c
  * par_utils.c
  * giop.h
  * get_niwa_iop.c
  * rayleigh.c
  * l1a_osmi.c
  * aerosol.c
  * get_es.c
  * l12_seawifs.c
  * l1_generic_write.c
  * vcal.c
  * init_l1.c
  * l1_oli.c
  * l1_czcs_hdf.c
  * get_par.c
  * l1_oli.h
  * l1_viirs_h5.c
  * las_iop.c
  * main_l1info.c
  * msl2val_struc.c
  * lonlat2pixline.c
  * convert_band.c
  * msl2val_struc.h
  * main_vcalmerge.c
  * l1_meris_N1.c
  * mumm.c
  * alloc_l1.c
  * CMakeLists.txt
  * get_Kd.c
  * init_l2.c
  * get_qaa.c
  * loadl1.c
  * l1_ocm2_hdf.c
  * get_smoke.c
  * gas_trans.c
  * soa_sma.c
  * convl12.c
  * main_l2binmatch.c
  * get_ndvi.c
  * carder.c
  * read_l3bin.cpp
  * convl21.c
  * l12_proto.h
  * owt.c
  * l1_aci_hdf.c
  * l1_goci.c
  * brdf.c
  * msl12_input.c
  * l1_meris_CC.c
  * get_l2prod_index.c
  * get_depth.c
  * l1_hdf_generic_write.c
  * get_chl.c
  * niwa_iop.c
  * l1_mos_hdf.c
  * target_io.c
  * cdom_morel.c
  * l1_hdf_generic_read.c
  * pml_iop_tables.c
  * main_l1mapgen.c
  * get_rhos.c
  * main_l2gen.c
  * atmocor1_land.c
  * l1a_seawifs.c
  * cdom_mannino.c
  * getl1rec.c
  * l1subpix.c
  * l2prod.h
  * l1_hico_h5.c
  * seawater_get.c
  * sst.c
  * seawater.c
  * main_l3gen.cpp
  * l1_imgscale.c
  * virtual_constellation.c
  * get_pml.c
  * l1_hmodis_hdf.c
  * l1_hmodis_hdf.h
  * main_lonlat2pixline.c
  * photic_depth.c
  * gsm.c
  * xcal.c
  * xcal.h
  * main_l1det2det.c
  * sstref.c
  * read_l3bin.h
  * setflags.c
  * l2_generic.c
  * myprod.c
  * filter.c
  * turbid.c
  * get_toa_refl.c
  * l1_ocmdb_hdf.c
  * l1_struc.h
  * l1_io.c
  * flags_iop.c
  * calcite.c
  * flags_iop.h
  * swim.c
  * l1_octs_hdf.c
  * atmocor1.c
  * cpl1rec.c
  * main_l1bgen.c
  * get_poc.c
  * getformat.c
  * l1_ocm_hdf.c
  * l1_xcal_hdf.c
  * get_owmc.c
  * l12_parms.h
  * fluorescence.c
  * calfile_utils.c


## 7.0.2 - 2014-09

Minor fixes for meta-data and formatting issues related to reprocessing. Added
support for OLI GEO file.  Added test for atmospheric correction MS to SS
root-finding algorithm, to prevent segfault on one VIIRS glint pixel.

## 7.0.1 - 2014-08

Assorted changes to defaults and minor algorithm mods in preparation for reprocessing.
Support for OLI GEO file.

## 7.0.0 - 2013 to 2014

Support for netcdf4 output.

## 6.7.1 - 2014-04-02

New common sensor table for PIC.  
New MLR coastal cdom products.
OCTS PAR support (PAR2)


## 6.7.0 - 2013-09-12

coincide with the release of SeaDAS 7.0.1

### Source Changed
  * l12_parms.h, version change.
  * l1_meris_N1.c, added ability to read pixel extracted N1 MERIS files.
  * swfinc/sensorDefs.h, swfinc/sensorInfo.h, filehandle.h, get_l2prod_index.c,
  getformat.c, l1_io.c, added landsat reader
  * atmocor2.c, get_rhown_nir.c, added red wavelength for landsat
  * l1_oli.c, l1_oli.h, landsat reader
  * bug fix for HMODIS lon/lat optimized read and added optimizatin to VIIRS, 
  GOCI, HICO, and OLI
        l1_goci.c
        l1_goci.h
        l1_hico_h5.c
        l1_hico_h5.h
        l1_io.c
        l1_oli.c
        l1_oli.h
        l1_viirs_h5.c
        l1_viirs_h5.h


## 6.6.8 - 2013-07-11
  * get_dem_height.f converted to c
  * CMakeLists.txt - removed all references to HDFEOS and added proj4
  * getformat.c - changed how GOCI L1B files are identified
  * l2_proto.h - added varable names to sunangs_ prototype
  * l1_goci.c - removed HDFEOS and use proj4, use solz/sola calculation 
        from libnav
  * sunangs.f, moved this file from libczcs to libnav
  * libgoci/goci.c - removed HDFEOS and use proj4
  * deleted tsm_clark, poc_clark, chl_clark products
  * l1_viirs_h5.c, changed print statements so lonlat2pixline is quiet
  * got rid of almost all compiler warnings.

  * inc/cdfinc/cdfdist.h
  * inc/cdfinc/cdflib.h
  * inc/cdfinc/windoz.h
  * inc/swfinc/config.h
  * auto_qc/fmt_check/chk_sds.c
  * geogen_hico/CMakeLists.txt
  * l2gen/get_dem_height.f
  * libswfnav/navtlm.f
  * libswfnav/read_analog.f
  * libswfnav/read_discrete.f
  * libswfnav/read_double.f
  * libswfnav/read_float.f
  * libswfnav/read_gps.f
  * libswfnav/read_long.f
  * libswfnav/read_short.f
  * libswfnav/scpar.f
  * orbgen_aquarius/read_float.f
  * orbgen_aquarius/read_long.f
  * orbgen_aquarius/sacdcomp.f
  * changed macro linux to LINUX to be conisistant with other code.


## 6.6.8 - 2013-06-24
  * l1_modis_hdf.c,  set NAVWARN when MODIS > 1 degree off nominal pointing
  * l1_hmodis_hdf.c, set NAVWARN when MODIS > 1 degree off nominal pointing


## 6.6.8 - 2013-04-26

### Source Changed
  * l12_parms.h, version change
  * l1_viirs_h5.c, Implement time-aggregated VIIRS granule processing and
  cability to use ellipsoidal geolocation (GMODO) properly.
  * main_l1info.c, enable it to work for VIIRS granules made using predicted 
  F-LUT tables.
  * get_l2prod_index.c, changed wavelength designation for a few products that
    will be written to the product XML file.
  * main_l2gen.c, product XML can now be output without specifying an ifile.
  * l1_goci.c, add GOCI support.
  * getformat.c, add GOCI support.
  * get_l2prod_index.c, add GOCI support.
  * get_chl.c, add GOCI support.
  * get_l2prod_index.c, add GOCI support.
  * l2prod.h, add GOCI support.
  * prodgen.c, add GOCI support.


## 6.6.7 - 2013-04-01

### Source Changed
  * input_struc.h, add picfile.
  * aerosol.c, comment-out loading messages.
  * calcite.c, replace input file discovery with parameter picfile.
  * msl12_input.c, add picfile.


## 6.6.7 2013-03-28 

### Source Changed
  * l12_parms.h, version change
  * l2prod.h, change chl_err to chl_owterr
  * prodgen.c, change chl_err to chl_owterr
  * get_l2prod_index.c, change chl_err to chl_owterr, change chl_error to chlor_a_owterr
  * owt.c, change chl_err to chl_owterr, change weighting
  * get_chl.c, add support for OCM2 in chl_oci.
  * main_vcalmerge.c, reduce max size of file reached before bailing.
  * main_l1ingo.c, fix handling of missing lon/lat values.


## 6.6.6 - 2013-03-21

### Source Changed
  * l12_parms.h, version change
  * l12_proto.h, proto for read_target_l3() and wrapper for C++.
  * read_l3bin.cpp, new function to read L3 bin files.
  * read_l3bin.h, interface for read_l3bin().
  * target_io.c, add read_target_l3() and support utilities.
  * getformat.c, fix segfault error when handling unknown filetypes.
  * main_l2gen.c, add option to id and read target data from L3 bin.
  Old flat binary capability retained (for now).  
  * l1_io.c, fix test for navfail.


## 6.6.5 - 2013-03-08 

### Source Changed
  * l12_parms.h, version change
  * convert_band.c, fix MERIS conversion for 560-555.
  * get_chl.c, informational message on chl_oci.
  * rdsensorinfo.c, convert to dynamic allocation, clean-up.
  * main_l1mapgen.c, main_l1brsgen.c, open file sooner, so bindex_set is 
  initialized by l1_io and calls to rdsensorinfo and bindex_set need not be 
  repeated here. 
  * get_Kd.c, new Kd_jamet product.
  * get_l2prod_index.c, new Kd_jamet product.
  * l2prod.h, new Kd_jamet product.
  * prodgen.c, new Kd_jamet product.
  * sssref.c, new sss climatology.
  * main_l1mapgen.c, fixed bug when searching for lat/lon extent of SeaWiFS
    file.  Changed to readl1_lonlat in search to speed things up.
  * l1_io.c, readl1_lonlat sets a few more variables in l1str so 
    scene_meta_put will work properly.

## 6.6.4 - 2013-02-29 

### Source Changed
  * l12_parms.h, version change, increment SENSOR_NUM, add HICO sensor
  * l1_io.c, add HICO L1B FMT
  * filehandle.h, add HICO L1B FMT
  * get_l2prod_index.c, add default chlor_a for HICO
  * get_chl.c, add default chlor_a for HICO
  * l1_hico_h5.h, header file for new HICO L1B I/O
  * l1_hico_h5.c, new HICO L1B I/O
  * getformat.c, addd HICO L1B
  * CMakeLists.txt, add l1_hico_h5.c, swim.c and elev.c to add_executable entries.
  Also added levmar, lapack abd blass libs to l2gen and l3gen fr swim. 
  * l1_struc.h, l2_struc.h, convl12.c, alloc_l1.c, init_l1.c, added elev for 
  elevation to structure.
  * msl12_input.c, input_struct.h, addedelevfile and elev_auxfile as input 
  parameters
  * get_l2prod_index.c, l2prod.h, prodgen.c, added products a_swim, bb_swim, 
  adg_swim, aph_swim, bbp_swim, and elev
  * loadl1.c, added loading of elevation to the l1 and l2 structures.
  * l12_proto.h, added functions for swim and elev products.


## 6.6.3 - 2013-02-18

### Source Changed
  * l1_meris_N1.c, fix bug that was causing FR ingest to puke.

## 6.6.3 - 2013-02-03 

### Source Changed
  * l12_parms.h, version change.
  * msl12_input.c, update GIOP defaults for giop_adg_opt and giop_rrs_opt 
  as per Werdell at al. 2013, add Jeremy's mgiop changes (multi-giop). Add iterate and rrsdiff_max
  options for mgiop.
  * get_l2prod_index.c, add mgiop prods, pic2, BSi.
  * giop.c, add functions to retrieve internal giop pointers for mgiop, add init reset for multiple calls..
  * giop.h, add functions to retrieve internal giop pointers for mgiop.
  * input_struc.h, add giop_iterate (?), and giop_rrsdiff_max.
  * prodgen.c, add mgiop products, pic2, BSi.
  * l12_proto.h, mgiop proto, get_bsi() proto.
  * l2prod.h, add mgiop products, pic2, BSi.
  * mgiop.c, new function to return multi-giop products.
  * calcite.c, add new pic2 product algorithm (combined alg. with lat-dep bbstar).
  * get_bsi.c, new biogenic silica algorithm.
  * get_par.c, add meris PAR support.
  * calc_par.f, add meris PAR support.


## 6.6.1 - 2013-01-02 

### Source Changed
  * l12_parms.h, version change.
  * qaa.c, new QAA V6 from Paul Martinolich.
  * get_qaa.c, name changes for V6.
  * qaa.h, name changes for V6.
  * l1subpix.c, fix indexing error in sw_a, sw_bb, sw_n subsetting.
  * xcal.h, increase XTNORDER to 6.
  * xcal.c, increase max noumber of coefficients to 6 (5th order).

## 6.6.0  2012-11-21

Adding MERIS and VIIRS support for PIC, VIIRS support for PAR, changing VIIRS
standard product suite to include pic, poc, par.  Improved handling of missing
values in VIIRS.  Improved handling of missing values in MODIS (swir
detectors). Support for HYCOM ancillary salinity.

### Source Changed
  * l12_parms.h, version change.
  * l1_viirs_h5.c, set missing values (e.g., bowtie deletion) to BAD_FLT.
  * calcite.c, add VIIRS and MERIS support.
  * calc_par.f, add VIIRS support.
  * get_par.c, add VIIRS functions.
  * get_rhos.c, set to BAD_FLT when input radiance is negative, else land blows-up.
  * l1_hmodis_hdf.c, re-enable averaging for hi-res interpolation with missing dets. Set
  missing or otherwise flagged radiances to BAD_FLT (instead of 1000).
  * l1_modis_hdf.c,  Set missing or otherwise flagged radiances to BAD_FLT (instead of 1000).
  * sssref.c, new get_hycom() function.
  * loadl1.c, don't apply xcal or vcal if Lt is negative or 1000.
  * sst.c, for cloud non-uniformity test, no need now to adjust threshold on Lt saturation
         based on vicarious gain.  Just use 999.  This changes quality levels in some 
         cases, as saturated radiances where getting into the test previously.


## 6.5.9 - 2012-09-19

Adding fields in l1/l2 rec to store band, salinity, and temperature-dependent
pure-seawater index of refraction, aw, and bbw.  Added products nw_nnn, aw_nnn,
bbw_nnn.  Add support for VIIRS recalibration and psuedo L1A.

### Source Changed
  * l12_parms.h, version change.
  * l1subpix.c, subset sw fields.
  * convl12.c, copy sw field pointers to l2 record.
  * seawater.c, new function for pure seawater values.
  * seawater_get.c, new functions to set and retrieve scan-specific seawater iops.
  * l2prod.h, add aw, bbw, nw prods.
  * alloc_l1.c, add space for sw fields.
  * l1_struc.h, add sw field pointers.
  * l2_struc.h, add sw field pointers.
  * prodgen.c, add aw, bbw, nw prods.
  * loadl1.c, compute and store nw, bbw, aw in l1rec and seawater_set().
  * main_l3gen.cpp, transfer sw fields.
  * get_l2prod_index.c, add aw, bbw, nw prods.
    use clo xml functions to escape illegal xml characters.
    changed products La, TLg, Lr, L_q, L_u and tLf to all wavelengths 
  * getl1rec.c, call seawater_set() to reset scan-specific iops.
  * init_l1.c, init sw fields.
  * l12_proto.h, add seawater() protos.
  * l1_viirs_h5.c,
  * l1_viirs_h5.h,
  * CMakeLists.txt, add seawater.c, seawater_get.c
  * get_Kd.c, update Zsd_gbr coeffs.


## 6.5.8 - 2012-08-26

Fixes for changes to support OCR-VC.  By-passed call to get_iops() was causing
KDLEE to fail.

### Source Changed
  * loadl1.c, don't call atmocor1 if ocrvc_opt set (no Rayleigh tables for OCRVC).
  * main_l3gen.cpp, remove call to force proc_ocean=0.
  * msl12_input.c, added atmocor to l1mapgen's help, changed absaer_opt default
    to -1.
  * main_l1mapgen.c, normal exit code is 0 instead of percent pixel coverage.
  * main_l2mapgen.c, got rid of the -1 exit value.


## 6.5.8 - 2012-08-10

Changes to support OCR-VC.

revisions
  * r6908, r6909, r6914

### Source Changed
  * l12_parms.h, version change, new sensor OCRVC.
  * input_struc.h, new inputs vcnnfile and ocrvc_opt.
  * l12_proto.h, new function virtual_constellation(). 
  * l2prod.h, new product entries Rrs_vc and chl_vc.
  * chl.h, new header to port capture standard chlorophyll algorithm defaults.
  * filehandle.h, add ocrvc sensor and ocrvc_opt field.
  * getformat.c, change sensor name to OCRVC for l3 input when ocrvc_opt is set.
  * get_chl.c, moved header info to chl.h, add default OC3 for OCRVC.
  * filehandle_init.c, init ocrvc_opt to 0.
  * main_l3gen.cpp, disable atmospheric correction (so no Rayleigh files are needed), 
  transfer ocrvc_opt, read Rrs_vc (if present) and operwrite Rrs field, copy sensor
  name for input filehandle (overwrite previous meta-data copy).
  * get_l2prod_index.c, ad Rrs_vc and chl_vc products, and chl default for OCRVC.
  * msl12_input.c, add support for ocrvc_opt and vcnnfile.
  * prodgen.c, add support for Rrs_vc and chl_vc.
  * virtual_constellation.c, new functions to convert Rrs to Rrs_vc and a chl_vc 
  algorithm that takes Rrs_vc as input.
  * CMakeLists.txt, add virtual_constellation.c.


## 6.5.7 - 2012-02-01

Standardized the separation between visible and nir/swir bands and IR bands, and 
extended vis range to 720nm.

### Source Changed
  * l12_parms.h, version change, add MAXWAVE_VIS and MAXWAVE_IR
  * prodlist.c, get nwaveVIS from rdsensorinfo().
  * get_qaa.c, get nbands from rdsensorinfo().
  * rdsensorinfo.c, add logic to return number of vis bands.
  * get_rhown_nir.c, use MAXWAVE_VIS to select red vs nir function.
  * convl21.c, get nwvis from rdsensorinfo().
  * gsm.c, get nbands from rdsensorinfo(). 
  * get_l2prod_index.c, use MINWAVE_IR to separate IR bands.  Changed CAT_Lt
    definition to all wavelengths and deleted CAT_Ltir.  Changed the product 
    definitions so "chlor_a" is written to the XML file.
  * atmocor2.c, get nwvis from rdsensorinfo().
  * prodgen.c, changed the definiton of product Lt so the XML will be
    correct.
  * msl12_input.c, added extra XML metadata.  Now use ifile and ofile CLO types.
    return 110 if n,s,e,w do not intersect file data.
  * main_vcaltarget.c, main_l2gen.c, main_l3gen.cpp, main_l1mapgen.c, 
    main_l1bgen.c, msl12_input.c, main_l1brsgen.c, l12_proto.h, added 
    selectOptionKeys to CLO to limit the help and XML option output.  
    Changed the usage calls.  Changed the input initalization.
  * aerosol.c, fixed a bug with the model selection which caused a seg fault.
  * prodlist.c, removed duplicate products.
  * main_l2gen.c,l12_parms.h, msl12_input.c, return 110 if n,s,e,w do not 
    intersect file data.

lib
  * clo.c,clo.h changed the XML output of the options.  Added ifile
    and ofile as parameter types.  Added selectOptionKeys to CLO to 
    limit the help and XML option output.


## 6.5.6 - 2012-02-01

### Source Changed
  * l12_parms.h, version change.
  * giop.c, use variance-covariance for giop uncertainties.
  * getformat.c, fix handling of viirs files in list input.
  * l1_viirs_h5.c, handle 47 scan granules
  * l2prod.h, add fsat2 support (test case, no floor on negatives)
  * fluorescence.c, add fsat2 support
  * prodgen.c, add fsat2 support.
  * get_l2prod_index.c, add fsat2 support.
  * l12_proto.h, add get_fsat2.
  * lonlat2pixline.c, lonlat2pixline.h, fixed resolution parameter for MODIS
    GEO files.  Reduced stack usage by putting instr on the heap.
  * main_l2gen.c, put a bunch of structures on the heap to reduce stack size.
  * setanc.c, made the error reporting more descriptive when there are
    missing anciliary files.
  * brightness.c, l1_viirs_h5.c, l1_meris_CC.c, added if statement to quite 
    output unless want_verbose is set.
  * getformat.c, reformatted file, fixed VIIRS file list and MERIS CC detection
  * main_l3gen.cpp, init path angles to zero, add approximation for solar zenith 
  assuming noon (fix for Kd_lee eval product), using sunangs.f.
  * photic_depth.c, use actual solar zenith in Lee algorithm (was zero).
  * l12_proto.h, add sunangs.f proto. 
  * setupflags.[ch] moved from many programs into libl2
  * added flaguse to l2brsgen
  * l2mapgen/l2mapgen_input.c, l2mapgen/main_l2mapgen.c, main_l1mapgen.c,
    msl12_input.c, l1mapgen.mk, added GeoTIFF output option

test/lonlat2pixline/
  * t37.txt, t38.txt, fixed the resolution parameter for MODIS GEO files.
  * t41.txt, added test for MERIS CC netCDF files
  * t42.txt, t43.txt, added test for VIIRS single file and list of files.

  * updated most of the test directory data files.

  * test/l2mapgen/t9, test/l1brsgen/A2002365234500.L1B_MAP.tiff - added tests 
    for GeoTIFF.


## 6.5.5 - 2011-10-14 

Add CoastColour file support. Rename chl_ci to chl_oci.

### Source Changed
  * l12_parms.h, version change.
  * filehandle.h, add FT_MERISCC format ID.
  * getformat.c, identify FT_MERISCC fiel type.
  * l1_io.c, send FT_MERISCC file type to new input routine.
  * l1_meris_CC.c, new MERIS CC input routines.
  * l1_meris_CC.h, new MERIS CC input routines.
  * l2prod.h, rename chl_ci to chl_oci.
  * get_chl.c, rename chl_ci to chl_oci.
  * prodgen.c, rename chl_ci to chl_oci.
  * get_l2prod_index.c, rename chl_ci to chl_oci.
  * msl12_input.c, add parameter and xml file output for l2gen command 
    line options

  * clo.c, clo.h added ability to dump options to XML or param file


## 6.5.4 - 2011-08-10

### Source Changed
  * l12_parms.h, version change.
  * l2_hdf_generic.c, l12_seawifs.c, l1_hdf_generis_write.c, ms12val_hdf.c
     mscal_struc.c, libhdfutils/hdf_utils.c, ice2hdf/ice2hdf.c, 
     inc/utils/hdf_utils.h, added standard name to the HDF metadata
  * HOWTO_Add_a_product.txt, modified to account for new l2prod structure
  * l2prod_struc.h, get_l2prod_index.c, changed structure for new product 
     definition and selection
  * main_l2gen.c, msl12_input.c, input_struc.h, l12_proto.h, 
     added prodxmlfile input parameter
  * added PNG creation to l2mapgen and l2brsgen
  * fixed a bug in lonlat2pixline that returned the wrong results when
     a L2 file was inside a requested box execpt the west side of the
     file is outside the box.
  * l2prod.h, no longer need MAX_CAT.
  * added PNG output to l1brsgen and l1mapgen
  * lonlat2pixline.c, fixed a bug for MODIS GEO files.

  * [libgenutils] inc/utils/genutils.h, moved get_product_table and 
     get_lut_file functions from l2gen directory into genutils library.
  * inc/utils/clo.h, added #ifdef to play nice with C++.  Added a function
     to output a param file from the current values of the options.
  * templates, simplified the templates and removed half of the files.

  * [lib3] updated to latest versions of hdf4 and hdf5.  Removed netCDF
    from HDF4.  Removed all of the bins and libs from svn.
  * added an architecture directory to the bin directory.

## 6.5.3 - 2011-07-19

### Source Changed
  * l12_parms.h, version change.
  * msl12_input.c, remove qaa_opt.
  * input_struc.h, remove qaa_opt.
  * main_l3gen.cpp, adapt to updated libbin++ interface.

## 6.5.2 - 2011-06-26

### Source Changed
  * l12_parms.h, version change.
  * get_chl.c, updated chl_hu algorithm, add chl_ci algorithm.
  * l2prod.h, add chl_ci.
  * prodgen.c, add chl_ci.
  * get_l2prod_index.c, add chl_ci.
  * convl21.c, include radcor delta.
  * photic_depth.c, update to Zsd_gbr coefficients.


## 6.5.1 - 2011-03-16

### Source Changed
  * l12_parms.h, version change.
  * msl12_input.c, error in reporting of gain_unc paramter setting.
  * getformat.c, l1_meris_N1.c, l1_meris_N1.h, lonlat2pixline.c, fixed a bug
    when using lonlat2pixline with l2gen for MERIS files.
  * l2prod.h, add fitpar_giop products.
  * giop.c, add fitpar_giop products.
  * prodgen.c, add fitpar_giop products.
  * get_l2prod_index.c, add fitpar_giop products.

### Data Added
  * updated cocclith coefficients for non-NASA sensors.


## 6.5.0 - 2011-03-16

Recompiled for production delivery.

### Source Changed
  * l12_parms.h, version change.


## 6.4.3 - 2011-03-15 

Remove all tabs in fortran code. Added test for '$' before passing demfile to
filenv function, as filenames with fully qualified paths (i.e. with expanded
environment variables) failed in gfortran on the Mac.  Repair h5 i/o seg fault
problem. Lots of changes to reduce compiler warnings. Modify l2gen to not treat
VIIRS bow tie delete as a stray light source. Fix the dateline bug in
lonlat2pixline. 

### Source Changed
  * aer_io.c
  * aerosol.c
  * alloc_l1.c
  * alloc_l2.c
  * atmcor_soa.f
  * bin_climatology.c
  * brightness.c
  * calc_par.f
  * carder.c
  * covariance_inversion.c
  * dem_s.fin
  * dtran_brdf.f
  * filter.c
  * filter.h
  * get_dem_height.f
  * getformat.c
  * getglint.f
  * getl1rec.c
  * get_l2prod_index.c
  * get_owmc.c
  * h5io.c
  * h5io.h
  * intpos.f
  * l12_proto.h
  * l1a_seawifs.c
  * l1_czcs_hdf.c
  * l1_hdf_generic_read.c
  * l1_hdf_generic_write.c
  * l1_hmodis_hdf.c
  * l1_io.c
  * l1_modis_hdf.c
  * l1_modis_hdf.h
  * l1_mos_hdf.c
  * l1_ocm2_hdf.c
  * l1_ocmdb_hdf.c
  * l1_ocm_hdf.c
  * l1_pci_hdf.c
  * l1subpix.c
  * l1_viirs_h5.c
  * l1_xcal_hdf.c
  * l2_flags.h
  * l2gen/l1a_seawifs.h
  * l2gen/l1_viirs_h5.c
  * l2gen/polcor.c
  * l2gen/smi_climatology.c
  * l2_hdf_generic.c
  * las_iop.c
  * loadl1.c
  * lonlat2pixline.c
  * main_l1bgen.c
  * main_l1info.c
  * main_l1mapgen.c
  * main_l2binmatch.c
  * main_l2gen.c
  * main_lonlat2pixline.c
  * main_vcalmerge.c
  * main_vcaltarget.c
  * mscal_struc.c
  * msl12_input.c
  * msl2val_hdf.c
  * owt.c
  * pml_iop_tables.c
  * polcor.c
  * rdsensorinfo.c
  * sst.c
  * target_io.c
  * viirs_utls.c


## 6.4.3 - 2010-12-10

VIIRSN changes for new RSR and band centers, support for thermal bands and SST.

### Source Changed
  * l12_parms.h, version change.
  * l1_viirsn_h5.c, add switching of calibration for VSWIR and THERM, and add use of
  reflectance for VSWIR.  Move CIRRUS band to rho_cirrus slot in l1rec.  Apply K&R
  formatting to source code.
  * nlw_outband.c, move 551 test to 550 for comfort.
  * prodlist.c, bug in handling of wavelengths > 9999 when translating iii.
  * brightness.c, extensive generalization of i/o for handling different sensors
  with varying dimensional information.
  * sst.c, reconize other sensors and treat as aqua, clean-up.
  * l1_viirs_h5.c
  * windex.c, expand wavelength range for bindex of thermal channels.
  * l1_io.c, move l1rec sensor dimensional information from loadl1 to here.
  * loadl1.c, see above.
  * msl12_input.c, add new parameter btfile.  Added a fix for MODIS resolution
    when using lonlat2pixline.
  * input_struc.h, add new parameter btfile.
  * lonlat2pixline.c, fixed the box retreval that crosses the dateline.

### Data Added
  * viirsn/rayleigh, new tables based on govbest RSRs
  * viirsn/aerosol, new tables based on Ahmad set03 and govbest RSRs
  * viirsn/msl12_sensor_info.dat, new bandpass info based on govbest RSRs
  * viirsn/msl12_defaults.par, new aerosol default, etc.
  * viirsn/cal/bt_viirsn.hdf, brightness tempeature conversion table.
  * *modis*/msl12_defaults.par, add btfile field.

  * [lib3] got rid of the shared lib for libz and proj4
  * [lib3] eliminated almost all compiler warnings on linux32, linux64, macintel
  * [lib3] consolidated a bunch of date functions into libgenutils

## 6.4.2 - 2010-12-08

OWT code and table update for cocco class.

### Source Changed
  * l12_parms.h, version change.  
  * fuzzy_func_v3.c, new version for 16 classes, 9 wts.
  * fuzzy_func.c, removed.
  * covariance_inversion.c, new code for covariance inversion.
  * owt.c, update for 16 classes, 9 wts.
  * ludcmp.c, needed for covariance inversion.

### Data Added
  * modis/class/owt16_modis_stats_101111.hdf, new fuzzy lut.
  * seawifs/class/owt16_seawifs_stats_101111.hdf, new fuzzy lut.

## 6.4.1 - 2010-11-09

GIOP changes to support uncertainty estimation.  Update to aph Bricaud to do
band averaging.

### Source Changed
  * l12_parms.h, version change.
  * init_l2.c, set Rrs and nLw uncertainties to 0.0.
  * l1subpix.c, handle new input fields, nobs and Lt_unc.
  * aph.c, add wrapper to call desired aph* function and compute band average.
  * convl12.c, transfer nobs and Lt_unc fileds to L2.
  * l2prod.h, add uncertainty fields and reorganize a bit.
  * alloc_l1.c, add nobs and Lt_unc fields.
  * l1_struc.h, add nobs and Lt_unc fields. 
  * l2_struc.h, add nobs and Lt_unc, Rrs_unc and nLw_unc fields.
  * giop.h, 
  * alloc_l2.c, add space for Lt_unc, Rrs_unc and nLw_unc fields.
  * prodgen.c, add support for input/output uncertainty fields.
  * loadl1.c, assign uncertainty to Lt (Lt_unc) based on input gain_unc, if Lt_unc was
  not set by the sensor i/o routine.
  * main_l3gen.cpp, read the number of observations and store in l1rec.
  * get_l2prod_index.c,  add support for input/output uncertainties and chisqr.
  * init_l1.c, init new Lt_unc field to BADFLT and nobs to 1. 
  * msl12_input.c, add support for vicarious gain uncertainties, gain_unc, and input Rrs 
  uncertainties, giop_rrs_unc.
  * atmocor2.c, propogate Lt_unc to nLw_unc and Rrs_unc uncertainties (first wag).
  * input_struc.h, add vicarious gain uncertainties, gain_unc, and input Rrs 
  uncertainties, giop_rrs_unc.
  * l12_proto.h, new get_aphstar() wrapper function.
  * giop.c, fix bug in iops_giop interface, add chisqr_giop and uncertainty products, 
  double the size of static arrays to store uncertainties, compute uncertainties 
  for LEVMARQ optimized parameters and derived products if given giop_rrs_unc input
  uncertainties, apply uncertainty as weights in LEVMARQ fit if  given giop_rrs_unc 
  input, standardize use of nbands and ipb vs ipb2, call new get_aphstar wrapper to 
  get band-averaged bricaud & ciotti, remove LU decomp code and option, pass foq 
  into model since it can change wavelengths between optimization and evaluation
  (this was a potential bug if 412 was dropped from fit), compute and return 
  chisquare in LEVMARQ, set default fit wavelengths to 0-671 range (exclude 678).

  * [lib3] moved the lib3 directories down a level or two to get rid of multiple
  -I and -L directives and elimiate the architecture directory.
  * [lib3] updated hdfeos, hdf5, hdfeos5, and sdtpk to the latest version.

## 6.4.0 - 2010-11-04 

Revised Hu chl.  Add Bricaud 1998 aph for giop.  Added lon,lat extraction to
l2gen.  Fixed l2gen help.  Added h4h5tools and proj4 to lib3.  Added run/bin3.
Moved spline and splint to libgenutils.

### Source Changed
  * l12_parms.h, version change.
  * aph.c, new wrapper for bricaud and caal to 1998 function.
  * get_chl.c, updated hu et al. coefficients.
  * main_l2gen.c, l2gen.mk, added lon,lat extraction to l2gen and fixed help
  * msl12_input.c, input_struc.h, added xbox,ybox to input structure, fixed help
  * main_l3gen.cpp, main_l1mapgen.c, main_l1brsgen.c, fixed help display
  * l2bin/l2bin_viirs.mk, makefile now includes libraries properly
  * bindepths/bindepths.c, bindepths/Makefile, libl2/readL2scan.c, 
l2brsgen/Makefile, nr.h, l2extract/Makefile, l2extract/main_l2extract.c, 
swfinc/readL2scan_viirs.h, use spline and splint from genutils

  * [lib3] handle optimization properly when building in debug mode 
  * [lib3] add ability to build without running the tests
  * [lib3] sdptk now links with the correct zlib library
  * [lib3] added and populated run/bin3 with useful binaries from lib3
  * [lib3] added proj4 and h4h5tools to lib3
  * [lib3] added -f to "#!/bin/csh" line to make the enviroment work properly
  * [lib3] added LIB3_CHECK and LIB3_BIN to OCSSW.env

## 6.3.9 - 2010-10-20 

Add support for recognizing and reading various VIIRS proxy data.

### Source Changed
  * getformat.c
  * l1_viirs_h5.c
  * h5io.c
  * h5io.h


## 6.3.9 - 2010-10-13

Add salinity reference climatology.  Add salinity and temperature dependence
options for backscatter of pure seawater.

### Source Changed
  * input_struc.h, add sssfile support, add giop_ts_opt.
  * l12_proto.h, add get_sssref function, add seawater_bb().
  * sssref.c, new code to read and interpolate salinity fields.
  * convl12.c, moved pointer transfer code to a function so l3gen could share it.
  * loadl1.c,  call get_sss() function.
  * main_l3gen.cpp, add new call to cpl1l2() to tranfer record pointers.
  * get_l2prod_index.c, add sssref prod.
  * msl12_input.c, add sssfile supporti, add giop_ts_opt.
  * giop.c, add call to seawater_bb().
  * giop.h, add ts_opt.
  * loadl1.c, clean-up.
  * seawater.c, new code containing equations of state for pure seawater.


### Data Added
  * common/sss_climatology_woa2009.hdf, new salinity climatology.
  * */msl12_defaults.par, add sssfile to all sensor defaults.

## 6.3.9 - 2010-10-06

Preparations for SSS i/o.  Fixing some valgrind complaints regarding demfile.
GIOP cleanup. Add chl_hu algorithm.

### Source Changed
  * l12_parms.h, version change.
  * convl12.c, add sssref.
  * l2prod.h, add sssref.
  * alloc_l1.c, add sssref.
  * l1_struc.h, add sssref.
  * l2_struc.h, add sssref.
  * prodgen.c, add sssref. 
  * get_l2prod_index.c, add sssref.
  * init_l1.c, add sssref.
  * msl12_input.c, remove GIOP characteristic wavelengths, init demfile to nulls, fix cases of
  sprintf using the variable for input and output (as this is an undefined behavior). 
  * get_dem_height.f, changed demfile string dimensions for consistency with calling routine.
  * l2prod.h, add chl_hu.
  * get_chl.c, add chl_hu.
  * prodgen.c, add chl_hu.
  * get_l2prod_index.c, add chl_hu.


## 6.3.8 - 2010-09-21

Change test for TOMS ozone to recognize new filenames of ozone files.  This
changes how interpolation is done. Fix handling of SeaWiFS cloud dilation.
GIOP fixes. CZCS navigation fix for dateline.  updated HDF4 to version 4.2.5.
Fixed the resolution=-1000 command for MODIS.

### Source Changed
  * l12_parms.h, version change.
  * main_vcaltarget.c, fix wavelength interpolation.
  * msl12_input.c, improve handling of flaguse.
  * l1_czcs_hdf.c, 
  * ll2vec.c, new code to interpolate lon/lat control points using vectors.
  * main_l1mapgen.c, 
  * giop.c, set spectral slopes to bad values when pixel is masked. Fix
  QAA bbp_s computation.
  * xcal.h, expand arrays to store up to NBANDS.
  * getformat.c, msl12_input.c, fixed the resolution=-1000 option
  * main_l1info.c, changed calling arguments for msl12_input_init
  * main_vcaltarget.c, main_l2gen.c, msl12_input, main_l1bgen.c, changed calling arguments for msl12_input, extra getFormat not needed anymore.
  * filehandle_init.c, need to init new member modis_subset_compat
  * main_l1mapgen.c, main_l1brsgen.c, changed calling arguments for msl12_option_input, extra getFormat not needed anymore.
input_struc.h, added modis_subset_compat member to instr
  * l12_proto.h, changed function prototypes to pass l1 filehandle in.

  * [libanc] getanc.c, change TOMS.OZONE to TOMS on string test.

### Data Added
  * seawifs/seawifs_gac_filter.dat, change dilate to stlight, as originally intended.

  * [lib3] hdf4, updated to version 4.2.5 with a fix in atom.c

## 6.3.7 - 2010-08-24

Modifications for VIIRSN simulation.

### Source Changed
  * l12_parms.h, version change.
  * aerosol.c, added test to verify that all model files are of the same dimensions.
  * filehandle.h, VIIRS to VIIRSN.
  * getformat.c, support for VIIRS file list.
  * h5io.c, 
  * l1_viirs_h5.c, 

### Data Added
  * virrsn/*, added additional bands and suite defaults.


## 6.3.6 - 2010-08-04 

Recompiled with updated lib3 suite.  Improved usage messages for GIOP, and some
simplification.

### Source Changed
  * l12_parms.h, version change.
  * msl12_input.c, usage messages.
  * giop.c, some reorg for simplification.
  * giop.h, change in nvec field names for clarity.
  * polcor.c, add secret option so Lu and Lq are computed but polcor is not,
  pol_opt=100.

### Data Added
  * seawifs/cal, new default cal table 2010.hdf
  * seawifs/msl12_defaults.par, new defaults gains.


## 6.3.5 - 2010-07-18

Add dynamic allocation of aerosol tables. GIOP support for multiple tabulated
vectors per component.  Fix GIOP for error introduced in 6.3.4 where input chl
was not updated.

### Source Changed
  * l12_parms.h, version change.
  * l12_proto.h, add table_io_wrapper.h.
  * aerosol.c, dynamic allocation for aerosol LUT i/o (from P.Martinolich).
  * loadl1.c, clean-up unused vars.
  * giop.c, add ability to handle multiple vectors per tabbed input file, save
  fitted params for full scanline. Fix error in input chl init.
  * giop.h, add fields for number of vectors per component.


## 6.3.4 - 2010-07-09 

Fix for Terra segfault. GIOP usage messages.

### Source Changed
  * l12_parms.h, version change.
  * l1_modis_hdf.c, change string handling to in get_hdfeos_meta to avoid overlap in strncpy. Not a problem
  but valgrind complains.
  * msl12_input.c, 
  * aerosol.c, trap divide by zero, sqrt of zero issues in rhoa_to_rhoas (due to possible issue with LUTs).
  * bin_climatology.c, error message fixes.
  * libgenutils/clo.c, error message fix.
  * src/templates/MakeDefine.xxx.tpl, fixed an if statement bug on some systems.
  * l1mapgen.mk, main_l2gen.c, main_l1mapgen.c, main_l1brsgen.c, msl12_input.c, input_struc.h, l12_proto.h, made changes 
  so l1mapgen and l1brsgen would have the same command line parameters as l2gen.
  * build/Makefile, src/Makefile, modifies the top Makefiles so they would bail when a compile error is encountered.


## 6.3.3 - 2010-07-01

### Source Changed
  * l12_parms.h, version change.
  * l2prod.h, new chl_cdomcorr_morel.
  * cdom_morel.c, add chl_cdomcorr_morel algorithm.
  * get_l2prod_index.c, add chl_cdomcorr_morel algorithm, add static identifier for default chl alg.
  * prodgen.c, add chl_cdomcorr_morel algorithm. 
  * msl12_input.c, fix typo in usage.

## 6.3.2 - 2010-06-21

### Source Changed
  * getformat, add support for generic L1B format for OCM2.
  * l1_hdf_generic_write, add start pixel meta-data for generic L1B.

lib
  * libseawifs, new caltable format for seawifs

### Data Added
  * seawifs/cal, new caltable 2009b (old cal, new format)
  * seawifs/msl12_defaults.par, new default for caltable

## 6.3.2 - 2010-05-24

Vcaltarget, Vcalmerge updates.

### Source Changed
  * l12_parms.h, version change.
  * main_vcaltarget.c, clean-up old code.
  * main_vcalmerge.c, minor clean-up.
  * msl12_input.c,  eliminate some un-used vcaltarget/vcalmerge inputs.  Saved input files from ilist into ifile# parameters.
  * input_struc.h, eliminate some un-used vcaltarget/vcalmerge inputs.
  * get_owmc.c, change Oliver et al. class function to float.
  * get_l2prod_index.c, change Oliver et al. class function to float.
  * l12_proto.h, change Oliver et al. class function to float.
  * lonlat2pixline.c, fixed a bug indexing the wrong element of an array. if the box requested is within the file, but too small to return any pixels, we now do a pixel search.
  * main_lonlat2pixline, use -x and -y to return a min # of pixels and lines when using box mode.

  * libgenutils/clo.c, print errors to stderror instead of stdout


## 6.3.1 - 2010-05-18

New Command line option library

### Source Changed
  * l2gen/l12_parms.h, version change.
  * libg2/Makefile, make mac work with gfortran
  * libgenutils/clo.c, libgenutils/Makefile, command line option library
  * l2gen/msl12_input.c,l2gen/l12_proto.h,l2gen/main_l1bgen.c,l2gen/main_vcaltarget.c, new command line arg lib used
  * src/libmap/get_l3m.c, Read data min/max for HDF5 files
  * libswfnav/proctim2.f, fix segfault in proctim2
  * libhdfutils/hdf5_utils.f, Use call exit(1) for 'not found' errors

  * utils/clo.h, command line option library header

### Data Added
  * mos/msl12_defaults_OC.par, osmi/msl12_defaults_OC.par,data/common/msl12_defaults.par, new data files needed for the command line lib changes.


## 6.3.0 - 2010-04-27

More GIOP enhancements.

### Source Changed
  * l12_parms.h, version change.
  * get_qaa.c, add function to retrieve QAA bbp for GIOP.
  * giop.c, add addition bbp options for fixed vector via QAA and LAS, and spectral slope derived
  from LAS.
  * giop.h, add BBPLASFIX, BBPQAAFIX, BBPSLAS.
  * prodgen.c, support for output bbp_s_las.
  * get_l2prod_index.c, support for output bbp_s_las.
  * msl12_input.c, new GIOP options.
  * las_iop.c, function to retrieve LAS bbp for GIOP, function to fit bbp spectral slope from LAS
  bbp and return it.
  * l12_proto.h, new functions for retrieving LAS bbp spectral slope and QAA bbp.
  * getformat.c, fixed typing of ocean subsetted MODIS files.
  * l1_modis_hdf.c, made the program terminate if you ask for high res with a low res MODIS file.
  * lonlat2pixline.c, made the program terminate if you ask for high res with a low res MODIS file.  Got rid of the reference to HQM and KQM files from the usage printout


## 6.2.9 - 2010-04-14

l2gen viirs unaggregated update. GIOP enhancements.

### Source Changed
  * input_struc.h, options to control deep-water mask.
  * msl12_input.c, options to control deep-water mask.
  * main_vcaltarget.c, options to control deep-water mask.
  * l1_struc.h, support for VIIRS superscan.
  * polcor.c, support for VIIRS superscan.
  * h5io.c, 
  * h5io.h, 
  * l1_viirs_h5.c
  * l1_viirs_h5.h
  * giop.h, add chl to control struc, add BBPLAS opts.
  * giop.c, replace chl_in local var with control struc, add BBPLAS opts.
  * photic_depth.c, fix Zsd_gbr coeffs.
  * main_l3gen.cpp, use copy_meta.
  * las_iop.c, add giop interface.
  * msl12_input.c, fix giop_bbp_opt options, add BBPLAS opts.
  * get_format.c, fix for subsetted modis bug.


## 6.2.9 - 2010-04-09 

### Source Changed
  * l12_parms.h, version change.
  * filehandle_init.c, initialize terrain_corrected field.
  * get_dem_height.f, trying to fix string handling complaint from valgrind.
  * l2prod.h, add rrsdiff_giop and mRrs_giop prods.
  * prodgen.c, add rrsdiff_giop and mRrs_giop prods.
  * giop.c, add rrsdiff_giop and mRrs_giop prods.
  * giop.h, remove giop_chl options.
  * main_l3gen.cpp, add computation of default chl and iops.
  * get_l2prod_index.c, add rrsdiff_giop and mRrs_giop prods.
  * flags_iop.h, add IOPF_RRSDIFF field.
  * msl12_input.c, remove giop_chl options, init xcal_nwave.
  * input_struc.h, remove giop_chl options.
  * alloc_l1.c, add want_verbose test to informational messages.
  * getformat.c, add want_verbose test to informational messages.
  * l1_hmodis_hdf.c, add want_verbose test to informational messages.
  * l1_meris_N1.c, add want_verbose test to informational messages.
  * lonlat2pixline.c, add want_verbose test to informational messages.
  * lonlat2pixline.h, add want_verbose test to informational messages.
  * main_lonlat2pixline.c, add want_verbose test to informational messages.
  * msl12_input.c, add want_verbose test to informational messages.
  * rdsensorinfo.c, add want_verbose test to informational messages.


## 6.2.8 - 2010-03-31

### Source Changed
  * l12_parms.h, version change.
  * photic_depth.c, add Zsd_gbr function.
  * l2prod.h, add Zsd_gbr product.
  * prodgen.c, add Zsd_gbr product.
  * get_l2prod_index.c, add Zsd_gbr product.

## 6.2.7 - 2010-03-19

### Source Changed
  * l12_parms.h, version change.
  * getformat.c, add OCM2 L1B support.
  * l1_io.c, add OCM2 L1B support.
  * filehandle.h, add OCM2 L1B support.
  * get_chl.c, add OCM2 chl algo.
  * get_l2prod_index.c, add OCM2 chl algo.
  * l1_ocm2_hdf.c, new OCM2 L1B i/o.
  * l1_ocm2_hdf.h, new OCM2 L1B i/o.
  * atmocor2.c, add 620 as a valid red band.
  * aerosol.c, fix RH switching for aer_opt=-1.


## 6.2.6 - 2010-03-17

### Source Changed
  * l1_modis_hdf.c, add orbit information to MODIS output
  * l1_hmodis_hdf.c, add orbit information to MODIS output
  * l1_hmodis_hdf.c, change LAC to [HQ]KM in filename only
  * get_rhown_nir.c, fix bbp extrapolation (this changes high-C pixels slightly)

## 6.2.6 -2010-03-15 

Update QAA, fix MODIS resolution problem.

### Source Changed
  * l12_parms.h, version number.
  * get_qaa.c, new get_qaa() from P. Martinolich, removed references to 640 band
  modified default params to specify the red wvl.
  * filehandle.h,
  * filehandle_init.c,
  * getformat.c,
  * l1_hmodis_hdf.c,
  * l1img_input.c,
  * main_l1bgen.c,
  * main_l1brsgen.c,
  * main_l1mapgen.c,
  * main_l2gen.c,
  * main_vcaltarget.c,
  * msl12_input.c,

### Data Added
  * */msl12_defaults.par, update QAA wavelengths


## 6.2.6 - 2010-02-26

### Source Changed
  * l12_parms.h, version number.
  * input_struc.h, minor mods to eliminate compile warnings.
  * l1_viirs_h5.c, minor mods to eliminate compile warnings.
  * main_l1info.c, minor mods to eliminate compile warnings.
  * xcal.c, correct type of norder, for PPC compatibility.

## 6.2.5 - 2010-02-12 

### Source Changed
  * l12_parms.h, version change.
  * getformat.c, auto-detect modis band-suite type (9 or 16-band), disable support for
  HKM/QKM input format.
  * input_struc.h, add tags for all input sst coeff files.
  * msl12_input.h, read input sst coeff file paths.
  * sst.c, use input sst coeff file paths (don't assume paths).
  * main_l2brsgen.c, revised calling args to readL2.
  * main_l2binmatch.c, revised calling args to readL2.
  * main_vcaltarget.c, revised calling args to readL2.
  * main_lonlat2pixline.c, revised calling args to readL2.
  * main_vcalmerge.c, revised calling args to readL2.

  * [libl2] add option to disable caching.

### Data Added
  * modisa, hmodisa, drop fqy as default product, add MODIGLINT masking to FLH suite, drop
  sst file inputs from hmodisa, rename sst file inputs for modisa.
  * modist, hmodist, drop fqy as default product, add MODIGLINT masking to FLH suite, drop
  sst file inputs from hmodisa, rename sst file inputs for modisa.

## 6.2.4 - 2010-01-28

### Source Changed
  * l12_parms.h, version change.
  * xcal.c, fix for alternate resolutions.
  * setflags.c, fix coccolith, atmwarn, and lowlw flags tests to use MODIS 547, not 555.

### Data Added
  * modisa, hmodisa, change defaults for reprocessing.
  * modisa, hmodisa, move polcor and xcal files from eval.

## 6.2.3 - 2010-01-24

### Source Changed
  * l12_parms.h, version change.
  * l2prod.h, add fdiff.
  * prodlist.c, fix iii and nnn options to append suffix.
  * fluorescence.c, return flh=0 if negative, add fdiff.
  * prodgen.c, add fdiff.
  * get_l2prod_index.c, add fdiff, add nflh alternate name.
  * l1_hmodis_hdf.c, fix bug in reading cirrus reflectance.

### Data Added
  * hmodisa/msl12_defaults, switch albedo to 0.027, since 869 is default cloud, change
  to read files from modisa where possible.
  * hmodisa/filter_file.dat, deleted (use modisa).
  * hmodisa/polcor files, deleted (use modisa).


## 6.2.2 - 2010-01-11

### Source Changed
  * l12_parms.h, version change.
  * filehandle.h, add terrain_corrected flag.
  * loadl1.c, call get_dem_height only when geolocation not already terrain-corrected.
  * l1_modis_hdf.c, read height from geolocation file
    read gflags to determine if geolocation has been terrain-corrected
  * l1_hmodis_hdf.c, read height from geolocation file
    read gflags to determine if geolocation has been terrain-corrected
    modis_geo_interp interpolates height
        add rho_cirrus code from l1_modis_hdf()
  * msl12_input.c, add cloud_wave, flh_offset.
  * input_struc.h, add cloud_wave, flh_offset.
  * setflags.c, use explicit cloud_wave for ocean cloud threshold test.
  * fluorescence.c, subtract input flh_offset from flh, fsat, add "typical angstrom" effect
  to aerosol adjustment, add limits on fqy to set PRODWARN.
  * get_l2prod_index.c, add scaling to int16 for flh products.


## 6.2.1 - 2009-12-15

### Source Changed
  * l12_parms.h, version change.
  * xcal.c, add band-dependent control of options and band-specific file format.
  * xcal.h, add band-dependent control of options.
  * polcor.c, support band-dependent xcal options.
  * loadl1.c, support band-dependent xcal options.
  * msl12_input.c, support band-dependent xcal options, add owtfile control..
  * input_struc.h, support band-dependent xcal options, add owtfile control.
  * owt.c, replaces optical_class.c.
  * get_l2prod_index.c, switch owtd product to float for now.

### Data Added
  * modisa/cal/xcal_modisa_nnn.hdf, new band-specific xcal files.

## 6.2.0 - 2009-12-14

### Source Changed
  * l12_parms.h, version change.
  * xcal.h, quick fix change to the dimensions of the xcal file.
  * pml_iop_tables.c, add byte swapping support to PML table i/o.
  * fuzzy_func.c, support for OWT fuzzy water classes.
  * gammln.c, support for OWT fuzzy water classes.
  * gcf.c, support for OWT fuzzy water classes.
  * gser.c, support for OWT fuzzy water classes.
  * nr.h, support for OWT fuzzy water classes.
  * nrutil.c, support for OWT fuzzy water classes.
  * nrutil.h, support for OWT fuzzy water classes.
  * optical_class.c, support for OWT fuzzy water classes.
  * sprsax.c, support for OWT fuzzy water classes.
  * sprsin.c, support for OWT fuzzy water classes.
  * get_l2prod_index.c, add OWT products.
  * l2prod.h, add OWT products.
  * prodgen.c, add OWT products.
  * giop.c, generalized above-to-below water code for OWT.
  * l12_proto.h, generalized above-to-below water code for OWT.
  * loadl1.c, support for new xcal_opt switches.
  * msl12_input.c, support for new xcal_opt switches.
  * polcor.c, support for new xcal_opt switches.
  * getformat.c, support LSU OCM data format.
  * main_l3gen.cpp, make noext=1 the default behavior.
  * l1img_input.c, bug fix for resolution control.
  * main_l1brsgen.c, bug fix for resolution control.


## 6.1.9 - 2009-11-30

### Source Changed
  * l12_parms.h, version number change.
  * l1_io.c, call to openl1_ocmdb_hdf was missing.
  * l1_ocmdb_hdf.c, don't set NAVFAIL for filled values, as this screws-up meta-data.

## 6.1.8 - 2009-11-18

### Source Changed
  * l12_parms.h, version number change, add AERRHMUMM and AERNULL aerosol options.
  * convert_band.c, new functions to shift band centers.
  * l12_proto.h, add conv_rrs_to_555().
  * l2gen.mk, add convert_band.c.
  * l3gen.mk, add convert_band.c.
  * get_poc.c, add generic band shift
  * cdom_morel.c, add generic band shift
  * aerosol.c, change dtran interp decision to use a float tolerance of 0.51 nm, add support for
  AERRHMUMM and AERRHNULL aerosol options, automatically switch non-RH options to RH options
  if model set suggests it is appropriate (30+ w/RH info).  Set angstrom to zero if AOT is 0.0.
  * msl12_input.c, revise aerosol option list.
  * atmocor2.c, add AERRHMUMM and AERRHNULL support, change NIR/SWIR switch to Wang 2008, record where SWIR is used via TURBIDW flag.
  * setflags.c, don't set TURBIDW if NIR/SWIR switching.

### Data Added
  * meris/msl12_defaults.par, switch to current R2009 defaults (AF aerosol models)
  * meris/msl12_defaults_OC.par, add meris OC defaults parameter file.
  * meris/aerosol/*, add AF aerosol models, remove GW aerosol models.


## 6.1.7 - 2009-11-09

### Source Changed
  * l12_parms.h, version number change.
  * get_poc.c, add 547 to 555 shift.
  * calc_par.f, disable zenith limit checks and error messages.  solar zenith limit of 83-deg
  is applied in get_par.c to all sensors.
  * atmocor2.c, add simple glint correction on glint_opt=2.
  * convl12.c, add simple glint correction on glint_opt=2.
  * msl12_input.c, add simple glint correction on glint_opt=2.
  * h5io.c, add grab_ds function.
  * h5io.h, add grab_ds function.
  * l1_viirs_h5.c, compute polarization frome rotation angle, alpha.

  * [libnav] ocorient.f, scrounged from OCl1bgen for VIIRS alpha angle computation.

## 6.1.6 - 2009-11-02 

### Source Changed
  * l12_parms.h, version number change.
  * get_rhown_nir.c, revised apg to chl relation, eta NaN fix, chl range 0.2-30.
  * atmocor2.c, change iteration control to use total Rrs670 on reset, only reset
  once, failed iteration now indigated by iter = max_iter + 1.
  * setflags.c, set AERITERMAX flag on > max_iter, also set ATMWARN flag.
  * qaa.c, bug in indexing fixed (idx555 not 555).
  * get_par.c, add MODIS PAR support.
  * calc_par.f, add MODIS PAR code.

## 6.1.5 - 2009-10-26

### Source Changed
  * l12_parms.h, version number change.
  * rayleigh.c, comment change.
  * l1subpix.c, add subsetting of rho_cirrus.
  * init_l1.c, add npix as argument.
  * l12_proto.h, init_l1() proto.
  * l1_io.c, pass npix to init_l1.c
  * ice_mask.c, fix missing initialization.
  * fluorescence.c, replace aerosol effect in flh with simple model.
  * l1_hmodis_hdf.c, read list of inputs to L1B from correct file


## 6.1.4 - 2009-10-19 

### Source Changed
  * l12_parms.h, version number change.
  * aph.c, remove normalization on aph* functions.
  * l2prod.h, add fqy, fqy2, fsat products, remove cfe, arp.
  * fluorescence.c, new fqy algorithms.
  * ipar.c, change i[ar to above-surface.
  * prodgen.c, add fqy, fqy2, fsat products, remove cfe, arp.
  * get_l2prod_index.c, add fqy, fqy2, fsat products, remove cfe, arp.
  * get_rhown_nir.c, fix error in apg6 relationship.
  * atmocor2.c, fix log(chl) error in iteration control for NIR algorithm.


## 6.1.3 - 2009-09-22

### Source Changed
  * l12_parms.h, version number change.
  * sst.c, switch some int to int32_t.
  * main_l2gen.c, if not specified, input resolution adopts file resolution (for MODIS).
  This allows l2gen to process MODIS at the resolution specified by the ifile.
  * polcor.c, just some indentation changes.
  * atmocor2.c, set Rrs to fill value whenever ATMFAIL occurs.  This will change the
  standard products.


## 6.1.2 - 2009-08-26

SeaWiFS reprocessing version.

### Source Changed
  * l12_parms.h, version number change, add subsensor identifiers for seawifs and MODII.
  * input_struc.h, add suite string and subsensorID.
  * msl12_input.c, usage messages on giop(), load suite defaults if specified, read
  subsensor defaults if required.
  * get_par.c, revised failure flagging.
  * get_l2prod_index.c, add aot_nnn support.
  * prodlist.c, fix null string test.
  * getformat.c, determine seawifs subsensor (LAC or GAC).
  * filehandle_init.c, add subsensorID.
  * filehandle.h, add sensorSub names.
  * calc_par.f, fix missing ozone correction for view path.
  * get_poc.c, change band 551 indexing to 550.


### Data Added
  * seawifs/msl12_defaults.par, replace l2prod spec with suite=OC
  * modisa/msl12_defaults.par, replace l2prod spec with suite=OC
  * hmodisa/msl12_defaults.par, replace l2prod spec with suite=OC
  * modist/msl12_defaults.par, replace l2prod spec with suite=OC
  * hmodist/msl12_defaults.par, replace l2prod spec with suite=OC
  * seawifs/msl12_defaults_OC.par, put l2prod list in new suite defaults.
  * hmodisa/msl12_defaults_SST4.par, new suite-specific defaults.
  * hmodisa/msl12_defaults_SST.par, new suite-specific defaults.
  * hmodist/msl12_defaults_OC.par, new suite-specific defaults.
  * hmodist/msl12_defaults_SST4.par, new suite-specific defaults.
  * hmodist/msl12_defaults_SST.par, new suite-specific defaults.
  * modisa/msl12_defaults_OC.par, new suite-specific defaults.
  * modisa/msl12_defaults_SST4.par, new suite-specific defaults.
  * modisa/msl12_defaults_SST.par, new suite-specific defaults.
  * modist/msl12_defaults_OC.par, new suite-specific defaults.
  * modist/msl12_defaults_SST4.par, new suite-specific defaults.
  * modist/msl12_defaults_SST.par, new suite-specific defaults.
  * seawifs/msl12_defaults_OC.par, new suite-specific defaults.
  * seawifs/msl12_defaults_seawifs_gac.par, new seawifs gac defaults.
  * seawifs/msl12_defaults_seawifs_lac.par, new seawifs lac defaults.


## 6.1.1 - 2009-08-14

### Source Changed
  * l12_parms.h, version number change.
  * l2_flags.h, dropped some flags, rectified some flag bit names, added
  PRODWARN.
  * l1a_osmi.c, updates for revised flag bit names.
  * l1a_seawifs.c, updates for revised flag bit names.
  * calcite.c, updates for revised flag bit names.
  * get_poc.c, updates for revised flag bit names.
  * get_tsm.c, updates for revised flag bit names.
  * setflags.c, updates for revised flags and flag bit names.
  * gsm.c, updates for revised flag bit names.
  * get_Kd.c, updates for revised flag bit names.
  * myprod.c, updates for revised flag bit names.
  * filter.c, updates for revised flag bit names.
  * get_tricho.c, updates for revised flag bit names.
  * photic_depth.c, updates for revised flag bit names.
  * get_par.c, updates for revised flag bit names.
  * l1_hdf_generic_read.c, updates for revised flag bit names.
  * aerosol.c, add explicit default aerosol model.
  * get_l2prod_index.c, change scaling on angstrom and pic.

### Data Added
  * */aerosol/aerosol_*_default.hdf, new tables.
  * */l2bin_defaults.par, rectify masking for all sensors.
  * seawifs/seawifs_gac_filter.dat, new filter file for GAC.
  * seawifs/seawifs_lac_filter.dat, new filter file for LAC.


## 6.1.0 - 2009-08-10

### Source Changed
  * convl12.c, set new fill values as default.
  * l2_hdf_generic.c, set new fill values as default.
  * l1a_seawifs.c, make new cal table format default.
  * l12_parms.h, delete old fill values, old unneeded eval switches.
  * whitecaps.c, make new model the default.
  * sst.c, set new fill values as default.
  * l1_io.c, init l1 record with new fill values.
  * get_chl.c, set new fill values as default.
  * setflags.c, set new fill values as default.
  * get_Kd.c, set new fill values as default.
  * aerosol.c, remove old AERTEST eval switch.
  * filter.c, set new fill values as default, add LTRreject filter.
  * filter.h, add LTRreject filter.
  * carder.c, set new fill values as default.
  * msl12_input.c, new defaults for coccolith and aermodels, new pversion parameter,
  remove rflag parameter.
  * atmocor2.c, set new fill values as default, new NIR Lw as default.
  * input_struc.h, new pversion.
  * get_l2prod_index.c, add scaling for new standard products.

### Data Added
  * seawifs/msl12_sensor_info.dat, new reprocessing values.
  * seawifs/msl12_defaults.par, new reprocessing values.
  * seawifs/rayleigh, new reprocessing tables.
  * modisa/msl12_sensor_info.dat, new reprocessing values.
  * modisa/msl12_defaults.par, new reprocessing values.
  * modisa/rayleigh, new reprocessing tables.
  * hmodisa/msl12_sensor_info.dat, new reprocessing values.
  * hmodisa/msl12_defaults.par, new reprocessing values.
  * hmodisa/rayleigh, new reprocessing tables.
  * modist/msl12_sensor_info.dat, new reprocessing values.
  * modist/msl12_defaults.par, new reprocessing values.
  * modist/rayleigh, new reprocessing tables.
  * hmodist/msl12_sensor_info.dat, new reprocessing values.
  * hmodist/msl12_defaults.par, new reprocessing values.
  * hmodist/rayleigh, new reprocessing tables.
  * modist/aerosol, new reprocessing aerosol tables.
  * hmodist/aerosol, new reprocessing aerosol tables.


## 6.0.2 - 2009-08-06

### Source Changed
  * l12_parms.h, version number.
  * calcite.c, read new calcite tables.
  * get_l2prod_index.c, add cdom_index.
  * aerosol.c, load and apply "default" aerosol model for low rhoa.

### Data Added
  * seawifs/iop/pic/seawifs_calcite_table.dat, new PIC table with expanded boundaries.
  * modis/iop/pic/modis_calcite_table.dat, new PIC table for modis
  * seawifs/iop/pic/seawifs_coccolith_table.dat, removed.
  * common/coccolith_table.dat, removed.
  * eval/seawifs/afrt, remove afrt files.
  * eval/modisa/afrt, remove afrt files.


## 6.0.1 - 2009-08-03

### Source Changed
  * l12_parms.h, version number.
  * aerosol.c, limit MS epsilon to 10 for RHAER options, to prevent segfaults at poles.
  set bad angstrom to BAD_FLT.
  * msl12_input.c, init cloud_eps.
  * l2_hdf_generic.c, add solar irradiance meta data.
  * calcite.c, rewrite to improve error handling.
  * l2prod.h, add logchl product.
  * prodgen.c, simplify calcite interface, add logchl.
  * get_chl.c, add logchl, don't set CHLFAIL, CHLRANGE (done in setflags).
  * get_l2prod_index.c, add logchl and scale Rrs.
  * l12_proto.h, revise calcite interface.


## 6.0.0  - 2009-07-24

NIR/SWIR switching with RH-based model selection. Add switch to force ocean
processing. Add eval version of NIR Lw iteration. SST mods for eval.

### Source Changed
  * l12_parms.h, add aer_opt AERRHSWIR and eval switch NEWNIRLW.
  * aerosol.c, add AERRHSWIR and AERRH tests, change sanity test on MS epsilon for RH aerosol
  methods to use reflectance, clean-up old AEROOBFIX stuff.
  * atmocor2.c, add AERRHSWIR and AERRH tests, add switch for get_rhown_eval().
  * msl12_input.c, usage for aer_opt expanded, expand proc_ocean usage.
  * loadl1.c, set all pixels to not LAND when proc_ocean is 2.
  * convl12.c, change proc_ocean tests for new switch.
  * get_rhown_nir.c, add get_rhown_eval() and associated functions.
  * l12_proto.h, add get_rhown_eval().
  * sst.c, for eval use night cubes for all cases and add 0.17 bias for skin.

## 5.9.9 - 2009-07-01

New SeaWIFS PIC table.

### Source Changed
  * aerosol.c, allow for extrpolation in FRRHNIR aerosol option.
  * atmocor2.c, disable NIR correction for Lw if itermax < 1.
  * calcite.c, use new seawifs table for seawifs.

### Data Added
  * seawifs/iop/pic/seawifs_coccolith_table.dat, new table.


## 5.9.8 2009-06-18

Add new option for specular cloud test, to unflag clouds that may actually be
turbid water. Add support for new aerosol model names with longer identifiers.
Add new KD2 Kd(490) algorithm with user-supplied coefficents.

### Source Changed
  * l12_parms.h, version number update, drop SPECLOUD eval switch, add NEWAERMODS eval switch, VIIRS.
  * input_struc.h, add cloud_eps input for new specular cloud test option, add KD2 algorithm
  options.
  * setflags.c, replace eval switch with cloud_eps input switch, add upper threshold on
  specular cloud test.
  * msl12_input.c, support cloud_eps input, add support for eval switch in input_init and
  add switching to new aerosol model suite names, add KD2 algorithm inpuuts.
  * aerosol.c, add better handling of failure conditions in rhoa_o_rhoas, fix message formats,
  modify get_angstrom to use 443 if valid band not provided. Allow for longer aerosol model
  names.
  * get_l2prod_index.c, recognize "angstrom" as 443 to 865, add KD2 Kd algorithm, remove OK2
  algorithm, make Kd_490 point to KD2, add ipar2, add poc and pic name variants,
  VIIRS default chl.
  * l12_proto.h, added eval switch to msl12_input_init().
  * input_struc.h, increase space for aerosol model names.
  * get_Kd.c, add KD2 algorithm, remove OK2 algorithm, reorganize.
  * l2prod.h, replace CAT_Kd_morel_ok2 with CAT_Kd_KD2, add CAT_ipar2.
  * prodgen, replace CAT_Kd_morel_ok2 with CAT_Kd_KD2, add ipar2 support.
  * ipar.c, rename function to get_ipar().
  * ipar_arp.c, rename old get_ipar() to get_ipar2().
  * l3gen.mk, add ipar.c, viirs support.
  * *.mk, add viirs files and HDF5 lib/inc.
  * getformat.c, add VIIRS support.
  * polcor.c, generalize for VIIRS.
  * get_chl.c, add VIIRS switch.
  * l1_viirs_h5.c, new VIIRS io code.
  * l1_viirs_h5.h, new VIIRS io code.
  * viirs_utls.c, new VIIRS io code.
  * l1_io.c, add VIIRS support.
  * filehandle.h, add VIIRS support.
  * h5io.c, new HDF5 code.
  * h5io.h, new HDF5 code.
  * cdom_morel.c, change failure test to chl > 0.

### Data Added
  * seawifs/aerosols/*r..f..v...hdf, new seawifs aerosol models.
  * */msl12_defaults.par, KD2 parameters.

## 5.9.7 2009-05-21

Support for new sea ice NRT data for sea ice flagging. Speed-up of
lonlat2pixline performance.  Continued preparations for evaluation of SST
reprocessing config.


### Source Changed
  * l12_parms.h, version number update.
  * sst.c, fix error in "eval" version of sst4 coeff read.
  * sstref.c, simplify title check of AVHRRSST to handle variations.
  * l2prod.h, added ice_frac product
  * get_l2prod_index.c, added ice_frac product
  * ice_mask.c, added ability to use raw NSIDC ice files, monthly and
  NRT HDF files.  added ice threshold input parameter.
  * input_struc.h, added ice threshold input parameter.
  * l12_proto.h, removed anc.h include file. Added ice_fraction functions.
  Changed parameters for ice_mask_init
  * l2gen.mk, added get_ice_frac.o product creation file.
  * l3gen.mk, added get_ice_frac.o product creation file.
  * loadl1.c, changed the calling parameters for ice_mask_init().  Only set
  the ice flag on if the pixel is not land.
  * msl12_input.c, added ice_threshold input parameter.
  * prodgen.c, added ice_frac product.
  * l1a_seawifs.h, l1a_seawifs.c, added lonlat read to speed up lonlat2pixline.
  * l1_modis_hdf.h, l1_modis_hdf.c, added lonlat read to speed up lonlat2pixline.
  * l1_hmodis_hdf.h, l1_modis_hdf.c, added lonlat read to speed up lonlat2pixline.
  * l1_meris_N1.h, added lonlat read to speed up lonlat2pixline.
  * l1_io.c, added lonlat read to speed up lonlat2pixline.
  * l12_proto.h, added lonlat read to speed up lonlat2pixline.
  * main_lonlat2pixline.c, added lonlat read to speed up lonlat2pixline.
  * aerosols.c, moved inline function out for solaris.

### Data Added
  * common/ice_fraction.hdf, NSIDC ice fraction climatology.



## 5.9.6 - 2009-04-22

Various GIOP enhancements. Add additional morel CDOM products. Add MS epsilon
product. Add AERTEST eval switch, eliminate MODIS cloud mask switch. Support
new daily ozone climatology.

### Source Changed
  * l12_parms.h, version number update.
  * giop.h, add BBPSCIOTTI, BBPSMM01, add BBPTAB and ADGTAB.  Add aph_tab, adg_tab,
  and bbp_tab fields (file, n, w, s) and drop taph fields.
  * giop.c, add support for ciotti and mm01 bbp basis functions, add tabbed input
  option for aph, adg, bbp read from external files, drop aphw and aphs inputs.
  Reorganize options.
  * msl12_input.c, add giop_aph_, adg_, and bbp_file inputs and set defaults, revise
  giop usage messages, drop aphw and aphs input support.
  * input_struc.h, drop giop_aphw, giop_aphs, and add giop_aph_file, giop_adg_file,
  and giop_bbp_file.
  * aph.c, don't extrapolate aph spectra.
  * brdf.c, remove messages.
  * l2prod.h, add morel CDOM products, add scattang, ms_epsilon.
  * cdom_morel.c, new function for morel CDOM products.
  * prodgen.c, add morel CDOM products, add scattang.
  * get_l2prod_index.c, add morel CDOM products, scattang, ms_epsilon.
  * convl12.c, copy scattang to l2rec.
  * alloc_l1.c, allocate scattang field.
  * l1_struc.h, add scattang field.
  * l2_struc.h, add scattang field.
  * loadl1.c, compute scattang.
  * setflags.c, comment-out MODIS cloudmask switch.
  * aerosol.c, add get_ms_epsilon function, AERTEST switch for epsilon.
  * l12_proto.h, add get_ms_epsilon().
  * setanc.c, add support for new daily ozone climatology.

### Data Added
  * common/ozone_climatology.hdf, new daily ozone climatology.
  * aph_default.txt, new default aph spectra.
  * adg_default.txt, new default adg spectra.
  * bbp_default.txt, new default bbp spectra.


## 5.9.5 - 2009-04-07

Various GIOP enhancements, fix for bad pixels in MERIS N1, adding cdom_index
product.  Incorporporate SST changes from Sue Walsh (RSMAS) for evaluation.

### Source Changed
  * l12_parms.h, version number update.
  * l1_meris_N1.c, flag bad pixels as NAVFAIL and set radiances to zero. Silently mark
  products missing in L2 files.  Fix bug where reading L1 flags even for a L2 file.
  Allow L2 pixels flagged as COSMETIC to be included.
  * giop.c, add support for Bricaud aph, iop to Rrs via f/Q, QAA-like adg_s option, and
  OBPG adg option. Add aph via Ciotti.
  * giop.h, add new option defines.
  * aph.c, new source to read aph (Bricaud, Ciotti) tables and derive aph.
  * l12_proto.h, add aph_bricaud() and morel_fqint(), cdom_morel().
  * msl12_input.c, add new GIOP options.
  * l2gen.mk, add aph.o, cdom_morel.o.
  * l3gen.mk, add aph.o, cdom_morel.o.
  * getglint.f, comment-out Ebuchi & Kizu.
  * cdom_morel.c, new function to derive morel CDOM index.
  * get_l2prod_index.c, add CDOM_index.
  * prodgen.c, check to see if requested MERIS product was found in the L2 file. Add CDOM_index.
  * brdf.c, add function to compute morel Q.
  * get_tsm.c, file for hmodis sensors.
  * l1a_seawifs.c, change int to int32.
  * l1_modis_hdf.c, add eval detector corrections for 4um bands.
  * l1_hmodis_hdf.c, add eval detector corrections for 4um bands.
  * sst.c, add eval support for coefficients by latitude band, eval support for
  rvs correction of sst.
  * sstref.c, revised filling scheme for 0.25-deg Reynolds.
  * main_vcaltarget.c, remove isleap().
  * main_l2binmatch.c, remove isleap().

  * [libgenutils] isleap.c, added to lib and inc/utils/genutils.h.

### Data Added
  * eval/modisa/cal/sst4_sses_modisa.hdf, updated hypercube.
  * eval/modisa/cal/sst4_modisa.dat, updated coeffs with latitude bands.
  * eval/modisa/cal/sst_sses_modisa.hdf, updated hypercube.
  * eval/modisa/cal/sst_modisa.dat, updated coeffs with latitude bands.
  * common/aph_bricaud_1995.txt
  * common/aph_ciotti_2002.txt
  * common/morel_cdm_index.dat
  * common/morel_chl_R490_R555.dat
  * common/morel_q0.dat
  * common/morel_sq.dat


## 5.9.4 - 2009-03-11

Added ability to initialize first guess aerosol model in NIR iteration.
Expanded available MERIS output products from input MERIS Level-2.  Fixed bug
in SOA code.  SeaWiFS caltable support for gain 3 correction. Bug fix for fixed
AOT option.

### Source Changed
  * l12_parms.h, version number update.
  * l12_proto.h, pass aer_opt into aerosol code directly.
  * aerosol.c, pass aer_opt into aerosol code directly.
  * atmocor2.c, added logic to init NIR iteration if aer_angstrom is set, enhanced
  logic to locate red and green bands.
  * get_rhown_nir.c, enhanced logic to locate red and green bands.
  * HOWTO_Add_a_product.txt, updated the procrdure.
  * atmcor_soa.f, fixed a bug in INBAND that was taking a log of a neg number.
  * convl12.c, copy in_prods pointer form l1rec into l2rec
  * get_l2prod_index.c, added catalog entry for chl1_meris, chl2_meris, adg_meris, tsm_meris,
  par_meris, angstrom_meris and Taua_meris.
  * l1_meris_N1.c, read in 7 products out of the L2 MERIS files.
  * l2_struc.h, added in_prods pointer.
  * l2_prod.h, added contants for chl1_meris, chl2_meris, adg_meris, tsm_meris, par_meris,
  angstrom_meris and Taua_meris.
  * prodgen.c, added case statements for chl1_meris, chl2_meris, adg_meris, tsm_meris, par_meris,
  angstrom_meris and Taua_meris.
  * aerosol.c, bug fix for fixed AOT option (iwtab not set when no interp), introduced in 5.9.3.

  * [libseawifs] get_cal_swf.c, add support for gain 3 correction.


## 5.9.3 - 2009-02-25

New ipar code (not yet active). Preparing to eliminate carder-based ipar, arp,
cfe. Fix l1mapimg.  Changes to support MERIS SMILE corrections. Updates for
RH-based aerosol model selection, including alternative RHAERFRNIR option (-13)
for selection in MS reflectance space. Updated qaa to version 5.

### Source Changed
  * l12_parms.h, version number, NOAEROOB switch, RHAERFRNIR aeosol option.
  * ipar.c, new ipar function.
  * whitecaps.c, returns foam reflectance instead of tLf, so rhof can be used in ipar.
  * atmocor1.c, compute tLf from foam reflectance.
  * brdf.c, modified fresnel functions to return reflectance instead of relative correction
  on request, so solar and sensor fresnel reflectances could be used in ipar.
  * alloc_l1.c, storage for rhof as well as csolz and csenz (to reduce recomputes),
  added space for radcor and pixdet, for MERIS.
  * l1_struc.h, storage for rhof as well as csolz and csenz, added radcor and pixdet.
  * l2_struc.h, pointer for rhof as well as csolz and csenz, added radcor and pixdet.
  * convl12.c, copy new fields from l1 to l2 rec, handle MERIS L2 data.
  * loadl1.c, compute and store csolz, csenz.
  * ipar_arp.c, debugging statements (commented-out).
  * l12_proto.h, changes to whitecaps(), add fresnel_coef(), fresnel_sen(), fresnel_sol(),
  gothic_R(), get_ipar2().
  * l2gen.mk, add ipar.c.
  * main_l1mapimg.c, force subsampling back to 1.
  * aerosol.c, fix for fixedaot() when LUT and sensor waves don't match, set aeroob to return 1.0
  for NOAEROOB. Fix bug in aerosol extrapolation above 80-deg solar zenith.  Add franzaer() and
  aerosol_select_franz() functions. Don't extrapolate RH in rhaer(), and explicitly handle
  bounding case and low aerosol reflectance case. Add switch for RHAERFRNIR through static aer_opt.
  * alloc_l1.c, added space for radcor and pixdet, for MERIS.
  * atmocor2.c, apply radcor offsets, add AERRHFRNIR recognition.
  * convl12.c, transfer radcor and pixdet.
  * init_l1.c, init MERIS fields.
  * l12_proto.h, proto for radcor().
  * l1_meris_N1.c, add functions to read tables and compute SMILE correction, use the new
  rad_opt switch in the l1_struct, read L2 files. Add incremental scan time adjustment.
  * l1subpix.c, transfer pixdet and radcor and fix bug in memmove length.
  * loadl1.c, add call to new radcor().
  * get_qaa.c, change call from qaaf_v4() to qaaf_v5()
  * qaa.c, change to version 5
  * qaa.h, change to version 5
  * getformat.c, better detection of L1 and L2 MERIS files.
  * input_struct.h, added rad_opt variable to structure.
  * msl12_input.c, added rad_opt.

  * [libseawifs] get_cal_swf.c, add support for exponential + linear cal format.

  * [libmeris] epr_msph.c, epr_product.c, add ability to read MER_FRS files.


## 5.9.2 - 2008-12-18 

Add code to select alternate sensor info file based on eval switch 8192.

### Source Changed
  * l12_parms.h, version number, eval define ALTSENSORINFO.
  * rdsensorinfo.c, add path switches.


## 5.9.1 - 2008-12-15

Add code to speed-up aerosol model computations (limit to only models that are
used).

### Source Changed
  * l12_parms.h, version number.
  * aerosol.c, add tests for last model to ss_to_ms_coef(), model_epsilon(),
  model_phase().


## 5.9.0 - 2008-12-10

Add access to MODIS 1380nm band for CIRRUS correction/detection.

### Source Changed
  * l12_parms.h, version number.
  * convl12.c, cp rho_cirrus pointer.
  * l1_modis_hdf.c, read EV_Band26 and stuff in rho_cirrus.
  * l2prod.h, add rho_cirrus product.
  * alloc_l1.c, add space for rho_cirrus and cirrus flag.
  * setflags.c, check rho_cirrus and flag cloud if thresholds provided.
  * l1_struc.h, add  rho_cirrus and cirrus flag.
  * l2_struc.h, add rho_cirrus.
  * prodgen.c, add rho_cirrus support.
  * get_l2prod_index.c, add rho_cirrus support.
  * init_l1.c, add rho_cirrus support.
  * msl12_input.c, read cirrus_thresh inputs.
  * input_struc.h, add cirrus_thresh inputs.


## 5.8.9 - 2008-12-04

Add NIWA IOP model. Improved memory usage for aerosol tables.

### Source Changed
  * l12_parms.h, version number, eval define for MODIS cirrus flagging, IOPNIWA.
  * l12_proto.h, new modis_cirrus_flag proto.
  * cloud_flag.c, new modis_cirrus_flag function.
  * setflags.c, eval switch to mask on MODIS cirrus flag.
  * convl12.c, add NIWA to iop selector.
  * l2prod.h, add NIWA a, bb, flags.
  * prodgen.c, add NIWA a, bb, flags.
  * get_l2prod_index.c, add NIWA a, bb, flags.
  * msl12_input.c, add NIWA iop selector to usage.
  * main_l1brsgen.c, fix sline eline spixl epixl indexing.
  * main_l1mapgen.c, fix sline eline spixl epixl indexing.
  * alloc_l1.c, add verification check for allocation size.
  * alloc_l2.c, add verification check for allocation size.
  * aerosol.c, switch from static to dynamic memory allocation for aerosol
  static table (more efficient and allows valgrind to work).
  * filter.c, in epsilon filtering don't use apparent epsilon if > 2.0.

## 5.8.8 - 2008-11-05

Add new RH-based aerosol model selection option.

### Source Changed
  * l12_parms.h, version number, expand number of aerosol models to 100.
  * l2prod.h, add flags_qaa product, bbps_las (incomplete).
  * giop.c, clean-up, remove BBPSPML.
  * prodgen.c, add flags_qaa product.
  * get_l2prod_index.c, add flags_qaa product.
  * atmocor2.c, add support for RH-based aerosol model selector.
  * msl12_input.c, make sure number of models does not exceed max.
  * aerosol.c, add RH-based aersol model selection scheme and modify
  wangaer and model select wang to support it.
  * whitecaps.c, add eval switch to Stramska model.
  * setflags.c, add eval switch for spectral cloud test.
  * xcal.c, added switch to handle alternative SDS name M11.

### Data Added
  * data/eval/modist/cal/xcal_modist_rvs_m12_m13_WITHglintNorm_tfit.hdf,
  normalized M11 variant.
  * data/eval/modist/cal/xcal/xcal_modist.hdf, set to normalized variant.

## 5.8.7 - 2008-09-18

Revamp global attribute scene meta-data handling to consolidate and standardize
across applications. Some reorganization in prep for l3gen development. Rewrite
l2_hdf_generic logic to move product generation to standalone function
prodgen().  New l3gen program. Also, change scaling on bb products and change
GSM to set BAD_FLT on failure. Add aph center wavelength specifier to support
Hoge & Lyon gaussian aph* formulation within giop. Try to fix all the missing
prototypes. Add flags_giop product and associated iop flag code. Change
linearization scheme in giop.c.

### Source Changed
  * l12_parms.h, version number change, expand MAXPIX to 10000.
  * scene_meta.c, new scene_meta_put(), scene_meta_get(), scene_meta_write() functions.
  * scene_meta.h, protos and struct for same.
  * l12_proto.h, added scene_meta.h, removed time_utils.h (also in genutils.h),
  removed parse_file_name() proto.
  * l1_io.c, pass each new read of l1rec to the scene_meta object.
  * convl12.c, remove scene meta-data computation code.
  * l2_hdf_generic.c, remove explicit scene meta writes and call scene_meta_write(),
  remove product extract/generation logic and call new prodgen() function.
  * l1_hdf_generic_write.c, ditto.
  * main_l1brsgen.c, ditto, also reindent to fix the mess.
  * main_l1mapgen.c, use scene_meta object to get map boundaries, also reindent.
  * cpl1rec.c, transfer scan and pix range info.
  * l1_struc.h, add scan and pix range boundaries, needed for bounding scene meta.
  * l2_struc.h, remove scene boundary meta-data fields from structure.
  * getformat.c, change calling sequence for GetFileDesc(), add recognition of Level-3
  bin format.
  * filehandle.h, new FT_L3BIN format.
  * parse_file_name.c, moved to libgenutils.
  * l1_mos_hdf.c, moved GetFileDesc() to libhdfutils.
  * l1_mos_hdf.h, removed proto for same.
  * get_rhos.c, clean-up.
  * loadl1.c, clean-up, change error returns to exits, don't set land flag for binned input.
  * prodgen.c, new function to replace product extract/generation logic in writel2_hdf().
  * get_l2prod_index.c, change slope from 0 to 1.0 for some integer products. A slope
  of 1.0 and offset of 0.0 now indicates no-scaling-required. Changed offset on scaling
  of bb products. Add flags_giop product.
  * prodlist.c, new function to consolidate l2 product list parsing.
  * input_struc.h, added sensorID and format fields.
  * msl12_input.c, capture sensorID and format, saves secondary call to getformat(),
  fixed lots of sprintf formatting errors.
  * main_l3gen.cpp, main file for new l3gen program.
  * l3gen.mk, makefile for l3gen.
  * Makefile, add l3gen.mk.
  * *.mk, add new scene_meta.c module, prodgen.c, prodlist.c.
  * gsm.c, set badval instead of -1.0 or fit result on error.
  * giop.c, clean-up, change matrix linearization.
  * flags_iop.c, new functions to set IOP flag.
  * flags_iop.h, proto and struct for above.
  * cloud_flag.c, fixed missing type on lastScan.
  * l1a_osmi.c, fixed missing type on  chindx, new call1a_proto_osmi.h include.
  * src/l2gen/soa_sma.c, protos explcitly added for fortran subs.
  * main_l2binmatch.c, added read_attrs proto.
  * main_vcaltarget.c, added read_attrs proto.
  * l1a_seawifs.c, add geonav_ proto.
  * etc. etc. etc.
  * gsm.c, set PRODFAIL for already masked cases.
  * las_iop.c, set solz=0 of MORELFOQ is enabled.
  * get_qaa.c, remove limitation that bb must be 1.05*bbw.

  * [libbin++] new library for bin file i/o.

  * [libhdfutils] hdf_utils.c, add GetFileDesc() function from l1_mos_hdf.c, and consolidate static data.

  * [libgenutils] parse_file_name.c, moved from l2gen dir.
  * [libgenutils] time_utils.c, extracted all routines and put them out as separate files.

  * [libsdsgen] parse_file_name.c, old version, removed.

  * hdf_utils.h, removed.
  * genutils.h. removed.

  * genutils.h, moved from swfinc and consolidated other includes.
  * hdf_utils.h, moved from swfinc, added SetF64GA() proto.
  * time_utils.h, added gmt_offset proto.

  * call1a_proto_osmi.h, was actually identical to call1a_proto.h, including ifndef ID.
  Changed calibrate_l1a() to calibrate_l1a_osmi().


## 5.8.6 - 2008-07-21

Prefix format identifiers with FMT for clarity. Add support for new 0.25-deg
AVHRR SST reference fields.  Add initial support for MERIS L2 and MERIS L1B.
Remove extrapolation from xcal application.

### Source Changed
  * l12_parms.h, version number change.
  * getformat.c, prefix format identifiers with FMT.
  * main_vcaltarget.c, prefix format identifiers with FMT, pass n_inprods to alloc_l1().
  * main_lonlat2pixline.c, prefix format identifiers with FMT, pass n_inprods to alloc_l1().
  * main_vcalmerge.c, prefix format identifiers with FMT, pass n_inprods to alloc_l1().
  * main_l2gen.c, prefix format identifiers with FMT, pass n_inprods to alloc_l1().
  * main_l1bgen.c, prefix format identifiers with FMT, pass n_inprods to alloc_l1().
  * main_l1info.c, pass n_inprods to alloc_l1(), pass n_inprods to alloc_l1().
  * l1_io.c, prefix format identifiers with FMT.
  * filehandle.h, prefix format identifiers with FMT.
  * l1_octs_hdf.c, prefix format identifiers with FMT.
  * sstref.c, add 0.25-deg AVHRR SST file support.
  * l1_meris_N1.c, new src reads MERIS L1B or L2 into L1 rec.
  * l12_proto.h, add meris api and pass n_inprods to alloc.
  * l2_struc.h, add pixdet field (detector number in pixel space).
  * l1_struc.h, add MERIS fields (pixdet, n_inprods, in_prods, in_flags).
  * getl1rec.c, pass n_inprods (number of L2 products to read).
  * l1subpix.c, copy new meris fields.
  * getformat.c, recognize MERIS formats FT_MERIS and FT_MERISL2.
  * cpl1rec.c, copy n_inprods header info.
  * alloc_l1.c, allocate MERIS fields.
  * l1_io.c, add MERIS L1B and L2 read support.
  * get_chl.c, add default chl algorithm for MERIS.
  * get_l2prod_index.c, add default chl algorithm for MERIS.
  * l1img_input.c, fix bug in sprintf statement.
  * l1img_input.h, change threshold field to float (was incorrectly long).
  * mscal_struc.c, add proto.
  * *.mk, add meris object, lib, and inc.
  * xcal.c, don't extrapolate.
  * get_pml.c, add aw to atot.
  * las_iop.c, bug fix in Loisel-Stramski implementation.

### Data Added
  * eval/modisa/msl12_sensor_info.dat, new out-of-band coefficients.
  * eval/seawifs/msl12_sensor_info.dat, new out-of-band coefficients.

## 5.8.5 - 2008-06-23

Major changes to incorporate additional IOP models. Fix for OCTS extracts.
Eliminated CZCS-specific aerosol option and replaced with generic fixed model
pair option with NIR iteration.  Added fixed angstrom options to support next
Terra calibration test. Added new aer_angstrom input param and aerbinfile
param, where the later is a bin file that can be used for extracting aerosol
fields (e.g., angstrom in lieu of using aer_angstrom parameter). Added more
aerosol option help. Deleted old czcsaer function from aerosols.  Updated octs
to properly handle extracts from l2extract. Add flexible cal table format and
code for SeaWiFS.

### Source Changed
  * l12_parms.h, version number change, added IOP model definitions, MERISL2 ID, added FIXANGSTROM,
  FIXANGSTROMNIR, FIXMODPAIRNIR aerosol options and removed CZCSAER..
  * l1_octs_hdf.c, comment-out old subscene support, add new l2extract subscene support.
  * amoeba.h, add meta field to pass auxilliary info to fit function.
  * convl12.c, add switches for new IOP models for internal a & bb computation.
  * filehandle.h, add MERIS files for NRL.
  * get_l2prod_index.c, add new IOP model products, change gsm01 to support gsm ID.
  * get_pml.c, new wrapper for PML IOP model.
  * get_qaa.c, NRL updates to QAA IOP model wrapper.
  * giop.c, new generic IOP model.
  * giop.h, new generic IOP model.
  * gsm.c, renamed gsm01.c and removed all 01 references.
  * input_struc.h, added input fields for new and updated IOP models, and aerbinfile, aer_angstrom.
  * l12_proto.h, prototypes for new IOP models.
  * l2gen.mk, added new IOP model objects.
  * l2_hdf_generic.c, added support for new IOP model products.
  * l2prod.h, added definitions for new IOP model products.
  * las_iop.c, new Loisel & Stramski IOP model.
  * msl12_input.c, added new and updated IOP model options, aerbinfile, aer_angstrom.
  * pml_bright.h, Plymouth Marine Lab IOP model.
  * pml.c, Plymouth Marine Lab IOP model.
  * pml.h, Plymouth Marine Lab IOP model.
  * pml_iop_calculate.c, Plymouth Marine Lab IOP model.
  * pml_iop_calculate.h, Plymouth Marine Lab IOP model.
  * pml_iop_config.c, Plymouth Marine Lab IOP model.
  * pml_iop_config.h, Plymouth Marine Lab IOP model.
  * pml_iop.h, Plymouth Marine Lab IOP model.
  * pml_iop_tables.c, Plymouth Marine Lab IOP model.
  * pml_iop_tables.h, Plymouth Marine Lab IOP model.
  * qaa.c, NRL updates to QAA IOP model.
  * qaa.h, NRL updates to QAA IOP model.
  * soa_sma.c, lose gsm01 reference.
  * water.c, add new function get_aw_bbw().
  * loadl1.c, don't add 360 to longitudes of -999.
  * atmocor2.c, fix for CZCS iteration control when 670 goes slight negative, add FIXMODPAIR option
  with NIR iteration, remove CZCSAER option, add FIXANGSTROM and FIXANGSTROMNIR.
  * aerosol.c, remove CZCSAER option, add FIXANGSTROM and FIXANGSTROMNIR, add function fixedmodpair(),
  add checks for negative input radiances for same.
  * bin_climatology.c, upgraded to allow bin file name as input.
  * l1_octs_hdf.c, fixes to properly index sample table and navigation when reading extracts.
  * l1a_seawifs.c, add new cal table i/o to eval switch.
  * time_utils.c, moved to libgenutils.
  * time_utils.h, moved to new inc/utils

  * [libseawifs] get_cal.c, _new routines renamed, and old old cal table support removed.
  * [libseawifs] get_cal_misc.c, o_new routines renamed, and old old cal table support removed.
  * [libseawifs] calibrate_l1a.c,  _new routines renamed, and old old cal table support removed.
  * [libseawifs] get_cal_swf.c, code to support new seawifs cal table format.

  * [libgenutils] time_utils.c, moved from l2gen.

  * time_utils.h, moved from l2gen so visible to new libseawifs.

  * call1a_proto.h, _new routines renamed, and old old cal table support removed.
  * getcal_proto.h, _new routines renamed, and old old cal table support removed.

### Data Added
  * seawifs/cal/seawifs_cal_200707.hdf, operational cal table in new format.


## 5.8.4 2008-04-05

OCM1 implementation, to support the preliminary version of the OCM1 and OCM2
HDF format and the HDF format developed by NRL for OCM1 direct broadcast users
(from terascan files). Add proto for xcal function. Add eval option to sst to
use cold-only test relative to reference.

### Source Changed
  * l12_parms.h, version number, replace OCM identifier with OCM1..
  * l12_proto.h, add get_xcal().
  * filehandle.h, add sensor ocm1 and formats OCML1B and OCML1BDB.
  * l1_ocm_hdf.c, new functions to read OCML1B.
  * l1_ocm_hdf.h, new functions to read OCML1B.
  * l1_ocmdb_hdf.c, new functions to read OCML1BDB.
  * l1_ocmdb_hdf.c, new functions to read OCML1BDB.
  * get_qaa.c, replace OCM with OCM1.
  * getformat.c, recognize OCML1B* format.
  * l1_io.c, support for OCML1B* format.
  * get_chl.c, default chl algorithm for OCM set to OC4.
  * get_l2prod_index.c, default chl algorithm for OCM set to OC4.
  * setanc.c, fix comment on ozone units.
  * atmocor2.c, FIXNIRTG test was treated as float.
  * sst.c, switch to cold reference test on eval.
  * *.mk, add l1_ocm_hdf.o and l1_ocmdb_hdf.o.

### Data Added
  * ocm/aerosol, using SeaWiFS tables.
  * ocm1/msl12_defaults.par, using SeaWiFS parameters.
  * ocm1/l1img_defaults.par, rgb defaults.
  * common/smigen_product_table.dat, add sstref support.

## 5.8.3 - 2008-03-14

### Source Changed
  * l12_parms.h, version number, add AVHRR, remove OCM2.
  * aerosols.c, add interpolation for diffuse transmittance if
  table wavelengths differ from sensor.
  * l1_hmodis_hdf.c, get_modis_calfile() was called incorrectly.
  * l1_modis_hdf.h, duplicate proto for get_modis_calfile().
  * rdsensorinfo.c, return -1 on error, instead of 0, to distinguish error
  from the case of zero vis/nir bands.
  * l1_io.c, check for -1 on return from rdsensorinfo(), allow 0, call
  avhrr functions as needed.
  * msl12_input.c, check for -1 on return from rdsensorinfo(), allow 0.
  * l2_hdf_generic.c, don't write vis/nir band meta-data if no vis/nir bands.
  * getformat.c, recognize AVHRR file.
  * l1_pci_hdf.c, pathfinder class AVHRR i/o functions.
  * l1_pci_hdf.h, prototypes for above.
  * l1_struc.h, add some fields for AVHRR support.
  * filehandle.h, add file definition for AVHRR type.
  * msl12_input.c, add degc switch for AVHRR (?)
  * input_struc.h, add degc switch for AVHRR (?)
  * runcal.h, AVHRR calibration definitions.
  * rawcal.h, AVHRR calibration definitions.

### Data Added
  * hmodisa/rayleigh, replace with Wang Rayleigh tables.
  * hmodisa/aerosol, replace with Wang modisa tables, extended with SWIR from
  Ahmad. VIS/NIR land bands will be interpolated.
  * hmodisa/msl12_defaults.par, revert to operational vicarious cal.
  * modist/msl12_sensor_info.dat, updated to Thullier, all bandpass quants
  recomputed for consistency.
  * hmodist/msl12_sensor_info.dat, updated to Thullier, all bandpass quants
  recomputed for consistency.
  * modist/rayleigh, tables recomputed using revised Taur.
  * hmodist/rayleigh, tables recomputed using revised Taur.
  * hmodist/aerosol, replaced with Wang tables, extended to SWIR via Ahmad.
  * modisa/aerosol/modis_aerosol_par.dat, new table for MODIS PAR
  * modist/msl12_defaults.par, new vicarious gains
  * hmodist/msl12_defaults.par, new vicarious gains
  * ocm1/*, new sensor info directory for OCM1
  * avhrr/*, new sensor info directory for AVHRR


## 5.8.2 - 2008-02-25

Fix Terra mside corrections, which were mis-indexed after the addition of more
thermal bands in 5.7.9.  This change and the adjustment to the Terra brightness
temperature tables was installed retroactively back to 5.7.9 (current
production version). The brightness temperature tables were recomputed using an
alternate integration technique to achieve consistency with Miami processing.
Modify no2 correction to use precomputed no2 fraction above 200m. Set PRODFAIL
on negative flh, rather than just setting flh to 0.0.

### Source Changed
  * l1_modis_hdf.c, fix ordering of Terra mside corrections, add 3.7um Terra detector
  corrections.
  * l1_hmodis_hdf.c, fix ordering of Terra mside corrections, add 3.7um Terra detector
  corrections, add calfile tracking. Fix error in 500m mrad and brad i/o used for
  subframe destriping.
  * aerosols.c, fix error in angstrom computation that effects atmospheric correction
  with fixed aerosol optical depth.
  * l2_hdf_generic.c, enable iii, vvv, nnn product specifiers to work for all cases,
  add new no2_frac product.
  * convl12.c, copy no2_frac.
  * setanc.c, read no2 fraction.
  * l2prod.h, added no2_frac product.
  * alloc_l1.c, allocate no2_frac field.
  * l1_struc.h, add no2_frac field.
  * l2_struc.h, add no2_frac field.
  * atmocor1.c, changed gaseous_transmittance arg list to support no2_frac.
  * get_l2prod_index.c, added no2_frac product.
  * gas_trans.c, changed gaseous_transmittance arg list to support no2_frac.
  * l12_proto.h, changed gaseous_transmittance arg list to support no2_frac.
  * time_utils.c, force UTC in gmttime for solaris port.
  * fluorescence.c, set PRODFAIL of flh < 0.0 or nLw(6xx) < 0.01.

### Data Added
  * modist/cal/bt_modist.hdf, recomputed using simpson integration.
  * hmodist/cal/bt_hmodist.hdf, recomputed using simpson integration.
  * trop_f_no2_200m.hdf, global tropospheric NO2 fraction above 200m.

## 5.8.1 - 2008-01-11

Add support for cross-calibration file.  Carry MODIS LUT info to output meta
data.  Fix bug in flag percentage reporting when multiple output files
requested. Remove some non-standard fortran code from SOA functions to allow
build on Macs.

### Source Changed
  * l12_parms.h, version change.
  * gsm01.c, require at least 3 positive Rrs for processing to proceed.
  * l1img_input.c, fix usage statement.
  * l1_modis_hdf.c, capture calibration information from L1B file, put in calfile string.
  * main_l2gen.c, copy calfile string to output file handle.
  * l2_hdf_generic.c, write new "Calibration Data" global attribute.
  * xcal.c, reading new full-mission XCAL hdf file, interpolating in time.
  * xcal.h, prototype and defines.
  * loadl1.c, apply xcal rvs.
  * polcor.c, apply xcal m12, m13.
  * l12_proto.h, remove old get_xcal proto.
  * l2_struc.h, remove percent flags (perc_flags) field.
  * filehandle.h, add percent flags (perc_flags) field.
  * l2_hdf_generic.c, use filehandle perc_flags.
  * alloc_l2.c, remove perc_flags initialization.
  * filehandle_init.c, add perc_flags initialization.
  * atmcor_soa.f, seadas portability fixes.
  * *.mk, seadas portability fixes.

### Data Added
  * eval/modist/cal/add xcal_modist.hdf, new xcal file.
  * eval/modist/cal/add xcal_modist_*.txt, remove old xcal files.


## 5.8.0  - 2007-12-07 

Install latest SOA code and add model parameter passing, etc.  Generalize GSM01
in the process. Change wavelength range for standard GSM01 to incude red bands.
Add additional Stramski POC variants. New endian test from Martinolich
(eliminates ifdefs).  New and improved L1 browse and mapping codes from Bailey,
including parameter file control. Drop polder. Drop flat binary L1B file
format. Change SeaWiFS NAVFAIL to NAVWARN for case of nflag[0] set, but check
all sensors for invalid geolocation. Change sign of polarization frame rotation
angle to conform with Gordon and Meister papers. Fix bug in aer_model_* output
that was causing corruption of first 3 pixels.

### Source Changed
  * l12_parms.h, version change, drop polder.
  * atmcor_soa.f, update from Kuchinke, made changes for parameter passing and multi-sensor
  input file location.
  * soa_sma_utils.f, update from Kuchinke.
  * get_qaa.c, changed QAA adg slope parameter name for consistency.
  * soa_sma.c, new SOA wrapper with model parameter passing, includes SMA (disabled).
  * soa.c, replaced with above.
  * convl12.c, removed SEADAS conditional compile, added SMA test.
  * l2_hdf_generic.c, add SMA products and one SOA product, fix badval scaling bug.
  * l2prod.h, add SMA products and one SOA product.
  * get_l2prod_index.c, add SMA products and one SOA product.
  * msl12_input.c, support for new GSM01 inputs adg_s, bbp_s, and aphw, aphs.
  * input_struc.h, support for new GSM01 inputs adg_s, bbp_s, and aphw, aphs.
  * gsm01.c, support input model params, extend wavelength range, clean-up.
  * l12_proto.h, atmcor_soa() expanded argument list.
  * get_poc.c, add 443/555 and 490/555 variations of Stramski.
  * b128_msk_get.c, replace ifdef with endianess().
  * sstref.c, replace ifdef with endianess().
  * main_l1mapgen.c, clean-up and add parameter file control.
  * main_l1brsgen.c, clean-up and add parameter file control.
  * l1_imgscale.c, consolidated RGB scaling code.
  * l1img_input.c, new input function for l1brsgen and l1mapgen.
  * main_l2binmatch.c, add nscenes field, remove bindx.
  * msl2val_struc.h, add nscenes field, remove bindx.
  * msl2val_hdf.c, add nscenes field, remove bindx.
  * calcite.c, remove bindx.
  * myprod.c, remove bindx.
  * atmocor1_land.c, remove bindx.
  * getformat.c, drop polder.
  * l1_io.c, drop polder, drop SIMBIOS generic L1B, add test for bad geolocation and
  set navfail flag accordingly..
  * get_chl.c, drop polder.
  * get_Kd.c, drop polder.
  * filehandle.h, drop polder, drop SIMBIOS generic L1B.
  * get_l2prod_index.c, drop polder.
  * msl12_input.c, drop polder.
  * l1a_seawifs.c, change nflag[0] condition from NAVFAIL to NAVWARN.
  * time_utils.c, add timegmt function to fix portability issue with extern timezone.
  * aerosol.c, fail fixed model case when input reflectance is negative.
  * main_l2binmatch.c, add nobs and nscene to validation output file.
  * l1_modis_hdf.c, change sign on alpha.
  * polcor.c, change sign on polarization equation.
  * time_utils.c, replace setenv with putenv for solaris compatibility.
  * main_l1mapgen.c, misplaced declaration.
  * l2_struc.h, switch aermodmin and aermodmax to long.
  * alloc_l2.c, switch aermodmin and aermodmax to long.
  * get_l2prod_index.c, switch aer_model_min and aer_model_max output to long.
  * l1info.mk, drop polder, drop l1_generic.
  * vcaltarget.mk, drop polder, drop l1_generic.
  * l1bgen.mk, drop polder, drop l1_generic.
  * lonlat2pixline.mk, drop polder, drop l1_generic.
  * l2gen.mk, replace soa.c with soa_sma.c, add libgenutils, drop polder, drop l1_generic.
  * l1brsgen.mk, l1_imgscale.o, l1img_input.o, libgenutils, drop polder, drop l1_generic.
  * l1mapgen.mk, l1_imgscale.o, l1img_input.o, libgenutils, drop polder, drop l1_generic.

  * [libgenutils] new library for general utilities
  * [libgenutils] endianess.c, returns 1 if big endian.

## defunct code
  * src/libpolder/*
  * inc/polinc/*
  * src/l2gen/l1_polder.h
  * src/l2gen/l1_polder.c
  * src/l2gen/l1_generic.c

### Data Added
  * */l1img_defaults.par, L1 browse and mapping defaults by sensor.
  * eval/seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200802



## 5.7.9 - 2007-10-24

Generalize handling of thermal bands to support varying number of bands.  Add
MODIS bands 20,27,28,29.  Add generalize BT output capability (by sensor
wavelength).  Add Lt output for thermal bands.  Add "_iii" support to indicate
all thermal bands. Also, did some clean-up of bindx and #defines.  Add
poc_stramski. Fixes from Martinolich.

### Source Changed
  * l12_parms.h, version change, increase NBANDSIR from 4 to 8. Eval switch for alternate
  sensor files reset to 32.
  * aerosol.c, protect MAXBAND macro with parenthesis.
  * cpl1rec.c, copy nbandsir.
  * filehandle.h, add nbandsir field.
  * filehandle_init.c, init new  nbandsir field.
  * get_l2prod_index.c, eliminate bindx, add BT_nnn and thermal Lt support, change lam to wave
  for some reason, separate numBands into thermal and visnir. Add poc_stramski and change units.
  * ipar_arp.c, eliminate bindx in call to get_l2prod_index and elsewhere.
  * l12_proto.h, get_l2prod_index() and bindex_get().
  * l1_hmodis_hdf.c, read thermal bands, eliminate bindx.
  * l1_io.c, call  bindex_set with nbands+nbandsir.
  * l1_modis_hdf.c, read thermal bands, eliminate bindx.
  * l1_struc.h, add nbandsir field.
  * l2_hdf_generic.c, add support for BT_nnn, Ltir, and  *_iii. Include all wavelengths in meta-data.
  Remove bindx. Add poc_stramski support.
  * l2struc.h, add nbandsir field.
  * l2prod.h, add BT_* product (phase-out old BT products, eventually). Add poc_stramski.
  * loadl1.c, remove bind, remove redundant setting windex (bindex_set)..
  * main_l2gen.c, add nbandsir (number of emmissive bands).
  * msl12_input.c, remove bindx.
  * msl2val_hdf.c, remove bindx from get_l2prod_index call.
  * rdsensorinfo.c, separate numBands into nbands (reflective) and nbandsir (emissive).
  * sst.c,  locate band index for BTs via bindex_get() instead of hard-coded, generalize
  code to handle sensors with different band combos (e.g., AVHRR).
  * water.c, add parens around #define NAWTAB to be safe.
  * windex.c, eliminate redundant windex_get, rename windex_set to bindex_set.
  * get_poc.c, add stramski algorithm, generalize clark, change units to mg/m^3.
  * mscal_struc.h, remove bindx from get_l2prod_index call.
  * alloc_aer.c, add memset proto via  <string.h>
  * alloc_l2.c, add memset proto via  <string.h>
  * alloc_target.c, add memset proto via  <string.h>
  * rayleigh.c, close SDS before exiting on error, remove bindx.
  * atmocor1.c, remove bindx.
  * setflags.c, remove bindx.
  * nlw_outband.c, remove bindx.
  * glint.c, remove bindx.
  * aerosol.c, remove bindx from arg list.
  * atmocor2.c, remove bindx from nlw_outband, glint_rad, and aerosol calls.
  * filter.c, remove bindx.

### Data Added
  * modisa/msl12_sensor_info.dat, add 8 thermal bands to wavelength list.
  * modist/msl12_sensor_info.dat, add 8 thermal bands to wavelength list.
  * hmodisa/msl12_sensor_info.dat, add 8 thermal bands to wavelength list.
  * hmodist/msl12_sensor_info.dat, add 8 thermal bands to wavelength list.
  * modisa/cal/bt_modisa.hdf, regenerated with 8 thermal bands.
  * modist/cal/bt_modist.hdf, regenerated with 8 thermal bands.
  * hmodisa/cal/bt_hmodisa.hdf, regenerated with 8 thermal bands.
  * hmodist/cal/bt_hmodist.hdf, regenerated with 8 thermal bands.
  * hmodisa/rayleigh/*, replaced old Zia tables with mix of Wang Aqua
  and Terra tables.


## 5.7.8 - 2007-10-02

Updates to vicarious calibration support. New capability to expand a products
specifier of the form ppp_nnn to product ppp for all sensor bands and ppp_vvv
to product ppp for all sensor visible bands. New vcalmerge to make use of l2gen
vcal capabilities (rather than reading the original and vicarious L1B).

### Source Changed
  * l12_parms.h, version number change.
  * target_struc.h, add mode field to indicate if record is enabled.
  * l2_hdf_generic.c, code to expand product list for ppp_nnn and ppp_vvv.
  * convl21.c, set in situ solz to sensor solz if in situ val is negative.
  * main_vcaltarget.c, lots of changes to support match-up with Level-3 input.
  * l2_struc.h, add pointer to calibration target record.
  * vcal.c, pass calibration target record.
  * main_l2gen.c, simplify handling of target file and don't require inverse mode.
  * get_l2prod_index.c, scale all transmittance products to INT16.
  * KDtree.c,
  * KDvector.h,
  * vcaltarget.mk, add Level-3 reader.
  * main_vcalmerge.c, read Level-2 file rather than two Level-1 files.
  * getformat.c, add identification of sensor ID for L2 files
  * filehandle.h, add names for OCM and MERIS.
  * mscal_struc.c,
  * mscal_struc.c,
  * vcalmerge.mk,

moved to archive
  * libcalfit, main_vcalgen.c, vcalgen.mk, mscal_func.c

## 5.7.7 - 2007-09-18

Add Kd_PAR_lee product. Initial changes to handle 17-band aerosol tables. QAA fix from
Martinolich (rrs vs Rrs). Changed most non-standard product algorithms to use PRODFAIL
rather than CHLFAIL to indicate algorithm failure.

### Source Changed
  * l12_parms.h, version number change.
  * l2_hdf_generic.c, add support for Kd_PAR_lee, add perc_flags logic (from convl12).
  * get_l2prod_index.c, add support for Kd_PAR_lee.
  * l2prod.h, add support for Kd_PAR_lee.
  * get_Kd.c, add support for Kd_PAR_lee, add PRODFAIL.
  * get_chl.c, add PRODFAIL, remove old Fo logic.
  * qaa.c, change rrs to Rrs.
  * aerosol.c, expand max bands to 17.
  * photic_depth.c, add special case of 1st optical depth.
  * convl12.c, remove perc_flags logic.
  * fluorescence.c, change CHLFAIL to PRODFAIL, change init value from 0.0 to BAD_FLT.
  * get_qaa.c, change CHLFAIL to PRODFAIL.
  * get_poc.c, change CHLFAIL to PRODFAIL.
  * get_tsm.c, change CHLFAIL to PRODFAIL.
  * myprod.c, change CHLFAIL to PRODFAIL.
  * turbid.c, change CHLFAIL to PRODFAIL.
  * ipar_arp.c, change CHLFAIL to PRODFAIL.
  * gsm01.c, change CHLFAIL to PRODFAIL.
  * carder.c, change CHLFAIL to PRODFAIL.

### Data Added
  * */l2bin_defaults.par, add PRODFAIL.


## 5.7.6 - 2007-09-05

Bug fix for handling timezone settings.  MSl1info was giving bad start dates in
som cases if TZ was set.

### Source Changed
  * l12_parms.h, version number change
  * time_utils.c, use timezone (extern variable) to get offset between GMT and local time.

### Data Added
  * hmodisa/rayleigh/*.hdf, new Rayleigh tables fix nadir singularity (zipper).

## 5.7.5 - 2007-08-16

QAA V4.0 updates from NRL. Uses 645 channel for MODIS when available. Updates to Lee
photic depth algorithm.

  * l12_parms.h, version number change
  * get_qaa.c, updates from NRL.
  * qaa.c, updates from NRL.
  * qaa.h, updates from NRL.
  * l2prod_struc.h, add badData value field to l2prod table..
  * get_l2prod_index.c, add badData value field, add function get_l2prod_badval() to retrieve
  badData value from catalog (not currently used).
  * l12_proto.h, get_l2prod_badval().
  * photic_depth.c, updates to Lee photic depth algorithm, change badData value.
  * l2_hdf_generic.c, add b_qaa and c_qaa products.
  * get_l2prod_index.c, change all Z products to scaled int, handle variable photic depth
  specifiers for Lee algorithm, add b and c qaa products.

## 5.7.4 - 2007-08-14

Initial implementation of new strategy for tracking invalid retrievals.

### Source Changed
  * l12_parms.h, version number change, NEWFILL eval bit (1), consolidate current suite of
  bad-value definitions, add new suite of type-specific bad values and fill values.
  * l12_proto.h, new init_l1() and init_l2() functions.
  * init_l1.c, initialize l1rec fields based on data type.
  * init_l2.c, initialize l2rec fields based on data type.
  * convl12.c, call init_l2() to clear the L2 record, if eval=1.
  * l1_io.c, call init_l1() to clear the L1 record, if eval=1.
  * l2_hdf_generic.c, add logic to insert bad_value and FillValue meta data into Level-2
  product attributes if eval=1 switch enabled.
  * atmocor2.c, replace all instances of CHL_BAD with badchl, and set badchl to BAD_FLT if
  eval=1. Same for BADVAL to badval.
  * sst.c, set missing sstref to BAD_FLT rather than SST_BAD.
  * get_chl.c, replace CHL_BAD with BAD_FLT if eval=1.
  * setflags.c, change CHLFAIL test to look for BAD_FLT rather than CHL_BAD, if eval=1.
  * alloc_l2.c, changed aermodmin and aermodmax from long to char, consistent with output.
  * l2_struc.h, changed aermodmin and aermodmax from long to char, consistent with output.
  * get_Kd.c, replaced all instances of KD_BAD with kdbad, and set kdbad to BAD_FLT when
  eval=1.
  * sstref.c, init sstref to BAD_FLT instead of SST_BAD.
  * loadl1.c, init sstref to BAD_FLT instead of SST_BAD.
  * filter.c, replaced all instances of BADVAL with badval and set to BAD_FLT if eval=1.
  * carder.c, set empirical chl to BAD_FLT on failure, if eval=1.
  * soa.c,  set badval to BAD_FLT is eval=1.
  * l2prod.h, added cat entries for new QAA prods b and c.
  * MSl1brsgen.c, limit pixel range from 0 to npix-1.
  * *.mk, add init_l1.o, init_l2.o.


## 5.7.3 - 2007-08-03

Bug fixes.

### Source Changed
  * l12_parms.h, version number change.
  * l1_czcs_hdf.c, fix problem that happens on SGIs with the log(negative).
  * l1a_seawifs.c, fix buffer overflow error reading inst_ana for new cal
  table format.


## 5.7.2 2007-07-11

Preliminary implementation of ZP Lee photic depth algorithm.  Add
stand-alone version of Carder empirical chlorophyll. Moved hdf
utilities to separate library. Added some code to free memory.

### Source Changed
  * l12_parms.h, version number change.
  * l2_hdf_generic.c, add chl_carder_empirical and Zeu_lee and Zsd_lee.
  * l2prod.h, add chl_carder_empirical and Zphotic.
  * get_l2prod_index.c, add chl_carder_empirical and Zeu_lee and Zsd_lee.
  * carder.c, add chl_carder_empirical().
  * photic_depth.c, add Zphotic_lee().
  * l12_proto.h, added chl_carder_empirical().
  * alloc_l1.c, comments.
  * alloc_l2.c, added free_l2().
  * getl1rec.c, generalize free_l1q() so it can be called externally.
  * MSl12.c, call free_l1(), free_l2(), and free_l1q() at end of program.
  * hdf_utils.c, moved to libhdfutils.a.
  * hdf_utils.h, moved to inc/swfinc.
  * *.mk, added linking to libhdfutils.a.


## 5.7.1 - 2007-07-11

New HMODISA gains and subframe striping corrections.

### Data Added
  * hmodisa/msl12_defaults.par, new gains.
  * hmodisa/cal/subframe*, revised indexing.

## 5.7.1 - 2007-07-09

Make new caltable the default

### Source Changed
  * l12_parms.h, version number change.
  * l1a_seawifs.c, switch eval tests to use new caltable 200707 as default.

### Data Added
  * data/seawifs/msl12_defaults.par, new 200707 caltable and gains.
  * data/seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200707, moved from eval path.

## 5.7.0

### Source Changed
  * l1_hmodis_hdf.c, add check for subsetted L1 file and exit with message.

## 5.7.0 2007-06-12

Various minor changes for SeaWiFS 5.2 Reprocessing. Switch to extended
Thuillier table for all sensors.  Make permanent a minor fix for handling
out-of-band water-vapor correction. Apply polarization correction when
computing tLw for vicarious calibration.  Add reference temperature to SeaWiFS
cal table. Discontinue tracking of cal-table fields in Level-2.

### Source Changed
  * l12_parms.h, version number change, remove eval bits AEROOBFIX, THUILLIEREXT, NEWPOLCAL, NEWOOB.
  * l2_flags.h, eliminate BADANC and TRICHO definitions.
  * setflags.c, moved isTricho() function to new tricho.c, disable setting of tricho flag.
  * loadl1.c, make extended Thullier F0 the default (slight effect on Rrs, Chl).
  * aerosol.c, apply out-of-band water-vapor correction to rhoa, not to log(rhoa), in SS to MS
  conversion, to be consistent with MS to SS (changes rhoa slightly).
  * convl21.c, apply polcor when computing tLw (changes MODIS vicarious cal slightly).
  * l1a_seawifs.c, add reference temperature to EVAL caltable.
  * tricho.c, new source file, not active (need to rewrite for multi-sensor).
  * l2_hdf_generic.c, drop vdata for SeaWiFS-specific calibration info, fix units on F0.
  * l12_seawifs.c, don't copy calibration-table fields to L2.
  * nlw_outband.c, remove "new" OOB switch.
  * l1_hdf_generic.c, drop vdata for SeaWiFS-specific calibration info.
  * brdf.c, fix error in dtran_brdf() call (no effect).

  * [libseawifs] calibrate_l1a.c, add reference temperature to caltable.
  * [libseawifs] get_cal.c, add reference temperature to caltable.
  * [libseawifs] get_cal_misc.c, add reference temperature to caltable.

  * swfinc/getcal_proto.h, add reference temperature to caltable.
  * swfinc/get_cal.h, add reference temperature to caltable.

### Data Added
  * common/f0_table_thuillier.dat, defunct.
  * seawifs/l2bin_defaults.par, add HIGLINT.
  * seawifs/msl12_defaults.par, new default gains and cal table (commented-out).
  * seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200707, new cal table with reference temperature.


## 5.6.8 - 2007-06-04

Add aerosol option for NIR-based to SWIR-based switching.

### Source Changed
  * l12_parms.h, version number change, remove MORELTFLAG eval switch.
  * turbid.c, new source file to store turbidy index algorithms.
  * atmocor2.c, add support for NIR to SWIR switch (aer_opt=AERWANGSWIR=-9)
  * aerosol.c, add informational message for new aerosol option.
  * setflags.c, moved tindx_morel() function to new turbid.c.
  * l2prod.h, add support for tindx_shi product.
  * get_l2prod_index.c, add support for tindx_shi product.
  * l2_hdf_generic.c, add support for tindx_shi product.
  * l12_proto.h, add tindx_shi() function.
  * MSl12.mk, add turbid.c.
  * MSscanpixlimits.c, reset min scan & pix to -1 if region not found (6.06).

## 5.6.7 - 2007-05-02

Provide ability to specify Rrs of aerosol selection bands, principally for
calibration in reflective waters. Other minor fixes and updates.

### Source Changed
  * l12_parms.h, version number change, add eval macro FIXNIRTG=65536.
  * l2_hdf_generic.c, internal change to variable Fo, to improve readability.
  * euphotic_depth.c, updated Zsd coefficients from Morel.
  * get_Kd.c, look for 551 rather than 555,to handle 547 shift.
  * input_struc.h, add aer_rrs_short and aer_rrs_long fields.
  * msl12_input.c, read aer_rrs_short and aer_rrs_long parameters.
  * atmocor2.c, add support for NIR correction via input Rrs, remove swirLw option, remove bindx,
  fix NIR correction options to remove double accounting for gaseous transmittance (used eval
  switch 65536 for standard NIR correction).
  * aerosol.c, change informational message.

## 5.6.6 - 2007-04-16

Crazy number of changes required to pass eval switch to rdsensorinfo().  This
was done to allow testing of alternate band-pass and nominal band center
effects.

### Source Changed
  * l12_parms.h, version number change.  NEWSENSINFO eval switch (16384)
  * get_qaa.c, pass evalmask to rdsensorinfo().
  * rayleigh.c, ""
  * rdsensorinfo.c, ""
  * l2_hdf_generic.c, ""
  * l1_hmodis_hdf.c, ""
  * calcite.c, ""
  * whitecaps.c, ""
  * l1_modis_hdf.c, ""
  * l1_hdf_generic_write.c, ""
  * l1_io.c, ""
  * nlw_outband.c, ""
  * aerosol.c, ""
  * loadl1.c, ""
  * ipar_arp.c, ""
  * msl12_input.c, ""
  * gas_trans.c, ""
  * gsm01.c, ""
  * l1_hdf_generic_read.c, ""
  * atmocor1.c, pass evalmask to whitecaps() and gaseous_transmittance().
  * msl2val_struc.h, add evalmask to input structure.
  * MSl2validate.c, read evalmask param and pass to rdsensorinfo().
  * l12_proto.h, pass evalmask in rdsensorinfo(),whitecaps(),gaseous_transmittance().

### Data Added
  * eval/modisa/msl12_sensor_info.dat, eval sensor info with alternate nLw OOB coefficients.
  * eval/modist/msl12_sensor_info.dat, eval sensor info with alternate nLw OOB coefficients.
  * eval/hmodisa/msl12_sensor_info.dat, eval sensor info with alternate nLw OOB coefficients.
  * eval/hmodist/msl12_sensor_info.dat, eval sensor info with alternate nLw OOB coefficients.
  * eval/seawifs/msl12_sensor_info.dat, eval sensor info with alternate nLw OOB coefficients.


## 5.6.5 - 2007-04-05 

Fix for negative radiances in brightness temperature conversion, which are
occurring for Terra.

  * l12_parms.h, version change.
  * brightness.c, trap zero or negative radiances before taking log.

## 5.6.4 - 2007-04-02

Fix xcal (eval=8192) to handle negative mside in terra. Also, to be safe,
changed modis i/o to set navfail for negative mside and then set mside to ## 0.
The former was already true, as negative mside was also bad lon/lat. The latter
is to protect against any subsequent array indexing to a bogus mside.

### Source Changed
  * l12_parms.h, version change.
  * loadl1.c, check for navfail before using mside (terra has mside -1).
  * l1_modis_hdf.c, change negative mside to 0 and set navfail.
  * l1_hmodis_hdf.c, change negative mside to 0 and set navfail.
  * polcor.c, limit mside to 0 to 1 (may have caused problems).

libsdsgen
  * parse_file_name.c, disallow ~ (old change which got lost).


## 5.6.3 - 2007-03-29

Fix MODIS times to include milliseconds. This will cause minor changes in
standard products due to change in F0.  It also fixes the granule start/end
times in meta-data. Aerosol corrections were modified to support fixed table
indexing for the case of more table bands than sensor bands (e.g., 9-band MODIS
with 16-band tables). OCTS processing may have been effected by a bug in the
reference time used for temporal gain corrections.

### Source Changed
  * l12_parms.h, version change.
  * l1_modis_hdf.c, fixed truncation error in msec field.
  * l1_hmodis_hdf.c, fixed truncation error in msec field.
  * aerosol.c, enable interp when more table bands than sensor bands.
  * convl21.c, don't try to compute chl for inversion when calibrating to zero reflectances.
  * l1_octs_hdf.c, refjsec was not declared static, this should have caused problems.


## 5.6.2 - 2007-03-08 

Calcite table and code updates from Balch & Bowler. Reduces calcite by 3-4x.
Bug fixes for NO2 climatology and vcal_opt. Allow processing of ocean products
over land if maskland is disabled.

### Source Changed
  * l12_parms.h, version change. New eval trigger for MODIS xcal corrections.
  * calcite.c, change PIC-specific backscatter. Set minimum for values the retrieve 0 or negative,
  and don't set CHLRANGE flag.
  * setanc.c, fix day argument to no2conc() and add _ to SDS names.
  * msl12_input.c, revised vcal_opt determination.
  * convl12.c, apply ocean algorithms over land if maskland=0.
  * dtran_brdf.f, updates from H. Gordon to improve efficiency of diffuse transmittance.
  * brdf.c, changed calling sequence for dtran_brdf().
  * l12_proto.h, new get_xcal().
  * polcor.c, apply alternate xcal polarization coefficients if requested.
  * loadl1.c, apply xcal RVS correction if requested.
  * xcal.c, new functions to read xcal coefficients.
  * MSl1tcpcbox.mk, add xcal.c to makefiles.
  * MScalmerge.mk, add xcal.c to makefiles.
  * MSl1brsgen.mk, add xcal.c to makefiles.
  * MSl12.mk, add xcal.c to makefiles.

### Data Added
  * common/coccolith.dat, updated table from Balch.
  * eval/modist/cal/xcal_modist_rvs.txt, xcal-derived RVS adjustment.
  * eval/modist/cal/xcal_modist_m12.txt, xcal-derived polarization sensitivities.
  * eval/modist/cal/xcal_modist_m13.txt, xcal-derived polarization sensitivities.

## 5.6.1 - 2007-02-26

Adding ability to generate vicarious calibration gains and associated
components as Level-2 products (with capability to produce vicarious Level-1
format retained).  Add support for NO2 climatology file.  Update Morel Kd
algorithm.

### Source Changed
  * l12_parms.h, version change, gas_opt bit defines.
  * l12_proto.h, changes to convl21(), ocbrdf(), get_rhown_nir(), new vcal().
  * input_struc.h, add vcal_nLw, vcal_Lw, vcal_chl, vcal_solz, vcal_opt.
  * msl12_input.c, add support for vcal params and init no2file to climatology.
  * vcal.c, new code to compute vicarious gains and auxiliarry data.
  * vcal_struc.h, structure to store auxilliary cal data for one scan.
  * alloc_vcal.c, allocate vicarious calibration auxilliary record.
  * convl21.c, pass output Lt array rather than pointer to l1rec, eliminate few references
  to l1rec (fsol), disable use of nominal band F0 (no effect of outband_opt). Use input
  chl to compute brdf, if supplied.
  * brdf.c, add chl argument. If valid chl is passed, use it to compute f/Q rather than iterating
  on Rrs and the default chl algorithm.
  * get_rhown_nir.c, pass aw and bbw rather than l2rec.
  * atmocor2.c, changes in arguments to ocbrdf() and get_rhown_nir().
  * MSl12.c, change old convl21() call to pass pointer to Lt.
  * get_Kd.c, replace Kd_morel with new (490-specific) LOV+NOMAD algorithm.  Disable old Kd_morel
  and Kd_morel_nomad.
  * gas_trans.c, moved gas_opt bit defines to l12_parms.h.
  * setanc.c, add no2 climatology support, only read no2file if gas_opt no2 bit set.
  * l2prod.h, remove morel_nomad, add vgain, vLt, vLw, vnLw, vtLw, vbsat, vbtgt
  * l2_hdf_generic.c, support for above.
  * get_l2prod_index.c, support for above.
  * MSscanpixlimits.c, subtract extract pixel offset when reporting pixel range from extracted
  MODIS GEO file.

### Data Added
  * hmodisa/msl12_defaults.par, remove NO2 bit from default gas_opt.
  * hmodist/msl12_defaults.par, remove NO2 bit from default gas_opt.


## 5.6.0 - 2007-02-16

MUMM updates and fixes based on review by Ruddick et al. Also MSscanpixlimits updated.

### Source Changed
  * l12_parms.h, version change.
  * get_l2prod_index.c, fix bug in MUMM rhom (was returning rhos)
  * mumm.c, don't compute rhom if pixel is L1 masked.
  * atmocor1_land.c, initialize polcor to 1.0.
  * loadl1.c, initialize polcor to 1.0.
  * MSscanpixlimits.c, removed and added back case of geo and l1b supplied.  Usage message
  updates. Version 6.01.

### Data Added
  * */msl12_defaults.par, add defaults for mumm_alpha, mumm_gamma, mumm_epsilon.


## 5.5.9 - 2007-01-17

Adding support for Gordon's BRDF correction to diffuse transmittance. Adding
alternate Levenberg Marquart fit option to GSM01.  Fixing est/est meta-data
error for CZCS.

### Source Changed
  * l12_parms.h, version change, BRDF option for adding Gordon.
  * dtran_brdf.f, new module direct from Gordon.
  * brdf.c, add dtran_brdf() function, remove bindx from all funcs while at it..
  * convl21.c, add ip and remove bindx from ocbrdf argument list.
  * msl12_input.c, update usage for brdf_opt, add gsm01_fit.
  * atmocor2.c, add ip and remove bindx from ocbrdf argument list.
  * l12_proto.h, add ip and remove bindx from ocbrdf argument list.
  * setflags.c, minor tweak for MODIS cloudmask file option.
  * convl12.c, fix for CZCS easternmost, westernmost meta data.
  * input_struc.h, add gsm01_fit.
  * gsm01.c, add model functions for LM fitting.
  * MSl12.mk, add dtran_brdf.f, add linking to gsl library.

### Data Added
  * eval/common/dtran_brdf/, new files for BRDF correction to diffuse transmittance.

## 5.5.8 - 2006-11-29

Change Kd_morel table to use MED coefficients for 490 (was using NOMAD). Add
Kd_490_morel_nomad product to access NOMAD coefficients. Remove 2nd optical
depth Kd_PAR_2_morel product, but add Zhl_morel product (depth of heated
layer).  Limit all Kd products to range 0.016-6.4.  Add Morel turbidity index
(tindx_morel) and use to set TURBIDW flag as eval=4096 option.  Store Rrs in
l2rec within atmocor2().

### Source Changed
  * l12_parms.h, version change, eval option to switch to Morel for TURBIDW flag..
  * l2_hdf_generic.c, add Zhl_morel, Kd_490_morel_nomad, Kd_490_morel_ok2, tindx_morel.
  * l2prod.h, add Zhl_morel, Kd_490_morel_nomad, Kd_490_morel_ok2, tindx_morel.
  * setflags.c, add Morel turbid-water index and option to set flag.
  * get_Kd.c, add Zhl_morel, Kd_490_morel_nomad, Kd_490_morel_ok2. Impose range limits.
  * get_l2prod_index.c, add Zhl_morel, Kd_490_morel_nomad, Kd_490_morel_ok2, tindx_morel.
  * l12_proto.h, add tindx_morel().
  * convl12.c, eliminate call to get_rrs().
  * atmocor2.c, store Rrs result in l2rec.
  * MSl12.mk, remove get_rrs.c.
  * get_rrs.c, defunct.

### Data Added
  * eval/seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200703, 2 temperature epochs + gain drift (take 2).


## 5.5.7 - 2006-11-14

Disable flagging of MODIS thermal bands for navfail or hilt. Fix MSl1info to
handle fully navfailed granules. Add Morel's Kd(PAR) and photic depth products.
Another SeaWIFS cal table for eval.

### Source Changed
  * l12_parms.h, version change.
  * MSl1info.c, improved NAVFAIL handling.
  * l1_modis_hdf.c, just set flagged thermal radiances to zero.
  * l1_hmodis_hdf.c, just set flagged thermal radiances to zero.
  * MSl1scanpixlimits.c, handle MODIS geolocation files.
  * filehandle.h, add MODIS geolocation file identifier.
  * get_format.c, identify MODIS geolocation file.
  * l2_hdf_generic.c, support Kd_PAR, Zeu, Zsd products.
  * l2prod.h, support Kd_PAR, Zeu, Zsd products.
  * get_Kd.c, add Kd_PAR.
  * get_l2prod_index.c, support Kd_PAR, Zeu, Zsd products.
  * photic_depth.c, algorithms for Zeu and Zsd.
  * MSl12.mk, add photic_depth.o.
  * l12_proto.h, add get_photic_depth().

### Data Added
  * eval/seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200702, 2 temperature epochs + gain drift.

## 5.5.6 - 2006-11-06

Support new SeaWiFS cal table format for evaluation.  Add MUMM NIR correction.

### Source Changed
  * l12_parms.h, version change. Add AERMUMM aerosol option.
  * aerosol.c, accept AERMUMM option.
  * atmocor2.c, call MUMM rhown(NIR) on request.
  * input_struc.h, add MUMM controls (alpha, gamma, epsilon).
  * l12_proto.h, new MUMM functions get_rho_mumm() and get_rhown_mumm().
  * l2_hdf_generic.c, support new rhom MUMM reflectance product.
  * l2prod.h, add rhom product.
  * msl12_input.c, support MUMM controls.
  * mumm.c, new module with MUMM functions.
  * MSl12.mk, add mumm.c.

  * [libseawifs] l1a_calibrate.c, decouple fp and inst temps.

### Data Added
  * eval/seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200701, 2 temperature epochs.

## 5.5.5 - 2006-09-18

### Source Changed
  * l12_parms.h, version change.
  * l1_octs_hdf.c, new cal table format.

### Data Added
  * octs/cal/cal_octs_notiltref_timedepnir_2seg.hdf, new table with two segment
  NIR calibration (pre/post heating).
  * octs/cal/cal_octs_notiltref_timedepnir.hdf, modified table with two segment
  NIR calibration, but both segments identical.


## 5.5.4 - 2006-08-30

Fix eval=NEWOOB error for 678 MODIS. Fix bug in MODIS cloud mask eval option. Change
scale on nLw_645.

### Source Changed
  * l12_parms.h, version change.
  * nlw_outband.c, missing a1 coefficient for 678.
  * cloud_flags.c, error in string compare for file existance.
  * get_l2prod_index.c, nLw_645 scaling.


## 5.5.3 - 2006-08-21

Reorganize atmospheric correction steps: move tLf outside of gaseous
transmittance to simplify.  New Morel spectral Kd product, revised CZCS
calibration. Fix bug in aer_opt=-9.

### Source Changed
  * whitecaps.c, removed ozone transmittance adjustment from tLf.
  * atmocor1.c, no need to pass ozone optical thickness from gas trans to whitecaps().
  * atmocor2.c, tLf now subtracted after Lt is corrected for ozone.
  * gas_trans.c, no need to pass ozone optical thickness out.
  * calcite.c, tLf now subtracted after Lt is corrected for ozone.
  * soa.c, tLf now subtracted after Lt is corrected for ozone.
  * setflags.c, in absaer(), tLf added before transmittance through gases.
  * filter.c, various filter funcs, tLf now subtracted after Lt is corrected for ozone.
  * polcor.c, don't subtract tLf from L_x.
  * convl21.c, tLf added before transmittance through gases.
  * l2_prod.h, add definition for Kd_morel.
  * get_l2prod_index.c, add definition for Kd_morel.
  * l2_hdf_generic, add support for Kd_morel.
  * get_Kd.c, add Kd_morel function.
  * l1_czcs_hdf.c, switch back to single exponential form for temporal
  calibration.
  * aerosol.c, fix to epsilon comp to allow multiple calls with different
  long-wave.
  * fluorescence.c, clean-up.

### Data Added
  * common/Kd_morel.dat, coefficent table for Morel spectral Kd algorithm.
  * czcs/cal/cal_czcs_200608.hdf, revised czcs cal table (single exponential)
  * eval/seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200603, new table for eval.
  * common/palettes/default.pal, new 256-level default map/browse palette.

## 5.5.2 - 2006-07-28

Add evaluation support for new seawifs calibration table with revised
temperature corrections.  Update no2 transmittance. Improve scaling for browse
and map RGB images. Add eval switches to fail MODIS pixels associated with a
specific mirror side.

### Source Changed
  * l12_parms.h, version change, new eval switch for testing seawifs cal.
  new eval switches for MODIS mirror-side masking.
  * l1a_seawifs.c, support for new temperature corrections and associated
  cal table format.
  * gas_trans.c, revert to simple airmass function for no2.
  * rdsensorinfo.c, remove no2 airmass coefficients.
  * l1_modis_hdf.c, set all data for a given mside to navfail, if eval switch.
  * convl21.c, eval switch FIXPOLCAL was always on.
  * MSl1brsgen.c, land_min from 0.03 to 0.01.
  * MSl1tcpcbox.c, match brsgen scaling, disable cloud-specific scaling.

  * [libseawifs] get_cal.c, add get_cal_new().
  * [libseawifs] get_cal_misc.c, add read_parm_data_new().
  * [libseawifs] calibrate_l1a.c, add calibrate_l1a_new().

  * getcal_proto.h, add  get_cal_new(), read_parm_data_new().
  * call1a_proto.h, add calibrate_l1a_new().
  * get_cal.h, add TFACTOR_FLDS_NEW.

### Data Added
  * eval/seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200601, new cal table for testing.
  * *modis*/msl12_sensor_info.dat, remove no2 airmass coefficients.

## 5.5.1 - 2006-07-19

Finalizing NO2 correction option, including loading of no2 fields from
ancillary sources (stratospheric & tropospheric), options to output of those
fields, and application for computing no2 transmittance.  Adding new
quality-level 4 for SST to indicate cases of true failure (SSTFAIL flag now set
for QL=4). Added a mechanism to simplify switching of tables for evaluation
purposes.  Add evaluation rayleigh and aerosol tables developed by Zia Ahmad
for modis and seawifs. Revised CZCS temporal calibration (code & table).
Remove polfile from input structure; add eval switch to test alternate
polfiles. The ozone transmittance product is now the gaseous transmittance
product (t_oz_sol is now tg_sol, t_oz_sen is now tg_sen).

### Source Changed
  * l12_parms.h, version change, new eval switches for testing rayleigh, aerosol, and
  polarization tables.
  * l12_proto.h, add evalmask to calling args on rayleigh().
  * l2prod.h, add no2_strat & no2_tropo products.
  * l2_hdf_generic.c, add no2 products.
  * get_l2prod_index.c, add no2 products.
  * l1_struc.h, add no2 fields.
  * l2_struc.h, add no2 fields.
  * alloc_l1.c, add space for no2 fields in l1rec.
  * msl12_input.c, add no2file, add usage for chlorophyll algorithm coeffs, remove polfile.
  * input_struc.h, add no2file, remove polfile.
  * setanc.c, call no2conc() to read no2 concentration from ancillary.
  * rdsensorinfo.c, read new no2 airmass function coefficients.
  * gas_trans.c, apply no2 algorithm from Z. Ahmad.
  * sst.c, set QL to 4 for unprocessed cases (land, saturation), and bias & stdv to -999, set
  SSTFAIL for QL=4 (was 2) and still set SSTWARN for QL=1.
  * sstref.c, fix potential indexing error.
  * aerosol.c, add eval switch for alternate "eval" path to aerosol tables.
  * rayleigh.c, add eval switch for alternate "eval" path to rayleigh tables.
  * atmocor1.c, pass evalmask input to rayleigh().
  * water.c, change data type of firstCall.
  * convl12.c, transfer no2 field.
  * atmocor1_land.c, fix error in the band indexing for SeaWiFS water-vapor over land.
  * l1_czcs_hdf.c, new double exponential temporal calibration.
  * MSl1brsgen.c, fix indexing of start/end pixel & line (off by 1).
  * polcor.c, remove polfile input and just assume location based on sensor.
  * *.c, changed all references of t_oz_* to tg_*, including product names, e.g.: t_oz_sol_412
  is now tg_sol_412.

### Data Added
  * eval/seawifs/rayleigh/*, new rayleigh tables from Zia
  * eval/seawifs/aerosol/*, new aerosol tables from Zia, also moved the SOA models
  here so SeaDAS does not have to filter them out.
  * eval/modisa/rayleigh/*, new rayleigh tables from Zia
  * eval/modisa/aerosol/*, new aerosol tables from Zia
  * czcs/cal/cal_czcs.hdf, updated temporal calibration.


## 5.5	- 2006-06-02

Continued development of HIRES MODIS.  Added SST capability.  Added SWIR
atmospheric correction and/or NIR Lw correction, with new inputs to control
SWIR band selection..  Added new gaseous transmittance computation, user
selectable, including support for no2, co2, and water-vapor, in addition to
ozone.  Changes inspired rewrite of much of the front-end of the atmospheric
computions, (e.g, Rayleigh, ozone) from fortran to C.  Eliminated need for band
indexing, and removed it from external tables.

### Source Changed

  * l12_parms.h, version change, revised all aer_opt identifiers, add eval switch GASTRANS,  add eval switch AEROOBFIX. and eval switch NLWPOLCAL.
  * aerosol.c, fixes to allow for the fact that the longest aerosol table wavelength is not neccessarily the long NIR wavelength, and to allow that aerosols be computed for all sensor
  wavelengths, including those that exceed the longest wavelength used for aerosol model
  selection.  Also add support for CZCS model selection based on angstrom climatology.  Added
  checks to ensure proper wavelength specification relative to aerosol option specification.
  Moved aersol-specific informational messages from main() to here. Moved the out-of-band
  water-vapor correction inside the log within rhoas_to_rhoa(), as it should be, and added
  a check for rhoas=0 before taking the log (these changes were wrapped in evalmask=AEROOBFIX
  as they will alter standard production slightly).
  * atmocor1.c, new C function to replace fortan version.  call gaseous transmittance function
  rather than ozone-specific function.  uses new rayleigh function.
  * atmocor2.c, add support for SWIR-based NIR Lw correction, all new aer_opt identifiers
  * atmocor_soa.f, stick old gordon_o2.f inside.
  * brightness.c, adjust detector numbers based on resolution.
  * carder.c, better fix for missing 510nm band.
  * convl12.c, switch MSSOA selector to AERSOA.
  * convl21.c, EVAL switch to fix polcor correction of tLw.
  * get_chl.c, fail chl if MBR exceeds 30 (czcs oc3 was blowing-up), put upper limit of 1000.
  * get_qaa.c, recognize HMODIS for appropriate 640 estimator.
  * get_tsm.c , call appropriate algorithm for HMODIS.
  * input_struc.h, add aer_swir wavelength fields.
  * l12_proto.h, passing resolution to radiancebt(). removed defunct fortran protos.
  * l1_hmodis_hdf.c, add loading and interpolation of thermal bands to support SST.
  * l1_modis_hdf.c, passing resolution to radiancebt().
  * loadl1.c, call new C function for atmocor1().
  * MSl12.c, moved aersol-specific messages to aerosol.c.
  * msl12_input.c, load aer_swir_short and aer_swir_long, simplify usage of aer_opt.
  * polcor.c, removed backward compatibiility for 60.0-deg AOI.
  * rayleigh.c, rewrite of getrayleigh.f, no change intended.
  * setflags.c, use SWIR band at 2130 for cloud screening, if it exists. Replace aerindex routine.
  * sst.c, recognize HMODIS variants, remove pixel-specific edge masking.
  * windex.c, fix for band differences longer than 1000nm.
  * whitecaps.c, switch to zero-based band indexing, as no called from C.
  * bin_climatology.c, function to read and interpolate 9km bin files.
  * gas_trans.c, new gaseous transmittance function.

  * MScalmerge.mk, add gas_trans.o, rayleigh.o, remove getozone.o, gordon_o2.o, getrayleigh.o
  * MSl12.mk, add gas_trans.o, rayleigh.o, remove getozone.o, gordon_o2.o, getrayleigh.o
  * MSl1brsgen.mk, add gas_trans.o, rayleigh.o, remove getozone.o, gordon_o2.o, getrayleigh.o
  * MSl1tcpcbox.mk, add gas_trans.o, rayleigh.o, remove getozone.o, gordon_o2.o, getrayleigh.o

  * atmocor1.f, defunct
  * atmocor_init.f, defunct
  * getozone.f, defunct
  * gordon_o2.f, defunct
  * getrayleigh.f, defunct
  * sensor_cmn.fin, defunct

### Data Added
  * */msl12_defaults.par, add default gas_opt, add aer_swir_short and aer_swir_long for HMODIS.
  * */msl12_sensor_info.dat, add entries for co2, no2, and h2o sensitivities.

## 5.4.2 - 2006-04-19

Continued development of HIRES MODIS; support for L1A extract processing.

### Source Changed
  * l12_parms.h,  version change.
  * l1_hmodis_hdf.c, recognize when L1B was generated from L1A extract, set proper start
  pixel and line, reduce fields to extract region. Also "fix" singularity in senzor
  azimuth interpolation. Fix interpolation at scan edge and along-track boundaries. Allow
  LAC and 1KM variants in filename convention. Don't check for HILT on swir bands (1640
  was setting dets 32-40 always.
  * l1_modis_hdf.c, removed support for old MODAPS "cookie-cutter" extracts.
  * MScalmerge.c, remove previous "HCALFIX" hacks and set resolution.
  * MScaltarget.c, set resolution.
  * MSl1bgen.c, set resolution.
  * MSl1tcpcbox.c, set band numbers of HIRES RGB image.


## 5.4.1 - 2006-03-14

Corrections for 16-band HIRES MODIS and fixes for proper band indexing for some
algorithms.  Change emperical chlorophyll algorithms to allow user-specified
coefficients.  Eliminate octsc and ndpi algorithms. Fix corrupted files with
DOS carriage returns. Update OCTS calibration algorithm and file.

### Source Changed
  * l12_parms.h,  version change, add eval mask for new out-of-band coefficients (NEWOOB).
  * l1_hmodis_hdf.c, fixed intepolation of solar and sensor path geometry (they were swapped). Change data buffering to
  float for better interpolation.  Support subframe destriping.  Support extract processing.
  * carder.c, added specific band selection (was assuming band positions).
  * water.c, handle SWIR bands by using longest available wavelength from table.
  * gsm01.c, fixed mis-named adgstar in CB model argument (no effect).
  * l1_octs_hdf.c, updated calibration table and flagging.
  * get_chl.c, major changes: coefficients no longer specified in code, now specified in
  default parameter files.  Old octsc and ndpi eliminated. Make OC3 the default for
  OSMI.
  * brdf.c, pass l2rec in calling sequence, so it can be passed to get_chl_default.
  * soa.c, change get_chl_default() calling args.
  * atmocor2.c, change get_chl_default() calling args, pass evalmask to nlw_outband().
  * l12_proto.h, change get_chl_default() calling args, new nlw_outband() args.
  * msl12_input.c, support reading of chl algorithm coefficients.
  * l2prod.h, drop octsc and ndpi.
  * l2_hdf_generic.c, drop octsc and ndpi.
  * get_l2prod_index.c, drop octsc and ndpi, use OC3 as default for OSMI.
  * qaa.c, remove CRs.
  * qaa.h, remove CRs.
  * get_qaa.c, look for 551 rather than 555 (so 16 and 9-band versions give same result),
  remove CRs.
  * get_rhown_nir.c, remove CRs.
  * MSl1brsgen.c, change RGB selection for 16-band MODIS.
  * MScalmerge.c, remove quickfix for handling HIRES modis.
  * getformat.c, recognize CZCS L1B generic.
  * nlw_outband.c, install new coefficients for evaluation.
  * l1_io.c, set band indexing (previously not done until loadl1, but needed in l1_hmodis_hdf).
  * aerosol.c, modified MS6 & MS6NIR options (czcs) to select model based on Angstrom, when
  evalmask contains 16.
  * bin_climatology,c, new module to read rkm bin files as climatologies.
  * convl21, pass l2rec to ocbrdf.

### Data Added
  * octs/cal/cal_octs_notiltref_timedep865.hdf, new OCTS cal table.
  * octs/msl12_defaults.par, new vicarious calibration.
  * */msl12_defaults.par, add chl algorithm coefficients specific to each sensor.
  * common/S*MC.hdf, updated climatologies, added angstrom_510 to 412, chlor_a.

## 5.4 - 2006-02-14

Significant changes to support 16-band HIRES MODIS concept.  Generalization of
various algorithm components, including new input parameters to select bands
for aerosol selection process.

### Source Changed
  * l12_parms.h, version change, add eval mask for extended Thuillier.
  * l12_proto.h, changes to argument list of get_rhown_nir() and get_default_chl().
  * input_struc.h, add aer_wave_short and aer_wave_long parameters to specify wavelengths
  to use for aerosol model selection. Add resolution parameter for switch control.
  * msl12_input.c, as above.
  * atmocor2.c, many changes to generalize handling of NIR corrections and specification
  of wavelengths used for model selection.
  * windex.c, adjust direct access band and wavelength indexing tables to allow for
  overlapping band passes (e.g., 551, 555).
  * l1_struc.h, add aw & bbw fields for band-averaged, sensor-specific coefficients.
  * l2_struc.h, add aw & bbw fields.
  * l2prod.h, added chl_oc3cb.
  * cpl1rec.c, replicate aw & bbw fields.
  * loadl1.c, load aw & bbw fields, allow use of extended Thuillier spectrum by
  non-HMODIS via eval=4.
  * get_rhown_nir.c, generalized to pass sensor wavelengths and band indices at which
  NIR water-leaving reflectance is desired. Use sensor-specific aw & bbw from l2rec.
  * l12_proto.h, new get_rhown_nir() prototype.
  * convl12.c, fix for scan-overlap when assessing orbit node. Transfer aw, bbw fields.
  * polcor.c, fix for computing detector number and AOI for varying resolutions.
  * aerosol.c, fixed epsilon interpolation for case of more bands than table contains.
  * nLw_outband.c, look for 551 rather than 555, since HMODIS has both.
  * get_chl.c, pass Rrs rather than nLw and Fo to get_default_chl.c. Added chl_oc3cb product.
  * brdf.c, change get_default_chl.c call.
  * soa.c, change get_default_chl.c call.
  * alloc_2d.c, added version for short int.
  * alloc_2d.h, added version for short int.
  * l1_hmodis_hdf.c, nearly complete rewrite to incorporate full 16-band suite and allow for
  switching resolutions. Geometry is now fully interpolated to the specified resolution.
  Lower-resolution radiances are also interpolated.  Input file can be any one of the three
  MODIS L1B files (QKM, HKM, 1KM [LAC]). If 1KM or LAC, specifying resolution=1000 will switch
  from 9-band OC processing to 16-band processing. Access to resolutions of 250 and 500
  require that the filenames follow a standard convention in that QKM, HKM, 1KM [or LAC]
  must appear once and only once in the filenames.
  * l1_modis_hdf.c, changed indexing variables to something meaningful.
  * filehandle.h, add resolution field.
  * filehandle_init.c, initialize resolution field.
  * MSl12.c, copy requested resolution to file handle.
  * getformat.c, check resolution when determining format.
  * getl1rec.c, when loading the queue, limit the record number to be no less than the start
  line, as specified by the user.
  * get_l2prod_index.c, added chl_oc3cb.
  * l2_hdf_generic.c, added chl_oc3cb.

### Data Added
  * hmodisa/cal/polcor*, additional polarization tables.
  * hmodist/cal/polcor*, additional polarization tables.
  * hmodisa/rayleigh/*, additional rayleigh tables.
  * hmodisa/rayleigh/*, additional rayleigh tables.
  * hmodisa/msl12_sensor_info.dat, expand to 16 bands.
  * hmodist/msl12_sensor_info.dat, expand to 16 bands.
  * */msl12_defaults.par, modified defaults for all sensors: add aerosol
selection wavelengths,  expand hmodis gains, offsets to 16 bands.


## 5.3.7 - 2006-02-03

Changes for absorbing aerosol flag evaluation. Updates to SST QL flagging. CZCS nav update.
Add ability to switch GSM01 model to Chesapeake model coefficients.

### Source Changed
  * l12_parms.h, version change.
  * input_struc.h, absaer_opt and gsm01_opt parameters, also additional changes for 5.4.
  * msl12_input.c, absaer_opt parameter and gsm01_opt parameters.
  * setflags.c, add functions to support use of climatologies for expected chl and nlw_412.
  * setupflags.c, include with MSl12, rename var null to nul.
  * sst.c, quality-level changes: allow sst to be 1-deg warmer than Reynolds for daytime
  before reaching QL=1.
  * l1_czcs_hdf.c, changes to utilize revised navigation.
  * l1_czcs_hdf.h, changes to utilize revised navigation.
  * MSl12.c, added informational message re along-track detectors.
  * gsm01.c, new gsm01_cb_model() coefficients, switched via gsm01_opt=1.

### Data Added
  * common/S*_L3b_MC.hdf, binned monthly climatologies of nLw_412 and chlor_a from SeaWiFS.


## 5.3.6 - 2006-01-11

Add SSES error statistics for SST (source tables from RSMAS) to support GHRSST. New
products bias_sst, bias_sst4, stdv_sst, stdv_sst4.  Support new OCTS calibration table
with temporal correction.

### Source Changed
  * l12_parms.h, version change,
  * sst.c, add functions to read RSMAS hypercube data from hdf tables and assign bias
  and standard deviation to each SST retrieval. Mask scan-edge pixels (<8 or >1345),
  and set to QL=3. Add test for 11-12um SST, in which the QL is set to 3 if the 4um
  BT difference is large relative to the reference 4um BT difference.
  * get_l2prod_index.c, add products bias_sst, bias_sst4, stdv_sst, stdv_sst4.
  * l2prod.h, add products bias_sst, bias_sst4, stdv_sst, stdv_sst4.
  * l2_hdf_generic.c, add products bias_sst, bias_sst4, stdv_sst, stdv_sst4.
  * l12_proto.h, prototypes for bias and stdv functions.
  * l1_octs_hdf.c, extensive changes to support new externalized calibration
  table, including temporal correction capabilities.

  * liboctscal, defunct (functionality incorporated into l1_octs_hdf.c).

### Data Added
  * octs/msl12_defaults.par, new calibration.
  * octs/cal/octs_cal_nasda_v4.hdf, original NASDA calibration for OCTS, formatted to HDF.
  * octs/cal/octs_cal_notiltref_timedepnir.hdf, remove tilt reflectance correction,
  add NIR time dependence.
  * modisa/cal/sst_sses_modisa.hdf, error table for MODISA sst
  * modisa/cal/sst4_sses_modisa.hdf, error table for MODISA sst4
  * modist/cal/sst_sses_modist.hdf, error table for MODIST sst
  * modist/cal/sst4_sses_modist.hdf, error table for MODIST sst4


## 5.3.5 - 2006-01-03

Simplification of chlorophyll algorithms to use precomputed Rrs.  Fix some
unititialized vars in dem height. Change SST coefficient look-up to use the
previous year when no new coefficients are provided.

### Source Changed
  * l12_parms.h, version change,
  * aerosol.c, fix for bad angstrom computation when fixed AOT option results in models
  out of range.
  * get_chl.c, pass Rrs where appropriate.
  * get_dem_height.f, add some vars to save command.
  * sst.c, change SST coefficient look-up to use the previous year when no new
  coefficients are provided.

### Data Added
  * modist/cal/sst*, extend coefficient files into 2006.

## 5.3.4 - 2005-12-15

Update for SST.  Fixes for HIRES MODIS calibration.  Other minor fixes.

### Source Changed
  * l12_parms.h, version change,
  * sst.c, reduce Reynolds difference from 6 to 5 for QL=3, SSTFAIL.
  * msl2val_hdf.c, at prototyping include.
  * getformat.c, add recognition of generic L1B formats from HMODISA and HMODIST.
  * l1_hmodis_hdf.c, fix frame number computation.
  * l1_io.c, initialize alpha and ndets (not set by generic i/o functions).
  * MSl12.c, for inverse processing, allow for target file to be pre-subsetted by scan.
  * polcor.c, prevent divide-by-zero.
  * MScalmerge.c, added fix to allow scan subsetted and non-subsetted L1B files to be
  merged, but wrapped fix in #if HCALFIX and disabled.


## 5.3.3 - 2005-12-06

Bug fix for SST (J.Gales). Went directly into production.

### Source Changed
  * sst.c, fix case of data day less than 100 when reading coefficient files.

## 5.3.3 - 2005-11-21

Changes to more fully support the HIRES MODIS processing. Major rewrite of
polarization code to use a revised, band-specific table format, and to handle
different numbers of detectors in each band/table.  New tables were generated
for standard MODIS as well.  The band suite was also expanded to include the
SWIR bands.  Also reverted back to using the 4um SST as night reference for
MODIS 11-12um SST and for quality tests.  Added new absorbing aerosol flag.
Added ability to read and apply MODIS cloud mask (PGE35), for evaluation only.

### Source Changed
  * l12_parms.h, new version number. Increase NBANDS to 12.  Add eval identifiers for alternate
  cloud flagging (eval=2) and forcing polarization AOI to 60 from 60.5 (eval=1).
  * l12_proto.h, new polcor() proto. Add compute_alpha() proto. glint_rad(), modis_cloud_flag()..
  * l1_modis_hdf.c, store number of detectors per band (ndets)
  * l1_hmodis_hdf.c, store number of detectors per band (ndets). Read 2130nm band.
  * l1_struc.h, add number of detectors per band (ndets) and number of lines (nscans).
  * l2_struc.h, add number of detectors per band (ndets).
  * l1_io.c, transfer number of detectors per band (ndets) and nscans.
  * filehandle.h, add number of detectors per band (ndets).
  * filehandle_init.c, initialize ndets.
  * cpl1rec.c, transfer number of detectors per band (ndets) and nscans.
  * MSl12.c, transfer number of detectors per band (ndets) to l2rec.
  * get_f0.c, add function to read extended solar spectrum.
  * polcor.c, replaces code in atmocor1.f and old polcor.f to read new polarization table format
  and compute polarization corrections.
  * atmocor2.c, limit atmospheric corrections to longest NIR wavelength. Pass NIR indicies to
  glint_rad().
  * glint.c, pass NIR indicies.
  * aerosol.c, remove wave and nwave parameters from aeroob. Remove assumption in get_angstrom()
  that the last band is the NIR band.
  * loadl1.c, standardized indexing variable names.  Removed polcor & dpol from atmocor1 call
  and added call to new polcor() function. Use extended Thuillier spectra for hires MODIS
  (change all sensors someday). Remove option to call Neckel&Labs.
  * atmocor1.f, remove polarization computations (polcor,dpol)
  * atmocor1_land.c, clean-up.
  * get_l2prod_index.c, remove NBANDS dependency from nLw scaling.
  * sst.c, changes to reinstate relationships between SST and SST4, and 11-12um and 4um BTs.
  * setflags.c, added aerindex() and revised aborbing aerosol test.
  * msl12_input.c, change default threshold on absorbing aerosol test (absaer).
  * l2_hdf_generic.c, change processing time meta to GMT.
  * cloud_flag.c, new function to read and apply standard MODIS cloud flag.
  * MSl2validate.c, clean-up.
  * getformat.c, clean-up.
  * convl12.c, clean-up.
  * polcor.f, defunct.

### Data Added
  * modisa/cal, add new polcor_modisa_*.hdf files and remove modisa_polsen_2004085.hdf.
  * hmodisa/cal, add new polcor_hmodisa_*.hdf files and remove modisa_polsen_2004085.hdf.
  * modist/cal, add new polcor_modist_*.hdf files and remove modist_polsen_2004175.hdf.
  * hmodist/cal, add new polcor_hmodist_*.hdf files and remove modist_polsen_2004175.hdf.
  * hmodisa/msl12_sensor_info.dat, include SWIR bands.
  * hmodist/msl12_sensor_info.dat, include SWIR bands.
  * modisa/msl12_defaults.par, new default polarization file, new absaer threshold.
  * modist/msl12_defaults.par, new default polarization file, new absaer threshold.
  * hmodisa/msl12_defaults.par, new default polarization file, new absaer threshold.
  * hmodist/msl12_defaults.par, new default polarization file, new absaer threshold.
  * czcs/l12_defaults.par, new absaer threshold.
  * seawifs/l12_defaults.par, new absaer threshold.
  * octs/l12_defaults.par, new absaer threshold.
  * mos/l12_defaults.par, new absaer threshold.
  * polder/l12_defaults.par, new absaer threshold.
  * czcs/l2bin_defaults.par, remove ABSAER masking.
  * seawifs/l2bin_defaults.par, remove ABSAER masking.
  * common/Thuillier_F0.dat, extended solar irradiance model.

## 5.3.2 - 2005-10-27

Changes to SST quality tests to remove relationships between short-wave and
long-wave sst.  Use Reynolds for day and night reference (not sst4 for night
reference).  Add RSMAS corrections for Terra radiometry (fixes striping in
short-wave SST) and add filter for Terra to perform detector averaging. Also
add capability to extend RGBN using aggregated MODIS land bands.  For CZCS, add
support for scene extact processing. Added ability to use aggregated bands to
extend ocean bands, for standard L1B files.  Added new combined calcite
algorithm (2-band with switching to 3-band).  Also added new set of stub
functions for user defined products.

### Source Changed
  * l12_parms.h,  version number change.
  * l1_modis_hdf.c, added (hardcoded) RSMAS corrections for Terra radiometry. Only destriping
  of thermal bands 22 and 23 is enabled. Added ability to extend RGBN bands using aggregated
  hires bands.
  * sst.c, changes to remove relationships between SST and SST4, and 11-12um and 4um BTs.
  * filter.c, add filter function to allow detector averaging/replacement.
  * filter.h, add filter function to allow detector averaging/replacement.
  * l1_czcs_hdf.c, mods to handle czcs scene extract files.
  * calcite.c, added combined calcite algorithm. Changed conversion from backscatter
  to calcite concentration.
  * myprod.c, new function to provide stubs for user-defined algorithms.
  * l2_hdf_generic.c, add support for calcite and myprod stubs.
  * l2prod.h, add support for calcite and myprod stubs.
  * get_l2prod_index.c, add support for calcite and myprod stubs.
  * l12_proto.h, add support for calcite and myprod stubs.
  * msl12_input.c, new coccolith coefficients (entered but not enabled).
  * MSl12.mk, added myprod.c.
  * MSl1info.c, fix for daynight flagging.
  * MSl1tcpcbox.c, remove relative referencing of includes.
  * convl12.c, fixed scene meta-data node calculation.

### Data Added
  * modist/msl12_filter.dat, add filter specs to average detectors 4 and 8 of band 22 (NIR band 3).


## 5.3.1 - 2005-09-13

Switch to RSMAS "collection 5" SST quality tests.  Fix gsm01 to include above
to below-water correction.

### Source Changed
  * l12_parms.h,  version number change.
  * sst.c, collection 5 quality test changes.
  * gsm01.c, add translation from surface to subsurface.

## 5.3.0 - 2005-09-13

Major update to SST algorithms in preparation for transition of official MODIS SST
processing responsibility to OBPG. Also support for MODIS high resolution bands.

### Source Changed
  * l12_parms.h, version number change. Added DAYSCENE, NIGHTSCENE, MIXEDSCENE defines for granule
  daynight identification. Increase MAXPIX to 5500 for hires MODIS.
  * sst.c, reads new external, time-dependent algorithm coefficient files. Computes quality flags
  and quality levels.  Moved brightness temp and sstref function to separate source files.
  * brightness.c, new module to compute brightness temperatures.  Handles 4um table inversion
  in log space.
  * sstref.c, new module of functions to read reference SST product (oisst or pathfinder).
  * l2_struc.h, added daynight meta field.
  * l2_hdf_generic.c, write daynight meta.  Support for sst flag and quality level products and
  4um BT products. Change buffer space to dynamic allocation (eliminate MAXPIX dependency).
  Compute daynight meta. Don't include SOA calls if SeaDAS (for now).
  * l2prod.h, support for sst flag and quality level products and 4um BT products.
  * get_l2prod_index.c, support for sst flag and quality level products and 4um BT products.
  * convl21.c, eliminate goto on negative calibration target, no effect to processing.
  * convl12.c, process sst by line rather than by pixel.  Compute daynight meta. Don't include
  SOA calls if SeaDAS (for now).
  * get_chl.c, added OC2 variant for hires MODIS bands.
  * l1_modis_hdf.c, always compute brightness temps (no proc_sst test). Set radiance to zero for
  all cases of flagged L1B, including saturation (to simplify detection of saturated bands in
  downstream processing). Change get_hdfeos_meta() to return rather than exit on error, so it
  might work for IMAPP L1B.
  * getl1rec.c, allow l1que to be visible to other source files. Set minimum queue size to NQMIN=3,
  to ensure that sufficient line buffering is available for the SST homogeneity tests.
  * l12_proto.h, add prototype for retrieving SST products and associated flags.
  * l1_struc.h, eliminate level-1 sstwarn, sstfail flags, sst field.
  * alloc_l1.c, eliminate level-1 sstwarn, sstfail flags, sst field.
  * setflags.c, eliminate level-1 sstwarn, sstfail flags.
  * carder.c, check that sst pointer was assigned.
  * l1_io.c, support for HMODISL1B files, don't init sst field of l1rec.
  * soa_sma_utils.f, error checking.
  * getformat.c, support for hires MODIS files.
  * b128_msk_get.c, support for larger scene size of hires MODIS.
  * filehandle.h, support for HMODISA & HMODIST sensors & HMODISL1B file type.
  * getrayleigh.f, support for HMODISA & HMODIST sensors.
  * l1_hmodis_hdf.c, new module to read hires MODIS files.
  * l1_hmodis_hdf.h, new module to read hires MODIS files.
  * MSl1brsgen.c, allow rgb image from hires MODIS bands.
  * setanc.c, eliminate MAXPIX dependencies.
  * atmocor2.c, fixed bug in fixed-aerosol model case (bug introduced in 5.2.3).
  * qaa.c, fixed indexing errors, use fixed aw and bbw for efficiency.
  * get_qaa.c, fixed indexing errors.
  * qaa.h, init function passes aw, bbw.

### Data Added
  * hmodisa, new directory for HiRes MODIS/Aqua
  * hmodist, new directory for HiRes MODIS/Terra
  * modisa/cal/sst_modisa.dat, new sst algorithm coefficients
  * modisa/cal/sst4_modisa.dat, new sst4 algorithm coefficients
  * modisa/cal/bt_modisa.hdf, updated to log for 4um temperatures
  * modist/cal/sst_modist.dat, new sst algorithm coefficients
  * modist/cal/sst4_modist.dat, new sst4 algorithm coefficients
  * modist/cal/sst_mside_modist.dat, new sst mirror-side corrections
  * modist/cal/bt_modist.hdf, updated to log for 4um temperatures
  * octs/msl12_filter.dat, add straylight filter
  * */msl12_defaults.par, change from l2prod1 to l2prod, for consistency.



## 5.2.3 - 2005-08-10

Generalization changes from Paul Martinolich. Add new ATMFAIL condition for
case of negative rhoa(nir). Fix potential for negative AOT (set rhoa to at
least 1e-6 if near zero). Trap bad Angstrom. Remove view-angles-per-band
support. Fix bug for ctl_pt_incr=0 (added in 5.2.0)..

### Source Changed
  * l12_parms.h, version number change.
  * MSl1brsgen.c, updates for CZCS & OCTS.
  * MScaltarget.c, pass nbands to allocation functions.
  * MSl12.c, pass nbands to allocation functions.
  * alloc_aer.c, add nbands argument to allocation functions.
  * alloc_l2.c, add nbands argument to allocation functions.
  * alloc_target.c, add nbands argument to allocation functions.
  * l12_proto.h, add nbands argument to allocation functions. Add proto
  for get_rhown_nir().
  * carder.c, pass nbands for allocation of static arrays.
  * get_chl.c, pass nbands of band selection. Add OC3 variant for seawifs.
  * atmcor_soa.f, limit to 72 cols, fix discrepancy in common.
  * l2_flags.h, standardize some flag names (TURBIDW, HITAU).
  * setflags.c, standardize some flag names.
  * l1_czcs_hdf.c, add ringing mask.
  * aerosol.c, fail cases of very negative input aerosol reflectances.  Trap bad
  angstrom products as CHLFAIL.
  * atmocor2.c, add ability to disable ramping of NIR correction, and do so for CZCS.
  Call new C version of NIR correction.
  * get_rhown_nir.c, new C version of NIR correction with update for CZCS.
  * l1_modis_hdf.c, remove tsolz, tsenz, tsola, tsena, remove radcor support..
  * l1_xcal_hdf.c, remove tsolz, tsenz, tsola, tsena.
  * l1subpix.c, remove tsolz, tsenz, tsola, tsena.
  * convl12.c, , remove tsolz, tsenz, tsola, tsena.
  * l1_polder.c, remove tsolz, tsenz, tsola, tsena.
  * l1_octs_hdf.c, remove tsolz, tsenz, tsola, tsena.
  * l1_mos_hdf.c, remove tsolz, tsenz, tsola, tsena.
  * l1_hdf_generic_write.c, remove tsolz, tsenz, tsola, tsena, tdelphi, EXTENDED format.
  * l1_hdf_generic_read.c, remove tsolz, tsenz, tsola, tsena, delphi, EXTENDED format.
  * l1_generic.c, remove tsolz, tsenz, tsola, tsena.
  * l1a_seawifs.c, remove tsolz, tsenz, tsola, tsena.
  * l1a_osmi.c, remove tsolz, tsenz, tsola, tsena.
  * alloc_l1.c, remove tsolz, tsenz, tsola, tsena, tdelphi.
  * l1_struc.h, remove tsolz, tsenz, tsola, tsena, tdelphi.
  * l2_struc.h, remove tsolz, tsenz, tsola, tsena, tdelphi.
  * loadl1.c, remove tdelphi computation.
  * l1_io.c, remove EL1B & EL1HDF format support.
  * filehandle.h, remove EL1B & EL1HDF format identifiers.
  * MScalmerge.c, remove reference to EL1HDF.

## defunct src:
  * get_rhown_nir.f, replaced with get_rhown_nir.c

### Data Added
  * czcs/aerosol/*, updated diffuse transmittance to use actual CZCS RSRs.


## 5.2.2

Clean-up and generalization changes from Paul Martinolich.

### Source Changed
  * l12_parms.h,
  * aerosol.c, clean-up unused vars.
  * alloc_l1.c, add nbands, nbandsir to arg list.
  * atmocor1_land.c, clean-up unused vars.
  * atmocor2.c, code indentation.
  * brdf.c, code indentation.
  * calcite.c, code indentation.
  * carder.c, rename init_carder() to alloc_carder().
  * getl1rec.c, pass nbands, nbandsir to alloc_l1().
  * get_l2prod_index.c, ???
  * get_qaa.c, pass nbands to alloc_qaa and change doubles to floats.
  * l12_proto.h, new alloc_l1() proto.
  * l1_polder.c, pass nbands, nbandsir to alloc_l1().
  * MScalmerge, pass nbands, nbandsir to alloc_l1().
  * MScaltarget, pass nbands, nbandsir to alloc_l1().
  * MSl12.c, pass nbands, nbandsir to alloc_l1().
  * MSl1bgen.c, pass nbands, nbandsir to alloc_l1().
  * MSl1brsgen.c, pass nbands, nbandsir to alloc_l1().
  * MSl1info.c, pass nbands, nbandsir to alloc_l1().
  * MSl1tcpcbox.c, pass nbands, nbandsir to alloc_l1().
  * MSscanpixlimits.c, pass nbands, nbandsir to alloc_l1().


## 5.2.1 - 2005-06-07 

Implement spectral optimization.

### Source Changed
  * l12_parms.h, version change, add MSSOA & MSSMA identifiers.
  * convl12.c, call run_soa if aer_opt is MSSOA.
  * MSl12.c, add informational messages for MSSOA and MSSMA.
  * l12_proto.h, added prototypes run_soa, atmcor_soa_ and related utilities.
  * l2prod.h, added chl_soa, bbp_soa, adg_soa, w0_soa, v_soa. Moved struct definition to
  new l2prod_struc.h and proto to l12_proto.h, so O.K. to include in fortran code.
  * l2_hdf_generic.c, add soa products.
  * get_l2prod_index.c, add soa products.
  * atmcor_soa.f, new module containing SOA models and optimization code.
  * soa_sma_utils.f, new utilities to store asnd retrieve optimization results.
  * soa.c, new wrapper to feed SOA routine.
  * l2prod_struc.h, new file, struc defintion moved from l2prod.h

### Data Added
  * seawifs/aerosol/seawifs_hzc*.dat, new haze-C aerosol tables for SOA


## 5.2.0 - 2005-05-20 

Finish transition of SeaWiFS to Thuillier (forgot Fonom). SeaWiFS chl will
change slightly.  OBPG Kd now default K_490.  Add new aerosol options for CZCS.
Initial implementation of 4um sst, and some checks to avoid processing of
visible bands for night data.  Fix indexing error in NDVI and EVI functions.
Add a new direct-access wavelength indexing scheme to simplify algorithm to
band mapping.  NIR rhown fix to use proper sensor variants..

### Source Changed
  * l12_parms.h, version change, NBANDSIR increased to 4 for sst4 support, added MS6
  and MS6NIR aerosol options for CZCS, SOLZNIGHT solarzenith threshold for recognizing
  nighttime data, eval switch for reverting to Neckel&Labs Fonom.
  * loadl1.c, default to Thuillier nominal-band F0 (eval=1 for N&L). Initialize new
  direct access wavelength index. Check nav before calling get_dem_height().
  * get_f0.c, fix initialization problem.
  * get_Kd.c, generalized and updated the OBPG Kd algorithm.
  * get_chl.c, generalized band selection, added CZCS variant of OC3.
  * l12_parms.h, new MS6 and MS6NIR and SOLZNIGHT defined.
  * MSl12.c, add informational messages for MS6 and MS6NIR. Dropped OUTPUT_NLW option. Force
  pixel control point increment to one for small scenes.
  * atmocor2.c, add switches for MS6 and MS6NIR, change incrementation loops for NIR correction.
  to handle case of nir_l = nir_s.
  * aerosol.c, switched model_taua() to use longest sensor wavelength rather than longest
  table wavelength (should only change CZCS).  Add czcsaer aerosol scheme (still developing).
  * nlw_outband.c, use 443/555 ratio for CZCS.
  * windex.c, added new windex_set, windex_get, bindex_get functions to simplify sensor to
  algorithm wavelength selection.
  * l12_proto.h, added prototypes for windex and sst4.
  * l1_modis_hdf.c, read 4um bands 22 and 23 for sst4.  Skip night data for visible bands.
  * sst.c, add 4um sst functions.
  * l2prod.h, add sst4.
  * l2_hdf_generic.c, add sst4.
  * get_l2prod_index.c, add sst4.  Add defaullt chlor_a for CZCS. Switch K_490 to use OBPG Kd.
  * setflags.c, avoid cloud, glint tests for night data. change upper limit on NEGLW check
  to 600nm (was 570) for MOS.
  * convl12.c, don't do OC processing for night data.
  * get_ndvi.c, fix band indexing error.
  * filehandle.h, clean-up.
  * get_rhown_nir.f, add CZCS variant (not working yet).  Also, sensorID params were not
  defined, so SEAWIFS was always selected. Fix will change MODISA slightly.
  * MSl1brsgen.c, change scaling back to old format for non-MODIS sensors.
  * MSl1tcpcbox.c, new RGB mapper.
  * MSl1tcpcbox.mk, makefile for above.

### Data Added
  * czcs/msl12_sensor_info.dat, updated values and dropped 750 band.
  * czcs/msl12_defaults.par, updated values and dropped 750 band.
  * modisa/cal/bt_modisa.hdf, expanded to include 3.9 and 4um channels.
  * modisa/cal/bt_modist.hdf, expanded to include 3.9 and 4um channels.

## 5.1.3 - 2005-04-06

Add alternate brdf_opt (Morel Q) for Ken Voss test.  Also, added iter_gsm01
product (number of iterations fitting GSM01), and fixed error in flags_carder
product. For Kd_lee, if no iop model was selected, default to use QAA.

### Source Changed
  * l12_parms.h, version number change.  New QMOREL brdf_opt identifier. Define
  default IOP model, IOPDEFAULT, as IOPQAA.
  * brdf.c, added switch in foq_morel() to use observed solz in numerator and
  denominator, effectively cancelling f and yielding a Q correction.
  * carder.c, changed form of get_flags_carder() to fix segfault error.
  * gsm01.c, added get_iter_gsm01() function. Also changed to return results
  even for max iteration case and chl < 0 case (CHLFAIL set).
  * l12_proto.h, changed form of get_flags_carder(), added get_iter_gsm01().
  * get_l2prod_index.c, added iter_gsm01 product.
  * l2prod.h, added iter_gsm01 product.
  * msl12_input.c, expanded usage for brdf_opt.
  * l1_modis_hdf.c, improved error message.
  * l2_hdf_generic.c, added iter_gsm01 product, and changed call to
  get_flags_carder().
  * get_Kd.c, use default IOP model when IOPNONE is selected and Kd_lee is
  requested as output product.


## 5.1.2 - 2005-04-01

Add QAA iop algorithm, Kd_lee spectral diffuse attenuation, and base iop
switching option. Also add basic L2 processing support for czcs. Fixed some
divide-by-zero errors. Fix error in reporting NIR nLw model values.

### Source Changed
  * l12_parms.h, version number change.  Add IOP selection identifiers, QAA opts.
  * alloc_l2.c, add space for a and bb IOP fields.
  * carder.c, add iops_carder() function.
  * convl12.c, call get_iops() to load a and bb fields.
  * get_l2prod_index.c, add QAA model products, and a and bb products.
  * gsm01.c, separate run_gsm01() from get_gsm01() functions to support new
  iops_gsm01() function.  Fixed model to l2rec band indexing for case of
  sensors with less than 9 bands.
  * input_struc.h, add iop_opt, qaa_opt, and qaa_s param fields.
  * l12_proto.h, new funcs invbindx(), get_iops(), get_qaa(), iops_qaa(),
  ipos_gsm01(), iops_carder().
  * l2_hdf_generic.c, new QAA prods and a and bb prods.
  * l2prod.h, new QAA prods and a and bb prods, Kd_lee.
  * l2_struc.h, add a and bb IOP fields.
  * msl12_input.c, handle new iop_opt, qaa_opt, and qaa_s params.
  * windex.c, add inverse band indexing function.
  * get_qaa.c, MSl12 interface to QAA algorithm.
  * qaa.c, Quasi-Analytical Algorithm (QAA).
  * get_Kd.c, new consolidated source for Kd algorithms, including spectral Kd
  algorithm of Z.P.Lee.
  * l1_czcs_hdf.c, updates from WR.
  * nr_spline.c, comments.
  * getrayleigh.f, CZCS selection.
  * aerosol.c, change rhoa_to_rhoas to double precision, to fix DBZ error.
  * getrayleigh.f, change pressure correction to double precision, to fix DBZ error.
  * atmocor2.c, fix indexing error when saving NIR nLw model values.
  * MSl12.mk, add qaa.c, get_qaa.c, get_Kd_lee.c.

## defunct src
  * get_k490.c, absorbed into new get_Kd.c

### Data Added
  * czcs/rayleigh, new rayleigh tables (approximated)


## 5.1.1 - 2005-03-02 

Fix for Solaris, and work-around for bad MODIS geo times.

### Source Changed
  * l12_parms.h, version number change.
  * l1_modis_hdf.c, replace atan2f wit atan2 in polarization alpha comp, to fix
  error in Solaris.  Use last good geolocation frame time when -1.0 is encountered.

## 5.1 - 2005-02-04

Upgrade evaluation changes to standard, in preparation for MODIS/Aqua
reprocessing.

### Source Changed
  * l12_parms.h, version number change.
  * getrayleigh.f, new pressure correction is only pressure correction.
  * atmocor2.c, removed old transmittance pressure correction.
  * aerosols.c, new transmittance pressure correction now standard.
  * polcor.f, new 1/(1-f) correction form now standard.
  * setflags.c, always subtract glint before cloud test, all sensors.
  * l1_modis_hdf.c, fixed handling of 65521 flag for non-rescaled L1A case.

## 5.0.8 - 2005-01-24

Rudimentary CZCS support.

### Source Changed
  * l12_parms.h, version number change. Add CZCS sensor identifier.
  * l1_czcs_hdf.c, new CZCS i/o routines.
  * l1_czcs_hdf.h, include for above.
  * MSl1info.c, add CZCS as new sensor and orbit = ASCENDING.  Fix csolz.
  * getformat.c, add CZCS branch and recognition
  * filehandle.h, add codes for czcs
  * l1_io.c, add calls to openl1_czcs, readl1_czcs, closel1_czcs
  * msl12_input.c, changed usage message.
  * MSscanpixlimits.c, changed usage message.
  * MSl1brsgen.c, changed usage message.
  * MSl2validate.c, changed usage message.

  * [libczcs] new czcs cal/nav support library.

### Data Added
  * czcs, new czcs sensor data directory


## 5.0.7 - 2005-01-11

Work-around for extract problem.

### Source Changed
  * l12_parms.h, version number change.
  * l1_modis_hdf.c, work-around for processing extracts with incorrect
  meta and/or sentinel values.

## 5.0.6 - 2005-01-04 

PAR bug fix.

### Source Changed
  * l12_parms.h, version number change.
  * get_par.c, bug fix (broke on v5.0.0).


## 5.0.5 - 2004-12-20

Handle new rescaled L1B MODIS.  Minor fixes.

### Source Changed
  * l12_parms.h, version number change.
  * l1_modis_hdf.c, handle new flags UNRESCALED_HIGH_SI and RESCALED_L1B_SI.
  * brdf.c, fix for fresnel_sol() implementation.
  * getglint.f, set lower limit on windspeed to prevent divide by zero.

## 5.0.4 - 2004-12-15

Fixes for OSX and Solaris ports.  Fixes for various discrepancies and
non-standard code reported by Paul Martinolich (beta tester). Update for
windspeed-dependence in Fresnel BRDF correction. New straylight filter.

### Source Changed
  * l12_parms.h, version number change.
  * aerosol.c, changed scope of compalphaT for OSX port.  Fixes for alternate aerosol
  selection scheme (models indexed from 0).
  * getrayleigh.f, use iand for bitwise test (more standard).
  * l1_mos_hdf.c, remove redundant memcpy of image to Lt.
  * l1_osmi_hdf.h, fixed prototype names.
  * msl12_input.c, fixed numTauas count in taua input.
  * MSl1brsgen.c, support for MODIS 6,5,3 rgb.
  * calc_par.f, replace jint, jmod with int, mod intrnsic functions.
  * l1_generic.c, replace memalign with malloc.
  * parse_file_name.c, informational.
  * l1_modis_hdf.c, add support for rescaled land bands and L1A extracts.
  * MSscanpixlimits.c, fix for L2 file identification support.
  * filehandle.h, add L2 file identification support.
  * getformat.c, add L2 file identification support.
  * MSl1brsgen.c,
  * brdf.c, remove return value on void function foqint.
  * atmocor1.f, break long line for solaris.
  * get_dem_height.f, break long line for solaris.
  * glint.c, no functions defined within functions for solaris.
  * water.c, changed return to exit on error.
  * get_dem_height.f, solaris did not like FILENAME_MAX
  * convl12.c, remove confusing OUTPUT_NLW option.
  * alloc_target.c, force initialization of record to 0.
  * alloc_aer.c, force initialization of record to 0.
  * rdsensorinfo.c, incorrect usage of memset (presumably ineffective?).
  * MSl12.c, corrected reporting for aer_opt selection.
  * brdf.c, add windspeed dependence to fresnel_sol().
  * atmocor2.c, passpass windspeed to brdf function.
  * convl21.c, pass windspeed to brdf function.
  * filter.c, add stlight filter function.
  * filter.h, add stlight filter.
  * l12_proto.h, at windspeed to ocbrdf().
  * polcor.f, alternate computation of polcor correction.
  * l1_octs_hdf.c, rename defines to reduce warning messages.
  * l1_hdf_generic_write.c, rename defines to reduce warning messages.

## 5.0.3 - 2004-11-02

### Source Changed
  * l12_parms.h, version number change.
  * gsm01.c, set CHLFAIL when optimized chl < 0.0.
  * ipar_arp.c, change integration of ipar an arp to use hc/lambda.
  * getrayleigh.f, clean-up.
  * filter.c, die if filter file not found.
  * msl12_input.c, usage message update.

## 5.0.2 - 2004-10-25

### Source Changed
  * l12_parms.h, version number change.
  * gsm01.c, limit bands used in optimization to those below 600nm.
  * convl21.c, generalize selection of visible wavelengths.
  * atmocor2.c, generalize selection of visible wavelengths.
  * input_struct.h, added metafile field for optional output meta-data file.
  * msl12_input.c, init and populate metafile field.
  * l2_hdf_generic.c, open and write to metafile.

## 5.0.1 - 2004-10-20

### Source Changed
  * l12_parms.h, version number change.
  * aerosol.c, change in form of transmittance pressure correction (no impact
  to standard processing).
  * msl12_input.c, usage message updated, require nband gains, offsets, tauas.
  * MSl2validate.c, generalized product list.
  * msl2val_hdf.c, generalized product list.
  * msl2val_struc.h, generalized product list.
  * read9km_mask.c, modified to handle other L3 resolutions.

## 5.0.0 - 2004-10-13

This is a major rewrite of the aerosol correction code, to improve flexibility
and allow for a variable number of sensor bands.  Much of the old code has been
rewritten in C.  Many old source files are no longer in use, and many new
source files have been added.  New data formats for Rayleigh and aerosol tables
have also been introduced.  Aerosol models are now stored in individual HDF
files.  Many band-related quantities have been moved to the relevant external
msl12_sensor_info.dat file.

### Source Changed

 The following source files are defunct

   * aeroplus.f
   * atmocor2.f
   * atmocor_set.f
   * calc_aero_refl.f
   * diffuse_transmittance.f
   * ff.f
   * fourier_a_b_c.f
   * franzaer.f
   * funct_eps.f
   * getaer.f
   * getaermod.f
   * getaerosol.f
   * get_brdf.f
   * getdtran.f
   * getglint_rad.f
   * get_pigment_nn.f
   * getwcap.f
   * gordon_tau_a.f
   * levmarq.f
   * linear_a_b_c.f
   * load_aer.f
   * load_ss11.f
   * lspline.f
   * nlw_outband.f
   * rhoa_from_taua.f
   * rho_a_sub_quad.f
   * spl1d.f
   * wangaer.f
   * wtaer.f
   * ylgint.f

 The following source files have been added:

   *  aerosol.c, all functions relating to aerosol model selection and use.
   *  airmass.c, airmass computations, including spherical variants.
   *  atmocor2.c, replaces fortan version of Lt to nLw function.
   *  brdf.c, all function relating to selection and use of ocean BRDF algorithm.
   *  fluorescence.c, fluorescence line height and fluorescence efficiency.
   *  glint.c, glint radiance function (replaces getglint_rad.f).
   *  ipar_arp.c, instantaneous PAR and ARP algorithms.
   *  lspline.c, linear interpolation and spline functions.
   *  nlw_outband.c, out-of-band nLw corrections (replaces nlw_outband.f).
   *  water.c, reads and interpolates water absorption and backscatter tables.
   *  whitecaps.c, whitecap radiances (replaces getwcap.f).
   *  windex.c, finds nearest wavelength table index for input wavelength.

 The following source files have been changed:

   * amoeba.c, add control for max iterations.
   * atmocor1.f, generalized band indexing, elimination of support for band-
   dependent view angles, new whitecaps() function.
   * atmocor_init.f, reduced parameter set in fortran common.
   * calcite.c, major rewrite, including addition of 2-Band algorithm.
   * carder.c, add bbp product, set CHLFAIL and CHLWARN flags as appropriate.
   * convl12.c, new atmocor2 interface (C function), eliminate atmocor_set().
   * convl21.c, generalized band indexing, new BRDF interface.
   * cpl1rec.c, new iwave, fwave copied (see l1_struc.h).
   * filter.c, generalized band indexing, clean-up.
   * get_chl.c, set CHLFAIL and CHLWARN flags as appropriate.
   * get_es.c, clean-up.
   * getglint.f, conflict with PI define.
   * get_l2prod_index.c, new products (see l2prod.h), generalized band indexing.
   * getozone.f, added but disabled spherical airmass for ozone.
   * get_par.c, generalized band indexing.
   * get_poc.c, set CHLFAIL and CHLWARN flags as appropriate.
   * getrayleigh.f, use standardized HDF Rayleigh tables for all sensors, which
   includes use of Q and U components.  Generalized band indexing. Elimination
   of support for band-dependent view angles.
   * get_rhown_nir.f, pass sensorID to avoid need for common block.
   * get_tsm.c, added MODISA, MODIST algorithm.  Set CHLFAIL and CHLWARN flags
   as appropriate.
   * gsm01.c, major rework with switch to amoeba optimization, generalization for
   non-SeaWiFS band, and expansion of output parameter suite to include IOPs
   at all sensor wavelengths. Set CHLFAIL and CHLWARN flags as appropriate.
   * input_struc.h, new fields include list of aerosol models to use and number
   of models, aermodels and naermodels, and model pair and ratio for new
   fixed model-pair option, aermodmin, aermodmax, and aermodrat.  Former
   taua input field has been converted to per band input, and former
   tau_a_per_band has been eliminated. New evalmask input for algorithm tests.
   * l12_parms.h, removed alternate or defunt aerosol option identifiers.  Added
   new BRDF option identifiers and EVAL option identifiers.  Also defined PI.
   * l12_proto.h, many new prototypes added, old prototypes deleted.
   * l12_seawifs.c, generalized band count.
   * l1a_osmi.c, generalized band indexing.
   * l1a_seawifs.c, generalized band indexing.
   * l1_hdf_generic_write.c, clean-up.
   * l1_modis_hdf.c, generalized band indexing. Radcor application disabled.
   Also added functions for vicarious RVS/destriping corrections, but disabled.
   * l1_mos_hdf.c, generalized band indexing.
   * l1_octs_hdf.c, generalized band indexing.
   * l1_polder.c, clean-up.
   * l1_struc.h, replaced wavelength array field wl with int and float variants
   iwave and fwave.
   * l2_hdf_generic.c, added new products for calcite, flh, ipar, arp, cfe,
   and gsm01.  Generalized band indexing.
   * l2prod.h, added new products for calcite, flh, ipar, arp, cfe, and gsm01.
   * l2_struc.h, replaced wavelength array field wl with int and float variants
   iwave and fwave, and added number of scans fields nscan.
   * loadl1.c, reduced parameters to atmocor_init.  Generalized band indexing.
   * MSl12.c, informational changes, new aerosol selection identifiers.
   * msl12_input.c, added eval, tau_a_per_band now called taua (old taua
   parameter dropped), generalized band indexing, added support for
   aermodels list of aerosol models to use.
   * polcor.f, generalized band indexing.
   * rdsensorinfo.c, added support for out-of-band coefs, water a and bb.
   * sensor_cmn.fin, removed many defunct variables from common.
   * setflags.c, disabled aerosol index flag, generalized band indexing.
   * sst.c, moved path to brightness temperature tables to cal subdir.
   * time_utils.c, added zulu2unix() function.
   * time_utils.h, added zulu2unix() function.

 The following utility-program source files have also changed:

   * MScalmerge.c
   * mscal_struc.c
   * mscal_struc.h
   * MScaltarget.c
   * MSl1info.c

 The following libraries are no-longer used by MSl12

   * libai, aerosol index
   * libcloud, cloud optical thickness

### Data Added

 The following data files are now defunct:

   * data/foq/*
   * data/*/transmittance
   * data/*/rayleigh/*pol_ray_B.dat
   * data/seawifs/aerosol/coef_quad*.dat
   * data/seawifs/aerosol/seawifs_aerosol*.dat (except seawifs_aerosol_par.dat).
   * data/seawifs/aerosol/abs_aerosol_seawifs.dat
   * data/modisa/aerosol/modis_aerosol*.dat
   * data/modist/aerosol/modis_aerosol*.dat
   * data/modisa/rayleigh/*iqu3.hdf
   * data/modist/rayleigh/*iqu3.hdf
   * data/modisa/bt_modisa.hdf (moved to cal subdir)
   * data/modist/bt_modisa.hdf (moved to cal subdir)

 The following files are new:

   * data/*/aerosol/*.hdf
   * data/*/rayleigh/*iqu.hdf
   * data/modisa/cal/bt_modisa.hdf
   * data/modist/cal/bt_modisa.hdf
   * data/common/morel_fq.dat
   * data/common/water_spectra.dat

 The following files have changed:

   * data/*/msl12_defaults.par
   * data/*/msl12_sensor_info.dat


## 4.0.4 - 2004-07-26

Add calculation for degree-of-polarization (dpol) and high-polarization flags (HIPOL)
based on threshold of dpol.  Cleaned-up product suite to eliminate many unneeded,
non-standard products including log and pigment variants of chlorophyll algorithms,
and the poorly documented neural-net chlorophyll algorithm.  Renamed the gs97 code
to gsm01 to better reflect the current algorithm state.  Updated 3-band calcite code.

### Source Changed
  * l12_parms.h, version number.
  * l2prod.h, remove old logchl prods, chl_nn, and pigment products  Add dpol. Change
  gs97 to gsm01. Remove cloud optical thickness, tauc.
  * get_l2prod_index.c, remove old logchl prods, chl_nn, and pigment products.  Add dpol.
  Remove alternate name gs97 for gsm01. Remove tauc.
  * l2_hdf_generic.c, remove old logchl prods, chl_nn, and pigment products. Add dpol.
  Remove tauc.
  * get_chl.c, remove old logchl prods, chl_nn, and pigment products.
  * calcite.c, update for out-of-band effects.
  * gs97.c, defunct.
  * gsm01.c, renamed from gs97.c. dropped logchl.
  * atmocor1.f, add dpol computation and return.
  * loadl1.c, return dpol from atmocor1.
  * l12_proto.h, add dpol to atmocor1() args.
  * l1_struc.h, add dpol field.
  * l2_struc.h, add dpol field.
  * alloc_l1.c, add dpol allocation.
  * convl12.c, copy dpol pointer to l2 record.
  * l2_flags.h, add HIPOL flag.
  * setflags.c, set HIPOL flag. Only call aeroindex for SeaWiFS (don't trust it).
  * input_struct.h, add hipol threshold for setting HIPOL flag.
  * msl12_input.c, add hipol threshold for setting HIPOL flag.
  * cpl1rec.c, copy input pointer.
  * filter.h, moved prototypes to l12_proto.h and removed include of l1q_struc.h to fix
  a conflict which prevented the input structure from being included in l1_struc.

  * [libcloud] no longer required, sice tauc product was removed.


## 4.0.3 - 2004-07-19

Add Carder semi-analytic chlorophyll model.  Update 3-band calcite algorithm.

### Source Changed
  * l12_parms.h, version number.
  * carder.c, Carder model and routines to interface with MSl12, function to read and interp. NDT.
  * l2prod.h, add support for Carder model products.
  * get_l2prod_index.c, add support for Carder model products.
  * l2_hdf_generic.c, add support for Carder model products. Update calcite call.
  * water_spectra.c, new routines for interpolating oceanic optical properties (e.g., aw of Pope & Fry)
  to specific wavelengths (courtesy of NRL).
  * spectral.h, include file for water_spectra.c.
  * table.h, include file for water_spectra.c.
  * l1_io.c, initialize sst.
  * l12_proto.h, carder function prototypes added, calcite function updated.
  * l2struc.h, added Rrs field.
  * alloc_l2.c, added Rrs field.
  * get_rrs.c, expanded to do full scanline processing as an option.
  * convl12.c, always compute Rrs.
  * calcite.c, add MODIS support, upgrade C-wrapper to become full algorithm (replaces get_calcite.c,
  coccolith.f).
  * get_dem_height.f, replace type statements with print statements (suggestion of Paul M.).
  * l1subpix.c, replace memcpy with safer memmove (suggestion of Paul M.)
  * get_calcite.c, defunct.
  * coccolith.f, defunct.
  * [libswfnav] *.f, replace type statements with print statements.
  * [libai] inter.f, remove unneeded readonly for linux (suggestion of Paul M.)

### Data Added
  * common/carder.par, Carder model coefficients.
  * common/ndt.hdf, global map of nitrogen depletion temperatures.



## 4.0.2 - 2004-05-10

Add ability to use Reynolds OISST binary files as input to NLSST algorithm.

### Source Changed
  * l12_parms.h, version number.
  * sst.c, new functions to determine sst reference file type, get_sstref(), and to read
  and interpolate Reynolds oisst files, get_oisst().
  * loadl1.c, change interface for loading sstref field in l2rec.
  * l12_proto.h, added get_sstref() proto.
  * get_l2prod_index.c, change sstref scaling from byte to int16.

### Data Added
  * common/sst_climatology.hdf, updated to fill mixed land/water pixels with neighboring water pixels.
  * modisa/msl12_defaults.par, updated brdf_opt, polfile, gain, albedo to reflect reprocessing.


## 4.0.1.eval - 2004-05-04

New EVAL pressure correction for Rayleigh.

### Source Changed
  * getrayleigh.f, added alternate pressure correction (Wang, 2004) to EVAL variant.
  * setflags.c, commented-out the alternate cloud test in EVAL variant.

## 4.0.1 - 2004-05-04 
Add brightness temperature products and change sst product scaling.
### Source Changed
  * l12_parms.h, version number.
  * l2prod.h, add BT_11, BT_12.
  * get_l2prod_index.c, add BT_11, BT_12, and change sst from byte to int16.
  * l2_hdf_generic.c, add BT_11, BT_12.

## 4.0 - 2004-04-05 
Change how defaults are specified and loaded by MSl12 and various MS utilities.
Upgrade new seawifs caltable support from EVAL to standard (remove piecewise
exponential support). Support for IMAPP L1B.  Fixes for OSMI and POLDER
calibration mode.  Fix for MOS under Fedora. BRDF option 4 is now f/Q +
Fresnel.
### Source Changed
  * l12_parms.h, version number, FOQFRESN.
  * atmocor1.f, pass polarization file name to modis_polcor (no more environment vars).
  * atmocor2.f, add brdf option for f/Q with Fresnel.
  * filehandle.h, add sensorDir[] array of sensor-specific data directories.
  * get_brdf.f, add brdf option for f/Q with Fresnel.
  * get_dem_height.f, accept dem filename as input (no more hard-coded filename).
  * getformat.c, add recognition of IMAPP direct broadcast L1B.
  * getrayleigh.f, moved polcor() function to polcor.f, add EVAL pressure correction.
  * input_struc.h, added polfile element, renamed elevfile to demfile.
  * l12_proto.h, msl12_input() protos, get_dem_height().
  * l1a_osmi.c, mask OSMI edges as HILT rather than NAVFAIL (it screws the meta), remove
  CAL_PATH environment variable usage
  * l1a_seawifs.c, load pixnum array in l1rec, remove CAL_PATH environment variable usage,
  upgraded new cal-table format support from EVAL to standard.
  * l12_seawifs.c, upgraded new cal-table format support from EVAL to standard.
  * l1_hdf_generic_write.c, fix indexing for less than 8 bands, upgraded new cal-table
  format support from EVAL to standard.
  * l1_hdf_generic_read.c, fix indexing for less than 8 bands.
  * l1_io.c, init pixnum field of l1rec.
  * l1_modis_hdf.c, remove CAL_PATH environment variable usage, clean-up.
  * l1_mos_hdf.c, make a large array in bilinear() static, to fix Fedora problem.
  * l2_hdf_generic.c,upgraded new cal-table format support from EVAL to standard.
  * loadl1.c, upgrade new E-S distance to standard processing, pass demfile and polfile.
  * Makefile, no more MSl1bingen.
  * MSl12.c, status messages.
  * msl12_input.c, extensive mods to load defaults from external files.
  * MSl12.mk, add polcor.c.
  * MSl1brsgen.c, call new msl12_input_defaults() to set input defaults.
  * MSl1brsgen.mk, add polcor.c.
  * MSl1info.c, call new msl12_input_defaults() to set input defaults.
  * MSll2snpx.c, call new msl12_input_defaults() to set input defaults.
  * polcor.f, new source file to consolidate polarization corrections.
  * rdsensorinfo.c, changed name of sensor info file, don't read gains, offsets.
  * setflags.c, reinstate EVAL version of glint-subtracted cloud test.
  * [libswfl1a] calibrate_l1a.c, upgraded new cal-table format support from EVAL to standard.
  * [libswfcal] get_cal.c, upgraded new cal-table format support from EVAL to standard.
  * [libswfcal] get_cal_misc.c, upgraded new cal-table format support from EVAL to standard.

  * get_cal.h, upgraded new cal-table format support from EVAL to standard.
  * getcal_proto.h, upgraded new cal-table format support from EVAL to standard.

### Data Added
  the following are new data files
  * seawifs/msl12_defaults.par
  * seawifs/msl12_sensor_info.dat
  * octs/msl12_defaults.par
  * octs/msl12_sensor_info.dat
  * osmi/msl12_defaults.par
  * osmi/msl12_sensor_info.dat
  * mos/msl12_defaults.par
  * mos/msl12_sensor_info.dat
  * polder/msl12_defaults.par
  * polder/msl12_sensor_info.dat
  * modisa/msl12_defaults.par
  * modisa/msl12_sensor_info.dat
  * modist/msl12_defaults.par
  * modist/msl12_sensor_info.dat

  the following files were renamed
  * seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200501 to SEAWIFS_SENSOR_CAL.TBL-200404
  * seawifs/seawifs_filter.dat to seawifs/msl12_filter.dat
  * octs/octs_filter.dat to octs/msl12_filter.dat

  the following files were moved
  * modisa/rayleigh/*polsen* to modisa/cal
  * modisa/rayleigh/*pol_cor* to modisa/cal
  * modist/rayleigh/*polsen* to modist/cal
  * modist/rayleigh/*pol_cor* to modist/cal

### Defunct Data
  the following files are no longer used
  * seawifs/seawifs_def_l2prod.dat
  * seawifs/seawifs_table.dat
  * octs/octs_def_l2prod.dat
  * octs/octs_table.dat
  * osmi/osmi_def_l2prod.dat
  * osmi/osmi_table.dat
  * mos/mos_def_l2prod.dat
  * mos/mos_table.dat
  * polder/polder_def_l2prod.dat
  * polder/polder_table.dat
  * modisa/modisa_def_l2prod.dat
  * modisa/modisa_table.dat
  * modist/modist_def_l2prod.dat
  * modist/modist_table.dat
  * seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200206
  * seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200401

  the following environment variables are no longer used
  * MODISA_POL_PATH
  * MODIST_POL_PATH
  * MODISA_CAL_PATH
  * MODIST_CAL_PATH
  * CAL_HDF_PATH
  * OSMI_CAL_PATH


## 3.5.4.eval - 2004-03-10
EVAL support for new cal table format.
### Source Changed
  * l1a_seawifs.c, l12_seawifs.c, l1_hdf_generic_write.c, l2_hdf_generic.c
  * [libswfl1a] calibrate_l1a.c
  * [libswfcal] get_cal.c
  * [libswfcal] get_cal_misc.c

  * get_cal.h, getcal_proto.h


## 3.5.3 - 2004-03-02
Support L2 processing of XCAL format.  Add output for full polarization
radiances (controlled by pol_opt). Allow redirection of polarization tables via
environment variables MODSA_POL_PATH, MODIST_POL_PATH.  Add additional output
products detnum, pixnum, mside, and polarization frame rotation angle alpha.
### Source Changed
  * l12_parms.h, version number.
  * l1_io.c, support new XCAL input functions.
  * l2prod.h, changed Lr_q and Lr_u to L_q and L_u.  Added detnum, mside,
  pixnum, and alpha products.
  * get_l2prod_index.c, changed Lr_q and Lr_u to L_q and L_u.  Added mside,
  detnum, pixnum, alpha.
  * l1_struc.h, changed Lr_q and Lr_u to L_q and L_u.
  * alloc_l1.c, changed Lr_q and Lr_u to L_q and L_u.
  * loadl1.c, changed Lr_q and Lr_u to L_q and L_u.
  * l2_struc.h, changed Lr_q and Lr_u to L_q and L_u.
  * convl12.c, changed Lr_q and Lr_u to L_q and L_u. Copy mside, detnum to
  L2 record.  Initialize start and end node meta-data to prevent seg faults
  when all navigation is bad.
  * l2_hdf_generic.c, changed Lr_q and Lr_u to L_q and L_u. Added detnum, mside,
  pixnum, and alpha products.
  * atmocor1.f, changed Lr_q and Lr_u to L_q and L_u.
  * getrayleigh.f, control MODIS polarization tables via env vars.
  * convl12.c, copy mside, detnum to L2 record.
  * l1_xcal_hdf.c, new XCAL input routine.
  * l1_xcal_hdf.h, new XCAL input routine header.
  * MSl2validate.c, new program to generate cross-sensor validation.
  * msl2val_hdf.c, new program to generate cross-sensor validation.
  * msl2val_struc.h, new program to generate cross-sensor validation.

## 3.5.2 - 2003-12-11
Fix bug in meta-data for files less than 3 scanlines in length. Fix bigger bug
in cloud test for modis.
### Source Changed
  * l12_parms.h, version number.
  * l1_hdf_generic.c, changed some elseifs to ifs.
  * setflags.c, bad indexing of solar irradiance on glint subtraction.
  * convl12.c, fix for inverse calibration with nLw=0.


## 3.5.1 - 2003-11-18 
Fix bug in SST masking.
### Source Chnaged
  * l12_parms.h, version number.
  * convl12.c, skip SST processing for any L1 masking condition.
  * l1_modis_hdf.c, set HILT flag if 11 or 12um channels saturated or invalid.

## 3.5 - 2003-11-04
Update to fix very old bug in SeaWiFS calibration library, which sometimes
mis-selected segments in a segmented calibration table (not an issue for repro
#4). Also, the sst product was separated from the input sst reference product,
and the latter is now available as sstref. SST processing was enabled by
default for the MODIS sensors. Installed sensor-specific defaults for cloud
threshold.  Eliminated EVAL code for subtracting glint before testing clouds
(now MODIS default behavior, not used for other sensors until it can be
tested). For MODIS processing, transition to seawifs-like glint correction, add
glint polarization correction, use new polarization tables for aqua (aqua_1a).
Transition of MODIS Aqua from Neckel and Labs to Thuillier solar irradiance.
Updated earth-sun distance function (was 1998-specific) for all but SeaWiFS
(until approved for SeaWiFS).
### Source Changed
  * l12_parms.h, new version, new definitions for MScalmerge, MScaltarget.
  * l12_proto.h, modified atmocor1() arg list.
  * l2_struc.h, added sstref field to l2 structure.
  * l1_struc.h, added sstref field to l1 structure.
  * l2prod.h, added sstref product specifier.
  * sst.c, user sstref for first guess sst in NLSST algorithm.
  * msl12_input.c, enabled sst procesing by default for MODIS sensors. Set proper default
  for MODIS cloud threshold. Incorporated additional changes to support MScalmerge. Changed
  default for MODIS polarization to Rayleigh+glint (pol_opt=3).
  * alloc_l1.c, allocate space for sst reference field.
  * alloc_l2.c, clean-up.
  * l2_hdf_generic.c, support for sstref product output.
  * get_l2prod_index.c, support for sstref product output.
  * get_tricho.c, use sstref rather than sst field.
  * convl12.c, copy sstref pointer.  Remove continue statement on ATMFAIL.  This was causing
  some pixels to skip sst processing, leaving the reference value in the sst field, and
  skipping SST flagging.  Cleaned-up some dead code.
  * setflags.c, always subtract glint before cloud test for MODIS (removed EVAL test). Don't
  modify the windspeed for MODIS.
  * getglint.f, use Cox and Munk, unless glint_opt>=2 (removed EVAL test)
  * atmocor2.f, handle glint_opt (glintON) equal to 2 or 3 (Ebuchi and Kizu).  Apply seawifs
  glint correction to MODIS. Use Thuilier for nominal F0 for all MODIS.
  * input_struc.h, expanded ifile list and parms and files meta strings to handle larger input
  file number, to support MScalmerge.
  * filehandle.h, new mode to identify L1XCAL (MScalmerge) files.
  * MScalibrate.c, development update.
  * mscal_func.c, development update.
  * mscal_struc.c, development update.
  * mscal_struc.h, development update.
  * getformat.c, add recognition for match-up file format.
  * loadl1.c, sst reference (climatology) loaded into sstref fields of l1rec.
  Use new e-s distance function for EVAL or non-SEAWIFS cases. Use Thuilier
  for nominal F0 for all MODIS.
  * l1_modis_hdf.c, use new e-s distance function.
  * l1_polder.c, use new e-s distance function.
  * getwcap.f, apply seawifs white-cap correction to MODIS
  * get_f0.c, split into two function, get_f0_neckel and get_f0_thullier.
  * atmocor1.f, moved glint calculation here from setflags(), needed for polcor().
  * setflags.c, glint coef and radiance calc moved elsewhere. Glint subtraction for cloud mask
  test now applied to all sensors in EVAL code, modis-only in standard code (until it can
  be evaluated for seawifs).  The glint subtraction for the cloud test now includes the
  direct Rayleigh transmittance.
  * getglint.f, added new glint coefficient function including polarization.
  * l1subpix.c, added alpha, Bt, and Ltir subpixelization support for polarization and thermal.
  * MScalmerge.c, new error codes, ability to read and append to match-up files.
  * getrayleigh.f, switch to new aqua_1a polarization file for MODIS/Aqua.
  * parse_file_name.c, don't overwrite a string.
  * MScaltarget.c, new error code. No pixel averaging. Better handling of polar cases.
  * readL2scan.c, defunct.  Now linking to libl2.
  * readL2scan.h, defunct.  Included from swfinc.

  * [libswfcal]  get_cal_misc.c, fixed cal segment selection.

  * [libswfnav] esdist.f, new earth-sun distance correction, valid for all years.

  * [libl2] readL2scan.c, updates to support MScaltarget.

  * readL2scan.h, updates to support MScaltarget.

### Data Added
  * common/f0_table.dat, defunct
  * common/f0_table_neckel.dat, renamed f0_table.dat
  * common/f0_table_thuillier.dat, new solar irradiance from Thuillier 2003.
  * modisa/modisa_table.dat, new Thuillier band-averaged F0 using new aqua RRS.
  * modisa/rayleigh/modis_pol_corr_aqua_1a.hdf, new polarization table
  * modisa/rayleigh/new_modis_pol_corr5a.hdf, defunct


## 3.4 - 2003-10-16
Rework of calibration modes for MSL12.  Now write HDF file for vicarious L1B, including
the l2_flags SDS. Had to rewrite MSl1bgen in the process.  Generalized f/Q and Fresnel
options and products. Eliminated flat binary output option for L2.  Added new MScaltarget,
MScalmerge, and MScalibrate programs, and modified some of the i/o functions to support
the new utilities.
### Source Changed
  * l12_parms.h, new version. New aerosol selection ID FIXMODPAIR.
  * MSl12.c, now writes vicarious L1B in HDF format.  Many changes to handle new parameter
  specifiers: ofmt eliminated, mode added to indicate FORWARD or INVERSE modes.
  * MSl1brsgen.c, revamped to use standardized i/o function calls.
  * MSl1bingen.c, init input struct and link to l1file.
  * MSl1info.c, init input struct and link to l1file.
  * MSll2snpx.c, init input struct and link to l1file.
  * alloc_l1.c, added allocation for flags array.
  * alloc_l2.c, removed flags allocation (moved to l1rec). Consolidated t_f and foq
  allocation to brdf field.
  * atmocor2.f, replace t_f and foq variables with brdf, moved fresnel code to get_brdf.f,
  no fresnel applied if full f/Q is requested. Eliminated t_f from argument list.
  * atmocor_init.f, replaced foqOn with brdf_opt.
  * convl21.c, replace t_f and foq variables with brdf, revamped to use new processing
  mode parameter from filehandle to determine how to invert.
  * convl12.c, improved l2_flags handling: set l2rec flag pointer to l1rec flag array,
  call new setflagbits() function (pre-atmocor2 flags set in loadl1), moved isTricho(),
  isCoccolith(), isTurbid() functions to setflags.c.  Moved easternmost(), westernmost()
  to l1_io.c (multiple instances were consolidated).  Changed call to atmocor2().
  * getdtran.f, comment-out error message for theta out-of-range
  * filehandle.h, new format extended L1 HDF, new mode indicators FORWARD, INVERSE_*, CALFIT.
  * getformat.c, initialize eosmeta to avoid warnings.
  * input_struc.h, added additional parameters to support MScaltarget.  Removed ofmt field
  (and option to output flat binary L2 files).  Extended mode option to indicate FORWARD
  or various types of INVERSE processing for MSl12.  Replaced foq_opt with generalized
  brdf_opt.  Replaced field target_nlw_file with tgtfile, and eliminated field
  new_l1b_file (now stored in ofile[0].
  * l12_proto.h, various function prototypes updated.
  * l1_generic.c, make sure pointer is not NULL before closing (flat binary L1).
  * l1_hdf_generic_read.c, added reading of mside and detnum SDSes.
  * l1_modis_hdf.c, cleaned-up radcor input to only read the one epoch of interest.
  * l1_struc.h, added flags field.
  * l2_struc.h, replace t_f and foq with brdf field.
  * msl12_input.c, changes to support MScaltarget and new mode and brdf options (see
  input_struc.h).  Consolidated usage messages here.
  * loadl1.c,a dded call to the flag-bit setting function (moved from convl12()).
  * l1_io.c, added support for l1_hdf_generic_write().  Also added call to the flag-bit
  setting function in readl1(). Moved easternmost(), westernmost() here.
  * l1_hdf_generic_write.c, revamped to standard style of other i/o routines, add
  writing of detnum and mside.
  * l1_hdf_generic_read.c, read detnum and mside.
  * l2prod.h, removed t_f and foq products and added brdf product.
  * get_l2prod_index.c, removed t_f and foq products and added brdf product.
  * l2_hdf_generic.c, removed t_f and foq products and added brdf product.
  * sensor_cmn.fin, changed foqOn to brdf_opt.
  * setflags.c,  new setflagbits() routine to consolidate the l2_flag setting.  Moved
  l1_mask() and several L2 flagging functions to this source file for ease of maintenance.

  * MScaltarget.c, new calibration target generation code
  * readL2scan.c, new code to support MScaltarget
  * readL2scan.h, new code to support MScaltarget
  * setupflags.c, new code to support MScaltarget
  * get_brdf.f, new code consolidation of brdf functions

  * MScalibrate.c, new calibration fitting facility.
  * mscal_func.c, new code to support MScalibrate

  * get_foq.f, defunct (moved functions to get_brdf.f)
  * l1_hdf_generic_write.h, defunct (moved into l12_proto.h)
  * l1_mask.c, defunct (function moved into setflags.c)
  * l1_mask.h, defunct (function moved into l12_proto.c)
  * l1a_seawifs.c, save mirror side in l1rec. Use mean dark restore rather than
  dark restore per scan for testing HILT.
  * l2_generic.c, defunct (eliminated non-hdf output options for MSl12)
  * l2_io.c, defunct (eliminated non-hdf output options for MSl12)
  * l2_small_generic.c, defunct (eliminated non-hdf output options for MSl12)
  * l2_small_generic.h, defunct (eliminated non-hdf output options for MSl12)
  * msl1b_input.c, defunct (replaced with msl12_input())

  * [libcalfit] new set of functions to support MScalibrate
  * new include dir (msinc) to support MScalibrate


## 3.3.1  14 August 2003
Bug fix for calibration mode.
### Source Changed
  * l12_parms.h, new version.
  * l1_generic.c, check for NULL file pointer before closing l1 file.
  * convl21.c, missing parens.

## 3.3 - 2003-08-13
Add MODIS SST support. Also add f/Q in calibration mode.
### Source Changed
  * l12_parms.h, new version, NBANDSIR, new SST and BT flag values.
  * l1_modis_hdf.c, read 11, 12um SST bands, apply radcor.
  * l1_struc.h, added Bt field.
  * alloc_l1.c, allocate Bt field.
  * atmocor2.f, set frensel tranmittance to 1.0 when applying f/Q.
  * convl12.c, call sst function, if requested.
  * convl21.c, apply f/Q correction, if requested, and pass in input struct.
  * input_struc.h, added proc_sst control.
  * l12_proto.h, new sst funcs, and changed convl21 proto.
  * l2_flags.h, new SSTWARN and SSTFAIL flags.
  * l2_struc.h, new Bt and sst fields.
  * loadl1.c, always load sst climatology.
  * MSl12.c, pass input struc to convl21().
  * msl12_input.c, new proc_sst input.
  * sst.c, functions for MODIS 11um SST (Aqua/Terra) and brightness temperature.

### Data Added
  * modist/bt_modist.hdf, terra radiance to brightness temperature table.
  * modist/cal/radcor_modist_2003224.hdf, terra combined radcor
  * modist/cal/radcor_modist_nocor.hdf, non-correcting radcor for terra
  * modisa/bt_modisa.hdf, aqua radiance to brightness temperature table.
  * modisa/cal/radcor_modisa_2003224.hdf, aqua combined radcor
  * modisa/cal/radcor_modisa_nocor.hdf, non-correcting radcor for aqua

## 3.2.1 - 2003-07-18
Bug in application of radcor. Wrong mirror-side RVS correction applied to pixels 1327-1354
### Source Changed
  * l12_parms.h, new version.
  * l1_modis_hdf.c, fix to RVS correction on scan edge.

## 3.2 - 2003-07-03
Bug fix for seg fault on very small SeaWiFS L1A files, other clean-up to
reduce benign memory access errors.  Updates for MODIS nLw out-of-band and
NIR water-leaving radiance calculations. Also, for the EVAL version of
MSl12, replaced Cox & Munk with Ebuchi & Kizu glint distribution, and
changed cloud flag algorithm to subtract glint.

### Source Changed
  * l12_parms.h, version number change.
  * atmocor2.f, call rhown NIR correction wrapper.
  * get_rhown_nir.f, added version for MODIS bands and sensor selection wrapper.
  * getglint.f, for EVAL version, replace Cox & Munk probablility distribution
  with Ebuchi & Kizu.
  * l1a_seawifs.c, don't pass dark_rest to calibrate_l1a(), as it doesn't
  use it.  Also, compute HILT test with mean dark restore (slight impact to
  standard SeaWIFS processing). Free dark_rest as soon as dark_mean is
  computed.
  * msl12_input.c, new defaults outband_opt=2, aer_opt=-3 for MODIS.
  * nlw_outband.f, added MODIS function, cleaned-up useless options.
  * setflags.c, for EVAL version, subtract glint from cloud reflectance.

  * [libswfl1a] calibrate_l1a.c, removed dark_rest argument.
  * [libswfl1a] inc/swfinc/call1a_proto.h, removed dark_rest argument from calibrate_l1a()
  * [libanc] getanc.c, init some variables relating to TOMS ozone test, as they get
  tested for non-TOMS data.

## 3.1 - 2003-06-26
MODIS support (phase 1). Miscellaneous clean-up and generalization.
### Source Changed
  * MSl12.c, attach geolocation filename pointer to file handle. Changed
  name of usage function, to avoid conflicts.
  * MSl1bgen.c, attach geolocation filename pointer to file handle. Also,
  call msl12_input() rather than msl1b_input(), to consolidate input
  parameter initialization.
  * MSl1bingen.c, accept geolocation file as additional input.
  * MSl1brsgen.c, accept geolocation file as additional input.
  * MSl1info.c, accept geolocation file as additional input.
  * MSll2snpx.c, accept geolocation file as additional input. Also, the
  -f command option was removed, as the current code can not distinguish
  commandline options from negative lon/lat positional parameters.
  * alloc_l1.c, allocate space for polarization fields and pixnum. Eliminate
  allocation for Es.
  * atmocor1.f, added support for modis polarization correction.
  * atmocor2.f, disable fresnel correction for modis and call modis-
  specific glint correction, to match standard modis code. Clean-up
  fresnel code. Remove Es parameter and calculation, only done on request.
  * atmocor_init.f, added polOpt polarization option.
  * calc_aero_refl.f, call sensor-specifc (ff) water-vapor correction.
  * convl12.c, copy new l1_struc fields to l2_struc. Pass polcor
  to atmocor2(). Use precomputed F0 in isTurbid().
  * cpl1rec.c, support new Fo fields and detnum, mside.
  * diffuse_transmittance.f, locate MODIST and MODISA files.
  * ff.f, added MODIS-specific out-of-band water-vapor corrections.
  * filehandle.h, added modis L1B file ID. geofile field.
  * filehandle_init.c, support geofile.
  * funct_eps.f, don't interpolate for MODIS, since models now MODIS-specific.
  * get_angstrom.c, added values for MODIS wavelengths (splined from SeaWiFS).
  * getaermod.f, call sensor-specifc (ff) water-vapor correction.
  * getaerosol.f, moved data statement so compiler stops complaining.
  * get_chl.c, added OC3M modis chlorophyll and MOS OC4 variant. Change
  to make use of precomputed solar irradiances.
  * get_depth.c, use precomputed F0 for normalization.
  * get_l2prod_index.c, added new OC3M chlorophyll and pigment products,
  polarizion products polcor and Lr_q, Lr_u, tsm_clark, poc_clark.
  * get_poc.c, new functions for particular organic carbon.
  * get_rhos.c, added missing comment end.
  * get_rrs.c, use precomputed F0 for normalization.
  * get_tsm.c, new functions for total suspended matter.
  * getdtran.f, add error message.
  * getformat.c, added support for MODIS/Terra and MODIS/Aqua HDF-EOS
  format.
  * getglint_rad.f, added routine to mimick standard modis glint
  correction (though the implementation is questionable).
  * getrayleigh.f, added new rayleigh_iqu routine, only used for
  modis, to include full Rayleigh stokes vector.
  * getwcap.f, added modis-specific coefficients and scaling. No
  windspeed limit for modis.  These changes were required to
  match standard modis code, but differences need review.
  * gs97.c, use precomputed F0 for normalization.
  * hdf_utils.c, moved some general HDF routines from l1_octs_hdf.c.
  * hdf_utils.h, moved some general HDF routines from l1_octs_hdf.c.
  * input_struc.h, added geofile and pol_opt (poalrization option).
  * l12_parms.h, version number change.  Added new sensors MODIST and
  MODISA.
  * l12_proto.h, add rdsensorinfo(), remove rdsensortab(), rdsensor_wl().
  Added parameters to atmocor1, atmocor2, atmocor_init to support modis.
  * l12_seawifs.c, fixed bug in l1b meta-data transfer.
  * l1_hdf_generic_read.c, standarized to use rdsensorinfo().
  * l1_hdf_generic_write.c, standarized to use rdsensorinfo().
  * l1_io.c, call new MODIS i/o routines. Generalized determination of
  nbands, bindx (band indexing) using rdsensorinfo().
  * l1_modis_hdf.c, new MODIS input routines.
  * l1_modis_hdf.h, new MODIS input routines.
  * l1_mos_hdf.c, removed header settings (handled in getFormat, l1_io).
  * l1_octs_hdf.c, moved hdf utilities to hdf_utils.c.
  * l1_octs_hdf.h, moved hdf utilities to hdf_utils.h.
  * l1_struc.h, added detector number (detnum), mirror side (mside),
  pixel number (pixnum), polarization frame rotation angle (alpha),
  Rayleigh Q and U components (Lr_q, Lr_u), and mean and nominal
  solar irradiances (Fobar, Fonom) as fields in the L1 record.
  * l1subpix.c, include the new pixnum field when pixel sub-setting and
  sub-sampling.
  * l1a_seawifs.c,  change navwarn, navfail setting scheme (see also version 3.0.2).
  * l2_hdf_generic.c, use of rdsensorinfo() in place of rdsensortab() to
  standardize the reading of the sensor information table. Support for
  new products in l2prod.h.
  * l2_struc.h, added detector number, mirror side as fields in the
  L2 structure.  Also added pointers to new L1 fields for pixel
  number, alpha, Lr_q, Lr_u, Fobar, and Fonom. Removed Es field.
  * l2prod.h, added OC3M chlorophyll, and a MOS variant of OC4. Also
  added Q and U Rayleigh components and polcor, the polarization
  correction as a product. Also D. Clark's TSM and POC.
  * load_aer.f, read MODIS-specific aerosol tables.
  * loadl1.c, use of rdsensorinfo() in place of rdsensortab() to
  standardize the reading of the sensor information table. Do terrain
  height correction to geolocation whenever land processing is enabled,
  regardless of pixel type (land or water).  This should improve
  navigation for high-altitude lakes.
  * msl12_input.c, add modis support, including geo-location and radcor
  file inputs, and polOpt polarization option.  Also standardize code
  for obtaining default vicarious gains and offsets. Added ofile (no
  number) so the same routine can be used by MSl1bgen.  Added l2prod
  (no number) so ofile, l2prod will now work for MSl12.
  * msl1b_input.c, defunct.
  * rdsensorinfo.c, generalized routine to replace rdsensortab, rdsensor_wl.
  * rdsensortab.c, defunct.
  * rdsensor_wl.c, defunct.
  * rho_a_sub_quad.f, added calls to MODIS and SeaWiFS-specific water-vapor
  correction.  Removed large block of long-dead code associated with
  test and re-use of previous model pair.
  * rhoa_from_taua.f, added calls to MODIS and SeaWiFS-specific water-vapor
  correction.
  * sensor_cmn.fin, added polOpt polarization option parameter.
  * setanc.c, do terrain height correction to pressure whenever land
  processing is enabled, regardless of pixel type (land, water).  This
  should allow processing over high-altitude lakes.
  * setflags.c, for modis only, increae windspeed 48% before calling
  getglint() (cox & munk).  This is to match modis standard code, but
  we will likely remove it later.

## Data Added
  * new data/modist directory and associated files.
  * new data/modisa directory and associated files.
  * data/seawifs/seawifs_table.dat, added Nband, Bindx info.
  * data/octs/octs_tabe.dat, added Nband, Bindx info.
  * data/polder/polder_table.dat, added Nband, Bindx info.
  * data/mos/mos_table.dat, added Nband, Bindx info.
  * data/osmi/osmi_table.dat, added Nband, Bindx info.


## 3.0.1 - 2002-08-03
Fixes for solaris compatibility. New depth algorithm.
### Source Changed
  * l12_parms.h, version number change.
  * msl12_input.c, initialize filehandle structure. Fix calfile init for extraneaous
  open on l1 file.
  * atmocor1_land.c, replace cosf, sqrtf, logf, expf with generalized variants.
  * get_rhos.c, replace cosf, sqrtf, logf, expf with generalized variants.
  * MSl1brsgen.c, fixed -p option.
  * get_depth.c, new depth algorithm from R. Stumph.  Depth is now float.
  * get_l2prod_index.c,  depth is now float.
  * l12_proto.h,  depth returns float.
  * l2_hdf_generic.c,  depth returns float.


## 3.0 - 2002-06-03
EVAL version becomes standard version. Old code removed.
### Source Changed
  * l12_parms.h, version number change.
  * atmocor2.f, convl12.c, get_chl.c, get_rhown_nir.f, get_rrs.c,
  getaerosol.f, gs97.c, l12_seawifs.c, l1_hdf_generic_write.c,
  l1a_seawifs.c, l2_hdf_generic.c, msl12_input.c, nlw_outband.f,
  setflags.c, remove non-eval code and #ifdef EVAL wrappers.
  * setflags.c, msl12_input.c, remove old cloud algorithm.
  * atmocor2.f, msl12_input.c, get_rhown_nir.f, MSl12.c, l12_parms.h,
  getaerosol.f, atmocor_init.f, remove Siegel NIR code and option.
  * get_tlw_nir.f, defunct code module removed.
  * get_l2prod_index.c, change chlorophyll and IOP products to floats.
  * l12_parms.h, change chl range from 0.01-64 to 0.0-100, for use in
  setting CHLWARN flag.
  * rho_a_sub_quad.f, added kludge for aerosol model selection problem
  when T99 and C50 cross.
  * l2_struc.h, atmocor2.f, l12_proto.h, convl12.c, alloc_l2.c, convl21.c,
  add Fresnel transmittance correction to nLw normalization.
  * msl12_input.c, maskstlight=1 is now default. Filtering for SeaWiFS
  is off by default.
  * get_l2prod_index.c, l2prod.h, l2_hdf_generic.c, add t_f (Fresnel
  transmittance) as output product.
  * get_chl.c, added spectral test to failure conditions for OC4.
  * atmocor2.f, back-out Fresnel correction if using Morel f/Q.
  * [libswfcal] get_cal.c, get_cal_misc.c, remove non-eval code and #ifdef EVAL
  wrappers.
  * [libswfl1a] calibrate_l1a.c, gac_st.c, remove non-eval code and #ifdef EVAL
  wrappers.
  * cf.h, get_cal.h, getcal_proto.h, remove non-eval code and #ifdef
  EVAL wrappers.


## 2.9.5 - 2002-03-21
### Source Changed
  * l12_parms.h, version number change. New AeroPlus aerosol identifiers
  * get_tlw_nir.f, update to adg670 and ap670 terms in Arnone.
  * get_rhown_nir.f, new function to compute NIR reflectance rather than
  radiance (EVAL only).
  * atmocor2.f, compute tLw_NIR from rhown_NIR (EVAL only). Use actual
  transmittance terms from last iteration.
  * filehandle.h, new recalibration mode TARGET_ZERO
  * convl21.c, new recalibration mode TARGET_ZERO
  * MSl12.c, don't open recal target file if TARGET_ZERO
  * msl12_input.c, don't require target file to enable recal mode. Change
  aer_iter_max default from 50 to 10. Set aerosol option to USERTAUAS if
  tau_a_per_band is specified.

## 2.9.4 - 2002-03-13
### Source Changed
  * l12_parms.h, version number change. New AeroPlus aerosol identifiers
  * aeroplus.f, new function implementing AeroPlus algorithm of J. O'Reilly.
  * getaerosol.f, call aeroplus() if requested.
  * MSl12.c, informational message added for AeroPlus.
  * get_tlw_nir.f, added adg670 computation to Arnone algorithm.
  * atmocor2.f, compute and pass Rrs555 to get_tlw_nir_arnone().
  * l2prod.h, new logchl products.
  * get_l2prod_index.c,  new logchl products.
  * l2_hdf_generic.c, new logchl products.
  * get_chl.c, new logchl products.
  * gs97.c, support logchl for GSM01 model results.
  * l12_proto.h, added chl2log() proto.
  * rho_a_sub_quad.f, removed those annoying parabola messages
  * atmocor_init.f, recognize new aerosol options

### DataAdded
  * common/watermask.dat, new etopo2 bathymetry

## 2.9.3 -  2002-02-20
### Source Changed
  * l12_parms.h, version number change.
  * msl12_input.c, usage message changes.
  * [libai] residue.f, new thresholds from C. Hsu.  (this will change standard version l2_flags as well)


## 2.9.2 - 2002-02-04
### Source Changed
  * l12_parms.h, version number change.
  * get_foq.f, fixed errors in Rgoth, f/Q interpolation. Switched CHL
  interpolation to log-linear.
  * atmocor2.f, modified NIR iteration control for EVAL. Stop when chl change is
  less than 2%.  Add damping between iterations. Reset chl and Rrs670 on ANY
  iteration where chl retieval fails.
  * msl12_input.c, changed default aer_iter_max to 50.
  * l1a_seawifs.c, changed sl_pixl default for GAC to 4
  * setflags.c, if cloud albedo threshold is negative, use rhos(865) test.
  * MSl12.c, fix handling of spixl when processing single-pixel files.

## 2.9.1 - 2002-01-29
### Source Changed
  * l12_parms.h, version number change.
  * get_tlw_nir.f, split code into Siegel and Arnone variants.
  * atmocor2.f, include Arnone NIR option.
  * atmocor_init.f, recognize Arnone option.
  * msl12_input.c, change default rhoamin to 0.0001, satzen to 60.
  * MSl12.c, informational messages for Arnone NIR correction.
  * [libswfl1a] gac_st.c, EVAL version modified to correct and not flag two additional pixels.
  * cf.h, include STLT correction factors for additional along-track pixels.

## 2.9.0 - 2002-01-26
### Source Changed
  * l12_parms.h, version number change.
  * nlw_outband.f, new seawifs correction coefficients (as before, but this
  time they're in the right place).
  * get_foq.f, modified to read and apply the latest f/Q tables.
  * foq.fin, new f/Q include file.
  * input_struc.h, added input fields to allow user over-ride of ancillary
  met and ozone fields with a fixed value.
  * msl12_input.c, code to init and set met and ozone over-ride parameters.
  Also added additional parameters to the parameter meta-data list.
  * loadl1.c, over-ride met and ozone fields from input struct.
  * convl12.c, raise EVAL turbid threshold on Rrs670 to 0.0012.
  * l1a_seawifs.c, limit HILT test to NIR bands.
  * getaerosol.f, set models to zero on failure (effectively re-initialize
  between NIR iterations). ** This changes eps_78 field in current standard
  processing, but only where ATMFAIL was already set **

### Data Added
  * data/foq/fq.dat, new f/Q table
  * data/foq/rgoth.dat, new gothic R table
  * data/foq/q0.dat, new Q0 table (not required)
  * data/foq/f.dat, new f table (not required)

## 2.8.9 - 2002-01-09

### Source Changed
  * l12_parms.h, version number change. Added band width parameter for
  nLw out-of-band correction.
  * get_f0.c, allow specification of band width.
  * l12_proto.h, modified get_f0() proto.
  * convl12.c, add band width parameter to get_f0() call.
  * get_chl.c, add band width parameter to get_f0() call.
  * atmocor2.f, add band width parameter to get_f0() call.
  * get_rrs.c, add band width parameter to get_f0() call.
  * gs97.c, add band width parameter to get_f0() call.
  * setflags.c, fixed transmittance usage for cloud albedo comp.
  * nlw_outband.f, new seawifs correction coefficients.
  * filter.c, minor fix.

## 2.8.8 - 2002-01-04

### Source Changed
  * l12_parms.h, version number change.
  * get_l2prod_index.c, meta-data and product name change for gs97 products to gsm01, Garver-Siegel-Maritorena 2001.  Algorithm was already gsm01.
  * get_tlw_nir.f, new bbp function from S. Bailey.
  * input_struc.h, added rhoamin parameter.
  * msl12_input.c, set rhoamin parameter.
  * sensor_cmn.fin, added space for rhoamin.
  * atmocor_init.c, pass rhoamin to sensor_cmn
  * l12_proto.h, new arg to atmocor_init()
  * loadl1.c, pass rhoamin to atmocor_init()
  * getaerosol.f, use rhoamin from common.
  * convl12.c, changed isTurbid algorithm for EVAL version.
  * filter.c, changed epsmean filter to epsiqmean
  * filter.h, added epsiqmean() prototype.

## 2.8.7 - 2001-11-29

Fix OCTS scan times.  Install several changes into the evaluation version
in anticipation of SeaWiFS reprocessing.

### Source Changed
  * l12_parms.h, version number change.
  * getaerosol.f, if NIR goes negative, return success but set aerosols
  to a small value.
  * atmocor2.f, don't automatically fail on negative conditions. Update to
  vanishing Siegel parameters.
  * get_tlw_nir.f, new absorption coefficients for bio-optical model.
  * l1_struc.h, dropped aerSelect field (defunct).
  * filter.c, modified to allow user specification of minimum fill.
  * filter.h, prototype mods for passing of minfill.
  * msl12_input.c, recording of filter minfill values in parameter meta-data.
  * get_f0.c, new function to return F0 at a specific wavelength (no band averaging).
  * get_chl.c, eval code to use nominal F0 when nLw is out-of-band corrected.
  * gs97.c, eval code to use nominal F0 when nLw is out-of-band corrected.
  * get_rrs.c, eval code to use nominal F0 when nLw is out-of-band corrected.
  * l1_octs_hdf.c, adjust msec, year, and day to handle day roll-over like seawifs.
  * l1a_seawifs.c, cal_mod needs to be declared static for Linux.
  * MSl1brsgen.c, clean-up.
  * Makefile.*, add get_f0.o module to MSl12.

### Data Added
  * common/f0_table.dat, new Neckel and Labs solar irradiance on 1-nm interval.
  * octs/octs_filter.dat, add new minimum fill field.
  * seawifs/seawifs_filter.dat, eval filter file for eps smoothing, no min fill.
  * seawifs/cal/SEAWIFS_SENSOR_CAL.TBL-200102, latest eval cal table for seawifs.
  * osmi/*, update to latest calibration and tables.


## 2.8.6 - 2001-10-04

Fix OCTS browse metadata.

### Source Changed
  * MSl12.c, limit eline input option to nscan
  * MSl1brsgen.c, fix for scene center coordinates. Rewrite to compute most
  meta-data rather than reading from input L1 hdf file, as OCTS L1A defines
  start, center, and end coordinates differently.
  * filehandle.h, eliminate scene start, center, end time fields.
  * filehandle_init.c, eliminated scene start, center, end time fields.
  * l1_octs_hdf.c, don't read scene boundary times.
  * l1a_seawifs.c, enable out-of-sequence reading and get orbit-node meta-data
  for MSl1brsgen.
  * msl12_input.c, set default for proc_land to 0 (off).

## 2.8.6 - 2001-09-18

Correct OSMI viewing/solar geometry.

### Source Changed
  * l12_parms.h, version number change.
  * l1a_osmi.c, interpolate s/c position as function of scan position.
  * intpos.f, new function to interpolate position and velocity vectors.
  * msl12_input.c, default OSMI caltable location to $OSMI_CAL_PATH
  * msl1b_input.c, default OSMI caltable location to $OSMI_CAL_PATH


## 2.8.5 - 2001-10-04

### Source Changed
  * MSl1brsgen.c, fix for scene center coordinates. Rewrite to compute most
  meta-data rather than reading from input L1 hdf file, as OCTS L1A defines
  start, center, and end coordinates differently.


## 2.8.5 - 2001-07-30

Fix some meta-data problems for OCTS, improve flagging for OCTS, and
add some new filtering capabilities, including the diamond kernel and
the interquartile mean.

### Source Changed
  * l12_parms.h, version number change.
  * setflags.c, don't set HILT for missing OCTS bands.
  * l2_hdf_generic.c, bogus node start/end for OCTS.
  * l1_octs_hdf.c, fixed month index in time format conversion.
  * filter.h, changed prototypes for filter functions.
  * filter.c, added alternative (diamond) kernel capability.
  * MSl1brsgen.c, fix eastern, western, etc. meta data.
  * filehandle.h, clean-up.
  * filehandle_init.c, clean-up.
  * get_k490.c, updated coefficients for 490/565 band pair from S. Bailey.
  * get_dem_height.f, bug fix for linux.
  * get_chl.c, updated OCTS OC4 coeffs, then reverted back.

  * [liboctscal] occal.c, apply calibration to saturated counts, rather than returning zero.
  * [liboctscal] l1a_occal.c, apply calibration to saturated counts, rather than returning zero.

## 2.8.4 - 2001-07-16

Add flagging for filter failure conditions.

### Source Changed
  * l12_parms.h, version number change.
  * l2_flags.h, added FILTER flag.
  * l1_struc.h, added filter flag.
  * alloc_l1.c, add space for filter flag.
  * filter.h, set FILTER flag when pixels are failed for insufficient window fill.
  * convl12.c, copy filter flag in l2_flags array.
  * l1_mask.c, set filter flag as mask.


## 2.8.3 - 2001-07-11

Add additional OCTS meta-data on DAAC request.

### Source Changed
  * l12_parms.h, version number change.
  * MSl12.c, copy additional meta-data from filehandle.
  * MSl1brsgen.c, write additional meta-data.
  * filehandle.h, add storage for additional meta-data.
  * filehandle_init.c. initialize new mate-data fields.
  * l1_octs_hdf.c, read additioanl meta-data into filehandle.
  * l2_hdf_generic.c, write additional meta-data from filehandle.

## 2.8.2 - 2001-07-05

Update to available filtering algorithms. Updates for OCTS, land, PAR
processing.

### Source Changed
  * l12_parms.h, version number change.
  * filter.h, added new epsilon targeted smoothing.
  * filter.c, added new epsilon targeted smoothing. Changed Lt-Lr median
  filter to just modify Lt, leaving Lr alone.
  * calc_par.f, raise solar zenith limit to 90-deg (it will still be
  limited by the run-time solz limit).
  * l2_hdf_generic, added orbit_nod_lon for OCTS.
  * MSl1brsgen.c, added various meta-data and changed lon/lat scaling.
  * filehandle.h, added storage for orbit_nod_lon.
  * l1_octs_hdf.c, read orbit_nod_lon.
  * l1_mos_hdf.c, fix bilinear interpolation.
  * get_ndvi.c, added update to EVI algorithm from Jacques Descloitres.

  * [liboctscal] l1a_occal.f, use gain from following scan (it appears to be 1-frame late).

## 2.8.1 - 2001-05-01

Compilation of exponential cal table support, as well as other "evaluation"
algorithm changes, are now controlled by defining EVAL=1 in the compiler
commandline.  Set EVAL=0 for standard versions.  The old SWCALTBLFMT
variable has been superceded.

### Source Changed
  * l12_parms.h, version number change.
  * MSl12.c, informational message change.
  * atmocor_init.f, fix taua_per_band for less than 8 bands.
  * atmocor_set.f, fix taua_per_band for less than 8 bands.
  * l1a_osmi.c, clean-up.
  * rhoa_from_taua.c, fix for sensors with less than 8 bands.
  * atmocor2.f, added vanishing Siegel EVAL code.
  * l1a_seawifs.c, look for EVAL=1 to compile exponential cal table support.
  * l2_hdf_generic.c,look for EVAL=1 to compile exponential cal table support.
  * l1_hdf_generic_write.c,look for EVAL=1 to compile exp cal table support.
  * l12_seawifs.c, look for EVAL=1 to compile exponential cal table support.
  * get_dem_height.f, added test to fail if DEM file not found.

  * [libswfl1a] calibrate_l1a.c, look for EVAL=1 to compile exponential cal table support.

  * [libswfcal] get_cal.c, look for EVAL=1 to compile exponential cal table support.
  * [libswfcal] get_cal_misc.c, look for EVAL=1 to compile exponential cal table support.

  * [libosmical] calibrate_l1a_osmi.c, separate electronic gain from calibration table gain.

  * get_cal.h, look for EVAL=1 to compile exponential cal table support.
  * getcal_proto.h, look for EVAL=1 to compile exponential cal table support.
  * cal_l1a_osmi.h, added definition for number of quadrants.


## 2.8 - 2001-03-29

Improvements to handling of inverse (calibration) mode. New capability to
read a per-pixel aerosol model specification file. Updates for PAR.  Fixes
for OCTS, DEM file handling.

### Source Changed
  * l12_parms.h, version number change. Add FORWARD, INVERSE defines. Give
  alternate version number when building eval version.
  * input_struc.h, add mode, aerfile fields.
  * msl12_input.c, set mode field (FORWARD or INVERSE processing). Handle aerfile
  input option. Moved default locations for climatology and mask files to
  data/common directory.
  * MSl12.c, use mode setting from input structure. Allocate and read aerfile
  if supplied.
  * convl12.c, pass mode to atmocor2().  Accept aerec and call atmocor_set() if
  aerec is active.
  * atmocor2.f, don't exit on dark pixel when in INVERSE calibration mode
  * l12_proto.h, add mod parameter to atmocor2(). Add aerec to convl12(). Add
  new atmocor_set() function, aer i/o and allocation functions.
  * calc_par.f, PAR algorithm changes from R. Frouin.
  * get_par.c, added water-vapor in call to calc_par().
  * sensor_cmn.fin, add model pair and ratio to common block.
  * atmocor_init.f, initialize model pair and ratio in common block.
  * get_aerosol.f, new function to compute aerosol from fixed model pair and
  ratio.
  * rhoa_from_taua.c, need to compute angstrom for every pixel now, to support
  aerfile option.
  * aer_struc.h, new data structure to hold one record of aerosol correction
  information.
  * alloc_aer.c, allocates one aer_struc record.
  * aer_io.c, functions to open, read, and close aerfile.
  * atmocor_set.f, function to set aerosol control paramters in common block.
  * l1_generic.c, fixed memory allocation bug when reading hdf and writing L1B.
  * l1_octs_hdf.c, increased maximum lines for array allocation.
  * get_dem_height.f, fixes for dateline and large sensor zenith conditions.

### Data Added
  * mos/rayleigh/mos_pol_ray_B.dat, corrected band 6 from 650 to 685 nm.
  * landmask.dat, moved to "common" subdirectory
  * watermask.dat, moved to "common" subdirectory
  * S19461993_COADS_GEOS1.MET_general, moved to "common" subdirectory
  * S19461993_COADS_GEOS1.MET_noon, moved to "common" subdirectory
  * S19891991_TOMS.OZONE, moved to "common" subdirectory



## 2.7.3

The f/Q tables were updated with new coefficients from Andre Morel.
Support was added for new exponential calibration table format. These
changes only get compiled if EXPCALTBL is defined.

### Source Changed
  * l12_parms.h, version number change.
  * get_foq.f, update to read new f/Q tables from Morel.
  * l1a_seawifs.c, support for new cal table format.
  * l12_seawifs.c, support for new cal table format.
  * l1_hdf_generic_write.c, support for new cal table format.
  * l2_hdf_generic.c, support for new cal table format.
  * get_chl.c, set OCTS and POLDER default chl to oc4v4_octs.
  * l1_octs_hdf.c, fix for solz > 90 cases.
  * msl12_input.c, default values stated in usage message.

  * [libswfl1a] calibrate_l1a.c, support for new cal table format.

  * [libswfcal] get_cal.c, support for new cal table format.
  * [libswfcal] get_cal_misc.c, support for new cal table format.

  * swfinc/cal_l1a.h, added <math.h> for new caltable support.
  * swfinc/get_cal.h, added definitions for new caltable.
  * swfinc/getcal_proto.h, added prototypes for new caltable.
  * swfinc/l1a.h, removed old code.
  * swfinc/l1a_proto.h, removed old code.

### Data Added
  * foq/fq.dat, updated f/Q table.
  * foq/q0.dat, updated Q0 table.

## 2.7.2

New capability to read and process OSMI data.

### Source Changed
  * l12_parms.h, version number change.
  * MSl1bgen.c, enable writing of generic L1B files with less
  than 8 bands.
  * diffuse_transmittance.f, add support for OSMI tables.
  * filehandle.h, add identifier for OSMI level-1a.
  * getformat.c, recognize OSMI L1A and L1B.
  * getrayleigh.f, add support for OSMI tables.
  * get_chl.c, added oc4v4 variant for OCTS.
  * l1_generic.c, add OSMI output option.
  * l1_hdf_generic_read.c, fix to handle generic L1B files with less
  than 8 bands.
  * l1_hdf_generic_write.c, fix meta-data for OSMI ascending orbit.
  * l1_io.c, add OSMI i/o support.
  * l1_octs_hdf.c, extend max lines for OCTS GAC, improve sample table
  usage.
  * l1a_osmi.c, new OSMI input routine.
  * l1a_seawifs.c, clean-up.
  * msl12_input.c, recognize OSMI as valid sensor.
  * msl1b_input.c, recognize OSMI as valid sensor.
  * rdsensor_wl.c, locate OSMI tables.
  * rdsensortab.c, locate OSMI tables.
  * setanc.c, fixed typo in informational message.
  * setflags.c, set OCTS saturation flag whenever Lt <= 0.
  * MSl1bgen.mk, link with l1a_osmi.o and libosmical.a, include osmiinc.
  * MSl1bingen.mk, link with l1a_osmi.o and libosmical.a, include osmiinc.
  * MSl1brsgen.mk, link with l1a_osmi.o and libosmical.a, include osmiinc.
  * MSl1info.mk, link with l1a_osmi.o and libosmical.a, include osmiinc.
  * MSll2snpx.mk, link with l1a_osmi.o and libosmical.a, include osmiinc.
  * hdf2pgm.mk, link with l1a_osmi.o and libosmical.a, include osmiinc.
  * Makefile, link with l1a_osmi.o and libosmical.a, include osmiinc.

  * [libosmical] calday.c
  * [libosmical] calibrate_l1a_osmi.c
  * [libosmical] ctogd.c
  * [libosmical] ctotc.c
  * [libosmical] eanom.c
  * [libosmical] earth.h
  * [libosmical] get_cal_misc_osmi.c
  * [libosmical] get_cal_osmi.c
  * [libosmical] gmha.c
  * [libosmical] julian.c
  * [libosmical] locate.c
  * [libosmical] obliq.c
  * [libosmical] sunpos.c
  * [libosmical] tconv.c
  * [libosmical] time-utils.c

  * [osmiinc] InstStatData.h
  * [osmiinc] cal_l1a_osmi.h
  * [osmiinc] call1a_proto_osmi.h
  * [osmiinc] earth.h
  * [osmiinc] get_cal_osmi.h
  * [osmiinc] getcal_proto_osmi.h
  * [osmiinc] l1a.h
  * [osmiinc] l1a_proto.h
  * [osmiinc] lunsol.h
  * [osmiinc] orbit.h

  * [liboctscal] l1a_occal, improve handling of GAC sampling table.

## 2.7.1 - 2000-10-31

Support for OCTS L1A GAC. Fix for stray light.

### Source Changed
- gs97.c, switch to Levenberg-Marquardt optimizer.
- filehandle.h, new OCTSL1A identifier.
- getformat.c, identify files of type "OCTS Level-1A GAC Data"
- l1_io.c, added OCTSL1A support.
- l1_octs_hdf.c, generalized to read OCTS L1A GAC as well as L1B LAC.

** changes some L1B radiances **
- [libswfl1a] lac_st.c, fix for rare mis-flagging of stray light.
- [libswfl1a] gac_st.c, fix for rare mis-flagging of stray light.

- [liboctscal] l1a_occal.f, new OCTS L1A system calibration routine.
- [liboctscal] occal.f, old OCTS L0 system calibration routine.
- [liboctscal] ocloadcal.f, old OCTS calibration coefficients.

## 2.7 - 2000-10-16

Refinements to land processing.

### Source Changed
  * l12_parms.h, version number change.
  * loadl1.c, don't compute surface reflectance over unprocessed pixels.
  Init pre-computed atmospheric quantities when not processing pixels.
  * setflags.c, rough improvement to cloud flag over land, use threshold
  on surface reflectance at 412. Threshold defined as albedo/10.
  * l1a_seawifs.c, bug fix to correct the handling of scan modulation
  corrections for pre-extracted SeaWiFS L1A subscene files.
  * l12_seawifs.c, corrected tilt_flags meta-data when subsetting.
  * setanc.c, changed correction of pressure for altitude from 100mb
  per 1000 meters to exponential.
  * levmarq.f, new Levenburg-Marquardt optimization routine.
  * amoeba.c, efficiency enhancements.
  * amoeba.h, efficiency enhancements.
  * gs97.c, added new function for Levenburg-Marquardt fitting.
  * water_vapor.c, new function to compute water_vapor transmittance
  for SeaWiFS. Only used in first estimate of surface reflectance.
  * l1_struc.h, added t_h2o field, water vapor transmittance.
  * l2_struc.h, added t_h2o field.
  * alloc_l1.c, allocate space for t_h2o.
  * convl12.c, transfer t_h2o.
  * atmocor1.f, pass water-vapor concentration, compute t_h2o.
  * atmocor1_land.c, compute t_h2o, set t_o2 to 1.0 everywhere.
  * loadl1.c, pass t_h2o, wv to atmocor1().
  * l12_proto.h, modified atmocor1() proto. New water_vapor() proto.
  * get_rhos.c, don't correct for spherical albedo over oceans. Do correct for
  water vapor using pre-computed transmittance.
  * l2prod.h, added t_h2o product.
  * get_l2prod_index, added t_h2o product.
  * l2_hdf_generic.c, added t_h2o product.

  * [libswfl1a] calibrate_l1a.c, pass nsamp, nsta, ninc to get_cal.
  * [libswfl1a] call1a_proto.h, new calibrate_l1a proto.

  * [libswfcal] get_cal.c, pass npix, nsta, ninc to setup_scanmod.
  * [libswfcal] get_cal_misc.c, generalize scan modulation table to use npix, nsta, ninc
  rather than data-type (GAC, LAC) for sub-setting, sub-sampling cases.
  * [libswfcal] getcal_proto.h, new get_cal proto.


## 2.6.1 28 September 2000

Linux compatibility.

### Source Changed
 * *.f,*.c, *.h, *.fin, various changes for linux compatibility.

## 2.6 - 2000-09-13

Fix non-seawifs out-of-band corrections. Mask glint in PAR.

### Source Changed
  * l12_parms.h, version change.
  * get_par.c, mask glint.
  * get_dem_height.f, fix to improve scan-edge interpolation.
  * calc_aero_refl.f, ignore out-of-band option when not seawifs.
  * getaermod.f, ignore out-of-band option when not seawifs.
  * rhoa_from_taua.f, ignore out-of-band option when not seawifs.
  * msl12_input.c, improved usage message.

## 2.5 - 2000-08-22

Support for 1997 Garver-Siegel bio-optical products. Updates for PAR, aerosol index.

### Source Changed
  * l12_parms.h, version change.
  * smi_climatology.c, new general reader for SMI-like climatological files.
  * smi_climatology.h, new header for smi_climatology.c.
  * alloc_2d.c, new 2-D array allocation utility.
  * alloc_2d.h, new header for alloc_2d.c.
  * get_par.c, modified to use climatology for aerosol properties, and pass
  bandpass-dependent constants to calc_par().
  * calc_par.f, modified to use angstrom(510) rather than angstrom(670). Also
  modified to accept parameters for wavelength, F0, kO3, and Tau_r.
  * loadl1.c, switched sst reader to use generic smi_climatology() function.
  * msl12_input.c, handle new input filenames for aerosol climatologies.
  * input_struc.h, new input filenames for aerosol climatologies.
  * sst.c, replaced with smi_climatology.c.
  * get_rrs.c, fixed error with day-of-year correction.
  * l2prod.h, added gs97 products.
  * get_l2prod_index.c, added gs97 products.
  * l2_hdf_generic.c, added gs97 products.
  * l1_struc.h, added tracking of current scan number.
  * l2_struc.h, added tracking of current scan number.
  * l1_io.c, set current scan number in readl1().
  * convl12.c, transfer current scan number to l2rec.
  * cpl1rec.c, transfer current scan number to duplicate l1rec.
  * l12_proto.h, added proto for get_gs97().
  * amoeba.c, new model optimization function.
  * amoeba.h, new header for amoeba().
  * gs97.c, new Garver-Seigel 1997 model and calling function, get_gs97().
  * nlw_outband.f, fixed bug in which POLDER was calling OCTS function.
  * get_par.c, convert to Einsteins/Day/m^2.
  * get_l2prod_index.c, change units meta-data for par.
  * calc_par.f, fixed error in day-of-year indexing.2.5 22 August 2000
  * [libai] residue.f, add screen for river sediment.
  * [libai] aero_indx.f, eliminate shallow-water mask.

### Data Added
  * common/sst_climatology.hdf, modified to use standardized SDS names.
  * common/alpha510_climatology.hdf, new angstrom climatology.
  * common/taua865_climatology.hdf, new aerosol optical thickness climatology.


## 2.4 - 2000-07-21

Added capability to specify aerosol optical thickness per band and
from that derive aerosol model, reflectance per pixel. New tricho
product and supporting SST climatology i/o. Fixed minor error in
angstrom calculation. Fixed DEM file handling.

### Source Changed
  * l12_parms.h, version change. New USERTAUAS aerosol option identifier.
  * rhoa_from_taua.f, new routine to compute aerosol reflectance based
  on user input aerosol optical thicknesses per band.
  * input_struc.h, added taua_per_band input parameter, sstfile parameter.
  * msl12_input.c, added taua_per_band input parameter, sstfile parameter.
  * sensor_cmn.fin, added storage for taua_per_band.
  * atmocor_init.f, store taua_per_band in common.
  * loadl1.c, pass taua_per_band to atmocor_init(), set SST per pixel. Fixed
  error in terrain-height trigger. Fixed error status returns.
  * l12_proto.h, extend proto for atmocor_init().
  * getaerosol.f, added rhoa_from_taua option.
  * get_angstrom.c, fixed error in which the order of the C70 and C90 models
  were swapped.
  * get_dem_height.f, updated to improve negative height calculation.
  * sst.c, new functions sst_init() and get_sst() to load sea surface
  temperature from climatology.
  * l1_struc.h, added sst pointer.
  * l2_struc.h, added sst pointer.
  * alloc_l1.c, added space for sst.
  * convl12.c, transfer sst pointer from l1 to l2rec.
  * l2_prod.h, added sst identifier.
  * get_l2prod_index.c, added support for sst. Changed epsilon to return
  as float (no impact to eps_78).
  * l2_hdf_generic.c, added support for sst
  * MSl12.c, informational messages enhanced.
  * get_tricho.c, limit calculation to low wind and warm water.
  * l12_proto.h, added sst protos.
  * ice_mask.c, fixed month indexing. Improved error checking.

### Data Added
  * common/sst_climatology.hdf, new sea surface temperature climatology file.
  * seawifs/aerosol/seawifs_aerosol.dat, new aerosol table which consolidates
  the content of the seawifs_aerosol_*.dat files, including angstrom coeffs.


## 2.3 - 2000-07-13
Set navigation failure flag. Mask glint in cloud optical thickness product.

### Source Changed
  * l12_parms.h, version change. Added new OUTPUT_QNLW format for
  calibration mode processing.
  * convl21.c, added option to output quasi-nLw.
  * convl12.c, set NAVFAIL flag.
  * filehandle.h, Added new OUTPUT_QNLW format for calibration mode
  processing.
  * l2_flags.h, added NAVFAIL flag.
  * get_calcite.c, set masked values to 0.0 instead of -0.01. Disable glint
  correction.
  * get_tauc.c, mask glint.
  * l12_seawifs.c, fixed problem with seawifs meta data when input
  parameter for start line is set to negative value.
  * [libai] aero_indx.f, disabled solar zenith angle switch at 40-deg on ppo2.


## 2.2 - 2000-06-30
### Source Changed
New calcite, cloud optical thickness, and remote sensing reflectance
products. New ice mask i/o for PAR masking.

### Source Changed
  * l12_parms.h, version change.
  * l2prod.h, added calcite, tauc, Rrs.
  * l12_proto.h, new calcite, cocolith, tauc, rrs  protos.
  * l2_hdf_generic.c, added call to get_calcite(), get_tauc(), get_rrs().
  * get_l2prod_index.c, added resolution of calcite, tauc, Rrs options.
  * get_calcite.c, new wrapper for calcite.f.
  * coccolith.f, new coccolithofore backscatter algorithm from H. Gordon.
  * atmocor2.f, slight change to improve retrieval of Lw(765) when Siegel 6-8 algorithm is used.
  * ice_mask.c, new ice mask function.
  * msl12_input.c, added input control for ice mask file. Set new atmospheric correction switch.
  * input_struc.h, added icefile to input structure. Added atmocor switch.
  * loadl1.c, init and call ice_mask().
  * l1_struc.h, add pointer for ice mask.
  * alloc_l1.c, make space for ice mask.
  * convl12.c, transfer ice mask to l2_flags. Enable atmocor toggle.
  * l2_flags.h, add ICE identifier.
  * get_par.c, don't compute par over land or ice.
  * get_tauc.c, new function to compute cloud optical thickness, seawifs GAC only, 865 only.
  * get_rrs.c, new function to compute remote sensing reflectance per band.
  * MSl12.c, added some informational messages.
  * [libcloud.a] new library of routines for cloud optical thickness

### Data Added
  * common/ice_mask.hdf, new sea ice climatology file.
  * seawifs/cloud/*, new data files for seawifs cloud optical thickness.

## 2.1 - 2000-06-21
Major changes to support land processing, including introduction of
the digital_elevation_map. New wind-dependent Rayleigh tables for
POLDER, MOS, OCTS.

### Source Changed 
  * alloc_l1.c, added storage for terrain height, t_o2, t_sol, t_sen, rhos.
  * alloc_l2.c, removed storage for t_o2, t_sol, t_sen, rhos.  These are
  now pointers to the l1 record.
  * atmocor1.f, moved computation of O2 transmittance from atmocor2() to
  atmocor1() so it can be used in surface reflectance computations without
  running through full atm. correction. Also added simple calc. of
  Rayleigh diffuse transmittance for same purpose.
  * atmocor1_land.c, new function to process to rhos over land, including
  new Rayleigh function for land. Called by loadl1(). Replaces land_refl.c.
  * atmocor2.f, call nlw_outband() for all sensors. Moved coputation of O2
  transmittance from here to atmocor1(), consistent with land.
  * convl12.c, transfer terrain height. Only do land processing over land.
  Eliminate call to land_refl() (done in loadl1()).
  * dem_s.fin, new DEM file structure include.
  * filter.c, dilation filter now sets dilated flag. Also added masking of
  the queue ROI when a dilation flag is detected. Flags are now specified
  by standard L2_flag numbers. All smoothing filters were changed to set
  the mask bit (rather than zeroing the data) when the ROI is less then
  half filled. Added test filter for epsilon smoothing (incomplete).
  * filter.h, removed old flag definitions and names.
  * get_dem_height.f, new function to read and apply terrain elevation data,
  adjust lon/lat and view angles accordingly, return height for pressure
  correction.
  * get_l2prod_index.c, added support for rhot and height products. Changed
  scaling on rhos. Added units of "dimensionless" for all products that
  had no units specified.
  * getl1rec.c, don't need to reset masking after filtering. Removed filter
  control parameter, which is now included in input control structure.
  * getrayleigh.f, call wind-dependent Rayleigh tables for MOS, OCTS, and
  POLDER.
  * get_rhos.c, new function to compute surface reflectance.
  * get_toa_refl.c, new source to compute TOA reflectance, rhot.
  * input_struc.h, added elevation filename field. (** actually added on
  4/25). Added filter control structure.
  * l12_parms.h, updated version number.
  * l12_proto.h, added get_toa_refl() prototype. Fixed getrayleigh proto
  for windspeed (not called by C, so no problem). Fixed get_es proto for
  readability.  Added getheight() proto. Eliminated land_refl() proto.
  Added atmocor1_land(). Removed fctl from getl1rec() proto. Changed
  atmocor1 proto.
  * l1_io.c, removed setting of defunct field, l1rec->input.
  * l1_mos_hdf.c, fixed bug in last-scan interpolation for lon/lat.
  * l1_struc.h, added fields height, rhos, t_o2, t_sol, and t_sen.
  The associated transmittance fields in the l2 rec become pointers to
  the l1 rec.
  * l1q_struc.h, added include statement.
  * l2_hdf_generic.c, added support for rhot and height products.
  * l2_struc.h, added terrain height field. (** actually added on 4/19).
  * l2prod.h, added rhot and height product identifiers.
  * loadl1.c, compute land/bath flags. Separate land and ocean streams.
  Call atmocor1_land() instead of caling land_refl() in convl12(). Call
  get_dem_height() to adjust lon/lat for terrain elevation. Added call
  to get_rhos(), and passing of transmittances to atmocor1().
  * msl12_input.c, added elev file specification option. (** actually added
  on 4/25). Changed outband_opt default to 2 for OCTS, POLDER. Added
  initialization of filter control filed and reading of filter file. Added
  echo of filters used in input meta-data.
  * nlw_outband.f, replaces nlw_outband_seawifs.f.  Includes OCTS and
  POLDER corrections.
  * rho_a_sub_quad.f, apply out-of-band wv correction to seawifs only.
  * setanc.c, read terrain height and adjust pressure, if over land and
  land processing is requested. (*Moved dem reading to loadl1()).
  * setflags.c, moved land/bath flags up to loadl1(), so we can use to decide
  if we need to read terrain height in setanc().
  * MSl12.c, removed filter init and reading, now done in msl12_input().
  * MSl1brsgen.c, set pointer to input struc in l1file.
  * MSl1brsgen.mk, added atmocor1_land.c, filter.c.
  * Makefile, added get_toa_refl.c, atmocor1_land.c.
  * [include] swfinc/nav_cnst.fin, new Fortan include for get_dem_height.f.

### Data Added
  * common/digital_elevation_map.hdf, new terrain height map.
  * mos/rayleigh/mos_pol_ray_B.dat, new wind-dependent Rayleigh table.
  * octs/rayleigh/octs/pol_ray_B.dat, new wind-dependent Rayleigh table.
  * polder/rayleigh/polder_pol_ray_B.dat, new wind-dependent Rayleigh table.
  * polder/polder_table.dat, new polder calibration.

## 2.0.2 - 2000-10-05

Fix for scan-modulation correction when processing subscenes.

### Source Changed
  * l12_parms.h, updated version number.
  * l1a_seawifs.c, bug fix to correct the handling of scan modulation corrections for pre-extracted SeaWiFS L1A subscene files.
  * l12_seawifs.c, corrected tilt_flags meta-data when subsetting.
  * [libswfl1a] calibrate_l1a.c, pass nsamp, nsta, ninc to get_cal.
  * [libswfl1a] call1a_proto.h, new calibrate_l1a proto.
  * [libswfcal] get_cal.c, pass npix, nsta, ninc to setup_scanmod.
  * [libswfcal] get_cal_misc.c, generalize scan modulation table to use npix, nsta, ninc rather than data-type (GAC, LAC) for sub-setting, sub-sampling cases.
  * [libswfcal] getcal_proto.h, new get_cal proto.


## 2.0.1 - 2000-09-11
### Source Changed
  * l12_parms.h, updated version number.
  * get_angstrom.f, corrected mis-ordering of C70 and C90 models.
  * [libai] residue.f,
  * [libai] aero_indx.f,
  * [libai] *.f, assorted changes for linux compatibility

## 2.0 - 2000-04-28
### Source Changed
  * convl1.c, set ATMWARN on NEGLW if band is greater than 2.  This will
  make ATMWARN more useful in L3 binning to mask bad K490.
  * setflags.c, don't compute glint for land pixel.
  * [libai] residue.f, updated to improve behavior of aerosol index over high chl.

### Data Added
  * $MSL12_DATA/seawifs/aerosol/abs_aerosol_seawifs.dat updated.

## 2.0 2000-04-26
### Source Changed
  * convl1.c, don't skip flag tests after ATMWARN.
  * atmocor2.f, don't include band 1 in ATMFAIL test for negative Lt-Lr.
  * setflags.c, moved dark pixel test outside land/glint/absaer test.

## 2.0 - 2000-04-25
### Source Changed
  * input_struc.h, extended replacement flag storage length.
  * convl12.c, fixed bug in tricho flag handling of bathymetry condition.
  * setanc.c, added informational messages.
  * l1a_seawifs.c, added informational messages re: seawifs cal table.
  * msl12_input.c, added informational messages re: seawifs cal table.

## 2.0 - 2000-04-19
### Source Changed
  * l2_flags.h, define new flag to indicate cloud-free ocean.
  * convl12.c, set flag to indicate cloud-free ocean, copy mask pointer and
  set mask array on atmospheric correction failure.
  * get_k490.c, set explicit max/min and failure conditions.
  * l2_struc.h, added pointer to L1 mask array.

## 2.0 - 2000-04-05
### Source Changed
  * msl12_input.c, set default solzmax to 75. Added diagnostic messages.
  * get_chl.c, added minimum cut-off on band ratio for OC4 algorithm.
  * getformat.c, added diagnostic messages.
  * l1_io.c, reset l1rec sensor ID on each read (required for generic L1B).
  * MSl12.c, moved setting of l1rec nbands and bindx to readl1().
  * cpl1rec.c, copy bindx field (for POLDER).

## 2.0 - 2000-03-31
### Source Changed
  * getrayleigh.f, updated to higher resolution wind tables.
  * getaerosol.f, require La(7)/La(8) > 0.1 before attempting to compute aerosols.
  * get_chl.c, updated OC4 algorithm to v4 coefs.

### Data Added
  * seawifs/rayleigh/seawifs_pol_ray_B.dat, updated with higher wind resolution.


## 2.0 - 2000-03-27
### Source Changed
  * get_l2prod_index.c, add support for tau_nnn, fixed rank problem for fsol.  Also changed default chl and pig to OC4.
  * l2_hdf_generic.c, fixed rank dimension problem for fsol product. Added global attribute for fsol.
  * l1_hdf_generic_read.c, modified to read standard seawifs-like L1B in
  addition to OCTS L1B form OCl1bgen.
  * MSl12.c, enhanced to include std seawifs meta-data when processing L1B.
  * hdf_utils.c, added SetF64GA() function.
  * get_chl.c, switch default chl algorithm to chl_oc4.
  * getformat.c, added informational messages.
  * loadl1.c, added some comments.
  * Makefile, *.mk, added utilities for new l1_hdf_generic_read().

## 2.0 - 2000-03-16
### Source Changed
  * msl12_input.c, updated to default aer_opt=-3, outband_opt=2 for seawifs. Also glint_opt=1 for all. Added rflag option for DAAC replacement flag control. Changed to recognize glint_thresh for glint.
  * input_struc.h, Added rflag.
  * l2_hdf_generic.c, Added rflag.
  * l1_octs_hdf.c, increased MAXPIX to handle fwd tilt scan width.
  * convl21.c, add correction for fsol on qnlw reconstruction.
  * l2prod.h, added support for fsol and t_o2 output.
  * get_l2prod_indx.c, added support for fsol and t_o2 output.
  * l2_hdf_generic.c, added support for fsol and t_o2 output.
  * MSl12.c, added "Processing Completed" message.
  * MSl1bgen.c, added "Processing Completed" message.
  * l1_mos_hdf.c, changed to read scan coeffs from external file.
  * getrayleigh.f, added new table read function for seawifs, to include windspeed.
  * atmocor1.f, pass windspeed to getrayleigh().

###Data Added
  * add relgain file to data/mos/cal.
  * add seawifs_pol_ray_B.dat to data/seawifs/rayleigh. This effectively
  replaces the eight files sea*pol.dat.

## 2.0 - 2000-02-29
### Source Changed
  * l2_hdf_generic.c, removed include for palette.h, added support for windangle, zwind, mwind, cloud_albedo.
  * l2_prod.h, added windangle, zwind, mwind, cloud_albedo.
  * get_l2prod_index.c, added windangle, zwind, mwind, cloud_albedo.
  * convl12.c, copy pointers for new products to l2rec.
  * alloc_l1.c, added storage for wd, cloud_albedo.
  * l1_struc.h, added pointer for wd, cloud_albedo.
  * l2_struc.h, added pointer for wd, zw, mw, cloud_albedo.
  * setanc.c, compute and store wind direction, l1rec->wd[i].
  * setflags.c, save cloud_albedo in l1rec->cloud_albedo[i].
  * get_chl.c, added chl_oc4v3. Changed oc4 option to v3.
  * earthsundist.f, moved to libswfnav().

  * MSl1bgen.c, fixed handling of line subsampling.
  * l12_seawifs.c,fixed handling of line subsampling.
  * l12_proto.h, changed proto for l1b_seawifs().
  * l1_hdf_generic_write.c, added csol_z metadata output.

## 2.0 - 2000-02-24
### Source Changed
  * msl12_input.c, moved usage message from MSl12.c, to enhance maintainability. Also updated the usage message to include all the latest options, and changed the taua input option to tau_a (for SeaDAS compatibility)
  * MSl12.c, moved usage message to msl12_input.c
  * getglint_rad.f, fix for glint over-correction.

  * MSl1bgen.c, extended l2_flags to 32-bit.
  * l1_hdf_generic_write.c, added additional meta-data and Vgroups to MSl1bgen output for SeaDAS compatibility. Extended l2_flags to 32-bit.
  * l1_hdf_generic_write.h, extended l2_flags to 32-bit in prototype.

## 2.0 - 2000-02-15
### Source Changed  
  * l1a_seawifs.c, fix for non-MSl12 programs. Check that input struct is initialized
  before trying to use it.
  * atmocor2.f, add pressure to input args. Correct local taur and adjust diffuse
  transmittances for pressure.
  * convl12.c, pass pressure to atmocor2().
  * l12_proto.h, change prototype for atmocor2().
  * land_refl.c, add pressure correction for Tau_r.

## 2.0 - 2000-02-10
  * new aerosol option, input taua with fixed aerosol model.
  * new land products: rhos_*, evi, new ndvi.
  * new chl and pigment products, chl_oc4, pig_oc4, etc.
  * fixed LAC straylight calculation to avoid error at scene edge

### Added  
  * get_smoke.c, computes Vermote smoke index.
  * land_refl.c, computes surface reflectance.
  * calc_aero_refl.f, computes aerosol reflectance from input taua and fixed model.

### Source Changed
  * Makefile, added get_smoke.o, land_refl.o, calc_aero_refl.o.
  * MSl12.c, added informational print statement.
  * alloc_l2.c, added allocation for rhos (surface reflectance).
  * atmocor2.f, changed chl test to look for chl=-1 instead of chl=0 as indicator of chl algorithm failure.  Also modified foq call to pass actual number of visible bands, for POLDER support.
  * atmocor_init.f, load new input option, tauaInput, into sensor common.
  * convl12.c, added switch to run land processing (compute surface
  reflectance).
  * get_chl.c, added several new chl and pigment algorithms, changed
  functions to return -1 on error, 0 as minimum.
  * get_foq.f, modified to use only the active visible bands (for POLDER), and to process through geometries that exceed table boundaries.
  * get_l2prod_index.c, added support for evi, smoke, assorted pigment and chl products, surface reflectance.  Eliminated CZCS_pigment.
  * get_ndvi.c, changed to use precomputed surface reflectances.
  Added evi.
  * getaerosol.f, added support for user-specified taua and fixed model.
  * input_struc.h, added input parameter taua.
  * l12_parms.h, added defines for CHL_MAX=64 and CHL_BAD=-1.
  * l12_proto.h, added prototypes for get_evi(), get_smoke(), land_refl().  Changed prototype for atmocor_init().
  * l2_hdf_generic.c, added function calls for the new output products. Fixed sds scaling test for minimum value.  Deleted some dead code.
  * l2_struc.h, added surface reflectance, rhos, to structure.
  * l2prod.h, added catalog entries for various new output products, deleted CZCS_pigment.
  * loadl1.c, pass tauaInput to atmocor_init().
  * msl12_input.c, added support for taua input parameter.

  * lac_st.c, added division by sum of weights (libswfl1a)


## 1.7 - 1999-03-11

### Source Changed
  * support for multiple, user-specifiable output files.

## 1.2 - 1998-07-22
### Source Changed
  * linear_a_b_c.f, fixed isun bug, changed wavelength variable name to
  eliminate confusion with sensor_cmn.
  * fourier_a_b_c.f, fixed isun bug, changed wavelength variable name to
  eliminate confusion with sensor_cmn.
  * getaermod.f, eliminated ctlaerosol common block ... unneccessary.
  * diffuse_transmittance.f, eliminated old aerosol models.
  * rho_a_sub_quad.f, removed extraneous save statement.

## 1.1 - 1998-06-20
### Source Changed
  * getrayleigh.f: fixed typo in rayleigh file name for OCTS.
  * earthsundist.f: updated with more accurrate formulation.
  * getwcap.f: reduced whitecap radiance by 75%.
  * diffuse_transmittance.f: add oceanic model selection.
  * fourier_a_b_c.f: new routine: reads fourier version of aerosol tables.
  * getaermod.f: added switch for reading fourier model coefs.
  * wangaer.f: added switch for reading fourier model coefs.
  * l12_parms.h: added definitions for aerosol model table selection.
  
### Data Added
  * seawifs/data/aerosol: added fourier tables and new oceanic aerosol models.
  * seawifs/data/transmittance: added oceanic aerosol models.
  * mos    /data/transmittance: added oceanic aerosol models.
  * octs   /data/transmittance: added oceanic aerosol models.

## 1.0 - 1998-06-01
### Added
  * Initial development

  
[Unreleased]: http://link-to-diff-of-1.1-and-head
[2]: http://svn101.domain.sdps/OCSSW/changeset?reponame=&new=87%40tags%2FV2.0&old=66%40tags%2FV2.0
[0]: http://svn101.domain.sdps/OCSSW/log/tags/V1.0

[bfranz]: mailto:bryan.a.franz@nasa.gov
