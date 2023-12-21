

/*
*************************************************
l1cgen
version 1.0.0
Originally created by Martin Montes on 3/27/2021
PACE-SDS team
*/***********************************************

//This software can handle the processing of l1b files derived from legacy (OCI) and non-legacy (SPEXone, HARP, HARP2) sensors
//a future l1cgen version will also read/process l2 and l1c files
//l1c runs indepdently from l2gen

//l1cgen can perform three main tasks: L1C grid creation (l1c_pflag flag=1) and binning of geophysical variables using the sbs (l1c_pflag=2) and area-weighting (l1c_pflag=3) method.
//sbs stands for saddleback binary search algorithm

//--------------------------------------------------------------------------------------------------------------------------
//Example #1
//--------------------------------------------------------------------------------------------------------------------------
//L1C grid generation for simulated OCI l1b files associated to an specific orbit (half-orbit daytime swath):

input files:
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T011500.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T012000.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T012500.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T013000.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T013500.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T014000.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T014500.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T015000.L1B.V5.nc
/accounts/mamontes/images/OCIS/sean/PACE_OCI_SIM.20190321T015500.L1B.V5.nc
 
//  run  l1cgen ifile=/accounts/mamontes/images/OCIS/sean/ocis_swtd_1.txt selyear=2019 selmon=3 selday=21 l1c_pflag=1


output:/out/PACE_OCIS.2019321T00:00:00ZL1Cgrid.nc

//this processing requires to specify the acquisition date of the images as year,month and day (selyear,selmon and selday, respectively).
//Likewise, l1cgen running directory must have a subdirectory /out for output files 



NOTE:in version 1.1.0, the orbital information is included in the l1b files. However,orbital velocity and scan time will be also distributed in future versions as a separate ephemerides file


//--------------------------------------------------------------------------------------------------------------------------
//Example #2
//--------------------------------------------------------------------------------------------------------------------------
//sbs binning radiance measurements derived from OCIS l1b files at L1C common grid resolution of 5.2 km

//  run  l1cgen ifile=/accounts/mamontes/images/OCIS/sean/ocis_swtd_1.txt selyear=2019 selmon=3 selday=21 l1c_pflag=2 selgran=2


output:/out/PACE_OCI.2019321T00:00:00ZbinLt_sbs_1.nc

// in this case, selgran=2 indicates the processing of file #2 (sequencial order) in the list ocis_swtd_1.txt previously created
//selgran can accept more than 1 granule, e.g., selgran=[1,3,5]. No arguments in selgran will process ALL granules for that swath 
//--------------------------------------------------------------------------------------------------------------------------
