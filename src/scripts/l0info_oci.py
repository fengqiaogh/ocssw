#!/usr/bin/env python3


__version__ = '2.7.0_2023-02-27'

__dtypes__ = ['','','_DARK','_SOL','_SPCA','_LIN','_LUN','_DIAG','_STAT',
    '_SPEC','','_SNAP-X','_SNAP-I','','_LUN-ST']

import numpy as np  
from l0info_utils import read_packet, read_apid, get_anc_packet_time, is_bad_packet  

def get_band_dims(apacket):

# extract aggregation information and compute data dimensions from ancillary packet

#       Arguments
#
#       Name    Type            I/O      Description
#       ----    ----            ---      -----------
#    apacket    byte(*)         I    Ancillary packet array
#    ncp    int         O    Number of CCD band pixels
#    nbb    int         O    Number of blue CCD bands
#    nrb    int         O    Number of red CCD bands
#    nsp    int         O    Number of SWIR band pixels
#    ndc    int         O    Number of CCD dark collect pixels
#    nds    int         O    Number of SWIR dark collect pixels
#    btaps    int         O    Blue CCD tap enable flags
#    rtaps    int         O    Red CCD tap enable flags
#    itable  struct         O    Spatial aggregation table        

#    orig: Frederick S. Patt, SAIC, 1 August 2018

    # Create spatial aggregation table structure
    itable = np.zeros(10, dtype={'names':('dtype', 'iagg', 'lines'),'formats':('i4','i4','i4')})

    # Extract spatial aggregation table and compute numbers of pixels
    ioff = 36
    nagg = [1,2,4,8]
    ncp = 0
    nsp = 0
    ndc = 0
    nds = 0
    for i in range(0,10):
        itable['dtype'][i] = apacket[ioff+3]%16
        itable['iagg'][i] = apacket[ioff+2]%4
        itable['lines'][i] = apacket[ioff]*256 + apacket[ioff+1]
        ioff = ioff + 4
        if ((itable['dtype'][i] > 0) and (itable['dtype'][i] <= 12) and (itable['dtype'][i] != 10)):
            if (itable['dtype'][i] == 2):
                ndc = ndc + itable['lines'][i]/nagg[itable['iagg'][i]] 
                nds = nds + itable['lines'][i]/8
            ncp = ncp + itable['lines'][i]/nagg[itable['iagg'][i]]
            nsp = nsp + itable['lines'][i]/8
    if (ncp == 0):
        ncp = 1
        nsp = 1

    if (ndc == 0):
        ndc = 1
        nds = 1
  
    ioff = ioff + 4

    # Extract spectral aggregation and compute numbers of bands
    #  Tap enable flags
    btap = apacket[ioff+2]*256 + apacket[ioff+3] 
    rtap = apacket[ioff]*256 + apacket[ioff+1]
    btaps = np.zeros(16)
    rtaps = np.zeros(16)
    
    #  Tap aggregation factors
    bagg = int.from_bytes(apacket[ioff+8:ioff+12],'big')
    ragg = int.from_bytes(apacket[ioff+4:ioff+8],'big')
    
    baggs = np.zeros(16)
    raggs = np.zeros(16)
    
    #  Compute number of bands for enabled taps
    nbb = 0
    nrb = 0
    ken = 1
    kag = 3
    lag = 1
    
    for i in range(15,-1,-1):
        btaps[i] = np.bitwise_and(btap, ken)/ken
        if (btaps[i]):
            baggs[i] = nagg[int(np.bitwise_and(bagg, kag)/lag)]
            nbb = nbb + 32/baggs[i]
        
        rtaps[i] = np.bitwise_and(rtap, ken)/ken
        if (rtaps[i]):
            raggs[i] = nagg[int(np.bitwise_and(ragg, kag)/lag)]
            nrb = nrb + 32/raggs[i]
        
        ken = ken*2
        kag = kag*4
        lag = lag*4
        
    return ncp,nbb,nrb,nsp,ndc,nds,btaps,rtaps,itable,baggs,raggs

def anc_compare(apacket0,apacket):

# function to compare spatial and spectral data collection parameters
#  from two OCI ancillary data packets

#        Arguments
# 
#        Name        Type    I/O      Description
#        ----        ----    ---      -----------
#     apacket0(*)    byte     I    Ancillary packet from first spin
#     apacket(*)    byte     I    Ancillary packet from next spin

#  Returns 1 (TRUE) if packets agree, 0 if not
    
    anc_compare = 1
    strmsg = ""
    # Compare spatial data collection fields
    ioff = 36
    ilen = 40
    if apacket[ioff:ioff+ilen] != apacket0[ioff:ioff+ilen]:
        c_time = get_anc_packet_time(apacket)
        strmsg += "\nSpatial table change at %s" % c_time.strftime('%Y-%m-%dT%H:%M:%S.%f')
        anc_compare = 0

    # Compare spectral data collection fields
    joff = 80
    jlen = 12
    if apacket[joff:joff+jlen] != apacket0[joff:joff+jlen]:
        c_time = get_anc_packet_time(apacket)
        strmsg += "\nSpectral table change at %s" % c_time.strftime('%Y-%m-%dT%H:%M:%S.%f')
        anc_compare = 0

    return anc_compare, strmsg

def is_bad_apid(apid,fh):
    if apid == -1:
        return False
    if (apid < 550) or (apid > 749):
        print("Invalid packet header at byte %d."%fh.tell())
        return True
    else:
        return False

def get_oci_data_type(fpacket):
    # To determine the OCI data type from and ancillary packet
    # Returns: dtype (1 - 12) if valid ancillary packet, -1 otherwise

    # Check APID
    apid = read_apid(fpacket)
    if (apid != 636): return -1

    dtypes = np.zeros(10)-1
    ioff = 36
    for i in range(0,10):
        dtypes[i] = fpacket[ioff+3] % 16
        ioff = ioff + 4
    kd = np.argwhere((dtypes != 2) & (dtypes != 10)) # Exclude dark and no processing types
    dtype = np.max(dtypes[kd])
    
    return int(dtype)


def l0info_oci(args, fh, output, bDSB):
    # procedure to get start and end times from Level-0 packet files for OCI
    print("Running l0info_oci (version: %s) \n" % __version__)
    
    dtype = -1
    str_stime = ''
    str_etime = ''
    
    status = 0
    if args.verbose:
        print("Reading OCI science data file.")
        
    bSPW = False
    lpoint = 0
    if bDSB:
        # chead = fh.read(64)
        lpoint = 64
        bSPW = True
    fh.seek(lpoint)
    
    # Get first ancillary packet
    apid = 0
    packetLength = 8
    while (apid != 636) and (packetLength > 7):
        try:
            fpacket, packetLength = read_packet(fh, bSPW)
            apid = read_apid(fpacket)
            if is_bad_apid(apid,fh):  return 104
        except:            
            return 103
    if is_bad_packet(packetLength,fh):  return 104
    
    # If no ancillary packets, rewind and get times from HKT packets
    if fpacket is None: 
        print('No science packets found in file.')
        status = 110
        
        fh.seek(lpoint)
        apid = 0
        packetLength = 8
        ccsec = 0
        # while ((apid < 550) or (apid > 750) or (apid == 558)) and (packetLength > 8):
        while (ccsec < 1900000000) and (packetLength>7):
            fpacket, packetLength = read_packet(fh, bSPW)
            if fpacket:
                apid = read_apid(fpacket)
            else:
                continue
            ccsec = int.from_bytes(fpacket[6:10],'big')
        # if is_bad_packet(packetLength,fh):  return 104       
        if fpacket:
            mpacket = fpacket
            try:
                stime = get_anc_packet_time(mpacket[6:12],0,'msec')
                str_stime = stime.strftime('%Y-%m-%dT%H:%M:%S.%f')
                print("start_time=%s" % str_stime)
            except:
                return 120

            while fpacket and (packetLength>7):
                fpacket, packetLength = read_packet(fh, bSPW)
                if fpacket:
                    apid = read_apid(fpacket)
                else:
                    continue
                if (apid >= 550) and (apid < 750):  mpacket = fpacket
            
            if is_bad_packet(packetLength,fh):  return 104
            etime = get_anc_packet_time(mpacket[6:12],0,'msec')
        else:
            if output:
                output.write("datatype=OCI\n")
            return status
    
    else:
        # Get data type
        dtype = get_oci_data_type(fpacket)
        if (dtype<1) or (dtype>12):
            print("Invalid data type %d."%dtype)
            return 104
                    
        # Get start time
        try:
            stime = get_anc_packet_time(fpacket,28,'usec')
            str_stime = stime.strftime('%Y-%m-%dT%H:%M:%S.%f')
            print("start_time=%s" % str_stime)
            epacket = fpacket
            etime = get_anc_packet_time(epacket,28,'usec')
        except:
            return 120        
        
        # Read to end of file
        apacket = fpacket
        acomp = 0
        nagg = np.array([1,2,4,8])
        spnp = int.from_bytes(apacket[24:28],'big')
        
        while fpacket:
            etime = get_anc_packet_time(apacket,28,'usec')
            if args.verbose:
                # spnum = swap_endian(long(apacket(24:27),0))
                spnum = int.from_bytes(apacket[24:28],'big')
                tmpstr = ""
                if ((spnum-spnp) > 10):
                    tmpstr += "\nSpin number gap: %d, %d" % (spnp, spnum)
                spnp = spnum
                if (not acomp):
                    tmpstr += "\n\nInstrument configuration change at spin %d, %s" % (spnum,etime)
                    ncp,nbb,nrb,nsp,ndc,nds,btaps,rtaps,itable,baggs,raggs = get_band_dims(apacket)
                    iz = np.where((itable['dtype'] != 0) & (itable['dtype'] != 10))[0]
                    tmpstr += "\nData Type       %s" % str(itable['dtype'][iz])
                    iagg = nagg[itable['iagg'][iz]]
                    tmpstr += "\nNumber of pixels %s" % str(itable['lines'][iz]/iagg)
                    tmpstr += "\nAggregation     %s" % str(iagg)
                    tmpstr += "\nBlue bands: %d   Red bands: %d" %(nbb, nrb)
                    tmpstr += "\nBlue aggregation per tap: %s" % str(baggs)
                    tmpstr += "\nRed aggregation per tap: %s\n" % str(raggs)
                    print(tmpstr)
                    
                    apacket = fpacket
                    
            apid = 0
            while (apid != 636) and (packetLength>0):
                fpacket, packetLength = read_packet(fh, bSPW)     
                apid = read_apid(fpacket)
                if is_bad_apid(apid,fh):  return 104
            
            if is_bad_packet(packetLength,fh):  return 104
            if fpacket:  epacket = fpacket
        
        try:
            etime = get_anc_packet_time(epacket,28,'usec')
        except:
            return 120
    
    str_etime = etime.strftime('%Y-%m-%dT%H:%M:%S.%f')
    print("stop_time=%s" % str_etime)
    
    datatype_name = ''
    if dtype > 0:
        try:
            datatype_name = 'OCI' + __dtypes__[dtype]
            if datatype_name=="OCI_LUN-ST":
                print("Unexpected datatype [%s] is detected." % datatype_name)
                status = 130
        except:
            print("OCI data type name reading error.")
            status = 130
    
    if datatype_name=='':
        datatype_name = 'OCI'
    
    print("datatype=%s" % datatype_name)
    if output:
        output.write("datatype=%s\n" % datatype_name)
        if str_stime!='':
            output.write("start_time=%s\n" % str_stime)
        if str_etime!='':
            output.write("stop_time=%s\n" % str_etime)
    
    return status
