#! /usr/bin/env python3

import sys
import numpy as np
import datetime
import os.path
import argparse
import netCDF4

from tlm.PacketUtils import *
from tlm.pace.APID108 import *
from tlm.pace.APID128 import *
from tlm.pace.APID198 import *
from tlm.timestamp import *

ignored_apids = [
    636,  # OCI Ancillary packet
    700,  # OCI science packet with spin number but no time field
    720,  # OCI SWIR science packet
    751,  # HARP science packet
    848,  # SPEX science packet
]

max_SC_packet = 2048
max_OCI_packet = 1618
max_HARP2_packet = 40
max_SPEXone_packet = 298

def main():

    # Read command line options
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Convert S-band housekeeping telemetry to NetCDF4'
    )
    parser.add_argument('ifile', nargs='?', type=str,
                        help='path to S-band telemetry (HSK file) OR list of input files, one per line, in chronological order')
    parser.add_argument("-o", "--output", metavar="ofile", dest="ofile",
                        help="output NetCDF4 file; defaults to PACE.yyyymmddThhmmss.HKT.nc")
    parser.add_argument('--verbose', '-v', action='store_true',
                        default=False, help="print status messages")

    args = parser.parse_args()

    if args.ifile is None:
        parser.print_help()
        sys.exit(1)

    # Is input file tlm or list of files?
    filelist = []
    infile = os.path.expandvars(args.ifile)

    with open(infile, mode='rt') as flist:  # try to read as list of files
        try:
            input_list = True
            for ifile in flist:
                f = os.path.expandvars(ifile.rstrip())
                if os.path.isfile(f):
                    filelist.append(f)
        except UnicodeDecodeError:
            input_list = False # contains invalid byte - infile is binary

    if not len(filelist):  # if that didn't work, it's a single HSK file
        filelist.append(infile)

    # Initialize record and packet lists
    att_recordlist = []
    orb_recordlist = []
    tilt_recordlist = []
    alltimes = []
    SC_HKT_packets = []
    OCI_HKT_packets = []
    HARP2_HKT_packets = []
    SPEXone_HKT_packets = []
    oversize_packets = bytearray()

    # Step through all input files
    for filename in filelist:
        print('reading: ',filename)

        try:
            ifile = open(filename, mode='rb')
        except BaseException:
            print("Cannot open file \"%s\": exiting." % filename)
            return 100

        # Read any file header(s)
        filehdr = readFileHeader(ifile)
        if filehdr:
            if args.verbose:
                print(filehdr)
                print()

            # Is it the right kind of file?
            desired = (filehdr['subtype'] == 101 and
                       filehdr['length'] == 64 and
                       filehdr['SCID'] == b'PACE' and
                       filehdr['processorid'] in (1,2,30) )
            if not desired:
                return 110

        # Now read CCSDS packets
        while True:

            # read header and data
            header, data, rawhdr = readPacket(ifile)
            if header is None:
                break  # EOF

            # check for invalid header
            if header['length'] > 16378:
                print('Invalid CCSDS packet header:', rawhdr)
                continue

            # check for invalid timestamp
            if (header['secondary'] == 1
                and data[0] > 112   # after  2017-07-18T05:49:15Z
                and data[0] < 192): # before 2060-01-28T16:50:35Z
                ts = readTimestamp(data)
                alltimes.append(ts)
                header['timestamp'] = ts
                # keep going to save packet

            if args.verbose:
                print(header)

            if header['APID'] not in ignored_apids:
                tmpDict = {}
                packet = rawhdr+data

                if header['APID'] < 550:  # Spacecraft packet

                    if header['APID'] == 108:  # Attitude
                        myDict = APID108(data)
                        keys = ['timestamp', 'Est_Time', 'q_EciToBrf_Est', 'w_EciToBrf_Brf_Est',]
                        for k in keys:
                            tmpDict[k] = myDict[k]
                        att_recordlist.append(tmpDict)

                    elif header['APID'] == 128:  # Ephemeris
                        myDict = APID128(data)
                        keys = ['timestamp', 'FSWTime', 'scPosJ2000', 'scVelJ2000','DCM_ecef2eci',]
                        for k in keys:
                            tmpDict[k] = myDict[k]
                        orb_recordlist.append(tmpDict)

                    elif header['APID'] == 198:  # Tilt
                        myDict = APID198(data)
                        keys = ['timestamp', 'TILT_ENCPOS',]
                        for k in keys:
                            tmpDict[k] = myDict[k]
                        tilt_recordlist.append(tmpDict)

                    else:  # not yet implemented or don't care
                        pass

                    # append rawhdr, data to sc_packets
                    if len(packet) > max_SC_packet:
                        oversize_packets += packet
                    else:
                        packet = pad_packet(packet, max_SC_packet)
                        SC_HKT_packets.append(packet)

                    # endif (Spacecraft packet)

                elif header['APID'] < 750:  # OCI packet
                    if len(packet) > max_OCI_packet:
                        oversize_packets += packet
                    else:
                        packet = pad_packet(packet, max_OCI_packet)
                        OCI_HKT_packets.append(packet)

                elif header['APID'] < 800:  # HARP2 telemetry packet
                    if len(packet) > max_HARP2_packet:
                        oversize_packets += packet
                    else:
                        packet = pad_packet(packet, max_HARP2_packet)
                        HARP2_HKT_packets.append(packet)

                elif header['APID'] < 850:  # SPEXone telemetry packet
                    if len(packet) > max_SPEXone_packet:
                        oversize_packets += packet
                    else:
                        packet = pad_packet(packet, max_SPEXone_packet)
                        SPEXone_HKT_packets.append(packet)

                else:  # APID >= 850
                    pass

            # endif (not in ignored_apids)

        # close input file
        ifile.close()

    # endfor filename in filelist


    # calculate orbit parameters
    if len(orb_recordlist) > 0:
        orbitparams = derive_orbitparams(orb_recordlist)

    # get start and end times
    if len(alltimes) == 0:
        print('No input packets with valid times')
        return 110

    stime = tai58_as_datetime(alltimes[0])
    etime = tai58_as_datetime(alltimes[-1])
    daystart = stime.replace(hour=0, minute=0, second=0, microsecond=0)

    # construct product name
    prod_name = stime.strftime('PACE.%Y%m%dT%H%M%S.HKT.nc')
    if args.ofile is None:
        args.ofile = prod_name

    # create new netcdf file
    try:
        ofile = netCDF4.Dataset(args.ofile, 'w')
    except BaseException:
        print("Cannot write file \"%s\": exiting." % args.ofile)
        return 1

    # define dimensions
    ofile.createDimension('vector_elements', 3)
    ofile.createDimension('quaternion_elements', 4)

    # write raw housekeeping data to file
    group = ofile.createGroup('housekeeping_data')

    if len(SC_HKT_packets) > 0:
        ofile.createDimension('SC_hkt_pkts', None)
        ofile.createDimension('max_SC_packet', max_SC_packet)
        var = group.createVariable('SC_HKT_packets', 'u1', ('SC_hkt_pkts', 'max_SC_packet'))
        var.long_name = "Spacecraft housekeeping telemetry packets"
        var[:] = SC_HKT_packets

    if len(OCI_HKT_packets) > 0:
        ofile.createDimension('OCI_hkt_pkts', None)
        ofile.createDimension('max_OCI_packet', max_OCI_packet)
        var = group.createVariable('OCI_HKT_packets', 'u1', ('OCI_hkt_pkts', 'max_OCI_packet'))
        var.long_name = "OCI housekeeping telemetry packets"
        var[:] = OCI_HKT_packets

    if len(HARP2_HKT_packets) > 0:
        ofile.createDimension('HARP2_hkt_pkts', None)
        ofile.createDimension('max_HARP2_packet', max_HARP2_packet)
        var = group.createVariable('HARP2_HKT_packets', 'u1', ('HARP2_hkt_pkts', 'max_HARP2_packet'))
        var.long_name = "HARP2 housekeeping telemetry packets"
        var[:] = HARP2_HKT_packets

    if len(SPEXone_HKT_packets) > 0:
        ofile.createDimension('SPEXone_hkt_pkts', None)
        ofile.createDimension('max_SPEXone_packet', max_SPEXone_packet)
        var = group.createVariable('SPEXone_HKT_packets', 'u1', ('SPEXone_hkt_pkts', 'max_SPEXone_packet'))
        var.long_name = "SPEXone housekeeping telemetry packets"
        var[:] = SPEXone_HKT_packets

    if len(oversize_packets) > 0:
        ofile.createDimension('os_pkts', None)
        var = group.createVariable('oversize_packets', 'u1', ('os_pkts'))
        var.long_name = "Buffer for packets exceeding maximum size"
        var[:] = oversize_packets

    # write navigation data to file
    group = ofile.createGroup('navigation_data')

    if len(att_recordlist) > 0:     # ATTITUDE
        ofile.createDimension('att_records', None)

        var = group.createVariable(  # att_time
            'att_time', 'f8', ('att_records'), fill_value=-999.9)
        var.long_name = "Attitude sample time (seconds of day)"
        var.valid_min = 0
        var.valid_max = 86400.999999
        var.units = "seconds"
        var[:] = [seconds_since(rec['timestamp'], basetime=daystart) for rec in att_recordlist]

        var = group.createVariable(  # att_quat
            'att_quat', 'f4', ('att_records', 'quaternion_elements'), fill_value=-999.9)
        var.long_name = "Attitude quaternions (J2000 to spacecraft)"
        var.valid_min = -1
        var.valid_max = 1
        var.units = "seconds"
        var[:] = [rec['q_EciToBrf_Est'] for rec in att_recordlist]

        var = group.createVariable(  # att_rate
            'att_rate', 'f4', ('att_records', 'vector_elements'), fill_value=-999.9)
        var.long_name = "Attitude angular rates in spacecraft frame"
        var.valid_min = np.array((-0.004), 'f4')
        var.valid_max = np.array(( 0.004), 'f4')
        var.units = "radians/second"
        var[:] = [rec['w_EciToBrf_Brf_Est'] for rec in att_recordlist]

    if len(orb_recordlist) > 0:     # EPHEMERIS (orbit)
        ofile.createDimension('orb_records', None)

        var = group.createVariable(  # orb_time
            'orb_time', 'f8', ('orb_records'), fill_value=-999.9)
        var.long_name = "Orbit vector time (seconds of day)"
        var.valid_min = 0
        var.valid_max = 86400.999999
        var.units = "seconds"
        var[:] = [seconds_since(rec['timestamp'], basetime=daystart)
                                for rec in orb_recordlist]

        var = group.createVariable(  # orb_pos
            'orb_pos', 'f4', ('orb_records', 'vector_elements'), fill_value= -9999999)
        var.long_name = "Orbit position vectors (ECR)"
        var.valid_min = -7200000
        var.valid_max = 7200000
        var.units = "meters"
        var[:] = orbitparams['posr']

        var = group.createVariable(  # orb_vel
            'orb_vel', 'f4', ('orb_records', 'vector_elements'), fill_value= -9999999)
        var.long_name = "Orbit velocity vectors (ECR)"
        var.valid_min = -7600
        var.valid_max = 7600
        var.units = "meters/second"
        var[:] = orbitparams['velr']

        var = group.createVariable(  # orb_lon
            'orb_lon', 'f8', ('orb_records'), fill_value=-999.9)
        var.long_name = "Orbit longitude (degrees East)"
        var.valid_min = -180
        var.valid_max = 180
        var.units = "degrees"
        var[:] = orbitparams['lon']

        var = group.createVariable(  # orb_lat
            'orb_lat', 'f8', ('orb_records'), fill_value=-999.9)
        var.long_name = "Orbit latitude (degrees North)"
        var.valid_min = -90
        var.valid_max = 90
        var.units = "degrees"
        var[:] = orbitparams['lat']

        var = group.createVariable(  # orb_alt
            'orb_alt', 'f8', ('orb_records'), fill_value=-999.9)
        var.long_name = "Orbit altitude"
        var.valid_min = 670
        var.valid_max = 710
        var.units = "meters"
        var[:] = orbitparams['alt']

    if len(tilt_recordlist) > 0:     # TILT
        ofile.createDimension('tilt_records', None)

        var = group.createVariable(  # tilt_time
            'tilt_time', 'f8', ('tilt_records'), fill_value=-999.9)
        var.long_name = "Tilt sample time (seconds of day)"
        var.valid_min = 0
        var.valid_max = 86400.999999
        var.units = "seconds"
        var[:] = [seconds_since(rec['timestamp'], basetime=daystart)
                                 for rec in tilt_recordlist]

        var = group.createVariable(  # tilt
            'tilt', 'f4', ('tilt_records'), fill_value=-999.9)
        var.long_name = "Tilt angle"
        var.valid_min = np.array((-20.1), 'f4')
        var.valid_max = np.array(( 20.1), 'f4')
        var.units = "degrees"
        var[:] = [rec['TILT_ENCPOS'] for rec in tilt_recordlist]

    # write global metadata

    # static
    ofile.title = "PACE HKT Data"
    ofile.instrument = "Observatory"
    ofile.processing_version = "V1.0"
    ofile.Conventions = "CF-1.6"
    ofile.institution = "NASA Goddard Space Flight Center, Ocean Biology Processing Group"
    ofile.license = "http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
    ofile.naming_authority = "gov.nasa.gsfc.sci.oceancolor"
    ofile.keywords_vocabulary = "NASA Global Change Master Directory (GCMD) Science Keywords"
    ofile.stdname_vocabulary = "NetCDF Climate and Forecast (CF) Metadata Convention"
    ofile.creator_name = "NASA/GSFC"
    ofile.creator_email = "data@oceancolor.gsfc.nasa.gov"
    ofile.creator_url = "http://oceancolor.gsfc.nasa.gov"
    ofile.project = "PACE Project"
    ofile.publisher_name = "NASA/GSFC"
    ofile.publisher_email = "data@oceancolor.gsfc.nasa.gov"
    ofile.publisher_url = "http://oceancolor.gsfc.nasa.gov"
    ofile.processing_level = "1"
    ofile.cdm_data_type = ""
    ofile.CDL_version_date = "2022-04-08"
    # dynamic
    ofile.product_name = args.ofile
    ofile.time_coverage_start = datetime_repr(stime)
    ofile.time_coverage_end = datetime_repr(etime)
    ofile.history = ' '.join([v for v in sys.argv])
    if input_list:
        ofile.source = ','.join([os.path.basename(f) for f in filelist])
    ofile.date_created = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    # close file
    ofile.close()

    return 0

if __name__ == '__main__':
    sys.exit(main())
