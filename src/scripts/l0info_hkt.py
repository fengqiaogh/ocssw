#!/usr/bin/env python3

from l0info_utils import read_cfe_header, read_dsb_header

__version__ = '2.1.1_2022-12-23'

def l0info_hkt(args, fh, output):
    # procedure to get start and end times from PACE S-band HKT data file name
    print("Running l0info_hkt (version: %s) \n" % __version__)
    
    if args.verbose:
        print("Reading PACE S-band data file.")
    
    try:
        stime, dtype = read_cfe_header(fh, bVerbose = False)
        print("start_time=%s" % stime)
        if output:
            output.write("start_time=%s\n" % stime)
    except:
        return 120
    
    try:
        etime = read_dsb_header(fh, bVerbose = False)
        print("stop_time=%s" % etime)
        if output:
            output.write("stop_time=%s\n" % etime)
    except:
        return 120

    return 0
