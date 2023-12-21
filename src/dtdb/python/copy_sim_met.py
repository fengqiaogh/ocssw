'''
* NAME: copy_sim_met.py
*
* DESCRIPTION: Executable program to copy meteorological data from simulation 
* files to a group called meteorology in another file.

* USAGE: python copy_sim_met.py -l [L1A directory] -o [output directory] 

* Created on May 10, 2019

* Author: Samuel Anderson
'''

import os
import sys
import StringIO
import optparse as optparse
import logging
import numpy as np
import tables 

LOG = logging.getLogger('copy_sim_met')
         
###################################################
#                  Main Function                  #
###################################################

def main():

    args = sys.argv
    
    description = '''Copy ancillary data from simulation files to output.'''
    usage = "usage: python copy_sim_met.py  -a [aer file] \
    -c [chm file] -m [met file] -o [output file] [options]"
    version = "v1"

    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments
    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    parser.add_option_group(mandatoryGroup)

    mandatoryGroup.add_option('-a','--aer',
                      action="store",
                      dest="aer" ,
                      type="string",
                      help="The full path of the input aer file")

    mandatoryGroup.add_option('-c','--chm',
                      action="store",
                      dest="chm" ,
                      type="string",
                      help="The full path of the input chm file")

    mandatoryGroup.add_option('-m','--met',
                      action="store",
                      dest="met" ,
                      type="string",
                      help="The full path of the input met file")

    mandatoryGroup.add_option('-o','--out',
                      action="store",
                      dest="out" ,
                      type="string",
                      help="The output directory path")

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize behavior of this program.")

    parser.add_option('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help='each occurrence increases verbosity 1 level from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG')

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (options, args) = parser.parse_args()

    # Set up the logging levels
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[options.verbosity])

    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = []
    mand_errors = []
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            print m_err
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    # Check that the output directory actually exist
    options.out = os.path.expanduser(options.out)
        
    # Check that we don't have mutually exclusive options

    LOG.info('Copying ancillary data  %s' % (options.out) )

    ydim = 1721
    xdim = 1271
    
    try:
        df = tables.open_file(options.out, "w")
        mg = df.create_group("/", "Meteorology", 'Meteorology data')
        atom = tables.Float32Atom()
        
        met_file = tables.open_file(options.met, "r")
        QV10M_o = df.create_carray(mg, 'QV10M', atom, (ydim,xdim), "")
        QV10M_o[:,:] = met_file.get_node("/", "QV10M")[0,:,:]
        SLP_o = df.create_carray(mg, 'SLP', atom, (ydim,xdim), "")
        SLP_o[:,:] = met_file.get_node("/", "SLP")[0,:,:]
        T10M_o = df.create_carray(mg, 'T10M', atom, (ydim,xdim), "")
        T10M_o[:,:] = met_file.get_node("/", "T10M")[0,:,:]
        TQV_o = df.create_carray(mg, 'TQV', atom, (ydim,xdim), "")
        TQV_o[:,:] = met_file.get_node("/", "TQV")[0,:,:]
        PS_o = df.create_carray(mg, 'PS', atom, (ydim,xdim), "")
        PS_o[:,:] = met_file.get_node("/", "PS")[0,:,:]
        U10M_o = df.create_carray(mg, 'U10M', atom, (ydim,xdim), "")
        U10M_o[:,:] = met_file.get_node("/", "U10M")[0,:,:]
        V10M_o = df.create_carray(mg, 'V10M', atom, (ydim,xdim), "")
        V10M_o[:,:] = met_file.get_node("/", "V10M")[0,:,:]
        met_file.close()
        
        chm_file = tables.open_file(options.chm, "r")
        chm_file.close()
        
        aer_file = tables.open_file(options.aer, "r")
        aer_file.close()

        df.flush()
        df.close()

    except Exception, e:
        print e.message + " ... exiting"
        df.close()

    LOG.info('Exiting...')
    return 0

if __name__=='__main__':
    sys.exit(main())        
        