'''
* NAME: dtdb_viirs_l1a.py
*
* DESCRIPTION: Executable program to automate download of L1A files in specified date-time range.

* USAGE: dtdb_viirs_l1a.py -b [start date-time] -e [end date-time] -o [output directory] 

* Created on February 27, 2021

* Author: Samuel Anderson
'''

import os
import sys
from datetime import datetime
import viirs_files as vf
import optparse as optparse
import logging
LOG = logging.getLogger('dtdb_viirs_l1a')

#str_url ="http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/"
str_url ="https://oceandata.sci.gsfc.nasa.gov/ob/getfile/"
str_wget ="wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --content-disposition "
str_wdir ="-P "

def init_dir(inpath):
    files = []
    dircontents = os.listdir(inpath)
    dircontents.sort()
    for x in dircontents:           
        if x[0] == ".":
            continue                
        infilepath = inpath + x            
        if not os.path.isfile(infilepath):
            continue            
        files.append(infilepath)
    return files

def hour2min(hour, min):
    tot_min = hour*6000 + min
    return tot_min

def min2hour(tot_min):
    h = tot_min//6000
    day = h//24
    hour = h - day*24
    min = tot_min - h*6000
    return day, hour, min

def date2day(year, month, day):    
    day_of_year = date(year, month, day).timetuple().tm_yday
    return day_of_year

def day2date(year, day_of_year):    
    daystr = '{:03}'.format(day_of_year)
    dt = datetime.strptime(str(year)+"-"+daystr,"%Y-%j")
    return dt.year, dt.month, dt.day      

def getDtFilename(l1b_name):
    DtFilename = l1b_name.replace("L1B", "L2.AER_DT")
    return DtFilename
      
def getDbFilename(l1b_name):
    DbFilename = l1b_name.replace("L1B", "L2.AER_DB")
    return DbFilename
             
###################################################
#                  Main Function                  #
###################################################

def main():

    args = sys.argv
    
    description = '''Download L1A files in specified date-time range.'''
    usage = "usage: get_viirs_l1a.py -b [start date-time] -e [end date-time] -o [output directory]"
    version = "v1"

    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments
    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    parser.add_option_group(mandatoryGroup)

    mandatoryGroup.add_option('-b','--begin',
                      action="store",
                      dest="begin" ,
                      type="string",
                      help="Begin date-time")

    mandatoryGroup.add_option('-e','--end',
                      action="store",
                      dest="end" ,
                      type="string",
                      help="End date-time")

    mandatoryGroup.add_option('-z','--platform',
                      action="store",
                      dest="platform" ,
                      type="string",
                      help="VIIRS Platform (snpp, jpss1 or jpss2")

    mandatoryGroup.add_option('-d','--alg',
                      action="store",
                      dest="alg" ,
                      type="string",
                      help="Dark Target (dt) or Deep Blue (db)")

    mandatoryGroup.add_option('--l1_par',
                      action="store",
                      dest="l1_par" ,
                      type="string",
                      help="The full path of the L1 PCF")

    mandatoryGroup.add_option('--l2_par',
                      action="store",
                      dest="l2_par" ,
                      type="string",
                      help="The full path of the L2 PCF")

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize behavior of this program.")

    optionalGroup.add_option('--bowtie',
                      action="store_true",
                      dest="bowtie" ,
                      default=False,
                      help="Remove bowtie fill values (default no)")
    
    optionalGroup.add_option('--geo',
                      action="store_true",
                      dest="geo" ,
                      default=False,
                      help="Write geolocation group (default no)")
    
    optionalGroup.add_option('--anc',
                       action="store_true",
                      dest="anc" ,
                      default=False,
                      help="Write ancillary group (default no)")

    optionalGroup.add_option('--obs',
                      action="store_true",
                      dest="obs" ,
                      default=False,
                      help="Write observations group (default no)")

    optionalGroup.add_option('--stats',
                      action="store_true",
                      dest="stats" ,
                      default=False,
                      help="Write statistics group (default no)")

    optionalGroup.add_option('--glint',
                      action="store_true",
                      dest="glint" ,
                      default=False,
                      help="Do not mask glint (default no)")

    optionalGroup.add_option('--cloud',
                      action="store_true",
                      dest="cloud" ,
                      default=False,
                      help="Do not mask clouds (default no)")

    optionalGroup.add_option('--float',
                      action="store_true",
                      dest="float" ,
                      default=False,
                      help="Floating point rather than short integer format (default no)")

    optionalGroup.add_option('-o','--output_path',
                      action="store",
                      dest="odir" ,
                      type="string",
                      help="The full path of the target directory for L1A files")

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
            print (m_err)
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    # Check that the output directory actually exist
    options.outputPath = os.path.expanduser(options.odir)
        
    # Check that we don't have mutually exclusive options

    LOG.info('Downloading L1A in  %s' % (options.odir) )

    if not options.outputPath:
        dir = "L1B_" + datetime.now().isoformat()
        o_path = os.getcwd() + dir
        os.mkdir(o_path)            
    elif os.path.isdir(options.outputPath):
        o_path = options.outputPath
    elif os.path.isdir(os.path.dirname(options.outputPath)):
        o_path = options.outputPath
        os.mkdir(o_path)            
    else:
        print ("Output path invalid")
        return
    
    str_l1a_dir = o_path + "/l1a/"
    if not os.path.isdir(str_l1a_dir):
        os.mkdir(str_l1a_dir)            
    str_l1b_dir = o_path + "/l1b/"
    if not os.path.isdir(str_l1b_dir):
        os.mkdir(str_l1b_dir)            
    str_geo_dir = o_path + "/geo/"
    if not os.path.isdir(str_geo_dir):
        os.mkdir(str_geo_dir)            
    str_dt_dir = o_path + "/dt/"
    if not os.path.isdir(str_dt_dir):
        os.mkdir(str_dt_dir)            
    str_db_dir = o_path + "/db/"
    if not os.path.isdir(str_db_dir):
        os.mkdir(str_db_dir)            
    str_anc_dir = o_path + "/anc/"
    if not os.path.isdir(str_anc_dir):
        os.mkdir(str_anc_dir)            

    # Download L1A files
    
    byear = int(options.begin[0:4])
    bday = int(options.begin[4:7])
    bhour = int(options.begin[7:9])
    bmin = ((int(options.begin[9:13]) + 599)//600) * 600
    btmin = hour2min(bhour, bmin)
    eyear = int(options.end[0:4])
    eday = int(options.end[4:7])
    ehour = int(options.end[7:9])
    emin = (int(options.end[9:13])//600) * 600
    etmin = hour2min(ehour, emin)

    str_platform = ".nc"
    if options.platform == "snpp":
        str_platform = "SNPP_VIIRS."
    elif options.platform == "jpss1":
        str_platform = "JPSS1_VIIRS."
    elif options.platform == "jpss2":
        str_platform = "JPSS2_VIIRS."

    algs = []        
    if options.alg == "darktarget":
        algs = ["darktarget"]
    elif options.alg == "deepblue":
        algs = ["deepblue"]
    elif options.alg == "all":
        algs = ["darktarget", "deepblue"]
    else:
        print ("No algorithm specified.  Exiting ...")
        return
    
    if byear != eyear:
        print ("Beginning and end years must be the same.  Exiting ...")
        return
    elif eday >= bday:
        t = btmin
        t24 = hour2min(24,0)
        tal = etmin + t24*(eday - bday)
        while t < tal:
# get l1a
            d, h, s = min2hour(t)
            year, month, day = day2date(byear, bday+d)
            str_date_time0 = "%04d%02d%02dT%02d%04d" % (byear, month, day, h, s)
            str_date_time1 = "%04d%02d%02dT%02d%04d" % (byear, month, day, h, s+1)
            str_l1a_name0 = str_platform + str_date_time0 + ".L1A.nc"
            str_l1a_name1 = str_platform + str_date_time1 + ".L1A.nc"
            print ("Processing: " + str_l1a_name0)
            sfl1a = str_l1a_dir + str_l1a_name0
            l1a_files = init_dir( str_l1a_dir ) 
            l1a_file = vf.l1_file(sfl1a)
            str_geo_name0 = str_platform + str_date_time0 + ".GEO.nc"   
            str_geo_name1 = str_platform + str_date_time1 + ".GEO.nc"   
            gmod = str_geo_dir + str_geo_name0
            geo_files = init_dir( str_geo_dir )    
            lmod = str_l1b_dir + l1a_file.getModFilename()                
            l1b_files = init_dir( str_l1b_dir )    
            if (not sfl1a in l1a_files) and (not lmod in l1b_files):
                srpath = str_url + str_l1a_name0
                command = str_wget + str_wdir + str_l1a_dir + " " + srpath
                print (command)                             
                result = os.system(command)
                if result > 0:
                    srpath = str_url + str_l1a_name1
                    command = str_wget + str_wdir + str_l1a_dir + " " + srpath
                    print (command)                             
                    result = os.system(command)
                    if result > 0:
                        print ("L1A download failure: " + str_l1a_name0)
                        return
                    else:
                        sfl1a = str_l1a_dir + str_l1a_name1
                        

# get met1 if necessary
            year, month, day = day2date(byear, bday+d)
            str_date_hour = "%04d%02d%02dT%02d%04d" % (byear, month, day, h, 0)
            str_anc_name = "GMAO_MERRA2." + str_date_hour + ".MET.nc"
            met1 = str_anc_dir + str_anc_name
            anc_files = init_dir( str_anc_dir )    
            if not met1 in anc_files:
                arpath = str_url + str_anc_name
                command = str_wget + str_wdir + str_anc_dir + " " + arpath
                print (command)                             
                result = os.system(command)
# get met2 if necessary
            if (h==23):
                year, month, day = day2date(byear, bday+d+1)
                str_date_hour = "%04d%02d%02dT%02d%04d" % (byear, month, day, 0, 0)
            else:
                str_date_hour = "%04d%02d%02dT%02d%04d" % (byear, month, day, h+1, 0)
            str_anc_name = "GMAO_MERRA2." + str_date_hour + ".MET.nc"
            met2 = str_anc_dir + str_anc_name
            if not met2 in anc_files:
                arpath = str_url + str_anc_name
                command = str_wget + str_wdir + str_anc_dir + " " + arpath
                print (command)                             
                result = os.system(command)
# geolocate
            if not gmod in geo_files:                           
                grpath = str_url + str_geo_name0
                command = str_wget + str_wdir + str_geo_dir + " " + grpath
                print (command)                             
                result = os.system(command)
                if result > 0:
                    grpath = str_url + str_geo_name1
                    command = str_wget + str_wdir + str_geo_dir + " " + grpath
                    print (command)                             
                    result = os.system(command)
                    if result > 0:
                        print ("L1A download failure: " + str_l1a_name0)
                        return
                    else:
                        gmod = str_geo_dir + str_geo_name1
            print ("geo file: " + gmod)
# calibrate
            if not lmod in l1b_files:                           
                command = "calibrate_viirs " + options.l1_par + " ifile=" + sfl1a + " l1bfile_mod=" + lmod
                print (command)                             
                result = os.system(command)
            print ("l1b file: " + lmod)
# resample to fill bowtie pixels
            if options.bowtie:
                cgmod = "geo_resam_" + datetime.now().strftime("%f") + ".nc"
                command = "cp " + gmod + " " + cgmod
                result = os.system(command)
                clmod = "l1b_resam_" + datetime.now().strftime("%f") + ".nc" 
                command = "cp " + lmod + " " + clmod
                result = os.system(command)
                command = "resam_viirs " + "ifile=" + clmod + " geofile=" + cgmod
                print ("Resampling granule to fill bowtie regions")                             
                result = os.system(command)
            else:
                cgmod = gmod
                clmod = lmod
                
# run the aerosol detection algorithms
            for alg in algs:
                if alg == "deepblue":
                    alg_str = "deepblue "
                    opath = str_db_dir + l1a_file.getDbFilename()
                else:
                    alg_str = "darktarget "
                    opath = str_dt_dir + l1a_file.getDtFilename()
                l2par = " par=" + options.l2_par
                ialg = " alg=" + alg_str
                il1b = " ifile=" + clmod
                igeo = " geofile=" + cgmod
                ianc1 = " met1=" + met1
                ianc2 = " met2=" + met2
                opth = " ofile=" + opath
                ogeo = " geolocation" if options.geo else ""        
                oanc = " ancillary" if options.anc else ""        
                oobs = " observations" if options.obs else ""        
                ostat = " statistics" if options.stats else ""        
                oglnt = " maskglint=off" if options.glint else ""        
                ocld = " maskcloud=off" if options.cloud else ""        
                oflt= " short_format=off" if options.float else ""        
                lprw = " lines_per_rw=10" 
                command = "dtdb " + l2par + ialg + il1b + igeo + ianc1 + ianc2 + opth + ogeo + oanc + oobs + ostat + oglnt + ocld + oflt + lprw
                print(command)                             
                result = os.system(command)
#            command = "rm " + cgmod + " " + clmod
            
            t += 600

    else:
        print ("Invalid start and end dates specified.  Exiting ...")
        return

    print ("All L1A files in range processed")
    
    LOG.info('Exiting...')
    return 0

if __name__=='__main__':
    sys.exit(main())        
        
