# ! /usr/bin/env python3

import os
import sys
import subprocess
import seadasutils.anc_utils as ga
from seadasutils.setupenv import env
from seadasutils.ParamUtils import ParamProcessing
import seadasutils.ProcUtils as ProcUtils
import seadasutils.LutUtils as Lut
from seadasutils.MetaUtils import readMetadata


class modis_geo:

    def __init__(self, filename=None,
                 parfile=None,
                 geofile=None,
                 a1=None, a2=None, a3=None,
                 e1=None, e2=None, e3=None,
                 download=True,
                 entrained=False,
                 terrain=False,
                 geothresh=95,
                 sensor=None,
                 anc_file=None,
                 ancdir=None,
                 curdir=False,
                 ancdb='ancillary_data.db',
                 refreshDB=False,
                 forcedl=False,
                 lutver=None,
                 lutdir=None,
                 log=False,
                 verbose=False,
                 timeout=10.0):

        # defaults
        self.filename = filename
        self.parfile = parfile
        self.geofile = geofile
        self.ancdir = ancdir
        self.ancdb = ancdb
        self.refreshDB = refreshDB
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3
        self.download = download
        self.forcedl = forcedl
        self.entrained = entrained
        self.terrain = terrain
        self.geothresh = geothresh
        self.lutversion = lutver
        self.lutdir = lutdir
        self.log = log
        self.proctype = 'modisGEO'
        self.curdir = curdir
        self.pcf_file = None
        self.verbose = verbose
        self.dirs = {}
        self.sensor = sensor
        self.sat_name = None
        self.start = None
        self.stop = None
        self.anc_file = anc_file
        self.timeout = timeout

        # version-specific variables
        self.collection_id = '061'
        self.pgeversion = '6.1.9'
        #        self.lutversion = '0'

        if self.parfile:
            print(self.parfile)
            p = ParamProcessing(parfile=self.parfile)
            p.parseParFile(prog='geogen')
            print(p.params)
            phash = p.params['geogen']
            for param in (list(phash.keys())):
                print(phash[param])
                if not self[param]:
                    self[param] = phash[param]

    def __setitem__(self, index, item):
        self.__dict__[index] = item

    def __getitem__(self, index):
        return self.__dict__[index]

    def chk(self):
        """
        Check parameters
        """

        if self.filename is None:
            print("ERROR: No MODIS_L1A_file was specified in either the parameter file or in the argument list. Exiting")
            sys.exit(1)
        if not os.path.exists(self.filename):
            print("ERROR: File '" + self.filename + "' does not exist. Exiting.")
            sys.exit(1)
        if self.a1 is not None and not os.path.exists(self.a1):
            print("ERROR: Attitude file '" + self.a1 + "' does not exist. Exiting.")
            sys.exit(1)
        if self.a2 is not None and not os.path.exists(self.a2):
            print("ERROR: Attitude file '" + self.a2 + "' does not exist. Exiting.")
            sys.exit(99)
        if self.a3 is not None and not os.path.exists(self.a3):
            print("ERROR: Attitude file '" + self.a3 + "' does not exist. Exiting.")
            sys.exit(99)
        if self.e1 is not None and not os.path.exists(self.e1):
            print("ERROR: Ephemeris file '" + self.e1 + "' does not exist. Exiting.")
            sys.exit(1)
        if self.e2 is not None and not os.path.exists(self.e2):
            print("ERROR: Ephemeris file '" + self.e2 + "' does not exist. Exiting.")
            sys.exit(1)
        if self.e3 is not None and not os.path.exists(self.e3):
            print("ERROR: Ephemeris file '" + self.e3 + "' does not exist. Exiting.")
            sys.exit(1)
        if self.a1 is None and self.e1 is not None or self.a1 is not None and self.e1 is None:
            print("ERROR: User must specify attitude AND ephemeris files.")
            print("       Attitude/ephemeris files must ALL be specified or NONE specified. Exiting.")
            sys.exit(1)
        if self.terrain is True and not os.path.exists(self.dirs['dem']):
            print("WARNING: Could not locate MODIS digital elevation maps directory:")
            print("         '" + self.dirs['dem'] + "/'.")
            print("")
            print("*TERRAIN CORRECTION DISABLED*")
            self.terrain = False

    def utcleap(self):
        """
        Check date of utcpole.dat and leapsec.dat.
        Download if older than 14 days.
        """

        lut = Lut.LutUtils(verbose=self.verbose,
                    mission=self.sat_name)
       # (verbose=self.verbose, mission=self.sat_name)
        utcpole = os.path.join(self.dirs['var'], 'modis', 'utcpole.dat')
        leapsec = os.path.join(self.dirs['var'], 'modis', 'leapsec.dat')

        if not (os.path.exists(utcpole) and os.path.exists(leapsec)):
            if self.verbose:
                print("** Files utcpole.dat/leapsec.dat are not present on hard disk.")
                print("** Running update_luts.py to download the missing files...")
            lut.get_luts()

        elif (ProcUtils.mtime(utcpole) > 14) or (ProcUtils.mtime(leapsec) > 14):
            if self.verbose:
                print("** Files utcpole.dat/leapsec.dat are more than 2 weeks old.")
                print("** Running update_luts.py to update files...")
            lut.get_luts()

    def atteph(self):
        """
        Determine and retrieve required ATTEPH files
        """

        self.attdir1  = self.attdir2  = self.attdir3  = "NULL"
        self.attfile1 = self.attfile2 = self.attfile3 = "NULL"
        self.ephdir1  = self.ephdir2  = self.ephdir3  = "NULL"
        self.ephfile1 = self.ephfile2 = self.ephfile3 = "NULL"

        # Check for user specified atteph files
        if self.a1 is not None:
            self.atteph_type = "user_provided"
            self.kinematic_state = "SDP Toolkit"
            self.attfile1 = os.path.basename(self.a1)
            self.attdir1 = os.path.abspath(os.path.dirname(self.a1))
            self.ephfile1 = os.path.basename(self.e1)
            self.ephdir1 = os.path.abspath(os.path.dirname(self.e1))

            if self.a2 is not None:
                self.attfile2 = os.path.basename(self.a2)
                self.attdir2 = os.path.abspath(os.path.dirname(self.a2))
            if self.a3 is not None:
                self.attfile3 = os.path.basename(self.a3)
                self.attdir3 = os.path.abspath(os.path.dirname(self.a3))

            if self.e2 is not None:
                self.ephfile2 = os.path.basename(self.e2)
                self.ephdir2 = os.path.abspath(os.path.dirname(self.e2))
            if self.e3 is not None:
                self.ephfile3 = os.path.basename(self.e3)
                self.ephdir3 = os.path.abspath(os.path.dirname(self.e3))

            if self.verbose:
                print("Using specified attitude and ephemeris files.")
                print("")
                print("att_file1:", os.path.join(self.attdir1, self.attfile1))
                if self.attfile2 == "NULL":
                    print("att_file2: NULL")
                else:
                    print("att_file2:", os.path.join(self.attdir2, self.attfile2))
                if self.attfile3 == "NULL":
                    print("att_file3: NULL")
                else:
                    print("att_file3:", os.path.join(self.attdir3, self.attfile3))
                print("eph_file1:", os.path.join(self.ephdir1, self.ephfile1))
                if self.ephfile2 == "NULL":
                    print("eph_file2: NULL")
                else:
                    print("eph_file2:", os.path.join(self.ephdir2, self.ephfile2))
                if self.ephfile3 == "NULL":
                    print("eph_file3: NULL")
                else:
                    print("eph_file3:", os.path.join(self.ephdir3, self.ephfile3))
        else:
            if self.verbose:
                print("Determining required attitude and ephemeris files...")
            get = ga.getanc(filename=self.filename,
                            atteph=True,
                            ancdb=self.ancdb,
                            ancdir=self.ancdir,
                            curdir=self.curdir,
                            refreshDB=self.refreshDB,
                            sensor=self.sensor,
                            start=self.start,
                            stop=self.stop,
                            download=self.download,
                            verbose=self.verbose,
                            timeout=self.timeout)

            # quiet down a bit...
            resetVerbose = 0
            if get.verbose:
                resetVerbose = 1
                get.verbose = False

            env(get)
            if resetVerbose:
                get.verbose = True
            get.chk()
            if self.filename and get.finddb():
                get.setup()
            else:
                get.setup()
                get.findweb()

            get.locate(forcedl=self.forcedl)
            get.cleanup()

            self.db_status = get.db_status
            # DB return status bitwise values:
            # 0 - all is well in the world
            # 1 - predicted attitude selected
            # 2 - predicted ephemeris selected
            # 4 - no attitude found
            # 8 - no ephemeris found
            # 16 - invalid mission
            if self.sat_name == "terra" and self.db_status & 15:
                self.kinematic_state = "MODIS Packet"
            elif self.db_status & 12:
                if self.db_status & 4:
                    print("Missing attitude files!")
                if self.db_status & 8:
                    print("Missing ephemeris files!")
                sys.exit(31)
            else:
                self.kinematic_state = "SDP Toolkit"
                if 'att1' in get.files:
                    self.attfile1 = os.path.basename(get.files['att1'])
                    self.attdir1 = os.path.dirname(get.files['att1'])
                else:
                    print("Missing attitude files!")
                    sys.exit(31)
                if 'eph1' in get.files:
                    self.ephfile1 = os.path.basename(get.files['eph1'])
                    self.ephdir1 = os.path.dirname(get.files['eph1'])
                else:
                    print("Missing ephemeris files!")
                    sys.exit(31)
                if 'att2' in get.files:
                    self.attfile2 = os.path.basename(get.files['att2'])
                    self.attdir2 = os.path.dirname(get.files['att2'])
                if 'att3' in get.files:
                    self.attfile3 = os.path.basename(get.files['att3'])
                    self.attdir3 = os.path.dirname(get.files['att3'])
                if 'eph2' in get.files:
                    self.ephfile2 = os.path.basename(get.files['eph2'])
                    self.ephdir2 = os.path.dirname(get.files['eph2'])
                if 'eph3' in get.files:
                    self.ephfile3 = os.path.basename(get.files['eph3'])
                    self.ephdir3 = os.path.dirname(get.files['eph3'])

    def geochk(self):
        """Examine a MODIS geolocation file for percent missing data
        Returns an error if percent is greater than a threshold"""

        thresh = float(self.geothresh)
        if not os.path.exists(self.geofile):
            print("*** ERROR: geogen_modis failed to produce a geolocation file.")
            print("*** Validation test failed for geolocation file:", os.path.basename(self.geofile))
            sys.exit(1)

        metadata = readMetadata(self.geofile)
        if metadata:
            if 'QAPERCENTMISSINGDATA' in metadata:
                pctmissing = metadata['QAPERCENTMISSINGDATA']
                if pctmissing is not None:
                    pctvalid = 100 - pctmissing
                    if pctvalid < thresh:
                        print("Percent valid data (%.2f) is less than threshold (%.2f)" % (pctvalid, thresh))
                        ProcUtils.remove(self.geofile)
                        sys.exit(1)
                    else:
                        if self.verbose:
                            print("Percentage of pixels with missing geolocation: %.2f" % pctmissing)
                            print("Validation test passed for geolocation file %s" % self.geofile)
            else:
                print("Problem reading geolocation file: %s" % self.geofile)
                sys.exit(2)
        else:
            print("Problem reading geolocation file: %s" % self.geofile)
            sys.exit(2)

    def run(self):
        """
        Run geogen_modis (MOD_PR03)
        """

        if self.verbose:
            print("")
            print("Creating MODIS geolocation file...")
        geogen = os.path.join(self.dirs['bin'], 'geogen_modis')
        status = subprocess.run(geogen, shell=False).returncode
        if self.verbose:
            print("geogen_modis returned with exit status: " + str(status))

        try:
            self.geochk()
        except SystemExit as e:
            if self.verbose:
                print("Validation test returned with error code: ", e)
            print("ERROR: MODIS geolocation processing failed.")
            self.log = True
            raise
        else:
            if self.verbose:
                print("geogen_modis created %s successfully!" % self.geofile)
                print("MODIS geolocation processing complete.")

        finally:
            if self.verbose:
                print("")

            ProcUtils.remove(os.path.join(self.dirs['run'], "GetAttr.temp"))
            ProcUtils.remove(os.path.join(self.dirs['run'], "ShmMem"))
            ProcUtils.remove('.'.join([self.filename, 'met']))
            ProcUtils.remove('.'.join([self.geofile, 'met']))
            if self.log is False:
                ProcUtils.remove(self.pcf_file)
                base = os.path.basename(self.geofile)
                ProcUtils.remove(os.path.join(self.dirs['run'], ('.'.join(['LogReport', base]))))
                ProcUtils.remove(os.path.join(self.dirs['run'], ('.'.join(['LogStatus', base]))))
                ProcUtils.remove(os.path.join(self.dirs['run'], ('.'.join(['LogUser', base]))))
