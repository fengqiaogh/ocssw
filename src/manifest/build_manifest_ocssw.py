#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:56:47 EST 2019

@author: dshea
"""

import sys
import argparse
import subprocess
import multiprocessing as mp

machineInfo = [
        # CentOS 7
        {
            "name" : "linux_64", 
            "login" : "seadasdev201",
            "initStr" : "source .bash_profile && source /opt/rh/devtoolset-7/enable",
            "exitStr" : ""
        },
        # RHEL 8
#        {
#            "name" : "linux_64", 
#            "login" : "seadasdev301",
#            "initStr" : "source .bash_profile",
#            "exitStr" : ""
#        },
        # Ubuntu 20.04
        {
            "name" : "odps", 
            "login" : "analysis701",
            "initStr" : "source .profile",
            "exitStr" : ""
        },
        # macOS
        {
            "name" : "macosx_intel", 
            "login" : "seadas4",
            "initStr" : "source .bash_profile",
            "exitStr" : " && fix_mac_rpath.py"
        },
]

firstRunStr = " && if [ ! -d ocssw ]; then mkdir ocssw; fi && if [ ! -d ocssw/share ]; then mkdir ocssw/share; fi && if [ ! -d ocssw/testdata ]; then mkdir ocssw/testdata; fi && if [ ! -d ocssw/var ]; then mkdir ocssw/var; fi && if [ ! -d ocssw/opt ]; then mkdir ocssw/opt; fi"
saveOcsswStr = " && rm -rf saveOcssw && mkdir saveOcssw && cd ocssw && mv share testdata var ../saveOcssw"
restoreOcsswStr = " && rm -rf share/modis && mv ../saveOcssw/share ../saveOcssw/testdata ../saveOcssw/var ."
saveOptStr = " && mv opt ../saveOcssw"
restoreOptStr = " && mv ../saveOcssw/opt ."
getOcsswStr = " && cd && rm -rf ocssw && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/obpg/ocssw.git && cd ocssw && source OCSSW_bash.env"
getSubmodulesStr = " && git submodule init && git submodule update"
buildOptStr = " && ./get_lib3_src.sh && cd opt/src && ./BuildIt.py && cd ../.."
buildOcsswStr = " && mkdir build && cd build && cmake .. -DBUILD_ALL=1 && make -j 20 install"
getViirsStr = " && cd && rm -rf viirs_l1 && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/viirs/viirs_l1.git && cd viirs_l1"
buildViirsStr = " && mkdir build && cd build && cmake .. && make -j 20 install"
getFocsStr = " && cd && rm -rf focs && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/obpg/focs.git && cd focs"
buildFocsStr = " && mkdir build && cd build && cmake .. && make -j 20 install"


def doIt(cmd, logFilename, out):
    logFile = open(logFilename, "w")
    result = subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT, stdout=logFile)
    logFile.close()
    out.put(result.returncode)


def run():
    processes = []

    parser = argparse.ArgumentParser(description="Build OCSSW on all of our architectures")
    parser.add_argument("-t", "--tag", default=None,
                        help="git tag or branch that you want to build")
    parser.add_argument("-a", "--arch", default=None,
                        help="comma separated list of architectures that you want to build (linux_64,odps,macosx_intel)")
    parser.add_argument("--build_opt", default=False, action="store_true", 
                        help="build opt (lib3) first")


    args = parser.parse_args()

    # make sure the arch list is valid
    if args.arch:
        machineList = []
        for info in machineInfo:
            machineList.append(info["name"])
        archList = args.arch.split(',')
        for arch in archList:
            if arch not in machineList:
                print("Architecture", arch, "is not in the list of supported architectures")
                sys.exit(1)

    runList = []

    # Define an output queue
    output = mp.Queue()

    archFound = False
    for info in machineInfo:
        if args.arch:
            if info["name"] not in args.arch:
                continue
        archFound = True
        
        runList.append(info["name"])
        print("\n---making", info["name"])
        
        commandLine = "ssh " + info["login"] + ' "' + info["initStr"] + firstRunStr + saveOcsswStr
        if not args.build_opt:
            commandLine += saveOptStr
        commandLine += getOcsswStr
        if args.tag:
            commandLine += " && git checkout " + args.tag
        commandLine += getSubmodulesStr
        if args.build_opt:
            commandLine += buildOptStr
        else:
            commandLine += restoreOptStr
        commandLine += restoreOcsswStr
        commandLine += " && rm -rf ../saveOcssw"

        # build VIIRS
        commandLine += getViirsStr
        if args.tag:
            commandLine += " && git checkout " + args.tag
        commandLine += getSubmodulesStr
        commandLine += buildViirsStr

        # build FOCS
        commandLine += getFocsStr
        if args.tag:
            commandLine += " && git checkout " + args.tag
        commandLine += getSubmodulesStr
        commandLine += buildFocsStr

        # build OCSSW
        commandLine += " && cd && cd ocssw"
        commandLine += buildOcsswStr


        commandLine += info["exitStr"]
        commandLine += '"'

        logFilename = "build." + info["name"] + ".log"

        processes.append(mp.Process(target=doIt, args=(commandLine, logFilename, output)))

    if not archFound:
        print("Architecture", args.arch, "is not valid")
        sys.exit(1)

    print("Starting Processes")
    
    # Run processes
    for p in processes:
        p.start()

    print("Waiting for Processes")

    # Exit the completed processes
    for p in processes:
        p.join()

    print("Checking results")
   
    # check the return result
    result = True
    count = 0
    for p in processes:
        tmp = output.get()
        if tmp != 0:
            result = False
            print(runList[count], "return code =", tmp)
        else:
            print(runList[count], "Success")
            
        count += 1            
   
    if result:
        print("Everyting completed successfully.")
        return 0
    
    return 1


if __name__ == "__main__":
    sys.exit(run())
