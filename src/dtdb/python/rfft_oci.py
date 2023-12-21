'''
* NAME: rfft_oci.py
*
* DESCRIPTION: Executable program to automate generation of VIIRS Deep Blue 
* Aerosol products from a directory of L1A files specified in the command line.

* USAGE: python mkern_viirs.py -i [input directory] -o [output directory]  

* Created on July 8, 2019

* Author: Samuel Anderson
'''

import os
import sys
import logging
import time
import optparse
import numpy as np
from numpy import linalg as la
from scipy.stats import lognorm
from scipy.integrate import quad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image as im
import tables 
import xarray as xr
from netCDF4 import num2date
import bhmie as bh
import mie_kern as mkrn
LOG = logging.getLogger('mkern_viirs')
from itertools import permutations
from random import sample
from scipy.linalg import dft
           
###################################################
#                  Main Function                  #
###################################################

def main():

    args = sys.argv
    
    description = '''Generate mie kernal transforms.'''
    usage = "usage: python mkern_viirs.py -i [input directory] -o [output directory]"
    version = "v1"

    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments
    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    mandatoryGroup.add_option('-i','--in',
                      action="store",
                      dest="idir" ,
                      type="string",
                      help="The full path of the input directory of L1B files to be processed")

    mandatoryGroup.add_option('-o','--out',
                      action="store",
                      dest="odir" ,
                      type="string",
                      help="The output directory path")

    parser.add_option_group(mandatoryGroup)

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize behavior of this program.")

    optionalGroup.add_option('-k','--kern',
                      action="store",
                      dest="kern" ,
                      type="string",
                      help="The full path of an input kernel file")

    parser.add_option('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help='each occurrence increases verbosity 1 level from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG')

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (opt, args) = parser.parse_args()
    # Set up the logging levels
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[opt.verbosity])
    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = []
    mand_errors = []
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not opt.__dict__[m]:
            isMissingMand = True
            print(m_err)
    if isMissingMand :
        print("Incomplete mandatory arguments, aborting...")
        return 
    
    # Check that the output directory actually exists
    opt.odir = os.path.expanduser(opt.odir)
    LOG.info('Processing files in  %s' % (opt.idir) )

    if not opt.odir:
        dir = alg + datetime.now().isoformat()
        o_path = os.getcwd() + dir
        os.mkdir(o_path)            
    elif os.path.isdir(opt.odir):
        o_path = opt.odir
    elif os.path.isdir(os.path.dirname(opt.odir)):
        o_path = opt.odir
        os.mkdir(o_path)            
    else:
        print("Output path invalid")
        return

    color_lvl = 8
    rgb = list(permutations(range(0,256,color_lvl),3))
    colors = np.array(sample(rgb,120))
    colors = colors/256.0
    
    start = 200
    stop = 10000
    m = 1.5+1.0E-10j 
    """
    for re in np.linspace(1.30, 1.60, 10): 
        m = re+0j 
        name = str(int(re*100)) + "_scaled_" 
        start = 100/(re-1)
        stop = 10000/(re-1)
        mk = mkrn.mie_kern(name, 'oci', m, start, stop, 50)
        kname = mk.save_all(opt.odir)
        mk.generate()        
        kname = mk.save_all(opt.odir)    
    """
    if not opt.kern:
        r = np.real(m)
        name = str(int(r*100))  
        mk = mkrn.mie_kern(name, 'oci', m, start, stop, 160)
        mk.generate()        
        kname = mk.save_all(opt.odir)
    else:
        mk = mkrn.mie_kern()
        mk.load(opt.kern)
       
    content_oci = os.listdir(opt.idir)
    content_oci.sort()
    for x in content_oci: 
        if x[0] == ".":
            continue                
        inpath = opt.idir + "/" + x            
        if not os.path.isfile(inpath):
            continue 
        # open oci file and read reflectances into matrix refl
        fin = tables.open_file(inpath, 'a')
        refl = np.append(fin.get_node("/", "observation_data/Lt_blue")[8:59], \
                        fin.get_node("/", "observation_data/Lt_red")[1:,:,:],axis=0)
        shp = refl.shape
               
        rayl = (mk.wl[0]/mk.wl[:])**4
        irayl = 1/rayl
        unit = np.ones((mk.nwl,mk.nwl))
        test = (unit.T*irayl).T
        DFT = dft(mk.nwl)
        RDFT = DFT*irayl
        IRDFT = la.inv(RDFT)
        
        toc0 = time.perf_counter()        
        frefl = np.tensordot(RDFT,refl,(((1),(0))))
        toc1 = time.perf_counter()
        tpp = (toc1-toc0)/shp[1]/shp[2]
        print("compute inversion", (toc1-toc0), "sec per oci granule,", tpp, "sec per pixel", flush=True)
        
        irefl1 = np.absolute(np.tensordot(IRDFT,frefl,(((1),(0)))))
        toc2 = time.perf_counter()
        tpp = (toc2-toc1)/shp[1]/shp[2]
        print("reconstruct reflectance", (toc2-toc1), "sec per oci granule,", tpp, "sec per pixel", flush=True)
        
        err = la.norm(refl-irefl1, ord=2, axis=0)
        toc3 = time.perf_counter()
        tpp = (toc3-toc2)/shp[1]/shp[2]
        print("compute rms error", (toc3-toc2), "sec per oci granule,", tpp, "sec per pixel", flush=True)
        
        frefl[0,:,:] = 0
        irefl2 = np.absolute(np.tensordot(IRDFT,frefl,(((1),(0)))))
        toc4 = time.perf_counter()
        tpp = (toc4-toc3)/shp[1]/shp[2]
        print("reconstruct reflectance", (toc2-toc1), "sec per oci granule,", tpp, "sec per pixel", flush=True)
        
        fname = opt.odir + "/RFFT_" + x + mk.name + ".png"      
        print( fname )
        
        plt.clf()
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(50, 30))
        fig.suptitle('Rayleigh adapted fourier transform: ' )
        axs[0,0].set(xlabel='Wavelength (Nm)', ylabel='Reflectance',
               title='Measured Reflectance as a function of Wavelength')
        axs[0,0].grid(b=True, which='both', axis='both')
        axs[1,0].set(xlabel='wavelength', ylabel='Rayleigh-corrected Reflectance',
                       title='Rayleigh adapted Fourier Transform')
        axs[1,0].grid(b=True, which='both', axis='both')
        axs[1,1].set(xlabel='Root sum of squares error ', ylabel='Frequency in Granule',
               title='Histogram of errors in reconstructed reflectance')
        axs[1,1].grid(b=True, which='both', axis='both')
        axs[0,1].set(title='Rayleigh-Corrected True color image')

        ro = irefl2[np.where(mk.wl==670)]
        go = irefl2[np.where(mk.wl==550)]
        bo = irefl2[np.where(mk.wl==485)]
        maxo_all = max(np.max(ro), np.max(go), np.max(bo)) + 0.001      
        ro /= maxo_all
        go /= maxo_all
        bo /= maxo_all
        gamma = 0.5 
        ga = np.ones_like(ro[0])*gamma   
        rg = np.power(ro[0], ga)
        gg = np.power(go[0], ga)
        bg = np.power(bo[0], ga)
        
        scale = 255 
        im_red = im.fromarray(np.uint8(np.clip(rg,0,1.0)*scale))
        im_green = im.fromarray(np.uint8(np.clip(gg,0,1.0)*scale))
        im_blue = im.fromarray(np.uint8(np.clip(bg,0,1.0)*scale))
        
        im_rgb = im.merge("RGB", (im_red, im_green, im_blue))
        axs[0,1].imshow(im_rgb)
        freq = range(mk.nwl)
        column = 1050
        for i in range(0,1700,100):
            axs[0,0].plot(mk.wl,refl[:,i,column])
            axs[0,0].plot(mk.wl,irefl1[:,i,column], color=colors[int(7*i)%120])    
            axs[1,0].plot(mk.wl,irefl2[:,i,column], color=colors[int(7*i)%120])
            axs[0,1].scatter(column,i,marker=".",color=colors[int(7*i)%120])
        ehist,be = np.histogram(err)
        axs[1,1].semilogy(be[:-1],ehist)    
    #    rhist,be = np.histogram(ratio, bins=20)
    #    axs[1,1].semilogy(be[:-1],rhist)    
    #    plt.show(block=True)
        plt.savefig(fname)
        plt.close(fig)
        fin.close()
        
    LOG.info('Exiting...')
    
if __name__=='__main__':
    sys.exit(main())        
    