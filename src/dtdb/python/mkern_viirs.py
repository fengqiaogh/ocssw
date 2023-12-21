'''
* NAME: mkern_viirs.py
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
import mie_kern
LOG = logging.getLogger('mkern_viirs')
from itertools import permutations
from random import sample

def normalize(a, order=2, axis=-1 ):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

def mie(sz,wl,m):
    s1,s2,qx,qs,qb,gs = bh.bhmie(2*np.pi*sz/wl,m,1)    
    return qs

           
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

    parser.add_option_group(mandatoryGroup)

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
        parser.error("Incomplete mandatory arguments, aborting...")
    # Check that the output directory actually exists
    opt.odir = os.path.expanduser(opt.odir)
    LOG.info('Processing VIIRS files in  %s' % (opt.idir) )

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
    stop = 4000
    awl = [ 412, 488, 550, 670, 865, 1240, 1610, 2250 ]
    wl = np.copy(awl)
    nwl = len(wl)
    lsz = np.linspace(np.log(start),np.log(stop),nwl)
    sz = np.exp(lsz)
    nsz = len(sz)
    #sz[nsz-1] = 10000
    x = np.arange(start,stop,1, )
    nx = len(x)
    # regularization matrices
    I = np.identity(nsz, dtype=np.float64)
    # first difference matrix
    H0 = -np.identity(nsz, dtype=np.float64)
    H0[0,0] = 0
    for i in range(1,nsz):
        for j in range(nsz):
            H0[i,j] = 1 if (j==i-1) else H0[i,j]
    # sum of first differences squared
    H1 = np.dot(H0.T,H0)
    # sum of second differences squared
    H2 = H1.copy()
    H2[0,:] = 0
    H2 = np.dot(H2.T,H2) 
    # sum of third differences squared
    H3 = np.zeros([nsz,nwl])  
    for i in range(3,nsz):
        for j in range(nwl):
            if (j==i):
                 H3[i,j] = -1  
            elif (j==i-1):
                H3[i,j] = 3
            elif (j==i-2):
                H3[i,j] = -3
            elif (j==i-3):
                H3[i,j] = 1
    H3 = np.dot(H3.T,H3) 

    # The standard deviations of the particle size distributions for each size bin     
    lsig = np.zeros(nsz)
    sig = np.zeros(nsz)
    for i in range(1,nsz-1):
        lsig[i] = (lsz[i+1]-lsz[i-1])/5.0
        sig[i] = (sz[i+1]-sz[i-1])/4.5
    lsig[0] = np.log((sz[1]-sz[0])/2)
    sig[0] = (sz[1]-sz[0])/2
    lsig[nsz-1] = (lsz[nsz-1]-lsz[nsz-2])/2
    sig[nsz-1] = (sz[nsz-1]-sz[nsz-2])/2
    
    # Index of refraction    
    m = np.zeros((nwl), dtype=np.complex64)
    m[0] = 1.33+6.0E-04j
    m[1] = 1.33+6.0E-04j
    m[2] = 1.33+2.0E-03j
    m[3] = 1.33+1.0E-05j
    m[4] = 1.33+1.0E-06j
    m[5] = 1.33+1.0E-03j
    m[6] = 1.33+8.0E-02j
    m[7] = 1.33+8.0E-02j
    
    #mk = mie_kern()
    kbh = np.zeros((nwl, nsz), dtype=np.float64)
    errs = np.zeros((nwl, nsz), dtype=np.float64)
    # compute kernel based on Mie scattering
    lpdf = lambda y, s, scale: lognorm.pdf(y, s=s, loc=0, scale=scale)
    for iwl in range(nwl):
        for isz in range(1,nsz):
            lmie = lambda y: lpdf(y,lsig[isz],sz[isz])*mie(y,wl[iwl],m[iwl])
            kbh[iwl,isz],errs[iwl,isz] = quad(lmie, sz[isz]-3*sig[isz], sz[isz]\
                                         + 3*sig[isz], epsabs=0.001,limit=100)
     #       kbh[iwl,isz] = mie(y,sz[isz],m)
            
# Add Rayleigh scattering component to the first column  
    kbh[:,0] = 1*(wl[0]/wl[:])**4 
# Normalize the kernal
    kbh = normalize(kbh, axis=0)
    mpl.rcParams["font.size"] = 18    
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    wlg, szg = np.meshgrid(wl, sz)
    ax.plot_wireframe(wlg, szg, kbh.T)
    plt.show()
    # Inverse kernel
    rp0 = 1e-1
    rp1 = 1e-5
    rp2 = 1e-5
    rp3 = 1e-2
    rikern = la.inv(np.dot(kbh.T,kbh) + rp0*I + rp1*H1 + rp2*H2 + rp3*H3)
       
    # Process and append files and generate db files
    content_viirs = os.listdir(opt.idir)
    content_viirs.sort()
    for x in content_viirs: 
        if x[0] == ".":
            continue                
        inpath = opt.idir + "/" + x            
        if not os.path.isfile(inpath):
            continue 
        # open viirs file and read reflectances into matrix refl
        fin = tables.open_file(inpath, 'a')
        refl = np.append(fin.get_node("/", "reflectance/TOA")[:6,:,:], \
                        fin.get_node("/", "reflectance/TOA")[7:,:,:],axis=0)
               
        psd = np.tensordot(np.dot(rikern,kbh.T),refl,(((1),(0)))) 
        irefl = np.tensordot(kbh,psd,(((1),(0))))
        err = la.norm(refl-irefl, ord=2, axis=0)/(la.norm(refl, ord=2, axis=0)+.001)
        ratio = la.norm(psd, ord=2, axis=0)/(la.norm(refl, ord=2, axis=0)+.001)
        
        print( inpath )
        fname = opt.odir + "/SZ_" + x + ".png"      
        
        plt.clf()
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(50, 30))
        fig.suptitle('Reflectance as a function of wavelength and particle size: ' + \
                     "{: .0e}".format(rp0) + " " + "{: .0e}".format(rp1) + " " + \
                     "{: .0e}".format(rp2) + " " + "{: .0e}".format(rp3))
        axs[0,0].set(xlabel='Wavelength (Nm)', ylabel='Reflectance',
               title='Reflectance as a function of Wavelength')
        axs[0,0].grid(b=True, which='both', axis='both')
        axs[1,0].set(xlabel='Particle Size (Nm)', ylabel='Inversion',
                       title='Representation as a function of Particle Size')
        axs[1,0].grid(b=True, which='both', axis='both')
        axs[1,1].set(xlabel='Root sum of squares error ', ylabel='Frequency in Granule',
               title='Histogram of errors in reconstructed reflectance')
        axs[1,1].grid(b=True, which='both', axis='both')
        axs[0,1].set(title='True color image')
    
        maxo_all = max(np.max(refl[3]), np.max(refl[2]), np.max(refl[1])) + 0.001      
        ro = refl[3]/maxo_all
        go = refl[2]/maxo_all
        bo = refl[1]/maxo_all
    
        gamma = 0.5 
        ga = np.ones_like(ro)*gamma   
        rg = np.power(ro, ga)
        gg = np.power(go, ga)
        bg = np.power(bo, ga)
        
        scale = 255 
        im_red = im.fromarray(np.uint8(np.clip(rg,0,1.0)*scale))
        im_green = im.fromarray(np.uint8(np.clip(gg,0,1.0)*scale))
        im_blue = im.fromarray(np.uint8(np.clip(bg,0,1.0)*scale))
        
        im_rgb = im.merge("RGB", (im_red, im_green, im_blue))
        axs[0,1].imshow(im_rgb)
        column = 1800
        for i in range(0,3200,10):
            axs[0,0].plot(wl,refl[:,i,column])
            axs[0,0].plot(wl,irefl[:,i,column], color=colors[int(7*i)%120])    
            axs[1,0].semilogx(sz,psd[:,i,column], color=colors[int(7*i)%120])
            axs[0,1].scatter(column,i,marker=".",color=colors[int(7*i)%120])
    #                    axs[0,1].scatter(i, err[500,i], marker='o', color=colors[int(13*i)%120])    
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
    