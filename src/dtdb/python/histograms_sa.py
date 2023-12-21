#!/usr/bin/env python
# encoding: utf-8

import os 
import sys
import datetime
import math
import argparse
from PIL import Image as im
from PIL import ImageDraw
from PIL import ImageFilter
from PIL import ImageMath
from PIL import ImageEnhance
from PIL import ImageOps as imo
from PIL import ImageColor as imc

import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

import viirs_files as vf
from dtdb import input 
from dtocean import dtocean
from dbocean import dbocean

D490 = 0
D550 = 1
D670 = 2
D860 = 3
D124 = 4
D164 = 5
D213 = 6
W410 = 0
W445 = 1
W490 = 2
W550 = 3
W670 = 4
W865 = 5
W1240 = 6
W1610 = 7
W2250 = 8
NWL =  9
NDWL = 7
NWS = 4
NSZA = 11
NVZA = 16
NRAA = 16
NAOT = 6
NBIG = 5
NSMALL = 4
DBNSZA = 22
DBNVZA = 20
DBNRAA = 21
DBNAOT = 14
DBNAOTM = 7
DBNFMF = 14
DBNFMFM = 8
DBNWS = 6
DBNCHL = 4

D2R = np.pi/180.0

wlstr = ["410","445","490","550","670","865","1240","1610","2250"]

Kb = 0.114
Kr = 0.299

RGB2YCbCr = np.array([[0.257, 0.504, 0.098],[-0.148, -0.291, 0.439],[0.439, -0.368, -0.071]])
YCbCr2RGB = np.array([[1.0, 0.0, 1.596],[1.0, -0.392, -0.813],[1.0, 2.017, 0.0]])

def rgb2ycbcr(r,g,b):
    y = RGB2YCbCr[0,0]*r + RGB2YCbCr[0,1]*g + RGB2YCbCr[0,2]*b
    u = RGB2YCbCr[1,0]*r + RGB2YCbCr[1,1]*g + RGB2YCbCr[1,2]*b
    v = RGB2YCbCr[2,0]*r + RGB2YCbCr[2,1]*g + RGB2YCbCr[2,2]*b
    return y,u,v

def ycbcr2rgb(y,cb,cr):
    r = YCbCr2RGB[0,0]*y + YCbCr2RGB[0,1]*cb + YCbCr2RGB[0,2]*cr
    g = YCbCr2RGB[1,0]*y + YCbCr2RGB[1,1]*cb + YCbCr2RGB[1,2]*cr
    b = YCbCr2RGB[2,0]*y + YCbCr2RGB[2,1]*cb + YCbCr2RGB[2,2]*cr
    return r,g,b


RGB2YUV = np.array([[0.299, 0.587, 0.114],[-0.14713, -0.28886, 0.436],[0.615, -0.51499, -0.10001]])
YUV2RGB = np.array([[1.0, 0.0, 1.13983],[1.0, -0.39465, -0.58060],[1.0, 2.03211, 0.0]])

def rgb2yuv(r,g,b):
    y = (RGB2YUV[0,0]*r + RGB2YUV[0,1]*g + RGB2YUV[0,2]*b)
    u = (RGB2YUV[1,0]*r + RGB2YUV[1,1]*g + RGB2YUV[1,2]*b)
    v = (RGB2YUV[2,0]*r + RGB2YUV[2,1]*g + RGB2YUV[2,2]*b)
    return y,u,v

def yuv2rgb(y,u,v):
    r = (YUV2RGB[0,0]*y + YUV2RGB[0,1]*u + YUV2RGB[0,2]*v)
    g = (YUV2RGB[1,0]*y + YUV2RGB[1,1]*u + YUV2RGB[1,2]*v)
    b = (YUV2RGB[2,0]*y + YUV2RGB[2,1]*u + YUV2RGB[2,2]*v)
    return r,g,b


def rgb2ycbcr2(r,g,b):
    y = Kr*r + (1-Kr-Kb)*g + Kb*b
    pb = (b-y)/(1-Kb)/2
    pr = (r-y)/(1-Kr)/2
    return y,pb,pr

def ycbcr2rgb2(y,pb,pr):
    r = 2*pr*(1-Kr) + y
    b = 2*pb*(1-Kb) + y
    g = (y - Kr*r - Kb*b)/(1-Kr-Kb)
    return r,g,b


def histeq(im,nbr_bins=256):
   imhist,bins = np.histogram(im.flatten(),nbr_bins)
   cdf = imhist.cumsum() #cumulative distribution function
   cdf = cdf/cdf[-1]
   #use linear interpolation of cdf to find new pixel values
   im2 = np.interp(im.flatten(),bins[:-1],cdf)   
   plt.plot(bins[:-1],imhist)
   plt.plot(bins[:-1],cdf)
   #plt.show()
   plt.clf()   
   plt.plot(np.divide(im2[0:10000],im.flatten()[0:10000]))
   #plt.show()
   plt.clf()   
   return im2.reshape(im.shape), cdf

def plot_scalar_array(ds_str, title_str, mina, maxa):  
    fig = plt.figure(facecolor='k')
    ds = xr.load_dataset(filepath,group='/geophysical_data',mask_and_scale=True)[ds_str].values
    #np.clip(ds, mina, maxa, ds)
    ml = mpl.cm.ScalarMappable(norm=None, cmap='nipy_spectral')
    ml.set_clim(mina, maxa)
    ml.set_array(ds)
    aec = ml.to_rgba(ds, alpha=None)        
    fig.figimage(aec,resize=True)
    loc = plt.LinearLocator(11)
    fig.set_figwidth(fig.get_figwidth()*1.1)
    fig.set_figheight(fig.get_figheight()*1.04)
    axa = fig.add_axes([0.94, 0.1, 0.012, 0.8], frameon=True)                       
    axa.tick_params(colors='c', labelsize=10, length=20, width=2, pad=5, right=True)
    axa.yaxis.set_major_locator(loc)
    tvals = loc.tick_values(mina, maxa)
    fig.colorbar(ml, cax=axa, ticks=tvals )
    fig.suptitle(x + title_str, color='c', size=10, y=0.99) 
#            plt.show()      
    plt.savefig(outpath, facecolor=fig.get_facecolor()) 
    plt.close(fig)
    fig.clf() 

def plot_rgb(filepath):   
    blue = xr.load_dataset(filepath,group='/observations',mask_and_scale=True)['rhot_490'].values
    green = xr.load_dataset(filepath,group='/observations',mask_and_scale=True)['rhot_550'].values
    red = xr.load_dataset(filepath,group='/observations',mask_and_scale=True)['rhot_670'].values
    ro = np.clip(np.nan_to_num(np.rot90(red,2)),0,1.0)
    go = np.clip(np.nan_to_num(np.rot90(green,2)),0,1.0)
    bo = np.clip(np.nan_to_num(np.rot90(blue,2)),0,1.0)
    ao = np.ones_like(ro)
    
    y, u, v = rgb2yuv(ro, go, bo)
    y += 0.001
    y = y/np.max(y)
    yh, Y_cdf = histeq(y)
    maxo_all = max(np.max(ro), np.max(go), np.max(bo)) + 0.001      
    ro /= maxo_all
    go /= maxo_all
    bo /= maxo_all
    
    gamma = 1.0 
    ga = np.ones_like(ro)*gamma   
    rg = np.power(ro, ga)
    gg = np.power(go, ga)
    bg = np.power(bo, ga)
    yg = np.power(yh, ga)
    
    ys = np.divide(yg,y)
    rg = np.multiply(ys,ro)
    gg = np.multiply(ys,go)
    bg = np.multiply(ys,bo)
    
    scale = 255 
    im_red = im.fromarray(np.uint8(np.clip(rg,0,1.0)*scale))
    im_green = im.fromarray(np.uint8(np.clip(gg,0,1.0)*scale))
    im_blue = im.fromarray(np.uint8(np.clip(bg,0,1.0)*scale))
    im_alpha = im.fromarray(np.uint8(np.clip(ao,0,1.0)*scale))
    
    im_rgb = im.merge("RGBA", (im_red, im_green, im_blue, im_alpha)) 
    return np.asarray(im_rgb)/255.0
      
# Start program execution

parser = argparse.ArgumentParser()
parser.add_argument("-id", "--idir", 
                    help="input file", required=True)
parser.add_argument("-dtlf", "--dtlut", type=argparse.FileType('r'), 
                    help="dt lookup table file", required=True)
parser.add_argument("-dblf", "--dblut", type=argparse.FileType('r'),
                    help="db lookup table file", required=True)
parser.add_argument("-al", "--alg", 
                    help="algorithm name", required=True)
parser.add_argument("-l", "--lut", 
                    help="show luts", required=False)
parser.add_argument("-s", "--start", 
                    help="start date and time", required=False)
parser.add_argument("-t", "--side", 
                    help="side -- left, right, or both", required=False)
parser.add_argument("-st", "--step", 
                    help="mask angle step size", required=False)
args = parser.parse_args()

if args.lut:
    lscy = np.zeros((DBNAOT,NDWL,DBNFMF,DBNCHL,DBNWS,DBNSZA,DBNVZA,DBNRAA))
    lscym = np.zeros((DBNAOTM,NDWL,DBNFMFM,DBNCHL,DBNWS,DBNSZA,DBNVZA,DBNRAA))
    lrx = np.zeros((DBNAOT,DBNFMF,DBNCHL,DBNWS,DBNSZA,DBNVZA,DBNRAA))
    lrxm = np.zeros((DBNAOTM,DBNFMFM,DBNCHL,DBNWS,DBNSZA,DBNVZA,DBNRAA))
    lscsy = np.zeros((NAOT,NDWL,NSMALL,NWS,NSZA,NVZA,NRAA))
    lscby = np.zeros((NAOT,NDWL,NBIG,NWS,NSZA,NVZA,NRAA))
    lrsx = np.zeros((NAOT,NSMALL,NWS,NSZA,NVZA,NRAA))
    lrbx = np.zeros((NAOT,NBIG,NWS,NSZA,NVZA,NRAA))
    if (args.alg == "deepblue"):
        dtdb = dbocean(args.dblut.name, 0)
        toa = np.transpose(np.clip(dtdb.lut,0.0001,0.999),(6,7,5,0,1,4,3,2))
        toam = np.transpose(np.clip(dtdb.lutm,0.0001,0.999),(6,7,5,0,1,4,3,2))
        for isz in range(DBNSZA):
            toa[:,:,:,:,:,isz,:,:] *= np.pi/np.cos(dtdb.sza_pts[isz]*D2R)
            toam[:,:,:,:,:,isz,:,:] *= np.pi/np.cos(dtdb.sza_pts[isz]*D2R)
        lrx = np.log10(np.clip(toa[:,D860,:,:,:,:,:,:],0.0001,0.999))
        lrxm = np.clip(np.log10(np.clip(toam[:,D860,:,:,:,:,:,:],0.0001,0.999)),-2.5,-0.1)
        for iw in range(NDWL):
            lscy[:,iw,:,:,:,:,:,:] = np.log10(toa[:,D490,:,:,:,:,:,:]/toa[:,iw,:,:,:,:,:,:])
            lscym[:,iw,:,:,:,:,:,:] = np.clip(np.log10(toam[:,D490,:,:,:,:,:,:]/toam[:,iw,:,:,:,:,:,:]),-0.1,4.0)
    elif (args.alg == "darktarget"):
        dtdb = dtocean(args.dtlut.name, 0)
        toas = np.transpose(np.clip(dtdb.toas,0.0001,0.999),(6,4,5,0,1,2,3))
        toab = np.transpose(np.clip(dtdb.toab,0.0001,0.999),(6,4,5,0,1,2,3))
        toar = np.transpose(np.clip(dtdb.toar,0.0001,0.999),(4,0,1,2,3))
        for ism in range(NSMALL):
            toas[0,:,ism,:,:,:,:] = toar
        for ibm in range(NBIG):
            toab[0,:,ibm,:,:,:,:] = toar
        for isz in range(NSZA):
            units = np.pi/np.cos(dtdb.sza_pts[isz]*D2R)
            toas[:,:,:,:,isz,:,:] *= np.pi/np.cos(dtdb.sza_pts[isz]*D2R)
            toab[:,:,:,:,isz,:,:] *= np.pi/np.cos(dtdb.sza_pts[isz]*D2R)
        lrsx = np.log10(np.clip(toas[:,D860,:,:,:,:,:],0.0001,0.999))
        lrbx = np.log10(np.clip(toab[:,D860,:,:,:,:,:],0.0001,0.999))
        for iw in range(NDWL):
            lscsy[:,iw,:,:,:,:,:] = np.log10(toas[:,D490,:,:,:,:,:]/toas[:,iw,:,:,:,:,:])
            lscby[:,iw,:,:,:,:,:] = np.log10(toab[:,D490,:,:,:,:,:]/toab[:,iw,:,:,:,:,:])
    else:
        print ("invalid algorithm")
        sys.exit()  
else:
    if (args.alg == "deepblue"):
        dtdb = dbocean(args.dblut.name, 0)
    elif (args.alg == "darktarget"):
        dtdb = dtocean(args.dtlut.name, 0)
 

dtdb_dirpath = args.idir
dtdb_dircontents = os.listdir(dtdb_dirpath)
dtdb_dircontents.sort()

for x in dtdb_dircontents:

    if x[0] == ".":
        continue
    
    sstr = ""
    if x.find("DT"): 
        sstr = "_dt"
    elif x.find("DB"): 
        sstr = "_db"
        
    filepath = dtdb_dirpath + "/" + x

    if not os.path.isfile(filepath):
        continue
    
    st = x.find("20")
    if st < 0:
        continue
    if x[st+8] != "T":
        continue
    if x[st:st+15] < args.start:
        continue
    ct = [3200,3200]
    st = [0,0]    
    vin = input(filepath, st, ct)
    
    D2R = np.pi/180.0
    R2D = 180.0/np.pi
    sa = 180.0 - np.arccos(np.cos(vin.sza*D2R)*np.cos(vin.vza*D2R) \
						+ np.sin(vin.sza*D2R)*np.sin(vin.vza*D2R) \
						* np.cos((vin.raa)*D2R))*R2D;
    szalo=0.0
    szahi=80.0
    raalo = -180
    raahi = 180
    maxz = np.max(vin.sza)
    minz = np.min(vin.sza)
    if maxz < szalo:
        continue
    if minz > szahi:
        continue
    if maxz < szahi:
        szahi = maxz
    if minz > szalo:
        szalo = minz
    condszalo = vin.sza < szalo
    condszahi = vin.sza > szahi
    condraalo = vin.raa < raalo
    condraahi = vin.raa > raahi
    
    image_path = dtdb_dirpath + "/IMAGES" 
    if not os.path.isdir(image_path):
        os.mkdir(image_path)
    if args.side == "left":
        image_path1 = image_path + "/histograms_sa_left" 
    elif args.side == "right":
        image_path1 = image_path + "/histograms_sa_right" 
    else:
        image_path1 = image_path + "/histograms_sa" 
    if (args.alg == "deepblue"):
        outfilename1 = x + "_hist_" 
    if (args.alg == "darktarget"):
        outfilename1 = x + "_hist_" 
    if not os.path.isdir(image_path1):
        os.mkdir(image_path1)
    outpath1 = image_path1 + "/" + outfilename1
    image_path2 = image_path + "/maps_sa" 
    outfilename2 = x + "_map_" 
    if not os.path.isdir(image_path2):
        os.mkdir(image_path2)
    outpath2 = image_path2 + "/" + outfilename2
    cmap = cm.get_cmap('nipy_spectral')

    nbins = 100
    hmin = 0
    hmax = 10
    max0 = 0.0
    min0 = -2.25
    max1 = 1.5
    min1 = -1.0
    fmf = 0
    
    hist = np.zeros((NWL, nbins, nbins))
    g0 = np.ones_like(vin.sza)*10
    cgeo = np.ones_like(vin.sza)
    
    lrfl865 = np.log10(np.clip(vin.rfl[:,:,W865], a_min=0.0001, a_max=1.0))    
    lrfl490 = np.log10(np.clip(vin.rfl[:,:,W490]*cgeo, a_min=0.0001, a_max=1.0))

    step = float(args.step)  
    for vlo in range(int(180.0/step)):
        salo = vlo*step
        condsalo = sa < salo
        sahi = salo+step
        condsahi = sa > sahi
        if np.max(sa) < salo or np.min(sa) > sahi:
            continue
        print ("Processing scattering angle range: " + str(salo))        
        fig1 = plt.figure(figsize=(24,15))
        fig2 = plt.figure(figsize=(16,10))
        
        irgb = plot_rgb(filepath)
        ax2 = fig2.add_subplot(2,4,1) 
        ax2.set_title(x)
        plt.xticks([])
        plt.yticks([])
        plt.imshow(irgb)
        
        ax2 = fig2.add_subplot(2,4,5) 
        ax2.set_title("log10(RHOT_"+wlstr[W865]+")")
        ml = mpl.cm.ScalarMappable(norm=None, cmap=cmap)
        ml.set_clim(min0, max0)
        lrflx=np.rot90(lrfl865,2)
        ml.set_array(lrflx)
        aec = ml.to_rgba(lrflx, alpha=None) 
        plt.xticks([])
        plt.yticks([])
        plt.imshow(aec)       
        bwl2 = True
      
        for wl in range(3, NWL):
            
            wl1 = wl-2
            if (wl<W1240): 
                wl2 = wl-1
            else: 
                wl2 = wl
                
 #           print ("Processing wavelength: " + str(wl1))        
    
            ax1 = fig1.add_subplot(2,3,wl1)           
            thrsh = 4*np.log10(float(wlstr[wl])/float(wlstr[W490]))
            thrsh0 = 0.09
            k = (thrsh-thrsh0)/thrsh0 + 0.0001
            r = 2.05
            c = 1.7
            t = 1.4
            h = -thrsh0
            ax1.set_ylim([min1,max1])
            F = c*thrsh/(1.0 + k*np.exp(r*lrfl865+t)) + h
            if args.lut:
                Flut = c*thrsh/(1.0 + k*np.exp(r*lrxm+t)) + h
            
    #       cgeo = dtdb.interp_geom(vin.rfl,vin.sza,g0,g0,vin.wnd,wl1,fmf,0)/dtdb.interp_geom(vin.rfl,vin.sza,vin.vza,vin.raa,vin.wnd,wl1,fmf,0)
            lrflwl = np.log10(np.clip(vin.rfl[:,:,wl]*cgeo, a_min=0.0001, a_max=0.999))
            
            lrflwlc = lrfl490 - lrflwl - F
      
            lmwlc = np.ma.filled(np.ma.masked_where(condsalo, lrflwlc))
            lmwlc = np.ma.filled(np.ma.masked_where(condsahi, lmwlc))
            lm865 = np.ma.filled(np.ma.masked_where(condsalo, lrfl865))
            lm865 = np.ma.filled(np.ma.masked_where(condsahi, lm865))

            try:
    #            hist[wl], edges = np.histogramdd((np.ma.notmasked_contiguous(lm865),np.ma.notmasked_contiguous(lmwlc)), range=[[min0,max0],[min1,max1]],bins=nbins, density=True) 
                if args.side == "both":           
                    hist[wl], edges = np.histogramdd((lm865.flatten(),lmwlc.flatten()), range=[[min0,max0],[min1,max1]],bins=nbins, density=True) 
                if args.side == "left":           
                    hist[wl], edges = np.histogramdd((lm865[:,:1600].flatten(),lmwlc[:,:1600].flatten()), range=[[min0,max0],[min1,max1]],bins=nbins, density=True) 
                if args.side == "right":           
                    hist[wl], edges = np.histogramdd((lm865[:,-1600:].flatten(),lmwlc[:,-1600:].flatten()), range=[[min0,max0],[min1,max1]],bins=nbins, density=True) 
            except Exception as inst:
                print(type(inst))
                print(inst) 
                print ("Numpy.histogramdd failure ... exiting")
                sys.exit()
    
            xpos0, ypos0 = np.meshgrid(edges[0][:-1]+edges[0][1:], edges[1][:-1]+edges[1][1:])
            cp1 = ax1.contourf(xpos0/2, ypos0/2, np.clip(hist[wl].T, a_min=hmin, a_max=hmax), levels=100, cmap=cmap)
            
            thrsh = 4*np.log10(float(wlstr[wl])/float(wlstr[W490]))
            thrsh0 = 0.09
            k = (thrsh-thrsh0)/thrsh0 + 0.0001
            r = 2.05
            c = 1.7
            t = 0.7
            h = -thrsh0
            ax1.set_ylim([min1,max1])
    #        f = c*thrsh/(1.0 + k*np.exp(r*(xpos0[0]/2+t))) + h;
            f = np.zeros(nbins)
            ax1.plot(xpos0[0]/2,f,color='white')        
    
            if args.lut and wl1>2:
                if (args.alg == "deepblue"):
                    colors = [cm.rainbow(x) for x in np.linspace(0, 1, DBNSZA)]
                    for ic in range(DBNSZA):
    #                    ax1.scatter(lrxm[0,0,0,0,ic,0:12,:],lscym[0,wl1,0,0,0,ic,0:12,:],marker='.',color=colors[ic],s=4.0,label=dtdb.sza_pts[ic])
    #                    ax1.scatter(lrx[0:6:5,13,0,0,ic,0:12,:],lscy[0:6:5,wl1,13,0,0,ic,0:12,:],marker='.',color=colors[ic],s=4.0,alpha=0.5,label=dtdb.sza_pts[ic])
    #                    ax1.scatter(lrx[0:6:5,0,0,0,ic,0:12,:],lscy[0:6:5,wl1,0,0,0,ic,0:12,:],marker='.',color=colors[ic],s=4.0,alpha=0.5,label=dtdb.sza_pts[ic])
                        ax1.scatter(lrxm[0,7,0,0,ic,:,:],lscym[0,wl1,7,0,0,ic,:,:]-Flut[0,7,0,0,ic,:,:],marker='.',color=colors[ic],s=4.0,alpha=0.5)
    #                    ax1.scatter(lrxm[0:2,0,0,0,ic,0:12,:],lscym[0:2,wl1,0,0,0,ic,0:12,:],marker='.',color=colors[ic],s=4.0,alpha=0.5)
    #                ax1.legend(markerscale=10,title="SZA")
                elif (args.alg == "darktarget"):
                    colors = [cm.rainbow(x) for x in np.linspace(0, 1, NSZA)]
                    for ic in range(NSZA):
                        ax1.scatter(lrsx[1:3:2,:,0,ic,0:8,:],lscsy[1:3:2,wl1,:,0,ic,0:8,:],marker='.',color=colors[ic],s=2.0,alpha=0.5)
    #                    ax1.scatter(lrbx[1:3:2,:,0,ic,0:8,:],lscby[1:3:2,wl1,:,0,ic,0:8,:],marker='.',color=colors[ic],s=2.0,alpha=1.0)
                        ax1.scatter(lrsx[0,0,0,ic,0:8,:],lscsy[0,wl1,0,0,ic,0:8,:],marker='.',color=colors[ic],s=2.0,alpha=1.0,label=dtdb.sza_pts[ic])
                    ax1.legend(markerscale=10,title="SZA")
    
    #        if thrsh<2.5:
    #            l1 = ax1.axhline(y=thrsh,color='white', ls=':') 
    #        l1 = ax1.axhline(y=0,color='white', ls=':') 
    #        l1.set_label('Rayleigh Threshold')
            ax1.set_title("log10(RHOT_"+wlstr[W490]+" / RHOT_"+wlstr[wl]+")")
            if wl1>3:
                ax1.set_xlabel("log10(RHOT_"+wlstr[W865]+")")
            if wl1==1 or wl1==4:
                ax1.set_ylabel("log10(RHOT_"+wlstr[W490]+" / RHOT_"+wlstr[wl]+")")
                bticks = True
            else:
                bticks = False
                
            ax1.tick_params(labelleft=bticks, left=bticks)
            ax1.grid(visible=True, which='major', color='#666666', linestyle='-')
            ax1.minorticks_on()
            ax1.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            
            if (bwl2):
                ax2 = fig2.add_subplot(2,4,wl2)           
                ax2.set_title("log10(RHOT_"+wlstr[W490]+" / RHOT_"+wlstr[wl]+")")
    #            cp2 = ax2.contourf(lrflwlc, levels=np.arange(-1.0,1.0,0.1), cmap=cmap)
                ml.set_clim(min1, max1)
                lrflx=np.rot90(lmwlc,2)
                ml.set_array(lrflx)
                aec = ml.to_rgba(lrflx, alpha=None) 
                plt.xticks([])
                plt.yticks([])
                plt.imshow(aec) 
                  
        loc = plt.LinearLocator(7)
        tvals = loc.tick_values(min1, max1)
        ml = mpl.cm.ScalarMappable(norm=None, cmap='nipy_spectral')
        ml.set_clim(min1, max1)
        fig2.subplots_adjust(bottom=0.07, top=0.93, left=0.05, right=0.92, wspace=0.02, hspace=0.0)
        cb_ax = fig2.add_axes([0.93, 0.15, 0.02, 0.7])
        fig2.colorbar(ml, cax=cb_ax, ticks=tvals)
        title = x + "   Scattering Angle:< " + str(salo) + ", " + str(sahi) + " >  " + args.side
        fig2.suptitle(title) 
    #        fig2.tight_layout()
    #        fig2.show()
        plt.savefig(outpath2+"_"+str(salo)+".png", facecolor=fig2.get_facecolor()) 
        plt.close(fig2)
        fig2.clf() 
        
    #        fig1.colorbar(cp1)
        loc = plt.LinearLocator(11)
        tvals = loc.tick_values(hmin, hmax)
        ml = mpl.cm.ScalarMappable(norm=None, cmap='nipy_spectral')
        ml.set_clim(hmin, hmax)
        fig1.subplots_adjust(bottom=0.07, top=0.92, left=0.05, right=0.92, wspace=0.05, hspace=0.15)
        cb_ax = fig1.add_axes([0.93, 0.15, 0.02, 0.7])
        fig1.colorbar(ml, cax=cb_ax, ticks=tvals)
        fig1.suptitle(title)
    #        fig1.tight_layout()
    #        fig1.show()
        plt.savefig(outpath1+"_"+str(salo)+".png", facecolor=fig1.get_facecolor()) 
        plt.close(fig1)
        fig1.clf() 
    print (x)
    continue
        
   
