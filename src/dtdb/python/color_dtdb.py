#!/usr/bin/env python
# encoding: utf-8

import os 
import sys
import datetime
import math
from PIL import Image as im
from PIL import ImageDraw
from PIL import ImageFilter
from PIL import ImageMath
from PIL import ImageEnhance
from PIL import ImageOps as imo
from PIL import ImageColor as imc

import viirs_files as vf
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

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

   #get image histogram
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


def pseudocolor(val, minval, maxval):
    # Scale val to be in the range [0, 1]
    val = (val - minval) / (maxval - minval)
    # Return RGBA tuple from nipy_spectral colormap
    return cm.nipy_spectral(val)

def save_scalar_array(ds_str, title_str, mina, maxa):  
    fig = plt.figure(facecolor='k')
    ds = xr.load_dataset(filepath,group='/geophysical_data',mask_and_scale=True)[ds_str].values
    ds = np.nan_to_num(np.rot90(ds,2))
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

def scalar_array(ds_str, title_str, mina, maxa):  
    ds = xr.load_dataset(filepath,group='/geophysical_data',mask_and_scale=True)[ds_str].values
    ds = np.nan_to_num(np.rot90(ds,2))
    ml = mpl.cm.ScalarMappable(norm=None, cmap='nipy_spectral')
    ml.set_clim(mina, maxa)
    ml.set_array(ds)
    aec = ml.to_rgba(ds, alpha=None)        
    return np.asarray(aec)

def plot_rgb():   
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

args = sys.argv

dtdb_dirpath = args[1]
dataset = args[2]

colors = "123"
if len(args) > 3:
    colors = args[3]
    bc = np.int(args[3][0])-1
    gc = np.int(args[3][1])-1
    rc = np.int(args[3][2])-1
else:
    bc = 1
    gc = 2
    rc = 3
    
logfile = ""
if len(args) > 4:
    logfile = args[4]
    command = "date > " + logfile
    result = os.system(command)

dtdb_dircontents = os.listdir(dtdb_dirpath)
dtdb_dircontents.sort()

for x in dtdb_dircontents:

    if x[0] == ".":
        continue
    
    sstr = ""
    if x.find("DT") >= 0: 
        sstr = "_dt"
    elif x.find("DB") >= 0: 
        sstr = "_db"
        
    filepath = dtdb_dirpath + "/" + x

    if not os.path.isfile(filepath):
        continue
    
    if (dataset == "retrievals"):
        image_path = dtdb_dirpath + "/IMAGES" 
        if not os.path.isdir(image_path):
            os.mkdir(image_path)
        image_path = image_path + "/geophysical_data" 
              
        outfilename = x + ".png" 
        path = image_path
        if not os.path.isdir(path):
            os.mkdir(path)
        outpath = path + "/" + outfilename
        
        fig = plt.figure(figsize=(12,13))
        title_str = ' Aerosol Optical Depth at 550 nm'
        ds_str = 'aot_550'
        maxa = 1.0
        mina = 0.0
        ax = fig.add_subplot(2,2,1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title_str) 
        sa = scalar_array(ds_str+sstr, title_str, mina, maxa)
        plt.imshow(sa)
        title_str = ' Fine Mode Fraction'
        ds_str = 'fmf_550'
        maxa = 1.0
        mina = 0.0
        ax = fig.add_subplot(2,2,2) 
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title_str) 
        sa = scalar_array(ds_str+sstr, title_str, mina, maxa)
        plt.imshow(sa)
        title_str = ' Angstrom Exponent'
        ds_str = 'angstrom'
        maxa = 2.0
        mina = 0.0
        ax = fig.add_subplot(2,2,3) 
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title_str) 
        sa = scalar_array(ds_str+sstr, title_str, mina, maxa)
        plt.imshow(sa)
        title_str = '- Log base 10 of Normalized Sum of Squares Error'
        ds_str = 'quality'
        maxa = 3.0
        mina = 0.0
        ax = fig.add_subplot(2,2,4) 
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title_str) 
        sa = scalar_array(ds_str+sstr, title_str, mina, maxa)
        plt.imshow(sa)
        fig.subplots_adjust(bottom=0.07, top=0.93, left=0.05, right=0.92, wspace=0.02, hspace=0.1)
        fig.suptitle(x) 
        plt.savefig(outpath, facecolor=fig.get_facecolor()) 
        plt.close(fig)
        fig.clf() 
          
        print (x)
        continue
       
    try:
        if (dataset == "observations"):
#            refl = xr.load_dataset(filepath,group='/observations')['rhot']
            bluec = xr.load_dataset(filepath,group='/observations',mask_and_scale=True)['rhot_490'].values
            greenc = xr.load_dataset(filepath,group='/observations',mask_and_scale=True)['rhot_550'].values
            redc = xr.load_dataset(filepath,group='/observations',mask_and_scale=True)['rhot_670'].values
        elif (dataset == "aot"):
            ocean = xr.load_dataset(filepath,group='/geophysical_data')['AOT_ocean']
            land = xr.load_dataset(filepath,group='/geophysical_data')['AOT_land']
            land = xr.load_dataset(filepath,group='/geolocation')['land_water']
            bluec = ocean[bc,:,:]
            greenc = ocean[gc,:,:]
            redc = ocean[rc,:,:]
            np.putmask(bluec, land>0, land[0,:,:])
            np.putmask(greenc, land>0, land[1,:,:])
            np.putmask(redc, land>0, land[2,:,:])
        elif (dataset == "cloudmask"):
            cldmsk = xr.load_dataset(filepath,group='/meteorology')['cloud_mask']
            cldtst = xr.load_dataset(filepath,group='/quality')['cloud_test']
            bluec = np.zeros_like(cldmsk, dtype=np.int16)
            greenc = np.zeros_like(cldmsk, dtype=np.int16)
            redc = np.zeros_like(cldmsk, dtype=np.int16)
            nomask = np.zeros_like(cldmsk, dtype=np.int16)
            bluec = cldmsk[:,:]
            np.putmask(greenc[:,:], cldtst[:,:]<0, -cldtst[:,:])
            np.putmask(bluec[:,:], cldtst[:,:]<0, nomask[:,:])
            np.putmask(redc[:,:], cldtst[:,:]>0, cldtst[:,:])
            np.putmask(bluec[:,:], cldtst[:,:]>0, nomask[:,:])
            colors = "rgb"
    except Exception as e:
        print (" Error loading dataset ... exiting")
        continue       

    if (dataset == "observations" or dataset == "aot" or dataset == "cloudmask"): 
        ro = np.clip(np.nan_to_num(redc),0,10.0)
        go = np.clip(np.nan_to_num(greenc),0,10.0)
        bo = np.clip(np.nan_to_num(bluec),0,10.0)
    
        if np.isnan(np.sum(ro)):
            print (np.max(ro), np.max(go), np.max(bo)) 
            continue
        
        y, u, v = rgb2yuv(ro, go, bo)
        y += 0.001
        y = y/np.max(y)
        yh, Y_cdf = histeq(y)
#        u *= 1
#        v *= 1
#        rh, gh, bh = yuv2rgb(yh, u, v)
        
        mino_all = min(np.min(ro), np.min(go), np.min(bo))
    #    ro -= mino_all
    #    go -= mino_all
    #    bo -= mino_all    
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
        
        im_rgb = im.merge("RGB", (im_red, im_green, im_blue))
        
    #    im_rgb.show()   
        set = dataset.partition("/")
        outfilename = x + "_" + set[2] + "_" + colors + ".png" 
        
        image_path = dtdb_dirpath + "/IMAGES" 
        if not os.path.isdir(image_path):
            os.mkdir(image_path)
        image_path = image_path + "/" + set[0]
        if not os.path.isdir(image_path):
            os.mkdir(image_path)
        image_path = image_path + set[1] + set[2] + "_" + colors
        if not os.path.isdir(image_path):
            os.mkdir(image_path)
        
        outpath = image_path + "/" + outfilename
        im_rgb.save( outpath ) 

    print (x)
   
