#!/usr/bin/env python
# encoding: utf-8

import os 
import sys
import math
import time
import argparse
import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib import *
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime

D490 = 0
D550 = 1
D670 = 2
D860 = 3
D124 = 4
D164 = 5
D213 = 6
W860 = D860+2
NDWL = 7
NWS = 4
NSZA = 11
NVZA = 16
NRAA = 16
NAOT = 6
NBIG = 5
NSMALL = 4
D2R = np.pi/180.0

class dtocean(object):
   
    def __init__(self, lut_filepath, mode):
        print ("Reading Darktarget LUT: " + lut_filepath)        
        self.toas = np.zeros((NDWL,NSMALL,NAOT,NWS,NSZA,NVZA,NRAA))
        self.toab = np.zeros((NDWL,NBIG,NAOT,NWS,NSZA,NVZA,NRAA))
        self.toar = np.zeros((NDWL,NWS,NSZA,NVZA,NRAA))
        self.taus = np.zeros((NDWL,NSMALL,NAOT))
        self.taub = np.zeros((NDWL,NBIG,NAOT))
        try:
            self.toas=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['AINTS'].values
            self.toab=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['AINTB'].values 
            self.toar=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['REF_RAYALL'].values 
            self.taus=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['TAUAS'].values 
            self.taub=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['TAUAB'].values 
            self.wl_pts=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['WAVE'].values
            self.raa_pts=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['PHC'].values
            self.vza_pts=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['THET'].values
            self.sza_pts=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL')['THET0'].values 
            self.taus = np.transpose(self.taus,(0,1,2))
            self.taub = np.transpose(self.taub,(0,1,2))
            self.toas = np.transpose(self.toas,(3,4,5,6,0,1,2))
            self.toab = np.transpose(self.toab,(3,4,5,6,0,1,2))
            self.toar = np.transpose(self.toar,(1,2,3,4,0))
            self.wnd_pts = np.array([0.0, 6.0, 10.0, 16.0])
            self.tau_pts = (self.taus[D550,0,:])                                
            self.pts = (self.wnd_pts, self.sza_pts, self.vza_pts, self.raa_pts)
        except Exception as inst:
            print(type(inst)) 
            print(inst)  
            print ("Unable to read LUT file ... exiting")
            sys.exit()

        self.tslut = np.zeros((NAOT,NWS,NSZA,NVZA,NRAA,NSMALL))
        self.tblut = np.zeros((NAOT,NWS,NSZA,NVZA,NRAA,NBIG))
        self.rfl_pts = np.linspace(0.0,0.2,NAOT)
        self.tpts = (self.rfl_pts, self.wnd_pts, self.sza_pts, self.vza_pts, self.raa_pts)

        try:
            tic = time.perf_counter()
            maxs = np.max(self.toas)
            maxb = np.max(self.toab)
                           
            for iwn in range(NWS):
                for isz in range(NSZA):
                    for ith in range(NVZA):
                        for iph in range(NRAA):
                            for ism in range(NSMALL):
                                self.tslut[:,iwn,isz,ith,iph,ism] = \
                                np.interp(self.rfl_pts,self.toas[iwn,isz,ith,iph,D860,ism,:],self.tau_pts)
                            for ism in range(NBIG):
                                self.tblut[:,iwn,isz,ith,iph,ism] = \
                                np.interp(self.rfl_pts,self.toab[iwn,isz,ith,iph,D860,ism,:],self.tau_pts)
    
            toc = time.perf_counter()
            
        except Exception as inst:
            print(type(inst)) 
            print(inst)  

    def interp_geom(self,rfl,sza,vza,raa,wnd,wl,fmf,im):
        xit = np.stack((wnd, sza, vza, raa),axis=-1)
        units = np.expand_dims(np.pi/np.cos(sza*D2R),axis=(2,3)) 
        unitsr = np.expand_dims(np.pi/np.cos(sza*D2R),axis=2) 
        shp = sza.shape
        trfl = np.zeros_like(sza)
        if fmf==1:
            mrfl = np.multiply(interp.interpn(self.pts, self.toas[:,:,:,:,:,im,:], xit, fill_value=None, bounds_error=False),units)
            mrfl[:,:,:,0] = np.multiply(interp.interpn(self.pts, self.toar, xit, fill_value=None, bounds_error=False),unitsr)
            for iy in range(shp[0]):
                for ix in range(shp[1]):
                    trfl[iy,ix] = np.interp(rfl[iy,ix,W860], mrfl[iy,ix,D860], mrfl[iy,ix,wl])
        else:
            mrfl = np.multiply(interp.interpn(self.pts, self.toab[:,:,:,:,:,im,:], xit, fill_value=None, bounds_error=False),units)
            mrfl[:,:,:,0] = np.multiply(interp.interpn(self.pts, self.toar, xit, fill_value=None, bounds_error=False),unitsr)
            for iy in range(shp[0]):
                for ix in range(shp[1]):
                    trfl[iy,ix] = np.interp(rfl[iy,ix,W860], mrfl[iy,ix,D860], mrfl[iy,ix,wl])
        return trfl
    
    def proc_gran(self,rfl,sza,vza,raa,wnd):
        xit = np.stack((wnd, sza, vza, raa),axis=-1)
        units = np.pi/np.cos(sza*D2R) 
        shp = sza.shape
        
        tstart = time.perf_counter()
        try:
            tic = time.perf_counter()
            mrfls = np.multiply(interp.interpn(self.pts, self.toas, xit),np.expand_dims(units,axis=(2,3,4)))
            mrflb = np.multiply(interp.interpn(self.pts, self.toab, xit),np.expand_dims(units,axis=(2,3,4)))
            mrflr = np.multiply(interp.interpn(self.pts, self.toar, xit),np.expand_dims(units,axis=2))
            toc = time.perf_counter()
        except Exception as inst:
            print(type(inst))
            print(inst) 
        print('\nGeometry interpolation: ', toc-tic, 'sec, ', (toc-tic)/shp[0]/shp[1], "per pixel")
            
        tauxs = np.zeros((shp[0],shp[1],NSMALL))
        tauxb = np.zeros((shp[0],shp[1],NBIG))
        trfls = np.zeros((shp[0],shp[1],NDWL,NSMALL))
        trflb = np.zeros((shp[0],shp[1],NDWL,NBIG))
        AA = np.zeros((shp[0],shp[1],NSMALL,NBIG,7,2))
        tic = time.perf_counter()
        for ism in range(NSMALL):
            mrfls[:,:,:,ism,0] = mrflr
        for ibm in range(NBIG):
            mrflb[:,:,:,ibm,0] = mrflr
        for iy in range(shp[0]):
            for ix in range(shp[1]):

                for ism in range(NSMALL):
                    tauxs[iy,ix,ism] = np.interp(rfl[iy,ix,W860], mrfls[iy,ix,D860,ism], self.taus[D550,ism])
                for ibm in range(NBIG):                    
                    tauxb[iy,ix,ibm] = np.interp(rfl[iy,ix,W860], mrflb[iy,ix,D860,ibm], self.taub[D550,ibm])
                        
                for iwl in range(NDWL):
                    for ism in range(NSMALL):
                        trfls[iy,ix,iwl,ism] = np.interp(tauxs[iy,ix,ism],self.taus[D550,ism],mrfls[iy,ix,iwl,ism])
                    for ibm in range(NBIG):
                        trflb[iy,ix,iwl,ibm] = np.interp(tauxb[iy,ix,ibm],self.taub[D550,ibm],mrflb[iy,ix,iwl,ibm])
                                
                for ism in range(NSMALL):
                    for ibm in range(NBIG):
                        AA[iy,ix,ism,ibm] = np.vstack((trfls[iy,ix,:,ism],trflb[iy,ix,:,ibm])).T

        toc = time.perf_counter()
        print('Per-pixel interpolation: ', toc-tic, 'sec, ', (toc-tic)/shp[0]/shp[1], "per pixel")

        tic = time.perf_counter()
        rfle = np.expand_dims(rfl[2:],axis=(2,3,5))
        try:
            xm = np.clip(np.matmul(np.linalg.pinv(AA),rfle),0,1)[:,:,:,:,0]
            mrfl = np.multiply(xm[:,:,:,:],AA[:,:,:,:,:,0]) + np.multiply((1-xm[:,:,:,:]),AA[:,:,:,:,:,1])
            err = np.squeeze(rfle,axis=5)-mrfl
            sse = np.linalg.norm(err[:,:,:,:,1:],axis=4)**2/(NDWL-2)
        except Exception as inst:
            print(type(inst))
            print(inst) 
        toc = time.perf_counter()
        print('FMF calculation: ', toc-tic, 'sec, ', (toc-tic)/shp[0]/shp[1], "per pixel")
        tstop = time.perf_counter()
        print('Total processing time: ', shp[0], 'lines, ', tstop-tstart, 'sec, ', (tstop-tstart)/shp[0]/shp[1], "per pixel")
        
        rmin = sse.reshape(shp[0],shp[1],NSMALL*NBIG).argmin(axis=2)
        im = divmod(rmin,NBIG)
        self.type = (NBIG-1)*im[0] + im[1]
        Y,X = np.ogrid[:shp[0],:shp[1]]
        self.sse = sse.reshape(shp[0],shp[1],NSMALL*NBIG)[Y,X,rmin]
        self.fmf = xm.reshape(shp[0],shp[1],NSMALL*NBIG)[Y,X,rmin]
        self.aot = np.multiply(self.fmf,tauxs[Y,X,im[0]]) + np.multiply((1.0-self.fmf),tauxb[Y,X,im[1]])

        return self.fmf, self.aot, self.sse, self.type
    
    def process(self,rfl,sza,vza,raa,wnd):
        xit = np.stack((wnd, sza, vza, raa))
        units = np.pi/np.cos(sza*D2R)
        rfld = rfl[2:]
        mrfls = interp.interpn(self.pts, self.toas, xit, bounds_error=False, fill_value=None)[0]*units
        mrflb = interp.interpn(self.pts, self.toab, xit, bounds_error=False, fill_value=None)[0]*units
        mrflr = interp.interpn(self.pts, self.toar, xit, bounds_error=False, fill_value=None)[0]*units
        
        tauxs = np.zeros(NSMALL)
        tauxb = np.zeros(NBIG)
        for ism in range(NSMALL):
            mrfls[:,ism,0] = mrflr
            tauxs[ism] = interp.interpn((mrfls[D860,ism,:],), self.taus[D550,ism,:], [rfld[D860]], bounds_error=False, fill_value=None)
        for ibm in range(NBIG):
            mrflb[:,ibm,0] = mrflr
            tauxb[ibm] = interp.interpn((mrflb[D860,ibm,:],), self.taub[D550,ibm,:], [rfld[D860]], bounds_error=False, fill_value=None)
                
        trfls = np.zeros((NDWL,NSMALL))
        trflb = np.zeros((NDWL,NBIG))
        for iwl in range(NDWL):
            for ism in range(NSMALL):
                trfls[iwl,ism] = interp.interpn((self.taus[D550,ism],), mrfls[iwl,ism], [tauxs[ism]], bounds_error=False, fill_value=None)
            for ibm in range(NBIG):
                trflb[iwl,ibm] = interp.interpn((self.taub[D550,ibm],), mrflb[iwl,ibm], [tauxb[ibm]], bounds_error=False, fill_value=None)
                        
        AA = np.zeros((NSMALL,NBIG,7,2))
        for ism in range(NSMALL):
            for ibm in range(NBIG):
                AA[ism,ibm] = np.vstack((trfls[:,ism],trflb[:,ibm])).T
        
        try:
            xm = np.clip(np.dot(np.linalg.pinv(AA),rfld)[:,:,0],0,1)
            xmd = np.expand_dims(xm,axis=2)
            mrfl = xmd*AA[:,:,:,0] + (1-xmd)*AA[:,:,:,1]
            err = rfld-mrfl
            sse = np.linalg.norm(err[:,:,1:],axis=2)**2/(NDWL-2)
        except Exception as inst:
            print(type(inst))
            print(inst) 
                
        im = divmod(sse.argmin(),NBIG)
        self.type = (NBIG-1)*im[0] + im[1]
        self.sse = sse[im[0],im[1]]
        self.fmf = xm[im[0],im[1]]
        self.aot = self.fmf*tauxs[im[0]] + (1.0-self.fmf)*tauxb[im[1]]
        self.trfls = trfls
        self.trflb = trflb
        self.srfl = trfls[:,im[0]]
        self.brfl = trflb[:,im[1]]
        self.mrfl = self.fmf*trfls[:,im[0]] + (1.0-self.fmf)*trflb[:,im[1]]
        self.rsd = (rfld - self.mrfl)/rfld
        self.rfl = rfld
        self.sse = np.dot(self.rsd,self.rsd)/(NDWL-2)
        
        return self.fmf, self.aot, wnd, self.sse, self.type
    
    def plot(self, iy, ix):    
        plt.clf()
        plt.grid(True)
        plt.plot(self.wl_pts, self.srfl, marker='.', color='b', label='fine mode')
        plt.plot(self.wl_pts, self.brfl, marker='.', color='k', label='coarse mode')
        plt.plot(self.wl_pts, self.trfls, marker='', color='b', label='fine mode', linewidth=1)
        plt.plot(self.wl_pts, self.trflb, marker='', color='k', label='coarse mode', linewidth=1)
        plt.plot(self.wl_pts, self.rfl, marker='.', color='r', label='measured')
        plt.plot(self.wl_pts, self.mrfl, marker='.', color='g', label='modeled')
#        plt.plot(self.wl_pts, self.rsd, marker='.', color='r', label='residual')
        plt.xlabel('wavelength (nm)')
        plt.ylabel('reflectance')
        tstr = "dtocean -- y={3:}, x={4:}  aot: {0:.3f}  fmf: {1:.3f}  sse: {2:.3}"
        plt.title(tstr.format(self.aot, self.fmf, self.sse, iy, ix))
        plt.legend(loc='upper right')
        plt.show()
         
    def plot_points(self):    
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        small = [0,1,2,3]        
        for s in small:
            xs = toast[D550,s,:]
            ys = toast[D670,s,:]
            zs = toast[D860,s,:]
            ax.plot(xs, ys, zs, c="b", label="small model")
            ax.scatter(xs, ys, zs, c="b", marker="o")            
        big = [0,1,2,3,5]        
        for b in big:
            xb = toabt[D550,b,:]
            yb = toabt[D670,b,:]
            zb = toabt[D860,b,:]
            ax.plot(xb, yb, zb, c="r", label="big model")
            ax.scatter(xb, yb, zb, c="r", marker="o")        
        for s in small:
            for b in big:
                xp = bspair[:,D550,s,b]
                yp = bspair[:,D670,s,b]
                zp = bspair[:,D860,s,b]
                ax.plot(xp, yp, zp, c="g", label="big/small continuum")

        xs = trfls[D550,:]
        ys = trfls[D670,:]
        zs = trfls[D860,:]
        ax.scatter(xs, ys, zs, c="b", marker="o", s=50)        
        xb = trflb[D550,:]
        yb = trflb[D670,:]
        zb = trflb[D860,:]
        ax.scatter(xb, yb, zb, c="r", marker="o", s=50)        
        xm = rfl[W550,row,col]
        ym = rfl[W670,row,col]
        zm = rfl[W860,row,col]
        ax.scatter(xm, ym, zm, c="k", marker="o", s=50)        
        ax.set_xlabel('D550')
        ax.set_ylabel('D670')
        ax.set_zlabel('D860')        
        #plt.show()
                
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')        
        ws = np.arange(NDWL)
        ts = np.arange(NAOT+1)
        ws, ts = np.meshgrid(ws, ts)
        ls = np.zeros((NDWL,NAOT+1))
        lb = np.zeros((NDWL,NAOT+1))        
        ls[:,:-1] = toast[:,s,:]
        ls[:,NAOT] = ls[:,NAOT-1]
        ax.plot_wireframe(ws, ts, ls, color='b', label="small model")            
        lb[:,:-1] = toabt[:,b,:]
        lb[:,NAOT] = lb[:,NAOT-1]
        ax.plot_wireframe(ws, ts, lb, color='r', label="big model")            
        plt.show()

