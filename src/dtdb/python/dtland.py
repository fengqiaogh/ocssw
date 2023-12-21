#!/usr/bin/env python
# encoding: utf-8

import os 
import sys
import math
import argparse
import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib import *
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
import Tasmanian

W470 = 0
W550 = 1
W659 = 2
W860 = 3
W124 = 4
W164 = 5
W213 = 6

D213 = 3
slope_644 = 0.559
yint_644 = 0.0
slope_466 = 0.645
yint_466 = 0.0

ISMALL = 0
IBIG = 1
NSEASONS = 4
NLATS = 180
NLONS = 360
NLWL = 4
NWL =7
NSZA = 11
NVZA = 15
NRAA = 16
NTAU = 7
NAOT = 6
NTAB = 5
PRESSURE_P0 = 1013.0
D2R = np.pi/180.0

class dtland(object):
   
    def __init__(self, lut_filepath, mode):
        print ("Reading Darktarget LUT: " + lut_filepath)        
        self.all = np.zeros((NSEASONS,NLATS,NLONS))
        self.int = np.zeros((NTAB,NLWL,NTAU,NSZA,NVZA,NRAA))
        self.t = np.zeros((NTAB,NLWL,NTAU,NSZA,NVZA))
        self.fd = np.zeros((NTAB,NLWL,NTAU,NSZA))
        self.sbar = np.zeros((NTAB,NLWL,NTAU,NSZA))
        self.opth = np.zeros((NTAB,NLWL,NTAU))
        self.massc = np.zeros((NTAB,NLWL,NTAU))
        self.extn = np.zeros((NTAB,NLWL,NTAU))
        self.ssa = np.zeros((NTAB,NLWL,NTAU))
        self.qext = np.zeros((NTAB,NLWL,NTAU))
        self.bext = np.zeros((NTAB,NLWL,NTAU))
        self.vext = np.zeros((NTAB,NLWL,NTAU))
        self.mext = np.zeros((NTAB,NLWL,NTAU))
        self.opth = np.zeros((NTAB,NLWL,NTAU))
        try:
            self.wl_pts=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['WAV_NL'].values 
            self.raa_pts=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['PHI_NL'].values 
            self.vza_pts=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['THE_NL'].values 
            self.sza_pts=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['THET0_NL'].values 
            self.mu0_pts=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['MU0_NL'].values 
            self.all=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['AEROSOL_ALL'].values 
            self.opth=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['OPTH_NL0'].values 
            self.int=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['INT_NL0'].values 
            self.t=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['T_NL0'].values 
            self.fd=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['Fd_NL0'].values 
            self.sbar=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['SBAR_NL0'].values 
            self.massc=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['MASSCOEF_NL0'].values 
            self.extn=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['EXTNORM_NL0'].values 
            self.ssa=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['SSA_NL0'].values 
            self.qext=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['QEXT_NL0'].values 
            self.bext=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['BEXT_NL0'].values 
            self.vext=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['VEXT_NL0'].values 
            self.mext=xr.load_dataset(lut_filepath,group='/LAND AEROSOL')['MEXT_NL0'].values 
        except:
            print ("Unable to read LUT file ... exiting")
            sys.exit()
        
        self.int = np.transpose(self.int, (3,4,5,1,0,2))
        self.t = np.transpose(self.t, (3,4,1,0,2))
        self.fd = np.transpose(self.fd, (3,1,0,2))
        self.sbar = np.transpose(self.sbar, (3,1,0,2))
        self.pts3 = (self.sza_pts, self.vza_pts, self.raa_pts)
        self.pts2 = (self.sza_pts, self.vza_pts)
        self.pts1 = (self.sza_pts,)

    def process(self,rfl,sza,vza,raa,height,mtab):
        xit3 = np.stack((sza, vza, raa))
        xit2 = np.stack((sza, vza))
        xit1 = [[sza]]
# interpolate angles
        int_nl = interp.interpn(self.pts3, self.int, xit3)[0]
        t_nl = interp.interpn(self.pts2, self.t, xit2)[0]
        fd_nl = interp.interpn(self.pts1, self.fd, xit1)[0]
        sbar_nl = interp.interpn(self.pts1, self.sbar, xit1)[0]
        fdt_nl = np.multiply(t_nl,fd_nl)
# interpolate elevation        
        eqwav_nl = np.zeros((NLWL))
        pressure = np.exp(-(height/7.5))
        for iw in range(3):
            wl = self.wl_pts[iw]
            expfactor = -4.15 + (0.2*wl)
            rod_pres = 0.0088*pressure*np.power(wl, expfactor)
            lambda1 = 0.1
            lambda2 = 2.0
            diff0 = 99.0
            while diff0 > 0.00001:
                lambda0 = (lambda1 + lambda2) / 2.0
                ftau0 = 0.0088*np.power(lambda0,-4.15+0.2*lambda0)
                ftau1 = 0.0088*np.power(lambda1,-4.15+0.2*lambda1)
                ftau2 = 0.0088*np.power(lambda2,-4.15+0.2*lambda2)
                if ((ftau1 > rod_pres) and (ftau2 < rod_pres)):
                    if (ftau0 > rod_pres): 
                        lambda1 = (lambda1 + lambda2)/2.0
                    else:
                        lambda2 = (lambda1 + lambda2)/2.0
                diff0 = np.abs(ftau0 - rod_pres)
            eqwav_nl[iw] = np.log(lambda0);
        logwl = (np.log(self.wl_pts),) 
        int_nl[:-1] = np.exp(interp.interpn(logwl, np.log(int_nl), \
                 (eqwav_nl[:-1],), bounds_error=False, fill_value=None))    
        t_nl[:-1] = np.exp(interp.interpn(logwl, np.log(t_nl), \
                 (eqwav_nl[:-1],), bounds_error=False, fill_value=None))    
        fd_nl[:-1] = np.exp(interp.interpn(logwl, np.log(fd_nl), \
                 (eqwav_nl[:-1],), bounds_error=False, fill_value=None))    
        sbar_nl[:-1] = np.exp(interp.interpn(logwl, np.log(sbar_nl), \
                 (eqwav_nl[:-1],), bounds_error=False, fill_value=None))    
        fdt_nl[:-1] = np.exp(interp.interpn(logwl, np.log(fdt_nl), \
                 (eqwav_nl[:-1],), bounds_error=False, fill_value=None))    
# simulate toa
        rho = np.zeros((NLWL,NTAB,NTAU))
        rm213 = np.ones((NTAB,NTAU))*rfl[W213]
        rho[D213] = np.divide(int_nl[D213]-rm213, \
                   np.multiply(sbar_nl[D213],(int_nl[D213]-rm213))-fdt_nl[D213])
        rho[W659] = slope_644*rho[D213] + yint_644
        rho[W470] = slope_466*rho[W659] + yint_466
        rho[W550].fill(0)
        rho_star = int_nl + np.divide(np.multiply(fdt_nl,rho), \
                           (np.ones_like(rho)-np.multiply(sbar_nl,rho)))

        tauxs = np.zeros(NSMALL)
        tauxb = np.zeros(NBIG)
        for ism in range(NSMALL):
            tauxs[ism] = np.interp(rfl[W860], mrfls[W860,ism], self.taus[W550,ism])
        for ibm in range(NBIG):
            tauxb[ibm] = np.interp(rfl[W860], mrflb[W860,ibm], self.taub[W550,ibm])
                
        trfls = np.zeros((NWL,NSMALL))
        trflb = np.zeros((NWL,NBIG))
        for iwl in range(NWL):
            for ism in range(NSMALL):
                trfls[iwl,ism] = np.interp(tauxs[ism],self.taus[W550,ism],mrfls[iwl,ism])
            for ibm in range(NBIG):
                trflb[iwl,ibm] = np.interp(tauxb[ibm],self.taub[W550,ibm],mrflb[iwl,ibm])
        
        aot = np.zeros((NSMALL,NBIG))
        fmf = np.zeros((NSMALL,NBIG))
        sse = np.zeros((NSMALL,NBIG))
        for ism in range(NSMALL):
            for ibm in range(NBIG):
                denom = rfl - mrflr + 0.01
                rbdif = np.divide((rfl - trflb[:,ibm]),denom)
                sbdif = np.divide((trfls[:,ism] - trflb[:,ibm]),denom)
                xm = np.dot(rbdif,sbdif)/np.dot(sbdif,sbdif)
                xm = np.max((np.min((xm,1.0)),0.0))
                mrfl = xm*trfls[:,ism] + (1.0-xm)*trflb[:,ibm]
                err = np.divide((rfl - mrfl),denom)
                sse[ism,ibm] = np.dot(err,err)
                aot[ism,ibm] = xm*tauxs[ism] + (1.0-xm)*tauxb[ibm]
                fmf[ism,ibm] = xm
        im = divmod(sse.argmin(), sse.shape[1])
        self.sse = sse[im[0],im[1]]
        self.aot = aot[im[0],im[1]]
        self.fmf = fmf[im[0],im[1]]
        self.rfl = rfl
        self.mrfl = self.fmf*trfls[:,im[0]] + (1.0-self.fmf)*trflb[:,im[1]]
        self.rsd = self.mrfl - self.rfl
        self.sse = np.dot(self.rsd,self.rsd)/(NWL-2)
        self.rayl = mrflr
        
        return self.fmf, self.aot, self.sse
    
    def plot(self, iy, ix):    
        plt.clf()
        plt.grid(True)
        plt.plot(self.wl_pts, self.rfl, marker='.', color='b', label='measured')
        plt.plot(self.wl_pts, self.mrfl, marker='.', color='g', label='modeled')
        plt.plot(self.wl_pts, self.rsd, marker='.', color='r', label='residual')
        plt.xlabel('wavelength (nm)')
        plt.ylabel('reflectance')
        tstr = "dtland -- y={3:}, x={4:}  aot: {0:.3f}  fmf: {1:.3f}  sse: {2:.3}"
        plt.title(tstr.format(self.aot, self.fmf, self.sse, iy, ix))
        plt.legend(loc='upper right')
        plt.show()
        
    def plot_points(self):    
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')        
        for s in small:
            xs = toast[W550,s,:]
            ys = toast[W659,s,:]
            zs = toast[W860,s,:]
            ax.plot(xs, ys, zs, c="b", label="small model")
            ax.scatter(xs, ys, zs, c="b", marker="o")            
        for b in big:
            xb = toabt[W550,b,:]
            yb = toabt[W659,b,:]
            zb = toabt[W860,b,:]
            ax.plot(xb, yb, zb, c="r", label="big model")
            ax.scatter(xb, yb, zb, c="r", marker="o")        
        for s in small:
            for b in big:
                xp = bspair[:,W550,s,b]
                yp = bspair[:,W659,s,b]
                zp = bspair[:,W860,s,b]
                ax.plot(xp, yp, zp, c="g", label="big/small continuum")

        xs = trfls[W550,:]
        ys = trfls[W659,:]
        zs = trfls[W860,:]
        ax.scatter(xs, ys, zs, c="b", marker="o", s=50)        
        xb = trflb[W550,:]
        yb = trflb[W659,:]
        zb = trflb[W860,:]
        ax.scatter(xb, yb, zb, c="r", marker="o", s=50)        
        xm = rfl[W550,row,col]
        ym = rfl[W659,row,col]
        zm = rfl[W860,row,col]
        ax.scatter(xm, ym, zm, c="k", marker="o", s=50)        
        ax.set_xlabel('W550')
        ax.set_ylabel('W659')
        ax.set_zlabel('W860')        
        #plt.show()
                
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')        
        ws = np.arange(NWL)
        ts = np.arange(NAOT+1)
        ws, ts = np.meshgrid(ws, ts)
        ls = np.zeros((NWL,NAOT+1))
        lb = np.zeros((NWL,NAOT+1))        
        ls[:,:-1] = toast[:,s,:]
        ls[:,NAOT] = ls[:,NAOT-1]
        ax.plot_wireframe(ws, ts, ls, color='b', label="small model")            
        lb[:,:-1] = toabt[:,b,:]
        lb[:,NAOT] = lb[:,NAOT-1]
        ax.plot_wireframe(ws, ts, lb, color='r', label="big model")            
        plt.show()

class input(object):
   
    def __init__(self, l1b_filepath):
        self.ifile = l1b_filepath
        print ("Reading VIIRS Data: " + self.ifile)
        try:
            self.rfl = np.append(xr.load_dataset(l1b_filepath,group='/reflectance')['toa_reflectance'][1:6,:,:].values, 
                            xr.load_dataset(l1b_filepath,group='/reflectance')['toa_reflectance'][7:,:,:].values,axis=0)
            self.sza = xr.load_dataset(l1b_filepath,group='/geolocation')['solar_zenith'].values
            self.vza = xr.load_dataset(l1b_filepath,group='/geolocation')['sensor_zenith'].values
            self.raa = xr.load_dataset(l1b_filepath,group='/geolocation')['relative_azimuth'].values
            self.height = xr.load_dataset(l1b_filepath,group='/geolocation')['elevation'].values
            self.lat = xr.load_dataset(l1b_filepath,group='/navigation_data')['latitude'].values
            self.lon = xr.load_dataset(l1b_filepath,group='/navigation_data')['longitude'].values
            self.cld = xr.load_dataset(l1b_filepath,group='/ancillary')['cloud_mask'].values
            self.wnd = xr.load_dataset(l1b_filepath,group='/ancillary')['wind_speed'].values
            self.rfl = np.transpose(self.rfl, (1,2,0))
        except:
            print ("Unable to read from input file ... exiting")
            sys.exit()


class output(object):
   
    def __init__(self, out_filepath, ydim, xdim):
        self.ofile = out_filepath
        self.ydim = ydim
        self.xdim = xdim
        try:
            self.rfl = xr.DataArray(np.zeros((NWL,ydim,xdim)),dims=('wl','y', 'x'))
            self.lat = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.lon = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.sza = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.vza = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.raa = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.wnd = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.aot = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.fmf = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.sse = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.ds = xr.Dataset({'rfl': self.rfl, 'lat': self.lat, 'lon': self.lon, 
                                  'sza': self.sza, 'vza': self.vza, 'raa': self.raa, 
                                  'wnd': self.wnd, 'aot': self.aot, 
                                  'fmf': self.fmf, 'sse': self.sse })
        except:
            print ("Unable to initialize output file ... exiting")
            sys.exit()

    def write(self):
        print ("Writing to file: " + self.ofile)
        self.ds.to_netcdf(self.ofile)
    
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ifile", type=argparse.FileType('r'), 
                        help="input file", required=True)
    parser.add_argument("-o", "--ofile", type=argparse.FileType('w'), 
                        help="output file", required=False)
    parser.add_argument("-l", "--lut", type=argparse.FileType('r'), 
                        help="lookup table file", required=True)
    parser.add_argument("-p", "--plot", type=bool, default=False,
                        help="plot pixel data", required=False)
    parser.add_argument("-m", "--mode", type=int, default=0,
                        help="mode option", required=False)
    parser.add_argument("y", type=int, help="start line")
    parser.add_argument("x", type=int, help="start pixel")
    parser.add_argument("z", type=int, help="square side ")
    
    args = parser.parse_args()
    vin = input(args.ifile.name)    
    vout = output(args.ofile.name, args.z, args.z)    
    dtl = dtland(args.lut.name, args.mode)
   
    print ("Processing mode {0}".format(args.mode))
    for iy in range(args.y, args.y+args.z):
        for ix in range(args.x, args.x+args.z):

            ilon = 180 + np.round(vin.lon + 0.5)
            ilat = 90 - np.round(vin.lat + 0.5) 
#            mtab = np.array([int(dtl.all[1,ilat,ilon]+1), 4])           
            mtab = np.array([1, 4])           

            if(vin.cld[iy,ix]):
                fmf,aot,sse = -999.9, -999.9, -999.9
            else:
                fmf,aot,sse = dtl.process(vin.rfl[iy,ix,:],vin.sza[iy,ix], 
                                          vin.vza[iy,ix],vin.raa[iy,ix], 
                                          vin.height[iy,ix]/1000.0,mtab)    
                if(args.plot): 
                    dtl.plot(iy,ix)

            vout.ds.rfl[:,iy-args.y,ix-args.x] = vin.rfl[iy,ix]
            vout.ds.lat[iy-args.y,ix-args.x] = vin.lat[iy,ix]
            vout.ds.lon[iy-args.y,ix-args.x] = vin.lon[iy,ix]
            vout.ds.sza[iy-args.y,ix-args.x] = vin.sza[iy,ix]
            vout.ds.vza[iy-args.y,ix-args.x] = vin.vza[iy,ix]
            vout.ds.raa[iy-args.y,ix-args.x] = vin.raa[iy,ix]
            vout.ds.wnd[iy-args.y,ix-args.x] = vin.wnd[iy,ix]
            vout.ds.aot[iy-args.y,ix-args.x] = aot
            vout.ds.fmf[iy-args.y,ix-args.x] = fmf
            vout.ds.sse[iy-args.y,ix-args.x] = sse
            
        print('.', end = '')
    print('  Done!')
    vout.write()

if __name__ == "__main__":

    main()