# cython: language_level=3

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

W470 = 0
W550 = 1
W659 = 2
W860 = 3
W124 = 4
W164 = 5
W213 = 6
NWL = 7
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
        self.toas = np.zeros((NWL,NSMALL,NAOT,NWS,NSZA,NVZA,NRAA))
        self.toab = np.zeros((NWL,NBIG,NAOT,NWS,NSZA,NVZA,NRAA))
        self.toar = np.zeros((NWL,NWS,NSZA,NVZA,NRAA))
        self.taus = np.zeros((NWL,NSMALL,NAOT))
        self.taub = np.zeros((NWL,NBIG,NAOT))
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
            self.tau_pts = (self.taus[W550,0,:])                                
            self.pts = (self.wnd_pts, self.sza_pts, self.vza_pts, self.raa_pts)
        except Exception as inst:
            print(type(inst)) 
            print(inst)  
            print ("Unable to read LUT file ... exiting")
            sys.exit()
        

    def interp_extrap( self, num, xin, x, y ):    
        if (xin <= x[0]):
            yout = y[0]+(xin-x[0])*(y[1]-y[0])/(x[1]-x[0])
        elif (xin >= x[num-1]):
            yout = y[num-2]+(xin-x[num-2])*(y[num-1]-y[num-2])/(x[num-1]-x[num-2])
        else:
            for i in range(num):
                if ((xin >= x[i]) and (xin <= x[i+1])):
                    yout = y[i]+(xin-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i])
                    break
        return yout;

    def proc_gran(self,rfl,sza,vza,raa,wnd):
        xit = np.stack((wnd, sza, vza, raa), axis=2)
        mrfls = np.zeros((NSMALL, NWL))
        mrflb = np.zeros((NBIG, NWL))
        mrflr = np.zeros((NWL))
        units = np.pi/np.cos(sza.astype(double)*D2R)
        
        max_rfl = np.max(self.toas)
        min_rfl = np.min(self.toas)
        shps = self.toas.shape
        shpb = self.toab.shape
        self.rfl_pts = np.linspace(0.0,0.5,shps[2])
        
        self.tslut = np.zeros_like(self.toas)
        self.tblut = np.zeros_like(self.toab)
        for iwl in range(shps[0]):
            for ism in range(shps[1]):
                for iwn in range(shps[3]):
                    for isz in range(shps[4]):
                        for ith in range(shps[5]):
                            for iph in range(shps[6]):
                                self.tslut[iwl,ism,:,iwn,isz,ith,iph] = \
                                np.interp(self.rfl_pts,self.toas[iwl,ism,:,iwn,isz,ith,iph],self.tau_pts)

        for iwl in range(shpb[0]):
            for ism in range(shpb[1]):
                for iwn in range(shpb[3]):
                    for isz in range(shpb[4]):
                        for ith in range(shpb[5]):
                            for iph in range(shpb[6]):
                                self.tblut[iwl,ism,:,iwn,isz,ith,iph] = \
                                np.interp(self.rfl_pts,self.toab[iwl,ism,:,iwn,isz,ith,iph],self.tau_pts)

        print("interpolate to reflectance for each wavelength and model")
        for iwl in range(NWL):
            for itau in range(NAOT):
                for ism in range(NSMALL):
                    mrfls[ism,iwl,itau,:] = units * \
                    interp.interpn(self.pts, self.toas[iwl,ism,itau,:,:,:,:], xit)
                for ibm in range(NBIG):
                    mrflb[ibm,iwl,itau,:] = units * \
                    interp.interpn(self.pts, toab[iwl,ibm,itau,:,:,:,:], xit)
            mrflr[iwl] = units * interp.interpn(self.r_pts, toab[iwl,ibm,:,:,:,:,:], xit)
        return 
    
    def process(self,rfl,sza,vza,raa,wnd):
        xit = np.stack((wnd, sza, vza, raa))
        units = np.pi/np.cos(sza*D2R)

        mrfls = interp.interpn(self.pts, self.toas, xit)[0]*units
        mrflb = interp.interpn(self.pts, self.toab, xit)[0]*units
        mrflr = interp.interpn(self.pts, self.toar, xit)[0]*units
        
        tauxs = np.zeros(NSMALL)
        tauxb = np.zeros(NBIG)
        for ism in range(NSMALL):
            mrfls[:,ism,0] = mrflr
            tauxs[ism] = np.interp(rfl[W860], mrfls[W860,ism], self.taus[W550,ism])
        for ibm in range(NBIG):
            mrflb[:,ibm,0] = mrflr
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
                denom = np.clip(rfl-mrflr,0,1.2) + 0.01
                rbdif = rfl - trflb[:,ibm]
                sbdif = trfls[:,ism] - trflb[:,ibm]
                sb2 = np.dot(sbdif[1:]/denom[1:],sbdif[1:]/denom[1:])
                rbsb = np.dot(rbdif[1:]/denom[1:],sbdif[1:]/denom[1:])
                xm = rbsb/sb2
                xm = np.max((np.min((xm,1.0)),0.0))
                mrfl = xm*trfls[:,ism] + (1.0-xm)*trflb[:,ibm]
#                err = np.divide((rfl - mrfl),denom)
                err = rfl - mrfl
                sse[ism,ibm] = np.dot(err[1:],err[1:])/(NWL-2)
                aot[ism,ibm] = xm*tauxs[ism] + (1.0-xm)*tauxb[ibm]
                fmf[ism,ibm] = xm
        im = divmod(sse.argmin(),NBIG)
        self.type = (NBIG-1)*im[1] + im[0]
        self.sse = sse[im[0],im[1]]
        self.aot = aot[im[0],im[1]]
        self.fmf = fmf[im[0],im[1]]
        self.rfl = rfl
        self.mrfl = self.fmf*trfls[:,im[0]] + (1.0-self.fmf)*trflb[:,im[1]]
        self.rsd = self.mrfl - self.rfl
        self.sse = np.dot(self.rsd,self.rsd)/(NWL-2)
        self.rayl = mrflr
        
        return self.fmf, self.aot, self.sse, self.type
    
    def plot(self, iy, ix):    
        plt.clf()
        plt.grid(True)
        plt.plot(self.wl_pts, self.rfl, marker='.', color='b', label='measured')
        plt.plot(self.wl_pts, self.mrfl, marker='.', color='g', label='modeled')
        plt.plot(self.wl_pts, self.rsd, marker='.', color='r', label='residual')
        plt.xlabel('wavelength (nm)')
        plt.ylabel('reflectance')
        tstr = "dtocean -- y={3:}, x={4:}  aot: {0:.3f}  fmf: {1:.3f}  sse: {2:.3}"
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
            self.rfl = np.append(xr.load_dataset(l1b_filepath,group='/reflectance')['toa_reflectance'][2:7,:,:].values, 
                            xr.load_dataset(l1b_filepath,group='/reflectance')['toa_reflectance'][8:,:,:].values,axis=0)
            self.sza = xr.load_dataset(l1b_filepath,group='/geolocation')['solar_zenith'].values
            self.vza = xr.load_dataset(l1b_filepath,group='/geolocation')['sensor_zenith'].values
            self.raa = xr.load_dataset(l1b_filepath,group='/geolocation')['relative_azimuth'].values
            self.lat = xr.load_dataset(l1b_filepath,group='/navigation_data')['latitude'].values
            self.lon = xr.load_dataset(l1b_filepath,group='/navigation_data')['longitude'].values
            self.cld = xr.load_dataset(l1b_filepath,group='/ancillary')['cloud_mask'].values
            self.wnd = xr.load_dataset(l1b_filepath,group='/ancillary')['wind_speed'].values
            self.rfl = np.transpose(self.rfl, (1,2,0))
        except Exception as inst:
            print(type(inst))   
            print(inst)          
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
            self.type = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.ds = xr.Dataset({'rfl': self.rfl, 'lat': self.lat, 'lon': self.lon, 
                                  'sza': self.sza, 'vza': self.vza, 'raa': self.raa, 
                                  'wnd': self.wnd, 'aot': self.aot, 'fmf': self.fmf, 
                                  'sse': self.sse, 'type': self.type })
        except Exception as inst:
            print(type(inst))   
            print(inst)          
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
    dto = dtocean(args.lut.name, args.mode)
   
    print ("Processing mode {0}".format(args.mode))
    for iy in range(args.y, args.y+args.z):
        for ix in range(args.x, args.x+args.z):
            
            if(vin.cld[iy,ix]):
                fmf,aot,sse,type = -999.9, -999.9, -999.9, -999
            else:
                fmf,aot,sse,type = dto.process(vin.rfl[iy,ix,:],vin.sza[iy,ix], 
                                          vin.vza[iy,ix],vin.raa[iy,ix], 
                                          vin.wnd[iy,ix])    
                if(args.plot): 
                    dto.plot(iy,ix)

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
            vout.ds.type[iy-args.y,ix-args.x] = type
            
        print('.', end = '')
    print('  Done!')
    vout.write()

if __name__ == "__main__":

    main()