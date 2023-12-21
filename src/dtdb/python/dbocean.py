#!/usr/bin/env python
# encoding: utf-8

import os 
import sys
import math
import time
import argparse
import numpy as np
import lmfit as lm
import xarray as xr
import matplotlib.pyplot as plt
import scipy.interpolate as trp

W470 = 0
W550 = 1
W659 = 2
W860 = 3
W124 = 4
W164 = 5
W213 = 6
NWL = 7
NSZA = 22
NVZA = 20
NRAA = 21
NAOT = 14
NAOTM = 7
NFMF = 14
NFMFM = 8
NF1 = 5
NF2 = 9
NWS = 6
NCHL = 4
D2R = np.pi/180.0

class dbocean(object):
   
    def __init__(self, lut_filepath, mode):
        print ("Reading Deepblue LUT: " + lut_filepath)
        self.mode = mode
        self.chlc = -2.0
        self.wndc = 1.0
        self.type = 0
        self.lut = np.zeros((NWL,NCHL,NWS,NFMF,NAOT,NRAA,NVZA,NSZA))
        self.lutm = np.zeros((NWL,NCHL,NWS,NFMFM,NAOTM,NRAA,NVZA,NSZA))
        self.fmf_pts = np.zeros((NFMF))
        try:
            self.lutm[W470][:,:,:,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MARITIME')['IoverF_m03'].values
            self.lutm[W550][:,:,:,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MARITIME')['IoverF_m04'].values
            self.lutm[W659][:,:,:,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MARITIME')['IoverF_m05'].values
            self.lutm[W860][:,:,:,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MARITIME')['IoverF_m07'].values
            self.lutm[W124][:,:,:,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MARITIME')['IoverF_m08'].values
            self.lutm[W164][:,:,:,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MARITIME')['IoverF_m10'].values
            self.lutm[W213][:,:,:,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MARITIME')['IoverF_m11'].values
            self.lutm = np.transpose(self.lutm, (1,2,5,6,7,3,4,0))
            self.lut[W470][:,:,0:NF1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['IoverF_m03'].values
            self.lut[W470][:,:,NF1-1:NF2-1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['IoverF_m03'].values
            self.lut[W470][:,:,NF2-2:NFMF,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['IoverF_m03'].values
            self.lut[W550][:,:,0:NF1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['IoverF_m04'].values
            self.lut[W550][:,:,NF1-1:NF2-1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['IoverF_m04'].values
            self.lut[W550][:,:,NF2-2:NFMF,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['IoverF_m04'].values
            self.lut[W659][:,:,0:NF1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['IoverF_m05'].values
            self.lut[W659][:,:,NF1-1:NF2-1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['IoverF_m05'].values
            self.lut[W659][:,:,NF2-2:NFMF,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['IoverF_m05'].values
            self.lut[W860][:,:,0:NF1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['IoverF_m07'].values
            self.lut[W860][:,:,NF1-1:NF2-1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['IoverF_m07'].values
            self.lut[W860][:,:,NF2-2:NFMF,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['IoverF_m07'].values
            self.lut[W124][:,:,0:NF1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['IoverF_m08'].values
            self.lut[W124][:,:,NF1-1:NF2-1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['IoverF_m08'].values
            self.lut[W124][:,:,NF2-2:NFMF,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['IoverF_m08'].values
            self.lut[W164][:,:,0:NF1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['IoverF_m10'].values
            self.lut[W164][:,:,NF1-1:NF2-1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['IoverF_m10'].values
            self.lut[W164][:,:,NF2-2:NFMF,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['IoverF_m10'].values
            self.lut[W213][:,:,0:NF1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['IoverF_m11'].values
            self.lut[W213][:,:,NF1-1:NF2-1,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['IoverF_m11'].values
            self.lut[W213][:,:,NF2-2:NFMF,:,:,:,:]=xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['IoverF_m11'].values
            self.chl_lut = xr.load_dataset(lut_filepath,group='/LOG_CHL')['LOG_CHL'].values
            self.raa_pts = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['Relative_Azimuth_Angle'].values
            self.vza_pts = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['View_Zenith_Angle'].values
            self.sza_pts = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['Solar_Zenith_Angle'].values
            self.wnd_pts = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['Wind_Speed'].values
            self.chl_pts = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['Chl_Conc'].values
            self.aot_pts = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['Aerosol_Optical_Depth_550'].values
            self.fmf_pts[0:NF1] = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_DUST')['Fine_Mode_Fraction_550'].values
            self.fmf_pts[NF1-1:NF2-1] = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['Fine_Mode_Fraction_550'].values
            self.fmf_pts[NF2-2:NFMF] = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_FINE')['Fine_Mode_Fraction_550'].values
            self.wl_pts = xr.load_dataset(lut_filepath,group='/OCEAN_AEROSOL_MIXED')['Band_Central_Wavelength'].values
            self.coef = np.array([0.06,0.06,0.04,0.04,0.07,0.06,0.1])
            self.pars = lm.Parameters()
            self.pars.add(name='fmf', value=0.5, min=0, max=1)
            self.pars.add(name='aot', value=1.0, min=0, max=100)
            if(self.mode==0 or self.mode==1):
                self.rpts = (self.chl_pts, self.wnd_pts, self.raa_pts, self.vza_pts, self.sza_pts)
                self.mpts = (self.fmf_pts, self.aot_pts)
                self.lut = np.transpose(self.lut, (1,2,5,6,7,3,4,0))
            elif(self.mode==2):
                self.rpts = (self.wnd_pts, self.raa_pts, self.vza_pts, self.sza_pts)
                self.mpts = (self.fmf_pts, self.aot_pts, self.chl_pts)
                self.lut = np.transpose(self.lut, (2,5,6,7,3,4,1,0))
                self.pars.add(name='chl', value=-2.0, min=-10.0, max=2.0)
            elif(self.mode==3):
                self.rpts = (self.chl_pts, self.raa_pts, self.vza_pts, self.sza_pts)
                self.mpts = (self.fmf_pts, self.aot_pts, self.wnd_pts)
                self.lut = np.transpose(self.lut, (1,5,6,7,3,4,2,0))
                self.pars.add(name='wnd', value=1.0, min=0, max=100)
            elif(self.mode==4):
                self.rpts = (self.raa_pts, self.vza_pts, self.sza_pts)
                self.mpts = (self.fmf_pts, self.aot_pts, self.chl_pts, self.wnd_pts)
                self.lut = np.transpose(self.lut, (5,6,7,3,4,1,2,0))
                self.pars.add(name='chl', value=-2.0, min=-10.0, max=2.0)
                self.pars.add(name='wnd', value=1.0, min=0, max=100)
        except Exception as inst:
            print(type(inst))
            print(inst) 
            print ("Unable to read LUT file ... exiting")
            sys.exit()

    def minfun(self, pars, data, scale, rlut):        
        if(self.mode==0 or self.mode==1):
            rxi = np.stack((pars['fmf'], pars['aot']))
        elif(self.mode==2):
            rxi = np.stack((pars['fmf'], pars['aot'], pars['chl']))
        elif(self.mode==3):
            rxi = np.stack((pars['fmf'], pars['aot'], pars['wnd']))
        elif(self.mode==4):
            rxi = np.stack((pars['fmf'], pars['aot'], pars['chl'], pars['wnd']))
        model = trp.interpn(self.mpts,rlut,rxi, 
                    bounds_error=False,fill_value=None )[0]
        return (model - data)*scale

    def process(self,rfl,sza,vza,raa,wnd,chl):
        cosz = np.cos(sza*D2R)
        rfl *= cosz/np.pi
        if(self.mode==0 or self.mode==1):
            if(self.mode==1):
                chl, wnd = self.chlc, self.wndc
            xi = np.stack((chl, wnd, raa, vza, sza))
        elif(self.mode==2):
            wnd = self.wndc
            xi = np.stack((wnd, raa, vza, sza))
        elif(self.mode==3):
            chl = self.chlc
            xi = np.stack((chl, raa, vza, sza))
        elif(self.mode==4):
            xi = np.stack((raa, vza, sza))
        tlut = trp.interpn(self.rpts, self.lut, xi)[0]
        scale = 1.0/(self.coef*(rfl + 0.00001))
        mfit = lm.minimize(self.minfun, self.pars, args=(rfl, scale, tlut))
        self.fmf = mfit.params['fmf'].value
        self.aot = mfit.params['aot'].value
        self.sse = mfit.redchi
        if(self.mode==0 or self.mode==1):
            self.chl = float(chl)
            self.wnd = float(wnd)
            mxi = np.stack((self.fmf,self.aot))
        elif(self.mode==2):
            self.wnd = float(wnd)
            self.chl = mfit.params['chl'].value
            mxi = np.stack((self.fmf,self.aot,self.chl))
        elif(self.mode==3):
            self.chl = float(chl)
            self.wnd = mfit.params['wnd'].value
            mxi = np.stack((self.fmf,self.aot,self.wnd))
        elif(self.mode==4):
            self.chl = mfit.params['chl'].value
            self.wnd = mfit.params['wnd'].value
            mxi = np.stack((self.fmf,self.aot,self.chl,self.wnd))
        mrfl = trp.interpn(self.mpts, tlut, mxi, 
                bounds_error=False, fill_value=None )[0]
# convert return values to unnormalized units
        self.rfl = rfl*np.pi/cosz
        self.mrfl = mrfl*np.pi/cosz
        self.rsd = self.mrfl - self.rfl
        self.sse = np.dot(self.rsd,self.rsd)/(NWL-2)
        return self.fmf, self.aot, self.chl, self.wnd, self.sse, self.type
    
    def plot(self, iy, ix):    
        plt.clf()
        plt.grid(True)
        plt.plot(self.wl_pts, self.rfl, marker='.', color='b', label='measured')
        plt.plot(self.wl_pts, self.mrfl, marker='.', color='g', label='modeled')
        plt.plot(self.wl_pts, self.rsd, marker='.', color='r', label='residual')
        plt.xlabel('wavelength (nm)')
        plt.ylabel('reflectance')
        tstr = "dbocean mode {3:d} -- y={4:}, x={5:}  aot: {0:.3f}  fmf: {1:.3f}  sse: {2:.3}"
        plt.title(tstr.format(self.aot, self.fmf, self.sse, self.mode, iy, ix))
        plt.legend(loc='upper right')
        plt.show()    


class input(object):
   
    def __init__(self, l1b_filepath):
        self.ifile = l1b_filepath
        print ("Reading sensor data: " + self.ifile)
        try:
            self.rfl[W470] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_490']
            self.rfl[W550] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_550']
            self.rfl[W659] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_670']
            self.rfl[W860] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_865']
            self.rfl[W124] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_1240']
            self.rfl[W164] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_1610']
            self.rfl[W213] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_2250']
            self.sza = xr.load_dataset(l1b_filepath,group='/geolocation')['solar_zenith'].values
            self.vza = xr.load_dataset(l1b_filepath,group='/geolocation')['sensor_zenith'].values
            self.raa = xr.load_dataset(l1b_filepath,group='/geolocation')['relative_azimuth'].values
            self.lat = xr.load_dataset(l1b_filepath,group='/navigation_data')['latitude'].values
            self.lon = xr.load_dataset(l1b_filepath,group='/navigation_data')['longitude'].values
            self.cld = xr.load_dataset(l1b_filepath,group='/ancillary')['cloud_mask'].values
            self.wnd = xr.load_dataset(l1b_filepath,group='/ancillary')['wind_speed'].values
            self.chl = xr.load_dataset(l1b_filepath,group='/ancillary')['chlorophyll'].values
            self.rfl = np.transpose(self.rfl, (1,2,0))
            self.shape = self.sza.shape
        except Exception as inst:
            print(type(inst))
            print(inst) 
            print ("Unable to read from file ... exiting")
            sys.exit()

class output(object):
   
    def __init__(self, out_filepath, ydim, xdim):
        self.ofile = out_filepath
        self.ydim = ydim
        self.xdim = xdim
        try:
            rfl = xr.DataArray(np.zeros((NWL,ydim,xdim)),dims=('wl','y', 'x'))
            lat = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            lon = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            sza = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            vza = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            raa = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            wnd = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            chl = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            aot = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            fmf = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            sse = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.ds = xr.Dataset({'rfl': rfl, 'lat': lat, 'lon': lon, 
                                  'sza': sza, 'vza': vza, 'raa': raa, 
                                  'wnd': wnd, 'chl': chl, 'aot': aot, 
                                  'fmf': fmf, 'sse': sse, })
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
    parser.add_argument("-if", "--ifile", type=argparse.FileType('r'), 
                        help="input file", required=True)
    parser.add_argument("-of", "--ofile", type=argparse.FileType('w'), 
                        help="output file", required=False)
    parser.add_argument("-lf", "--lut", type=argparse.FileType('r'), 
                        help="lookup table file", required=True)
    parser.add_argument("-pf", "--plot", type=bool, default=False,
                        help="plot pixel data", required=False)
    parser.add_argument("-mf", "--mode", type=int, default=0,
                        help="mode option", required=False)
    parser.add_argument("y", type=int, help="start line")
    parser.add_argument("x", type=int, help="start pixel")
    parser.add_argument("z", type=int, help="square side ")
    
    args = parser.parse_args()
    vin = input(args.ifile.name)
    dimx = dimy = args.z
    if (args.z <= 0):    
        args.y = args.x = 0
        dimy = vin.shape[0]
        dimx = vin.shape[1]
    
    vout = output(args.ofile.name, dimy, dimx)
 
    dbo = dbocean(args.lut.name, args.mode)
   
    print ("Processing mode {0}".format(args.mode))
    for iy in range(args.y, args.y+dimy):
        tic = time.perf_counter()
        for ix in range(args.x, args.x+dimx):
            
            if(vin.cld[iy,ix] or vin.chl[iy,ix]<-999):
                fmf,aot,chl,wnd,sse = -999.9, -999.9, -999.9, -999.9, -999.9
            else:
                try:
                    fmf,aot,chl,wnd,sse = dbo.process(vin.rfl[iy,ix],vin.sza[iy,ix], 
                                                vin.vza[iy,ix],vin.raa[iy,ix], 
                                                vin.wnd[iy,ix],vin.chl[iy,ix])  
                except Exception as inst:
                    print(type(inst))
                    print(inst) 
                    print ("processing error at pixel", iy, ix)
                
            if(args.plot): 
                dbo.plot(iy,ix)

            vout.ds.rfl[:,iy-args.y,ix-args.x] = vin.rfl[iy,ix]
            vout.ds.lat[iy-args.y,ix-args.x] = vin.lat[iy,ix]
            vout.ds.lon[iy-args.y,ix-args.x] = vin.lon[iy,ix]
            vout.ds.sza[iy-args.y,ix-args.x] = vin.sza[iy,ix]
            vout.ds.vza[iy-args.y,ix-args.x] = vin.vza[iy,ix]
            vout.ds.raa[iy-args.y,ix-args.x] = vin.raa[iy,ix]
            vout.ds.wnd[iy-args.y,ix-args.x] = wnd
            vout.ds.chl[iy-args.y,ix-args.x] = chl
            vout.ds.aot[iy-args.y,ix-args.x] = aot
            vout.ds.fmf[iy-args.y,ix-args.x] = fmf
            vout.ds.sse[iy-args.y,ix-args.x] = sse
            
        toc = time.perf_counter()
        tpp = (toc-tic)/dimx
#        print('.', end = '', flush=True)
        print(iy, tpp, flush=True)
    print('  Done!')
    vout.write()

if __name__ == "__main__":

    main()