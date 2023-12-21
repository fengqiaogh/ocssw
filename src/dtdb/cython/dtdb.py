#!/usr/bin/env python
# encoding: utf-8

import os 
import sys
import math
import time
import argparse
import xarray as xr
import numpy as np
from dtocean import dtocean
from dbocean import dbocean

NWL = 7

class input(object):
   
    def __init__(self, l1b_filepath):
        self.ifile = l1b_filepath
        print ("Reading sensor data: " + self.ifile)
        try:
            self.rfl = np.append(xr.load_dataset(l1b_filepath,group='/reflectance')['toa_reflectance'][2:7,:,:], 
                            xr.load_dataset(l1b_filepath,group='/reflectance')['toa_reflectance'][8:,:,:],axis=0)
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
    parser.add_argument("-af", "--alg", 
                        help="algorithm name", required=True)
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
 
    if (args.alg == "deepblue"):
        dtdb = dbocean(args.lut.name, args.mode)
    elif (args.alg == "darktarget"):
        dtdb = dtocean(args.lut.name, args.mode)
    else:
        print ("invalid algorithm")
        sys.exit()   
   
    print ("Processing mode {0}".format(args.mode))
    for iy in range(args.y, args.y+dimy):
        tic = time.perf_counter()
        tip =0
        for ix in range(args.x, args.x+dimx):
            
            if(vin.cld[iy,ix] or vin.chl[iy,ix]<-999):
                fmf,aot,chl,wnd,sse = -999.9, -999.9, -999.9, -999.9, -999.9
            else:
                try:
                    fmf,aot,chl,wnd,sse = dtdb.process(vin.rfl[iy,ix],vin.sza[iy,ix], 
                                                vin.vza[iy,ix],vin.raa[iy,ix], 
                                                vin.wnd[iy,ix],vin.chl[iy,ix])  
                except Exception as inst:
                    print(type(inst))
                    print(inst) 
                    print ("processing error at pixel", iy, ix)
                
            if(args.plot): 
                dbdt.plot(iy,ix)

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