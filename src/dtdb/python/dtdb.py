#!/usr/bin/env python
# encoding: utf-8

import os 
import sys
import math
import time
import argparse
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from dtocean import dtocean
from dbocean import dbocean

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

class input(object):
   
    def __init__(self, l1b_filepath, st, ct):
        self.ifile = l1b_filepath
        print ("Reading sensor data: " + self.ifile)
        try:
            self.rfl = np.zeros((NWL,ct[0],ct[1]))
            self.rfl[W410] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_410'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W445] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_445'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W490] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_490'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W550] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_550'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W670] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_670'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W865] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_865'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W1240] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_1240'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W1610] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_1610'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.rfl[W2250] = xr.load_dataset(l1b_filepath,group='/observations')['rhot_2250'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.sza = xr.load_dataset(l1b_filepath,group='/geolocation')['solz'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.vza = xr.load_dataset(l1b_filepath,group='/geolocation')['senz'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.raa = xr.load_dataset(l1b_filepath,group='/geolocation')['relaz'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.lat = xr.load_dataset(l1b_filepath,group='/navigation_data')['latitude'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.lon = xr.load_dataset(l1b_filepath,group='/navigation_data')['longitude'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.cld = xr.load_dataset(l1b_filepath,group='/ancillary')['cloud_mask'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
            self.wnd = xr.load_dataset(l1b_filepath,group='/ancillary')['windspeed'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
#            self.chl = xr.load_dataset(l1b_filepath,group='/ancillary')['chlorophyll'].values[st[0]:st[0]+ct[0],st[1]:st[1]+ct[1]]
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
            typ = xr.DataArray(np.zeros((ydim,xdim)),dims=('y', 'x'))
            self.ds = xr.Dataset({'rfl': rfl, 'lat': lat, 'lon': lon, 
                                  'sza': sza, 'vza': vza, 'raa': raa, 
                                  'wnd': wnd, 'chl': chl, 'aot': aot, 
                                  'fmf': fmf, 'sse': sse, 'typ': typ})
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
    parser.add_argument("-um", "--unmask", type=bool, default=False,
                        help="ignore masks", required=False)
    parser.add_argument("y", type=int, help="start line")
    parser.add_argument("x", type=int, help="start pixel")
    parser.add_argument("z", type=int, help="square side ")
    
    args = parser.parse_args()
    dimx = dimy = args.z
    if (args.z <= 0):    
        args.y = args.x = 0
        vin = input(args.ifile.name)
        dimy = vin.shape[0]
        dimx = vin.shape[1]
    else:
        vin = input(args.ifile.name, [args.y,args.x], [dimy,dimx])
    
    vout = output(args.ofile.name, dimy, dimx)
 
    if (args.alg == "deepblue"):
        dtdb = dbocean(args.lut.name, args.mode)
    elif (args.alg == "darktarget"):
        dtdb = dtocean(args.lut.name, args.mode)
    else:
        print ("invalid algorithm")
        sys.exit()   
    
    if (args.mode == 2):
        print ("\nProcessing statistics mode {0}".format(args.mode))

        refwl = W865
        nbins = 100
        pca10 = np.clip(vin.rfl[:,:,refwl].flatten(), a_min=0.0001, a_max=1.0)
        pca20 = np.log10(pca10)
        hist = np.zeros((NWL, nbins, nbins))
        cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
        fig = plt.figure()
        bticks = True
        for wl in range(NWL):
            
            pca1 = np.clip(vin.rfl[:,:,wl].flatten(), a_min=0.0001, a_max=1.0)
            pca2 = np.log10(pca1)

            try:
                hist[wl], edges = np.histogramdd((pca20,pca2-pca20), range=[[-2.0,0.0],[-1.0,1.0]],bins=nbins, density=True) 
            except Exception as inst:
                print(type(inst))
                print(inst) 
                print ("Numpy.histogramdd failure ... exiting")
                sys.exit()
            
#            ax2 = fig.add_subplot(121, projection='3d')           
            ax2 = fig.add_subplot(1,9,wl+1)           
            xpos0, ypos0 = np.meshgrid(edges[0][:-1]+edges[0][1:], edges[1][:-1]+edges[1][1:])           
            cp2 = ax2.contourf(xpos0/2, ypos0/2, hist[wl].T, levels=100, cmap=cmap)
            plt.tick_params(labelleft=bticks, left=bticks)
            plt.grid(b=True, which='major', color='#666666', linestyle='-')
            plt.minorticks_on()
            plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            ax2.set_title("RHOT_"+wlstr[wl])
            ax2.set_xlabel("log10 (RHOT_"+wlstr[refwl]+")")
            if bticks:
                ax2.set_ylabel("log10(RHOT_X) - log10(RHOT_"+wlstr[refwl]+")")
            bticks = False
        fig.colorbar(cp2)
        path, fname = os.path.split(args.ifile.name)
        plt.suptitle(fname)
        plt.show()

    elif (args.mode == 1):
        print ("\nProcessing tensor mode {0}".format(args.mode))
        
        sy = [args.y,args.y+dimy]
        sx = [args.x,args.x+dimx]
        
        remaining = dimy
        count = 0
        CHUNK = 10
        cy = [args.y,args.y+dimy]
        while (remaining > CHUNK):
            
            chnk = min(CHUNK,remaining)
            cy[0] += CHUNK
            cy[1] = cy[0] + chnk
            rfl = vin.rfl[cy[0]:cy[1],sx[0]:sx[1]]
            lat = vin.lat[cy[0]:cy[1],sx[0]:sx[1]]
            lon = vin.lon[cy[0]:cy[1],sx[0]:sx[1]]
            sza = vin.sza[cy[0]:cy[1],sx[0]:sx[1]]
            vza = vin.vza[cy[0]:cy[1],sx[0]:sx[1]]
            raa = vin.raa[cy[0]:cy[1],sx[0]:sx[1]]
            wnd = vin.wnd[cy[0]:cy[1],sx[0]:sx[1]]
#            chl = vin.chl[cy[0]:cy[1],sx[0]:sx[1]]
            
            fmf,aot,sse,typ = dtdb.proc_gran(rfl,sza,vza,raa,wnd)
            
            oy = [cy[0] - args.y, cy[1] - args.y]
            vout.ds.rfl[:,oy[0]:oy[1],sx[0]:sx[1]] = rfl.transpose((2,0,1))
            vout.ds.lat[oy[0]:oy[1],sx[0]:sx[1]] = lat
            vout.ds.lon[oy[0]:oy[1],sx[0]:sx[1]] = lon
            vout.ds.sza[oy[0]:oy[1],sx[0]:sx[1]] = sza
            vout.ds.vza[oy[0]:oy[1],sx[0]:sx[1]] = vza
            vout.ds.raa[oy[0]:oy[1],sx[0]:sx[1]] = raa
            vout.ds.wnd[oy[0]:oy[1],sx[0]:sx[1]] = wnd
#            vout.ds.chl[oy[0]:oy[1],sx[0]:sx[1]] = chl
            vout.ds.aot[oy[0]:oy[1],sx[0]:sx[1]] = aot
            vout.ds.fmf[oy[0]:oy[1],sx[0]:sx[1]] = fmf
            vout.ds.sse[oy[0]:oy[1],sx[0]:sx[1]] = sse
            vout.ds.typ[oy[0]:oy[1],sx[0]:sx[1]] = typ
            
            remaining -= CHUNK
            count += 1
    else:
        print ("\nProcessing pixel mode {0}".format(args.mode))
        for iy in range(dimy):
            tic = time.perf_counter()
            tip =0
            for ix in range(dimx):
                
                if(vin.cld[iy,ix] and not args.unmask):
                    fmf,aot,chl,wnd,sse = -999.9, -999.9, -999.9, -999.9, -999.9
                else:
                    try:
                        fmf,aot,wnd,sse,aero = dtdb.process(vin.rfl[iy,ix],vin.sza[iy,ix], 
                                                    vin.vza[iy,ix],vin.raa[iy,ix],vin.wnd[iy,ix])  
                    except Exception as inst:
                        print(type(inst))
                        print(inst) 
                        print ("processing error at pixel", iy, ix)
                    
                    if(args.plot): 
                        dtdb.plot(iy,ix)
#                        dtdb.plot_points()
    
                vout.ds.rfl[:,iy,ix] = vin.rfl[iy,ix]
                vout.ds.lat[iy,ix] = vin.lat[iy,ix]
                vout.ds.lon[iy,ix] = vin.lon[iy,ix]
                vout.ds.sza[iy,ix] = vin.sza[iy,ix]
                vout.ds.vza[iy,ix] = vin.vza[iy,ix]
                vout.ds.raa[iy,ix] = vin.raa[iy,ix]
                vout.ds.wnd[iy,ix] = wnd
#                vout.ds.chl[iy,ix] = chl
                vout.ds.aot[iy,ix] = aot
                vout.ds.fmf[iy,ix] = fmf
                vout.ds.sse[iy,ix] = sse
            
        toc = time.perf_counter()
        tpp = (toc-tic)/dimx
#        print('.', end = '', flush=True)
        print(iy, tpp, " sec per pixel", flush=True)
    print('  Done!')
    vout.write()

if __name__ == "__main__":

    main()