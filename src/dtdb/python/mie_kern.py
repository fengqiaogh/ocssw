#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import json
import numpy as np
from numpy import linalg as la
from scipy.stats import lognorm
from scipy.integrate import quad
from scipy.spatial import distance
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmath 
import bhmie as bh 

oci_wl = [310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385,
    390, 395, 400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465,
    470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545,
    550, 555, 560, 565, 570, 575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625, 
    630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 700, 705, 
    710, 715, 720, 725, 730, 735, 740, 745, 750, 755, 760, 765, 770, 775, 780, 785, 
    790, 795, 800, 805, 810, 815, 820, 825, 830, 835, 840, 845, 850, 855, 860, 865, 
    870, 875, 880, 885, 890, 895, 940, 1040, 1250, 1378, 1615, 2130, 2260 ]

viirs_wl = [ 412, 488, 550, 670, 865, 1240, 1610, 2250 ]

def normalize(a, order=2, axis=-1 ):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

def mie(sz,wl,m):
    s1,s2,qx,qs,qb,gs = bh.bhmie(2*np.pi*sz/wl,m,1)    
    return qs


class mie_kern(object):
   
    def __init__(self, name="", platform='oci', m=1.4+0j, start_sz=200, stop_sz=20000, points=50):
        self.name  = name
        if (platform == 'oci'):
            self.wl = np.copy(oci_wl[8:-7])
        elif (platform == 'viirs'):
            self.wl = np.copy(viirs_wl)
        self.nwl = len(self.wl)
        self.start = start_sz
        self.stop = stop_sz
        self.points = points
        self.lsz = np.linspace(np.log(self.start),np.log(self.stop),self.nwl)
        self.sz = np.exp(self.lsz)
        self.nsz = len(self.sz)
        self.m_real = np.zeros((self.nwl), dtype=np.float64) 
        self.m_imag = np.zeros((self.nwl), dtype=np.float64) 
        self.m_real[:] = np.real(m) 
        self.m_imag[:] = np.imag(m) 
        self.kernel = np.zeros((self.nwl, self.nsz), dtype=np.float64)
        # log-normal distribution standard deviation 
        self.lsig = np.zeros(self.nsz)
        self.lsig[0] = np.log((self.sz[1]-self.sz[0])/2)
        self.sig = np.zeros(self.nsz)
        self.sig[0] = (self.sz[1]-self.sz[0])/2
        for i in range(1,self.nsz-1):
            self.lsig[i] = (self.lsz[i+1]-self.lsz[i-1])/5.0
            self.sig[i] = (self.sz[i+1]-self.sz[i-1])/4.5
        self.lsig[0] = self.lsig[1]
        self.lsig[self.nsz-1] = self.lsig[self.nsz-2]
        self.sig[0] = self.sig[1]
        self.sig[self.nsz-1] = self.sig[self.nsz-2]
        
    def generate(self):
        print('Generating kernel: '+ self.name, end = ' ')
        # compute kernel based on Mie scattering
        m = self.m_real + self.m_imag*1j
        lpdf = lambda y, s, scale: lognorm.pdf(y, s=s, loc=0, scale=scale)
        for iwl in range(self.nwl):
            for isz in range(1,self.nsz):
        #        lmie = lambda y: lpdf(y,lsig[isz],sz[isz])*mie(y,wl[iwl],m)
                xs = np.linspace(self.sz[isz]-5*self.sig[isz], self.sz[isz]+ 5*self.sig[isz], self.points )
                ys = np.zeros(len(xs))
                for k in range(len(xs)):
                    ys[k] = lpdf(xs[k],self.lsig[isz],self.sz[isz])*mie(xs[k],self.wl[iwl],m[iwl])
                self.kernel[iwl,isz] = np.trapz(ys,xs)
        #        self.kernel[iwl,isz],errs[iwl,isz] = quad(lmie, self.sz[isz]-3*self.sig[isz], self.sz[isz] \
        #                                     + 3*self.sig[isz], epsabs=0.0001,limit=self.points)
        #        self.kernel[iwl,isz] = mie(self.sz[isz],self.wl[iwl],m[iwl])
            print('.', end = '')
        print('  Done!')
        self.kernel[:,0] = 4*(self.wl[0]/self.wl[:])**4 
#        self.kernel = normalize(self.kernel, axis=0)
        return self.kernel

    def inverse(self, rp0, rp1, rp2, rp3):
        # regularization matrices
        I = np.identity(self.nsz, dtype=np.float64)
        # first difference matrix
        H0 = -np.identity(self.nsz, dtype=np.float64)
        H0[0,0] = 0
        for i in range(1,self.nsz):
            for j in range(self.nsz):
                H0[i,j] = 1 if (j==i-1) else H0[i,j]
        # sum of first differences squared
        H1 = np.dot(H0.T,H0)
        # sum of second differences squared
        H2 = H1.copy()
        H2[0,:] = 0
        H2 = np.dot(H2.T,H2) 
        # sum of third differences squared
        H3 = np.zeros([self.nsz, self.nwl])  
        for i in range(3,self.nsz):
            for j in range(self.nwl):
                if (j==i):
                     H3[i,j] = -1  
                elif (j==i-1):
                    H3[i,j] = 3
                elif (j==i-2):
                    H3[i,j] = -3
                elif (j==i-3):
                    H3[i,j] = 1
        H3 = np.dot(H3.T,H3)
        
        self.ikernel = la.inv(np.dot(self.kernel.T,self.kernel) + rp0*I + rp1*H1 + rp2*H2 + rp3*H3)
        return self.ikernel
    
    def plot(self, what):
        mpl.rcParams["font.size"] = 18    
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        wlg, szg = np.meshgrid(self.wl, self.sz)
        if what == "ikern":
            ax.plot_wireframe(wlg, szg, self.ikernel.T)
        else:
            ax.plot_wireframe(wlg, szg, self.kernel.T)
        plt.show()

    def save(self, odir):
        ofilepath = odir + '/KERN_' + self.name + '.txt'
        with open(ofilepath, 'w') as outfile:
            list = self.kernel.tolist()
            json.dump(list, outfile, separators=(',', ':'), indent=4)
        return ofilepath

    def to_json(self):
        jd = json.dumps(self.__dict__, default=lambda o: o.tolist(), indent=4)
        return jd 
       
    @classmethod
    def from_json(cls, data):
        return cls(data)
        
    def save_all(self, odir):
        ofilepath = odir + '/KERN_' + self.name + '_all.txt'
        with open(ofilepath, 'w') as outfile:
            json_obj = self.to_json()
            outfile.write(json_obj)
        return ofilepath

    def load(self, ifilepath):
        with open(ifilepath, 'r', encoding='utf-8') as infile:
            obj_text = infile.read()
            dict = json.loads(obj_text)
            self.name  = dict['name']
            self.wl = np.array(dict['wl'])
            self.nwl = dict['nwl']
            self.start = dict['start']
            self.stop = dict['stop']
            self.points = dict['points']
            self.lsz = np.array(dict['lsz'])
            self.sz = np.array(dict['sz'])
            self.nsz = dict['nsz']
            self.m_real = np.array(dict['m_real']) 
            self.m_imag = np.array(dict['m_imag']) 
            self.kernel = np.array(dict['kernel'])
            self.lsig = np.array(dict['lsig'])
            self.sig = np.array(dict['sig'])
        return self.kernel

