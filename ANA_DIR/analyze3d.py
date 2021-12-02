import os, sys
import h5py
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import units as un
from skimage.measure import block_reduce
import plotly.graph_objects as go

import yt
from yt.funcs import mylog
mylog.setLevel(40)


class data:
    """
    Class to load, analyse and plot HDF5 data
    """
    def __init__(self, prim_files, units):
        self.prim_files = prim_files
        self.u = units
        
        self.rho = []
        self.press = []
        self.v = []
        self.v1 = []
        self.v2 = []
        self.v3 = []
        self.temp = []
        self.x1s = []
        self.x2s = []
        self.x3s = []
        self.x1p = []
        self.x2p = []
        self.x3p = []
        self.load()
        
        self.v = np.sqrt(np.array(self.v1)**2+np.array(self.v2)**2+np.array(self.v3)**2)
        
        print(np.shape(self.rho), np.shape(self.v))
        
        self.tonp()
        self.temp = self.press/self.u.presref / (self.rho/self.u.densref) * self.u.tempref
        
        self.phi = np.zeros(np.shape(self.rho[0,:,:]))
        self.vr = np.array([])
        self.vphi = np.array([])
        #self.polar_vel()
        
    
    def load(self):
        unit_base= {"length_unit":(self.u.abin,"AU"), "time_unit":(self.u.tbin,"s"), "mass_unit":(self.u.mbin,"Msun")}
        for fprim in self.prim_files:
            print(fprim)
            ds_prim = yt.load(fprim, units_override = unit_base)
            all_data_level_1_prim = ds_prim.covering_grid(level=1, \
                                            left_edge=ds_prim.domain_left_edge, dims=ds_prim.domain_dimensions * 2)
            rho = all_data_level_1_prim["rho"].in_units("g/cm**3").to_ndarray() # rho is a NumPy array
            press = all_data_level_1_prim["press"].in_units("g/(cm*s**2)").to_ndarray()
            v1 = all_data_level_1_prim["vel1"].in_units("km/s").to_ndarray()
            v2 = all_data_level_1_prim["vel2"].in_units("km/s").to_ndarray()
            v3 = all_data_level_1_prim["vel3"].in_units("km/s").to_ndarray()
            print(ds_prim.field_list)
            self.rho.append(rho)
            self.press.append(press)
            self.v1.append(v1)
            self.v2.append(v2)
            self.v3.append(v3)        
        
    def tonp(self):
        self.rho = np.array([r for r in self.rho])
        self.v = np.array([v for v in self.v])
        self.v1 = np.array([v for v in self.v1])
        self.v2 = np.array([v for v in self.v2])
        self.v3 = np.array([v for v in self.v3])
            
    def plot(self, simdata, fname, lab="", low=None, hig=None):
        if not low:
            low = np.min(simdata)
        if not hig:
            hig = np.max(simdata)
        for filename, data in zip(self.prim_files, simdata):
            #mean = np.mean(data)
            path = os.path.join(os.getcwd(), 'plots/'+ filename.replace(".athdf", fname))
            
            reduced = block_reduce(data, block_size=(8, 8, 8), func=np.mean)
            
            X, Y, Z = np.mgrid[-10:10:64j, -10:10:64j, -5:5:32j]

            fig = go.Figure(data=go.Volume(
                x=X.flatten(),
                y=Y.flatten(),
                z=Z.flatten(),
                value=reduced.flatten(),
                isomin=low,
                #isomax=1e-11,
                opacity=0.3, # needs to be small to see through all surfaces
                surface_count=15, # needs to be a large number for good volume rendering
                ))
            #fig.update_layout(scene_camera = dict(
            #    up=dict(x=1, y=0, z=0),
            #    center=dict(x=0, y=0, z=0),
            #    eye=dict(x=1.5, y=1.5, z=1.5)
            #))

            fig.write_image(path)
            print(filename)
    
    def polar_vel(self):
        #decompose velocity in radial and azimuthal parts
        
        #compute polar angle on mesh
        l = np.shape(self.phi)[0]
        X = np.array([[i-l/2+1./2 for i in range(l)] for j in range(l)])
        Y = X.T
        
        self.phi = np.arctan2(Y,X)
        
        #compute polar velocity
        self.vr   = self.v1*np.cos(self.phi) + self.v2*np.sin(self.phi)
        self.vphi = -self.v1*np.sin(self.phi) + self.v2*np.cos(self.phi)
    
    def column_dens(self, scale, low=None, hig=None):
        frame = self.u.abin * 10
        lim = self.u.abin * 10.*scale
        for filename, data in zip(self.prim_files, self.rho):
            column = block_reduce(data, block_size=(1, 1, 128*2), func=np.mean)[:,:,0]
            print(np.shape(column))
            if not low:
                lo = np.min(column)
            else:
                lo = low
            if not hig:
                hi = np.max(column)
            else:
                hi = hig
            
            path = os.path.join(os.getcwd(), 'plots/column/'+ filename.replace(".athdf", "_rho.png"))
            fig = plt.figure(figsize=(10,10))
            cmap="plasma"
            ax = fig.add_subplot(111)
            pos = ax.imshow(column, cmap=cmap, vmin=lo, vmax=hi, origin='lower', \
                            extent=[-frame,frame,-frame,frame])
            fig.colorbar(pos, ax=ax, label="$\Sigma ~[g/cm^3]$")

            ax.set_aspect('equal', 'box')
            ax.set_xlabel("$x_0~[AU]$")
            ax.set_ylabel("$x_1~[AU]$")
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)

            plt.savefig(path)
            plt.close()
            print(filename)
            
    def column_temp(self, scale, low=None, hig=None):
        frame = self.u.abin * 10
        lim = self.u.abin * 10.*scale
        for filename, data in zip(self.prim_files, self.temp):
            column = block_reduce(data, block_size=(1, 1, 128*2), func=np.max)[:,:,0]
            print(np.shape(column))
            if not low:
                lo = np.min(column)
            else:
                lo = low
            if not hig:
                hi = np.max(column)
            else:
                hi = hig
            
            path = os.path.join(os.getcwd(), 'plots/column/'+ filename.replace(".athdf", "_temp.png"))
            fig = plt.figure(figsize=(10,10))
            cmap="plasma"
            ax = fig.add_subplot(111)
            pos = ax.imshow(column, cmap=cmap, vmin=lo, vmax=hi, origin='lower', \
                            extent=[-frame,frame,-frame,frame])
            fig.colorbar(pos, ax=ax, label="$T ~[K]$")

            ax.set_aspect('equal', 'box')
            ax.set_xlabel("$x_0~[AU]$")
            ax.set_ylabel("$x_1~[AU]$")
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)

            plt.savefig(path)
            plt.close()
            print(filename)