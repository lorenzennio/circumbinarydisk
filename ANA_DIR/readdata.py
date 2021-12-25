import os, sys
import h5py
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import units as un
from astropy import constants as const
from skimage.measure import block_reduce
import plotly.graph_objects as go

import yt
from yt.funcs import mylog
mylog.setLevel(40)

class load:
    """
    define variables, units and coordinates
    """
    
    def __init__(self, prim_files, units, dims, level):
        self.prim_files = prim_files
        self.u = units

        #definitions
        self.rho = []
        self.press = []
        self.v1 = []
        self.v2 = []
        self.v3 = []

        #load data
        self.load(level)

        #convert to numpy array
        self.tonp()

        #remove levels form 2d
        if dims==2:
            self.rmlevels()
        elif dims ==3:
            pass
        else:
            raise ValueError('Can only process 2d or 3d.')

        #compute temperature
        self.temp = self.press/self.u.presref / (self.rho/self.u.densref) * self.u.tempref



    def load(self, level):
        unit_base= {"length_unit":(self.u.abin,"AU"), "time_unit":(self.u.tbin,"s"), "mass_unit":(self.u.mbin,"Msun")}
        for fprim in self.prim_files:
            print(fprim)
            ds_prim = yt.load(fprim, units_override = unit_base)
            # level = number of refinements
            # dims = total dimensions 2**(refinement level)
            all_data_level_1_prim = ds_prim.covering_grid(level=level, \
                                            left_edge=ds_prim.domain_left_edge, dims=ds_prim.domain_dimensions*2**level)
            rho = all_data_level_1_prim["rho"].in_units("g/cm**3").to_ndarray()
            press = all_data_level_1_prim["press"].in_units("g/(cm*s**2)").to_ndarray()
            v1 = all_data_level_1_prim["vel1"].in_units("km/s").to_ndarray()
            v2 = all_data_level_1_prim["vel2"].in_units("km/s").to_ndarray()
            v3 = all_data_level_1_prim["vel3"].in_units("km/s").to_ndarray()
            self.rho.append(rho)
            self.press.append(press)
            self.v1.append(v1)
            self.v2.append(v2)
            self.v3.append(v3) 

    def tonp(self):
        #self.rho = np.array([r for r in self.rho])
        self.rho = np.array([self.rho])
        self.press = np.array([self.press])
        self.v1 = np.array([self.v1])
        self.v2 = np.array([self.v2])
        self.v3 = np.array([self.v3])

    def rmlevels(self):
        self.rho = self.rho[:,:,:,0]
        self.press = self.press[:,:,:,0]
        self.v1 = self.v1[:,:,:,0]
        self.v2 = self.v2[:,:,:,0]
        self.v3 = self.v3[:,:,:,0]
            
#----------------------------------------------------------------              
class polarcoords:
    def __init__(self, dims):
        pass
                
#----------------------------------------------------------------       
class plot:
    def plot(self, data, path, low, hig, frame, lab, lim):
        cmap="plasma"
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        pos = ax.imshow(data, cmap=cmap, vmin=low, vmax=hig, origin='lower', \
                        extent=[-frame,frame,-frame,frame])
        fig.colorbar(pos, ax=ax, label=lab)

        ax.set_aspect('equal', 'box')
        ax.set_xlabel("$x_0~[AU]$")
        ax.set_ylabel("$x_1~[AU]$")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)

        plt.savefig(path)
        plt.close()

#----------------------------------------------------------------  

class data2d(load, plot):
    
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
    
    def draw(self, simdata, fname, lab="", scale = 1., low=None, hig=None):
        frame = self.u.abin * 10.
        lim = frame*scale
        cmap="plasma"
        for filename, data in zip(self.prim_files, simdata):
            
            path = os.path.join(os.getcwd(), 'plots/'+ filename.replace(".athdf", fname))
            
            if not low:
                lo = np.min(column)
            else:
                lo = low
            if not hig:
                hi = np.max(column)
            else:
                hi = hig
            
            self.plot(self, data, path, lo, hi, frame, lab, lim)
            print(filename)
    
class data3d(load, plot):
    
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
        
    def plot3d(self, simdata, fname, lab="", low=None, hig=None):
        if not low:
            low = np.min(simdata)
        if not hig:
            hig = np.max(simdata)
            
        for filename, data in zip(self.prim_files, simdata):
            #mean = np.mean(data)
            path = os.path.join(os.getcwd(), 'plots/'+ filename.replace(".athdf", fname))
            
            blocksize = np.shape(data)[0] / 64
            reduced = block_reduce(data, block_size=(blocksize, blocksize, blocksize), func=np.mean)
            
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
            #camera tilt
            #fig.update_layout(scene_camera = dict(
            #    up=dict(x=1, y=0, z=0),
            #    center=dict(x=0, y=0, z=0),
            #    eye=dict(x=1.5, y=1.5, z=1.5)
            #))

            fig.write_image(path)
            print(filename)
        
    def column_dens(self, simdata, fname, scale, low=None, hig=None):
        frame = self.u.abin*10.
        lim = frame*scale
        
        for filename, data in zip(self.prim_files, self.rho):
            path = os.path.join(os.getcwd(), 'plots/column/'+ filename.replace(".athdf", fname))
            
            zcolumn = np.shape(data)[-1]
            column = block_reduce(data, block_size=(1, 1, zcolumn), func=np.mean)[:,:,0]
            print(np.shape(column))
            
            if not low:
                lo = np.min(column)
            else:
                lo = low
            if not hig:
                hi = np.max(column)
            else:
                hi = hig
            
            self.plot(self, data, path, lo, hi, frame, lab, lim)
            
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
    

#----------------------------------------------------------------              