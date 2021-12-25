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

class data:
    """ 
    Class to load, analyse and plot HDF5 data for 2D simultions
    """
    def __init__(self, units, prim_files, temp_files = []):
        self.prim_files = prim_files
        self.temp_files = temp_files
        
        self.u = units
        
        self.rho = []
        self.press = []
        self.v = []
        self.v1 = []
        self.v2 = []
        self.v3 = []
        #self.temp = []
        self.x1s = []
        self.x2s = []
        self.x3s = []
        self.x1p = []
        self.x2p = []
        self.x3p = []
        self.load()
        
        self.v = np.sqrt(np.array(self.v1)**2+np.array(self.v2)**2+np.array(self.v3)**2)
        
        self.kepler = []
        self.cent = np.array([self.x1p,self.x2p]).T
        
        #Centered Kepler, or shifted center
        self.kepler_vel()
        #self.kepler_vel(center=self.cent[85:86], size=0.6)
        
        self.tonp()
        self.temp = self.press/self.u.presref / (self.rho/self.u.densref) * self.u.tempref
        
        self.phi = np.zeros(np.shape(self.rho[0,:,:]))
        self.vr = np.array([])
        self.vphi = np.array([])
        self.polar_vel()
        
    
    def load(self):
        unit_base= {"length_unit":(self.u.abin,"AU"), "time_unit":(self.u.tbin,"s"), "mass_unit":(self.u.mbin,"Msun")}
        """
        for fprim, ftemp in zip(self.prim_files, self.temp_files):
            ds_prim = yt.load(fprim, units_override = unit_base)
            all_data_level_1_prim = ds_prim.covering_grid(level=1, \
                                            left_edge=ds_prim.domain_left_edge, dims=ds_prim.domain_dimensions * 2)
            ds_temp = yt.load(ftemp, units_override = unit_base)
            all_data_level_1_temp = ds_temp.covering_grid(level=1, \
                                            left_edge=ds_temp.domain_left_edge, dims=ds_temp.domain_dimensions * 2)
            rho = all_data_level_1_prim["rho"].in_units("g/cm**3").to_ndarray() # rho is a NumPy array
            press = all_data_level_1_prim["press"].in_units("g/(cm*s**2)").to_ndarray()
            v1 = all_data_level_1_prim["vel1"].in_units("km/s").to_ndarray()
            v2 = all_data_level_1_prim["vel2"].in_units("km/s").to_ndarray()
            v3 = all_data_level_1_prim["vel3"].in_units("km/s").to_ndarray()
            print(ds_prim.field_list)
            temp = all_data_level_1_temp["user_out_var0"].to_ndarray()
            self.rho.append(rho)
            self.press.append(press)
            self.v1.append(v1)
            self.v2.append(v2)
            self.v3.append(v3)
            self.temp.append(temp)
        """
        for fprim in self.prim_files:
            ds_prim = yt.load(fprim, units_override = unit_base)
            # level = number of refinements
            # dims = total dimensions 2**(refinement level)
            all_data_level_1_prim = ds_prim.covering_grid(level=2, \
                                            left_edge=ds_prim.domain_left_edge, dims=ds_prim.domain_dimensions * 2**2)
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
        #self.load_sink()
        
    def load_sink(self):
        #sink paticle positions
        binary_orbit = np.loadtxt('/u/lgaertner/RUN_DIR/binary_orbit.tab')
        self.x1s = binary_orbit[:,4]
        self.x2s = binary_orbit[:,5]
        self.x3s = binary_orbit[:,6]
        self.x1p = binary_orbit[:,7]
        self.x2p = binary_orbit[:,8]
        self.x3p = binary_orbit[:,9]
        
        self.x1s = np.array(self.x1s)
        self.x2s = np.array(self.x2s)
        self.x3s = np.array(self.x3s)
        self.x1p = np.array(self.x1p)
        self.x2p = np.array(self.x2p)
        self.x3p = np.array(self.x3p)
        
    def tonp(self):
        self.rho = np.array(self.rho)[:,:,:,0]
        self.rho = np.array([r.T for r in self.rho])
        self.press = np.array(self.press)[:,:,:,0]
        self.press = np.array([r.T for r in self.press])
        self.v = np.array(self.v)[:,:,:,0]
        self.v = np.array([v.T for v in self.v])
        self.v1 = np.array(self.v1)[:,:,:,0]
        self.v1 = np.array([v.T for v in self.v1])
        self.v2 = np.array(self.v2)[:,:,:,0]
        self.v2 = np.array([v.T for v in self.v2])
        self.v3 = np.array(self.v3)[:,:,:,0]
        self.v3 = np.array([v.T for v in self.v3])
        #self.temp = np.array(self.temp)[:,:,:,0]
        #self.temp = np.array([v.T for v in self.v3])
        
        
        self.kepler = np.array(self.kepler)
            
    def plot(self, simdata, fname, lab="", scale = 1., height = []):
        if not height:
            height = np.array([np.min(simdata), np.max(simdata)])
        frame = self.u.abin * 10
        lim = self.u.abin * 10.*scale
        for filename, data in zip(self.prim_files, simdata):
            
            path = os.path.join(os.getcwd(), 'plots/'+ filename.replace(".athdf", fname))
            fig = plt.figure(figsize=(10,10))
            cmap="plasma"
            ax = fig.add_subplot(111)
            pos = ax.imshow(data, cmap=cmap, vmin=height[0], vmax=height[1], origin='lower', \
                            extent=[-frame,frame,-frame,frame])
            fig.colorbar(pos, ax=ax, label=lab)
            
            #add point on sink particle
            #ax.scatter(*(self.cent[85].T*256/5))
            
            ax.set_aspect('equal', 'box')
            ax.set_xlabel("$x_0~[AU]$")
            ax.set_ylabel("$x_1~[AU]$")
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            
            plt.savefig(path)
            plt.close()
        
            
    def kepler_vel(self, center = [], size=15):
        """
        Compute mesh of azimuthal Kepler velocity
        """
        if not np.any(center):
            center = np.array([[0.,0.] for i in range(len(self.rho))])
        
        l = np.shape(self.rho[0])[0]
        for c in center:
            kepler = np.zeros((l, l))
            for i in range(l):
                for j in range(l):
                    x = (20.*i/l-10.-c[0])
                    y = (20.*j/l-10.-c[1])
                    kepler[j,i] = 1/(x**2+y**2)**(1./4)
                    if 1/kepler[j,i]**2 > size: kepler[j,i] = 0.
            kepler[kepler==np.inf] = 0.
            #units
            kepler = kepler*np.sqrt(self.u.mbin* const.GM_sun / (self.u.abin *const.au))
            kepler = kepler.to('km/s')
            self.kepler.append(kepler.value)
    
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
        
        