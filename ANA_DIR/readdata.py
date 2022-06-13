import os, sys
import h5py
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import units as un
from astropy import constants as const
from skimage.measure import block_reduce
import plotly.graph_objects as go
import cv2

import yt
from yt.funcs import mylog
mylog.setLevel(40)

class load:
    """
    define variables and read in the data from the athena++ simulation output
    
    prim_files  = files with primitive data from athena++
    units       = units class instance
    level       = refinement level
    dims        = dimensions (2/3)
    """
    
    def __init__(self, prim_files, units, level, dims):
        print("Initializing data readin")
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

        #remove levels from 2d
        if dims==2:
            self.rmlevels()
        elif dims ==3:
            pass
        else:
            raise ValueError('Can only process 2d or 3d.')

        #compute temperature
        try:
            self.temp = self.press/self.u.presref / (self.rho/self.u.densref) * self.u.tempref
        except:
            pass
        
        #compute velocity magnitude
        self.v = np.sqrt(self.v1**2 + self.v2**2 + self.v2**2)

        
        print("Processed all data")



    def load(self, level):
        """
        function to load the data from the hdf5 files into arrays
        
        level = number of refinements
        """
        unit_base= {"length_unit":(self.u.abin,"AU"), "time_unit":(self.u.tbin,"s"), "mass_unit":(self.u.mbin,"Msun")}
        for fprim in self.prim_files:
            ds_prim = yt.load(fprim, units_override = unit_base)
            
            all_data_level_1_prim = ds_prim.covering_grid(level=level, \
                                            left_edge=ds_prim.domain_left_edge, \
                                            dims=ds_prim.domain_dimensions*2**level)
                                            #dims = total dimensions 2**(refinement level)

            
            rho = all_data_level_1_prim["rho"].in_units("g/cm**3").to_ndarray()
            try:
                press = all_data_level_1_prim["press"].in_units("g/(cm*s**2)").to_ndarray()
            except:
                pass
            v1 = all_data_level_1_prim["vel1"].in_units("km/s").to_ndarray()
            v2 = all_data_level_1_prim["vel2"].in_units("km/s").to_ndarray()
            v3 = all_data_level_1_prim["vel3"].in_units("km/s").to_ndarray()
            
            self.rho.append(rho)
            try:
                self.press.append(press)
            except:
                pass
            self.v1.append(v1)
            self.v2.append(v2)
            self.v3.append(v3) 
            
            print("Read " + fprim)
            

    def tonp(self):
        """
        function to convert to numpy arrays
        """
        self.rho = np.array(self.rho)
        try:
            self.press = np.array(self.press)
        except:
            pass
        self.v1 = np.array(self.v1)
        self.v2 = np.array(self.v2)
        self.v3 = np.array(self.v3)

    def rmlevels(self):
        """
        function to remove the outer refinement levels form the read data
        """
        self.rho = self.rho[:,:,:,0]
        try:
            self.press = self.press[:,:,:,0]
        except:
            pass
        self.v1 = self.v1[:,:,:,0]
        self.v2 = self.v2[:,:,:,0]
        self.v3 = self.v3[:,:,:,0]
            
#----------------------------------------------------------------

class cartesian:
    """
    class to define cartesian, cylindrical and polar coordinates on a meshgrid with the same size as the input data
    """
    def __init__(self, mesh, codelength):
        #cartesian
        xylen = len(mesh)
        self.x = np.array([[i-(xylen-1)/2 for i in range(xylen)] for j in range(xylen)])
        self.x *= (2*codelength/xylen)
        self.y = -self.x.T
        
        #3d - TODO
        if len(mesh.shape)==3:
            zlen = len(mesh[0,0])
            self.z = np.array([[i-(zlen-1)/2 for i in range(zlen)] for j in range(zlen)])
            self.z *= (2*codelength/xylen)

class cylindrical(cartesian):
    """
    class to define 2d cylidrical polar coordinates on a meshgrid with the same size as the input data
    """
    def __init__(self, mesh, codelength):
        super().__init__(mesh, codelength)
        self.r = np.sqrt(self.x**2+self.y**2)
        self.phi = np.arctan2(self.y,self.x)

        
class spherical(cylindrical):
    """
    class to define 3d spherical polar coordinates on a meshgrid with the same size as the input data
    """
    def __init__(self, mesh, codelength):
        super().__init__(mesh, codelength)

class coordinates(spherical, cylindrical, cartesian):
    def __init__(self, mesh, codelength):
        super().__init__(mesh, codelength)
#----------------------------------------------------------------

class plot:
    """
    class for specific plotting routines
    """
    
    
    def plot(self, data, path, low, hig, frame, lab, lim, log=False):
        cmap="plasma"
        fs = 16
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        if log:
            pos = ax.imshow(data, cmap=cmap, origin='lower', extent=[-frame,frame,-frame,frame], norm=matplotlib.colors.LogNorm(vmin=low, vmax=hig))
        else:
            pos = ax.imshow(data, cmap=cmap, origin='lower', vmin=low, vmax=hig, extent=[-frame,frame,-frame,frame])
        
        cbar = fig.colorbar(pos, ax=ax)
        
        cbar.ax.tick_params(labelsize=fs)
        cbar.set_label(lab, size=fs)
        
        ax.set_aspect('equal', 'box')
        ax.set_xlabel("$x_0~[AU]$", fontsize=fs)
        ax.set_ylabel("$x_1~[AU]$", fontsize=fs)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(fs)
        
        plt.savefig(path)
        plt.close()
        
    def fieldlines(self, coord, vx, vy, v, path, frame, lab, lim, step):
        cmap="plasma"
        lab = "$v~[km/s]$"

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        
        pos = ax.quiver(coord.x[::step, ::step],coord.y[::step, ::step],\
                  vx[::step, ::step], vy[::step, ::step], v[::step, ::step],\
                 cmap=cmap)
        fig.colorbar(pos, ax=ax, label=lab)

        ax.set_aspect('equal', 'box')
        ax.set_xlabel("$x_0~[AU]$")
        ax.set_ylabel("$x_1~[AU]$")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)

        plt.savefig(path)
        plt.close()
        
    def video(self, path):
        #VIDEO
        images = []
        os.chdir(path)
        files = np.sort([f for f in os.listdir(os.getcwd()) if f.endswith("rho.png")])
        for f in files[:1000]:
            images.append(cv2.imread(f))
        height,width,layers=np.shape(images[1])
        fourcc = cv2.VideoWriter_fourcc(*"XVID")
        #fourcc = cv2.VideoWriter_fourcc(*'DIVX')
        video=cv2.VideoWriter('evolution.avi',fourcc, 5,(width,height), True)
        for i in images:
            video.write(i)
        cv2.destroyAllWindows()
        video.release()
        os.chdir("..")
        print("Created video.")

#----------------------------------------------------------------  

class data2d(load, plot):
    """
    class for 2d specific data and plotting methods
    """
    
    def __init__(self, prim_files, units, level):
        dims = 2
        super().__init__(prim_files, units, level, dims)
    
    def polarvel(self):
        #decompose velocity in radial and azimuthal parts
        
        #compute polar angle on mesh
        l = np.shape(self.phi)[0]
        X = np.array([[i-l/2+1./2 for i in range(l)] for j in range(l)])
        Y = X.T
        
        self.phi = np.arctan2(Y,X)
        
        #compute polar velocity
        self.vr   = self.v1*np.cos(self.phi) + self.v2*np.sin(self.phi)
        self.vphi = -self.v1*np.sin(self.phi) + self.v2*np.cos(self.phi)
    
    def draw(self, simdata, fname, lab="", scale = 1., low=None, hig=None, log=False):
        frame = self.u.abin * 10.
        lim = frame*scale
        for filename, data in zip(self.prim_files, simdata):
            path = os.path.join(os.getcwd(), 'plots/'+ filename.replace(".athdf", fname))
            
            if not low:
                lo = np.min(data)
            else:
                lo = low
            if not hig:
                hi = np.max(data)
            else:
                hi = hig
            
            self.plot(data, path, lo, hi, frame, lab, lim, log)
            
            print("Drew - " + filename)
            
    def vectorfield(self, fname, lab="$v~[km/s]$", scale = 1., step=10):
        frame = self.u.abin * 10.
        lim = frame*scale
        
        coord = coordinates(self.v[0], frame)

        for filename, vx, vy, v in zip(self.prim_files, self.v1, self.v2, self.v):
            path = os.path.join(os.getcwd(), 'plots/'+ filename.replace(".athdf", fname))
            
            self.fieldlines(coord, vx, vy, v, path, frame, lab, lim, step)
            
            print("Drew field lines - " + filename)
#----------------------------------------------------------------  
    
class data3d(load, plot):
    """
    class for 3d specific data and plotting methods
    """
    def __init__(self, prim_files, units, level):
        dims = 3
        super().__init__(prim_files, units, level, dims)
    
    
    def polar_vel(self):
        """
        decompose velocity in radial and azimuthal parts
        
        TODO
        """
        pass
        
    def plot3d(self, simdata, fname, lab="", low=None, hig=None):
        """
        function to plot the data in 3d
        """
        if not low:
            low = np.min(simdata)
        if not hig:
            hig = np.max(simdata)
            
        for filename, data in zip(self.prim_files, simdata):
            #mean = np.mean(data)
            path = os.path.join(os.getcwd(), 'plots/'+ filename.replace(".athdf", fname))
            
            blocksize = int(np.shape(data)[0] / 64)
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
            print("Drew 3d: " + filename)
        
    def column(self, simdata, function, fname, lab="", scale = 1., low=None, hig=None, log=False):
        """
        function to plot the top view colums of the data
        """
        frame = self.u.abin*10.
        lim = frame*scale
        
        for filename, data in zip(self.prim_files, simdata):
            path = os.path.join(os.getcwd(), 'plots/column/'+ filename.replace(".athdf", fname))
            
            zcolumn = np.shape(data)[-1]
            column = block_reduce(data, block_size=(1, 1, zcolumn), func=function)[:,:,0]
            
            if not low:
                lo = np.min(column)
            else:
                lo = low
            if not hig:
                hi = np.max(column)
            else:
                hi = hig
            
            self.plot(column, path, lo, hi, frame, lab, lim, log)
            
            print("Drew " + filename)
            
#----------------------------------------------------------------  

def masses(dens, rmin=0, rmax=np.inf):
    """
    DEPRECATED - use simulation output instead!!!
    compute masses within radius for 2d simulations
    """
    masses = np.array([])
    
    xyl = len(dens[0])
    coords = cylindrical(dens[0], 10)
    rflat = coords.r.flatten()
    
    dx = 2*10./xyl * 50.*const.au.to('cm').value
    dy = 2*10./xyl * 50.*const.au.to('cm').value
    dz = 2*0.5 * 2.*const.au.to('cm').value
    
    for d in dens:
        m=0
        dflat = d.flatten()
        for rf, df in zip(rflat, dflat):
            if rf >= rmin and rf <= rmax:
                m+=df
        m *= dx*dy*dz
        masses = np.append(masses, m)
    return masses