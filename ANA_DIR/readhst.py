import os, sys
import h5py
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import units as un
import re

class readhst:
    """
    class to load and process the history file 
    """
    
    def __init__(self, file):
        self.keys = self.genkeys(file)
        self.data = dict(zip(self.keys, np.genfromtxt(file).T))
        self.orbits()
        self.accperorbit()
        
    def genkeys(self, file):
        f = open(file, 'r')
        lines = f.readlines()
        f.close()
        l = re.sub("# ", "", lines[1])
        l = re.sub(" +", " ", l)
        l = re.sub(r'\[[0-9]+\]=', "", l)
        l = re.sub(" +", " ", l)
        keys = np.array(l.split(" "))[:-1]
        return keys
    
    def orbits(self):
        self.data["orbits"] = self.data["time"]/0.62832
    
    def accperorbit(self):
        self.data["accm1_orbit"] = np.array([j - i for i, j in zip(self.data["accm1"][:-1], self.data["accm1"][1:])])
        self.data["accm2_orbit"] = np.array([j - i for i, j in zip(self.data["accm2"][:-1], self.data["accm2"][1:])])
        
    def accperyear(self, period):
        self.data["accm1_year"] = np.array([j - i for i, j in zip(self.data["accm1"][:-1], self.data["accm1"][1:])])/(period/31557600)
        self.data["accm2_year"] = np.array([j - i for i, j in zip(self.data["accm2"][:-1], self.data["accm2"][1:])])/(period/31557600)
        
    def years(self, period):
        self.data["years"] = self.data["orbits"]*period/31557600
        self.accperyear(period)
        
    def plot(self, mbin):
        
        fs = 16
        xval = self.data['orbits']
        xlab = 'Time [orbits]'
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5), constrained_layout=True)

        ax1.plot(xval, self.data['accm1']*mbin)
        ax1.plot(xval, self.data['accm2']*mbin)
        ax1.set_xlabel(xlab, fontsize=fs)
        ax1.set_xscale('log')
        ax1.set_ylabel('Cumulative accreted mass [M_solar]', fontsize=fs)
        ax1.tick_params(labelsize=fs)
        ax1.legend(["m1", "m2"], fontsize=fs)

        ax2.plot(xval[1:], self.data['accm1_orbit']*mbin)
        ax2.plot(xval[1:], self.data['accm2_orbit']*mbin)
        ax2.set_xlabel(xlab, fontsize=fs)
        ax2.set_xscale('log')
        ax2.set_ylabel('Accreted mass / orbit [M_solar/orbit]', fontsize=fs)
        ax2.tick_params(labelsize=fs)
        ax2.legend(["m1", "m2"], fontsize=fs)

        ax3.plot(xval, self.data['m1']*mbin)
        ax3.plot(xval, self.data['m2']*mbin)
        ax3.set_xlabel(xlab, fontsize=fs)
        ax3.set_xscale('log')
        ax3.set_ylabel('Stellar disk mass (25AU radius) [M_solar]', fontsize=fs)
        ax3.tick_params(labelsize=fs)
        ax3.legend(["m1", "m2"], fontsize=fs)
        
        path = os.path.join(os.getcwd(), 'plots/hst.png')
        plt.savefig(path)