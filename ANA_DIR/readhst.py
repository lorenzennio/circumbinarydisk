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
    
    def accperorbit(self):
        self.data["accm1_orbit"] = np.array([j - i for i, j in zip(self.data["accm1"][:-1], self.data["accm1"][1:])])
        self.data["accm2_orbit"] = np.array([j - i for i, j in zip(self.data["accm2"][:-1], self.data["accm2"][1:])])