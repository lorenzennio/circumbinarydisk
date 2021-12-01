import numpy as np

class units:
  def __init__(self, mbin, abin):
    #cgs constants
    self.msol = 1.98847e33 #g
    self.mh   = 1.673534e-24 #g
    self.mu   = 2.3
    self.gn   = 6.6743e-8 #cm^3 g^-1 s^-2
    self.kb   = 1.3807e-16 #cm^2 g s^-2 K^-1
    self.au   = 1.495978707e13 #cm
    self.mbin    = mbin
    self.abin    = abin
    self.tbin    = np.sqrt((self.abin*self.au)**3/(self.mbin*self.msol*self.gn))
    self.densref = self.mbin*self.msol/(self.abin*self.au)**3
    self.presref = self.mbin*self.msol/(self.tbin**2 * self.abin*self.au) #g/(s^2 cm)
    self.tempref = self.mh * self.mu /self.kb*(self.mbin*self.msol*self.gn/(self.abin*self.au))
