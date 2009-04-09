## 
# Dave Williams, 20090315
# This script executes to produce a large matrix of the xxCG model's
# probabilities of diffusing to a given location
# 
##

# FIXME: Find correct scaling factor for kT to be used with our pN forces and nM scales

from numpy import array, pi, sin, cos, tan, arctan2, sqrt, hypot
import numpy as np
import contour

class xxCG():
    def __init__(self):
        """Create the values we'll be referencing for the XB"""
        self.Cs = pi/3 # rest angle of converter domain
        self.Ck = 100  # torsional spring const of converter domain
        self.Cv = (pi/3, pi/3, 1.2*pi/3) # normal and rigor values of Cs
        self.Gs = 10.5 # rest length of globular domain
        self.Gk = 5    # spring constant of globular domain
        self.Gv = (5, 5, 5) # normal and rigor values of Gs
        self.Bd = 0.55 # dist at which binding becomes likely
        # Current state and identity of XB
        self.rest_head_loc() # Sets head_loc to rest location
        self.bound = False
        self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly
        # Diffusion related values
        self.T = 288             #the temperature (in K) that this runs at
        self.K = 1.381 * 10**-23 #Boltzman const (in J/K)
        self.kT = self.K * self.T * 10**21 # kT with pN/nM conversion
        self.Gz = sqrt((pi * self.kT) / (2 * self.Gk))
        self.Cz = sqrt((2 * pi * self.kT) / self.Ck)
    
    def bop(self):
        """Bop the xb to a new location, based on an exponential distribution 
        of energies for each independent segment of the crossbridge as 
        determined by Boltzmann's law. Return the new head location.
        
        This makes use of numpy's exponential distribution:
            numpy.random.exponential(scale=1.0, size=None)
        Where scale is B such that the probability density function of the 
        resulting distribution is: 
            f(x, B) = 1/B exp(-x/B)
        If we feed numpy's exponential distribution B=kT, and take x as 
        energy U, we have a distribution of 
            f(U) = 1/kT exp(-U/kT)
        Now, we will have to scale this output as we need our distributions 
        to have probability density functions of:
            f(U) = 1/Z exp(-U/kT) 
                Where Z is in the format of:
                Zc = sqrt(2 pi kT / Kc) or
                Zg = sqrt(pi kT / (2 Kg))
        This means we need to multiply the numpy output by kT/Z, which, in 
        these cases would work out to:
            kT/Zc = sqrt(Kc kT/ (2 pi)) or
            kT/Zg = sqrt(2 Kg kT / pi)
        This gives us the energy distributions we want, from which we can 
        backtrack to get the length or angle values for each component.
        """
        Cu = self.kT / self.Cz * np.random.exponential(scale=self.kT)
        C  = sqrt(2 * Cu / self.Ck) + self.Cs
        Gu = self.kT / self.Gz * np.random.exponential(scale=self.kT)
        G  = sqrt(2 * Gu / self.Gk) + self.Gs
        self.set_conv_and_head_from_segments(C, G)
        return self.head_loc
        
    def probability(self):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location"""
        G = self.glob_len()
        C = self.conv_ang()
        U = self.energy()
        pG = 1/self.Gz * np.exp(-U / self.kT)
        pC = 1/self.Cz * np.exp(-U / self.kT)
        return pG * pC
    
    def energy(self):
        """Return the energy stored in the XB, 
        given the current head_loc
        """
        G = self.glob_len()
        C = self.conv_ang()
        return (0.5 * self.Gk * (G-self.Gs)**2 + 
        0.5 * self.Ck * (C-self.Cs)**2)

    def glob_len(self):
        """Return the globular length at the current head_loc"""
        return hypot(self.head_loc[0], self.head_loc[1])
    
    def conv_ang(self):
        """Return the converter angle at the current head_loc"""
        return arctan2(self.head_loc[1], self.head_loc[0])

    def rest_head_loc(self):
        """Set the head loc to its rest location"""
        self.head_loc = (self.Gs * cos(self.Cs),
                         self.Gs * sin(self.Cs))
        
    def set_conv_and_head_from_segments(self, C, G):
        """Set the converter and head loc from passed segment values"""
        self.head_loc = (G * cos(C),
                         G * sin(C))

    

## Begin the script that will produce the matrix of stored probabilities
trials = 2000000
x_locs = np.arange(-15, 15, .2) 
y_locs = np.arange(-5, 25, .2)
probs = np.zeros((y_locs.size, x_locs.size))
hits = np.zeros((y_locs.size, x_locs.size))
# Instantiate the xb
xb = xxCG()
# Cycle through all iterations and collect the head locations
for i in range(trials):
    loc = xb.bop()
    x_ind = np.searchsorted(x_locs, loc[0]) - 1 #FIXME Check that there is not  
    y_ind = np.searchsorted(y_locs, loc[1]) - 1 #    an off by one error here
    hits[x_ind, y_ind] = hits[x_ind, y_ind] + 1
# Normalize the hit likelihood
min = np.min(hits)
max = np.max(hits)
hits = (hits - min)/(max-min)
contour.title = "Probability of an xxCG crossbridge being\n found at a given head locations"
contour.xlabel = "Location of XB head (nm)"
contour.ylabel = "Location of XB head (nm)"
contour.levels = [.1, .3, .5, .7, .9] 
contour.contour(x_locs, y_locs, hits)
print(hits)