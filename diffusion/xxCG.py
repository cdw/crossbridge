## 
# Dave Williams, 20090315
# This script executes to produce a large matrix of the xxCG model's
# potential energies and transition rates
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
        self.kT = self.K * self.T * 10**23 # kT with pN/nM conversion
        self.Gz = sqrt((pi * self.kT) / (2 * self.Gk))
        self.Cz = sqrt((2 * pi * self.kT) / self.Ck)

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
    

## Begin the script that will produce the matrix of stored energies
x_locs = np.arange(-3, 13, .2) 
y_locs = np.arange(0, 16, .2)
probs = np.zeros((y_locs.size, x_locs.size))
# Instantiate the xb
xb = xxCG()
# Cycle through and collect all the probabilities
n = [0,0]
for y in y_locs:
    for x in x_locs:
        xb.head_loc = (x,y)
        probs[n[0], n[1]] = xb.probability()
        n[1] = n[1] + 1
    n[0] = n[0] + 1
    n[1] = 0
# Normalize the probabilities
min = np.min(probs)
max = np.max(probs)
probs = (probs - min)/(max-min)
contour.title = "Probability of an xxCG crossbridge being\n found at a given head locations"
contour.xlabel = "Location of XB head (nm)"
contour.ylabel = "Location of XB head (nm)"
contour.levels = [.9, .95, .98, .99, .999] 
contour.contour(x_locs, y_locs, probs)
print(probs)