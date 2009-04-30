#!/usr/bin/env python
## This file defines the system for the xNxG crossbridge. This 
## crossbridge type has linear springs representing the neck and
## globular regions. These springs are angled such that they replicate
## a version of the one-dimensional model. There is some ambiguity as 
## to how we represent the one-dimensional model in a two-dimensional 
## space, with no principle present to guide use between the choices 
## of having a globular domain spring with no stiffness (allowing the 
## pseudo-one-dimensional model to range freely over all lattice 
## spacings) or having a globular domain with infinite stiffness 
## (fixing the head at one particular lattice spacing). I am leaning 
## towards the all-lattice spacings interpretation to indicate that the
## one-dimensional crossbridge is independent of lattice spacing, but the 
## fixed-spacing version does accuretly conveigh that the one-dimensional
## crossbridge is simply incapable of looking at lattice spacings.
##      H   - Head of the myosin
##      |  
##      G   - Globular domain, linear spring/fixed length
##      |  
##      C   - Converter region, fixed angle
##     /   
##    N     - Neck region, linear spring
##   /     
##==T====== - Thick filament, fixed angle


from numpy import array, pi, sin, cos, tan, arctan2, sqrt, hypot
import numpy as np
from scipy.optimize import fmin_powell as fmin
import time
import contour


class xNxG():
    def __init__(self):
        """Create the values we'll be referencing for the XB"""
        self.Ts = 0 # angle (rad) of the connection to the thick fil
        self.Tk = 100  # spring constant of connection to the thick fil
        self.Ns = 5    # rest length of neck region
        self.Nk = 5 #was 10    # spring constant of neck region
        self.Cs = pi/2+(pi-self.Ts) #rest angle (rad) of the converter domain
        self.Ck = 100  # torsional spring const of converter domain
        self.Cv = (pi/3, pi/3, 1.2*pi/3) # normal and rigor values of Cs
        self.Gs = 3    # rest length of globular domain
        self.Gk = .001 # spring constant of globular domain
        self.Gv = (5, 5, 5) # normal and rigor values of Gs
        self.Bd = 0.55 # dist at which binding becomes likely
        # Current state and identity of XB
        self.rest_conv_and_head_loc() # Set conv_loc and head_loc to rest locs
        self.bound = False
        self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly
        # Diffusion related values
        self.T = 288             #the temperature (in K) that this runs at
        self.K = 1.381 * 10**-23 #Boltzman const (in J/K)
        self.kT = self.K * self.T * 10**23 # kT with pN/nM conversion
        self.Tz = sqrt((2 * pi * self.kT) / self.Tk)
        self.Nz = sqrt((pi * self.kT) / (2 * self.Nk))
        self.Cz = sqrt((2 * pi * self.kT) / self.Ck)
        self.Gz = sqrt((pi * self.kT) / (2 * self.Gk))
        
    def __repr__(self):
        """Return a string representation of the XB"""
        # Angles and lengths
        T = self.thic_ang()
        N = self.neck_len()
        C = self.conv_ang()
        G = self.glob_len()
        # Energies
        Tu = 0.5 * self.Tk * (T-self.Ts)**2
        Nu = 0.5 * self.Nk * (N-self.Ns)**2 
        Cu = 0.5 * self.Ck * (C-self.Cs)**2
        Gu = 0.5 * self.Gk * (G-self.Gs)**2
        Total = Tu + Nu + Cu + Gu
        return ("Angles/Lengths and Energies\n" +
                "===========================\n" +
                "x = ang/len : energy\n" +
                "T = %02.3fpi : %02.3f fixed \n" %(T/pi, Tu) +
                "N = %02.3f   : %02.3f \n" %(N, Nu) + 
                "C = %02.3fpi : %02.3f fixed \n" %(C/pi, Cu) + 
                "G = %02.3f   : %02.3f \n" %(G, Gu) +
                "Tot energy  = %02.3f" %Total)
        
    def probability(self):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location"""
        U = self.minimize() #gives energy and sets C to lowest U loc
        N = self.neck_len()
        C = self.conv_ang()
        G = self.glob_len()
        pN = 1/self.Nz * np.exp(-U / self.kT)
        pC = 1/self.Cz * np.exp(-U / self.kT)
        pG = 1/self.Gz * np.exp(-U / self.kT)
        return pN * pC * pG
        
    def minimize(self):
        """Set the conv_loc to minimize the XB's energy 
        and return the newly located minimum energy
        """
        (x, y) = self.head_loc
        self.conv_loc = (x, self.conv_loc[1])
        return self.energy()
        
    def energy(self):
        """Return the energy of the xb without altering any positions"""
        T = self.thic_ang()
        N = self.neck_len()
        C = self.conv_ang()
        G = self.glob_len()
        return (0.5 * self.Tk * (T-self.Ts)**2 + 
                0.5 * self.Nk * (N-self.Ns)**2 +
                0.5 * self.Ck * (C-self.Cs)**2 + 
                0.5 * self.Gk * (G-self.Gs)**2)
        
    def glob_len(self):
        """Return the globular length at the current head_loc"""
        x = self.head_loc[0] - self.conv_loc[0]
        y = self.head_loc[1] - self.conv_loc[1]
        return hypot(x, y)
        
    def conv_ang(self):
        """Return the converter angle at the current head_loc"""
        x = self.head_loc[0] - self.conv_loc[0]
        y = self.head_loc[1] - self.conv_loc[1]
        return arctan2(y, x) + pi - self.thic_ang()
        
    def neck_len(self):
        """Return the neck length at the current conv_loc"""
        return hypot(self.conv_loc[0], self.conv_loc[1])
        
    def thic_ang(self):
        """Return the angle of the neck's attachment to the thick filament"""
        x = self.conv_loc[0]
        y = self.conv_loc[1]
        return arctan2(y, x)
        
    def rest_conv_and_head_loc(self):
        """Set the converter and head loc to their rest locations"""
        self.conv_loc = (self.Ns * cos(self.Ts),
                         self.Ns * sin(self.Ts))
        x = self.conv_loc[0] + self.Gs * cos(self.Cs + self.Ts - pi)
        y = self.conv_loc[1] + self.Gs * sin(self.Cs + self.Ts - pi)
        self.head_loc = (x, y)

## b = xNxG()
## window = graphXB.graphXB()
## window.draw_XB(b, free=(True, True, True, True))
## for i in np.arange(14,0,-.25):
##     b.head_loc = (12-i, i)
##     b.minimize()
##     print(b)
##     window.update_XB(b)
##     time.sleep(.01)
## for i in np.arange(14,0,-.25):
##     b.head_loc = (4, i)
##     b.minimize()
##     print(b)
##     window.update_XB(b)
##     time.sleep(.01)
## for i in np.arange(14,-3,-.25):
##     b.head_loc = (i, 5)
##     b.minimize()
##     print(b)
##     window.update_XB(b)
##     time.sleep(.01)

# Begin the script that will produce the matrix of stored probabilities
x_locs = np.arange(-3, 13, .1) 
y_locs = np.arange(0, 16, .1)
probs = np.zeros((y_locs.size, x_locs.size))
# Instantiate the xb
xb = xNxG()
pT = time.time()
cT = time.time()
# Cycle through and collect all the probabilities
for n,x in enumerate(x_locs):
    for m,y in enumerate(y_locs):
        xb.head_loc = (x,y)
        probs[m, n] = xb.probability()
    # Tell me how much time is left, about
    cT = time.time()
    rT = (cT-pT)*(x_locs.size - (n+1))
    print('On col %(c)04d of %(t)04d, about %(m)02d:%(s)02d left' \
          %{'c':n, 't':x_locs.size, 'm':rT//60, 's':rT%60})
    pT = cT
# Normalize the probabilities
min = np.min(probs)
max = np.max(probs)
probs = (probs - min)/(max-min)
contour.title = "Probability of an xBxG crossbridge being\n found at a given head locations"
contour.xlabel = "Location of XB head (nm)"
contour.ylabel = "Location of XB head (nm)"
contour.levels = [.9, .95, .98, .99] 
contour.contour(x_locs, y_locs, probs)