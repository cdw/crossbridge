#!/usr/bin/env python
## This file defines the system for the xNCG crossbridge. This 
## crossbridge type has linear springs representing the neck and
## globular regions, and a torsional spring representing the 
## converter domain.
##      H   - Head of the myosin
##      |  
##      G   - Globular domain, linear spring
##      |  
##      C   - Converter region, torsional spring
##     /   
##    N     - Neck region, linear spring
##   /     
##========= - Thick filament


from numpy import array, pi, sin, cos, tan, arctan2, sqrt, hypot
import numpy as np
from scipy.optimize import fmin_powell as fmin
import time
import contour
import graphXB


class xNCG():
    def __init__(self):
        """Create the values we'll be referencing for the XB"""
        self.Ts = pi/4 # angle (rad) of the connection to the thick fil
        self.Tk = 100  # spring constant of connection to the thick fil
        self.Ns = 5    # rest length of neck region
        self.Nk = 10    # spring constant of neck region
        self.Nv = (6, 6, 6) # normal and rigor values of Ns
        self.Cs = pi/3+(pi-self.Ts) #rest angle (rad) of the converter domain
        self.Ck = 100  # torsional spring const of converter domain
        self.Cv = (pi/3, pi/3, 1.2*pi/3) # normal and rigor values of Cs
        self.Gs = 3    # rest length of globular domain
        self.Gk = 5    # spring constant of globular domain
        self.Gv = (5, 5, 5) # normal and rigor values of Gs
        self.Fm = 0    # mean of forces exerted on myosin heads
        self.Fv = 10.0 # variance of forces exerted on myosin heads
        self.Bd = 0.55 # dist at which binding becomes likely
        # Current state and identity of XB
        self.rest_conv_loc() # Sets conv_loc to rest location
        self.rest_head_loc() # Sets head_loc to rest location
        self.bound = False
        self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly
        
    def __repr__(self):
        """Return a string representation of the XB"""
        # Angles and lengths
        T = self.thic_ang()
        N = self.neck_len()
        G = self.glob_len()
        C = self.conv_ang()
        # Energies
        Tu = 0.5 * self.Tk * (T-self.Ts)**2
        Nu = 0.5 * self.Nk * (N-self.Ns)**2 
        Cu = 0.5 * self.Ck * (C-self.Cs)**2
        Gu = 0.5 * self.Gk * (G-self.Gs)**2
        return ("Angles/Lengths and Energies\n" +
                "===========================\n" +
                "x = ang/len : energy\n" +
                "T = %02.3fpi : %02.3f (fixed)\n" %(T/pi, Tu) +
                "N = %02.3f   : %02.3f \n" %(N, Nu) + 
                "C = %02.3fpi : %02.3f \n" %(C/pi, Cu) + 
                "G = %02.3f   : %02.3f" %(G, Gu))
        
    def minimize(self):
        """Set the conv_loc to minimize the XB's energy 
        and return the newly located minimum energy
        """
        e = lambda (l): self.e_dep_conv(l)
        min_n_len = fmin(e, self.Ns, disp=0)
        return self.energy()
        
    def e_dep_conv(self, N_len):
        """Update the conv_loc and return the xb energy"""
        N_len = float(N_len) # Was a numpy array, weirding TK
        self.conv_loc = (N_len * cos(self.Ts),
                         N_len * sin(self.Ts))
        return self.energy()
    
    def energy(self):
        """Return the energy of the xb without altering any positions"""
        N = self.neck_len()
        G = self.glob_len()
        C = self.conv_ang()
        return (0.5 * self.Nk * (N-self.Ns)**2 +
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
    
    def rest_conv_loc(self):
        """Set the converter loc to its rest location"""
        self.conv_loc = (self.Ns * cos(self.Ts),
                         self.Ns * sin(self.Ts))
    
    def rest_head_loc(self):
        """Set the head loc to its rest location"""
        x = self.conv_loc[0] + self.Gs * cos(self.Cs + self.Ts - pi)
        y = self.conv_loc[1] + self.Gs * sin(self.Cs + self.Ts - pi)
        self.head_loc = (x, y)

## b = xNCG()
## window = graphXB.graphXB()
## window.draw_XB(b, free=(False, True, True, True))
## for i in np.arange(14,0,-.25):
##     b.head_loc = (3.7, i)
##     b.minimize()
##     print(b)
##     window.update_XB(b)
##     time.sleep(.1)

## Begin the script that will produce the matrix of stored energies
print("This might take a while")
x_locs = np.arange(-3, 13, .1) 
y_locs = np.arange(0, 16, .1)
energs = np.zeros((y_locs.size, x_locs.size))
# Instantiate the xb
xb = xNCG()
# Cycle through and collect all the energies
n = [0,0]
for y in y_locs:
    for x in x_locs:
        xb.head_loc = (x,y)
        energs[n[0], n[1]] = xb.minimize()
        n[1] = n[1] + 1
    n[0] = n[0] + 1
    n[1] = 0
contour.title = "Energy level of an xNCG crossbridge at different head locations"
contour.xlabel = "Location of XB head (nm)"
contour.ylabel = "Location of XB head (nm)"
contour.levels = [1, 4, 7, 10, 13]
contour.contour(x_locs, y_locs, energs)
