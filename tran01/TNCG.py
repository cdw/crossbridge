#!/usr/bin/env python
## This file defines the system for the TNCG crossbridge. This 
## crossbridge type has linear springs representing the neck and
## globular regions, and torsional springs representing the 
## converter domain and thick filament attachment site.
##      H   - Head of the myosin
##      |  
##      G   - Globular domain, linear spring
##      |  
##      C   - Converter region, torsional spring
##     /   
##    N     - Neck region, linear spring
##   /     
##==T====== - Thick filament, torsional spring

# FIXME: Find correct scaling factor for kT to be used with our pN forces and 
#        nM scales
# FIXME: Make sure that contour binning isn't off-by-one

from numpy import array, pi, sin, cos, tan, arctan2, sqrt, hypot
import numpy as np
from scipy.optimize import fmin_powell as fmin
import time
import cPickle as pkl
import contour
#import graphXB


class TNCG:
    def __init__(self):
        """Create the values we'll be referencing for the XB"""
        self.Ts = pi/4 # angle (rad) of the connection to the thick fil
        self.Tk = 100  # spring constant of connection to the thick fil
        self.Ns = 5    # rest length of neck region
        self.Nk = 10    # spring constant of neck region
        self.Cs = pi/3+(pi-self.Ts) #rest angle (rad) of the converter domain
        self.Ck = 100  # torsional spring const of converter domain
        self.Cv = (pi/3, pi/3, 1.2*pi/3) # normal and rigor values of Cs
        self.Gs = 3    # rest length of globular domain
        self.Gk = 5    # spring constant of globular domain
        self.Gv = (5, 5, 5) # normal and rigor values of Gs
        self.Bd = 0.55 # dist at which binding becomes likely
        # Current state and identity of XB
        self.rest_conv_and_head_loc() # Set conv_loc and head_loc to rest locs
        self.bound = False
        self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly
        # Diffusion related values
        self.T = 288             #the temperature (in K) that this runs at
        self.K = 1.381 * 10**-23 #Boltzman const (in J/K)
        self.kT = self.K * self.T * 10**21  # kT without pN/nM conversion
        self.Tz = sqrt(2 * pi * self.kT / self.Tk)
        self.Nz = sqrt(2 * pi * self.kT / self.Nk)
        self.Cz = sqrt(2 * pi * self.kT / self.Ck)
        self.Gz = sqrt(2 * pi * self.kT / self.Gk)
        self.Tsig = sqrt(self.kT / self.Tk)
        self.Nsig = sqrt(self.kT / self.Nk)
        self.Csig = sqrt(self.kT / self.Ck)
        self.Gsig = sqrt(self.kT / self.Gk)
        
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
                "T = %02.3fpi : %02.3f \n" %(T/pi, Tu) +
                "N = %02.3f   : %02.3f \n" %(N, Nu) + 
                "C = %02.3fpi : %02.3f \n" %(C/pi, Cu) + 
                "G = %02.3f   : %02.3f \n" %(G, Gu) +
                "Tot energy  = %02.3f" %Total)
    
    def tran01(self, bSite):
        """Given an (x,y) location of an open binding site, bind or not after
        bopping the cross-bridge head to a new location. Return a boolean, True
        for a binding event and False for no binding event.
        """
        # Bop
        hLoc = self.bop()
        # Dist
        dist = hypot(bSite[0]-hLoc[0], bSite[1]-hLoc[1])
        # Bind?
        bProb = np.exp(-dist) # Binding prob is dependent on the exp of dist
        binds = bProb > np.random.rand()
        # Return
        return binds
    
    def bop(self):
        """Bop the xb to a new location, based on an exponential distribution 
        of energies for each independent segment of the crossbridge as 
        determined by Boltzmann's law. Return the new head location.
            
        Justification of technique
        --------------------------
        We are dealing with the energy stored in a spring, this has a 
        dependence on the displacement of the spring thusly:
            U(x) = 1/2 k x^2
        Where k is the spring constant of that spring and x is really x-xs,
        where xs is the rest length of the spring. If we are looking for the 
        probability that the particle or myosin head on the end of the spring 
        will be in a given location at a given time, we can fall back onto 
        Boltzmann's law which tells us that:
            p(x)/p(xs) = exp(-U(x)/kT)
        Where p in this instance is the relative probability of finding the 
        particle at x (as opposed to xs), k is Boltzmann's constant, and T is 
        the system's temperature. This could also be represented as 
        P(x)/P(xs). Where P(x) is the absolute probability of being located at 
        position x. We would like to deal with the absolute probability and so 
        must find Z, a factor that allows us to normalize the distribution 
        such that 
            \int_{-\inf}^\inf 1/Z exp(-U(x)/(kT)) dx = 1
        This is calculated in Mathematica by taking the inverse of the  
        integral
            \int_{-\inf}^\inf 1/Z exp(-.5 k x^2 /(kT)) dx
        This evaluates to 
            Z = \sqrt{2 pi kT / k}
        Meaning we are presented with a probability density function of
            P(x) = \sqrt{k / (2 pi kT)} exp(-(k x^2)/(2 kT))
        Details of implementation
        -------------------------
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
        #Tu = self.kT / self.Tz * np.random.exponential(scale=self.kT)
        #T  = sqrt(2 * Tu / self.Tk) + self.Ts
        #Nu = self.kT / self.Nz * np.random.exponential(scale=self.kT)
        #N  = sqrt(2 * Nu / self.Nk) + self.Ns
        #Cu = self.kT / self.Cz * np.random.exponential(scale=self.kT)
        #C  = sqrt(2 * Cu / self.Ck) + self.Cs
        #Gu = self.kT / self.Gz * np.random.exponential(scale=self.kT)
        #G  = sqrt(2 * Gu / self.Gk) + self.Gs
        T = np.random.normal(self.Ts, self.Tsig)
        N = np.random.normal(self.Ns, self.Nsig)
        C = np.random.normal(self.Cs, self.Csig)
        G = np.random.normal(self.Gs, self.Gsig)
        self.set_conv_and_head_from_segments(T, N, C, G)
        return self.head_loc
    
    def probability(self):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location
        """
        U = self.minimize() #gives energy and sets C to lowest U loc
        T = self.thic_ang()
        N = self.neck_len()
        C = self.conv_ang()
        G = self.glob_len()
        pT = 1/self.Tz * np.exp(-U / self.kT)
        pN = 1/self.Nz * np.exp(-U / self.kT)
        pC = 1/self.Cz * np.exp(-U / self.kT)
        pG = 1/self.Gz * np.exp(-U / self.kT)
        return pT * pN * pC * pG
    
    def minimize(self):
        """Set the conv_loc to minimize the XB's energy 
        and return the newly located minimum energy
        """
        e = lambda (l): self.e_dep_conv(l)
        min_n_len = fmin(e, self.conv_loc, disp=0)
        return self.energy()
    
    def e_dep_conv(self, conv):
        """Update the conv_loc and return the xb energy"""
        self.conv_loc = (conv[0], conv[1])
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
    
    def set_conv_and_head_from_segments(self, T, N, C, G):
        """Set the converter and head loc from passed segment values"""
        self.conv_loc = (N * cos(T),
                         N * sin(T))
        x = self.conv_loc[0] + G * cos(C + T - pi)
        y = self.conv_loc[1] + G * sin(C + T - pi)
        self.head_loc = (x, y)
    



## Begin the script that will produce the matrix of stored rates
# Instantiate the xb
xb = TNCG()

# Bebop that xb head
trials = 1000000
x = np.zeros(trials)
y = np.zeros(trials)
for i in range(trials):
    x[i], y[i] = xb.bop()
# Save it
dataout = open('bop_TNCG.pkl', 'wb')
pkl.dump([xb, trials, x, y], dataout, 2)
dataout.close()
# Plot it all pretty like
contour.title = "Location histogram of a diffusing TNCG crossbridge head"
contour.xlabel = "X loc of head (nm)"
contour.ylabel = "Y loc of head (nm)"
contour.hexbin(x, y, [-2, 15, -2, 15])

# Find transition rates
trials = 2000 #per (x,y) location
x_locs = np.arange(-2, 15, .1)
y_locs = np.arange(-2, 15, .1)
rates = np.zeros((y_locs.size, x_locs.size))
# Cycle through head locations, collecting trans rates
for n,x in enumerate(x_locs):
    print('On col %(c)04d of %(t)04d' %{'c':n, 't':x_locs.size})
    for m,y in enumerate(y_locs):
        for i in range(trials):
            rates[m, n] += xb.tran01((x, y))
# Save it
dataout = open('trans01_TNCG.pkl', 'wb')
pkl.dump([xb, trials, x_locs, y_locs, rates], dataout, 2)
dataout.close()
# Plot the output
contour.title = ("0->1 transition rate of an TNCG crossbridge as a /n" \
                "function of head location")
contour.xlabel = "Location of XB head (nm)"
contour.ylabel = "Location of XB head (nm)"
contour.levels = np.array([.2, .4, .6, .8, .9])*rates.max()
contour.contour(x_locs, y_locs, rates)
# Display all the graphs produced in this script
contour.show()