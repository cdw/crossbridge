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

# FIXME: Find correct scaling factor for kT to be used with our pN forces and nM scales


from numpy import array, pi, sin, cos, tan, arctan2, sqrt, hypot
import numpy as np
from scipy.optimize import fmin_powell as fmin
import time
import contour
#import graphXB


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
        self.Bd = 0.55 # dist at which binding becomes likely
        # Current state and identity of XB
        self.rest_conv_loc() # Sets conv_loc to rest location
        self.rest_head_loc() # Sets head_loc to rest location
        self.bound = False
        self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly
        # Diffusion related values
        self.T = 288             #the temperature (in K) that this runs at
        self.K = 1.381 * 10**-23 #Boltzman const (in J/K)
        self.kT = self.K * self.T * 10**21 # kT with pN/nM conversion
        self.Tz = sqrt((2 * pi * self.kT) / self.Tk)
        self.Nz = sqrt((pi * self.kT) / (2 * self.Nk))
        self.Cz = sqrt((2 * pi * self.kT) / self.Ck)
        self.Gz = sqrt((pi * self.kT) / (2 * self.Gk))
        
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
        Nu = self.kT / self.Nz * np.random.exponential(scale=self.kT)
        N  = sqrt(2 * Nu / self.Nk) + self.Ns
        Cu = self.kT / self.Cz * np.random.exponential(scale=self.kT)
        C  = sqrt(2 * Cu / self.Ck) + self.Cs
        Gu = self.kT / self.Gz * np.random.exponential(scale=self.kT)
        G  = sqrt(2 * Gu / self.Gk) + self.Gs
        self.set_conv_and_head_from_segments(N, C, G)
        return self.head_loc
        
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
        
    def set_conv_and_head_from_segments(self, N, C, G):
        """Set the converter and head loc from passed segment values"""
        self.conv_loc = (N * cos(self.Ts),
                         N * sin(self.Ts))
        x = self.conv_loc[0] + G * cos(C + self.Ts - pi)
        y = self.conv_loc[1] + G * sin(C + self.Ts - pi)
        self.head_loc = (x, y)

## Begin the script that will produce the matrix of stored probabilities
print("This might take a while")
trials = 2000000
x_locs = np.arange(-5, 15, .2) 
y_locs = np.arange(-5, 15, .2)
probs = np.zeros((y_locs.size, x_locs.size))
hits = np.zeros((y_locs.size, x_locs.size))
# Instantiate the xb
xb = xNCG()
# Cycle through all iterations and collect the head locations
for i in range(trials):
    loc = xb.bop()
    x_ind = np.searchsorted(x_locs, loc[0]) - 1 #FIXME Check that there is not  
    y_ind = np.searchsorted(y_locs, loc[1]) - 1 #    an off by one error here
    hits[x_ind, y_ind] = hits[x_ind, y_ind] + 1
# Normalize the probabilities
min = np.min(hits)
max = np.max(hits)
hits = (hits - min)/(max-min)
contour.title = "Probability of an xNCG crossbridge being\n found at a given head locations"
contour.xlabel = "Location of XB head (nm)"
contour.ylabel = "Location of XB head (nm)"
contour.levels = [.1, .3, .5, .7, .9] 
contour.contour(x_locs, y_locs, hits)
