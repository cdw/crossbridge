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

## def xb_energy_func(xb):
    ## """If you give me a state of the XB segments,
    ## I'll give you a function that takes in their values
    ## and returns the total energy stored in them.
    ## """
    
    ## #Currently only worked on for TNCG case
    ## f = lambda T, N, C, G: (SPR_T*.5*(1/N)*pow(T-R_T,2) +
                            ## SPR_N*.5*pow(N-R_N,2) +
                            ## SPR_C*.5*(1/G)*pow(C-R_C,2) + 
                            ## SPR_G*.5*pow(G-R_G,2))
    ## return f

class xNCG():
    def __init__(self):
        """Create the values we'll be referencing for the XB"""
        self.Ts = pi/4 #rest angle (rad) of the connection to the thick fil
        self.Ns = 6    # rest length of neck region
        self.Nk = 5    # spring constant of neck region
        self.Nv = (5, 5, 5) # normal and rigor values of Ns
        self.Cs = pi/3+(pi-self.Ts) #rest angle (rad) of the converter domain
        self.Ck = 200  # torsional spring const of converter domain
        self.Cv = (pi/3, pi/3, 1.2*pi/3) # normal and rigor values of Cs
        self.Gs = 3    # rest length of globular domain
        self.Gk = 5    # spring constant of globular domain
        self.Gv = (5, 5, 5) # normal and rigor values of Gs
        self.Fm = 0    # mean of forces exerted on myosin heads
        self.Fv = 10.0 # variance of forces exerted on myosin heads
        self.Bd = 0.55 # dist at which binding becomes likely
        # Current state and identity of XB
        self.rest_head_loc() # Sets head_loc to rest location
        self.bound = False
        self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly
    
    def energy(self):
        pass
        
    def rest_head_loc(self):
        """Set the head loc to its rest location"""
        self.head_loc = (self.Gs * cos(self.Cs),
                         self.Gs * sin(self.Cs))
    