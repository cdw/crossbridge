## 
# Dave Williams, 20090315
# This script executes to produce a large matrix of the xxCG model's
# potential energies and transition rates
# 
##

from numpy import array, pi, sin, cos, tan, arctan2, sqrt, hypot
import numpy as np
import matplotlib as matp
import matplotlib.pyplot as plt

class xxCG():
    def __init__(self):
        """Create the values we'll be referencing for the XB"""
        self.Cs = pi/3 # rest angle of converter domain
        self.Ck = 200  # torsional spring const of converter domain
        self.Cv = (pi/3, pi/3, 1.2*pi/3) # normal and rigor values of Cs
        self.Gs = 10.5 # rest length of globular domain
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
        """Return the energy stored in the XB, 
        given the current head_loc
        """
        G = self.glob_len()
        C = self.conv_ang()
        return ((1/2) * self.Gk * (G-self.Gs)**2 + 
        1/(2*G) * self.Ck * (C-self.Cs)**2)

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
    
def contour(x,y,z):
    """Take 1-d x and y arrays, and a 2-d z array and plot out a contour"""
    # Set something?
    matp.rcParams['xtick.direction'] = 'out'
    matp.rcParams['ytick.direction'] = 'out'
    # Produce the coordinate matrices from two coordinate vectors
    X, Y = np.meshgrid(x, y)
    # Plot the contours with default labeling
    plt.figure()
    CS = plt.contour(X, Y, z)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Simplest default with labels')
    # Display our contour
    plt.show()
    

## Begin the script that will produce the matrix of stored energies
x_locs = np.arange(3, 8, .1) 
y_locs = np.arange(7, 11, .1)
energs = np.zeros((y_locs.size, x_locs.size))
# Instantiate the xb
xb = xxCG()
# Cycle through and collect all the energies
n = [0,0]
for y in y_locs:
    for x in x_locs:
        xb.head_loc = (x,y)
        energs[n[0], n[1]] = xb.energy()
        n[1] = n[1] + 1
    n[0] = n[0] + 1
    n[1] = 0
