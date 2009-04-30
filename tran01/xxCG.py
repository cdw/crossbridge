## 
# Dave Williams, 20090315
# This script executes to produce a large matrix of the xxCG model's
# probabilities of diffusing to a given location
# 
##

# FIXME: Find correct scaling factor for kT to be used with our pN forces and nM scales

from numpy import array, pi, sin, cos, tan, arctan2, sqrt, hypot
import numpy as np
import cPickle as pkl
import contour

class xxCG:
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
        self.Csig = sqrt(self.kT / self.Ck)
        self.Gsig = sqrt(self.kT / self.Gk)
    
    def tran01(self, bSite):
        """Given an (x,y) location of an open binding site, bind or not after
        bopping the cross-bridge head to a new location. Return a boolean, 
        True for a binding event and False for no binding event.
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
        where xs is the rest length of the spring. At x=xs the spring sees no 
        strain and U(x)=0. If we are looking for the probability that the 
        particle or myosin head on the end of the spring will be in a given 
        location at a given time, we can fall back onto Boltzmann's law which 
        tells us that:
            p(x)/p(xs) = exp(-U(x)/kT)/exp(-U(xs)/kT)
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
            \int_{-\inf}^\inf 1/Z \exp(-.5 k x^2 /(kT)) dx
        This evaluates to 
            Z = \sqrt{2 pi kT / k}
        Meaning we are presented with a probability density function of
            P(x) = \sqrt{k / (2 pi kT)} \exp(-(k x^2)/(2 kT))
        This looks very similar to the PDF of a normal distribution, 
        typically written as:
            NormDis(x) = 1/sqrt(2 pi sig^2) exp(-(x-mu)^2/(2 sig^2))
        Where sig is the standard deviation of the distribution and mu is its
        mean. This means that we can interpert our desired PDF as that of a 
        normal distribution having:
            mu = xs
            sig = sqrt(kT / k)
            
        Details of implementation
        -------------------------
        This makes use of numpy's normal distribution:
            numpy.random.normal(loc=0.0, scale=1.0, size=None)
        Where loc is the mean is mu and scale is the standard deviation is sig 
        such that the probability density function of the resulting 
        distribution is: 
            f(x) = 1/sqrt(2 pi sig^2) exp(-(x-mu)^2/(2 sig^2))
        If we feed numpy's exponential distribution mu=xs and sig=sqrt(kT/k), 
        we have a distribution of: 
            f(x) = 1/sqrt(s pi kT / k) exp(-(x-xs)^2/(2 sqrt(kT/k)^2))
        Or, simplifiying:
            f(x) = sqrt(k / (2 pi kT)) exp(-k (x-xs)^2/(2 kT))
        The PDF that we derived above (with an explicit xs term this time).
        We can customize this for each spring with which we deal with by 
        plugging in their own means and spring constants.
        """
        C = np.random.normal(self.Cs, self.Csig)
        G = np.random.normal(self.Gs, self.Gsig)
        self.set_conv_and_head_from_segments(C, G)
        return self.head_loc
        
    def probability(self):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location
        """
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
# Instantiate the xb
xb = xxCG()

# Bebop that xb head
trials = 5000000
x = np.zeros(trials)
y = np.zeros(trials)
for i in range(trials):
    x[i], y[i] = xb.bop()
# Save it
dataout = open('bop_xxCG.pkl', 'wb')
pkl.dump([xb, trials, x, y], dataout, 2)
dataout.close()
# Plot it all pretty like
contour.title = "Location histogram of a diffusing \n xxCG crossbridge head"
contour.xlabel = "X loc of head (nm)"
contour.ylabel = "Y loc of head (nm)"
contour.hexbin(x, y)

# Find transition rates
trials = 5000 #per (x,y) location
x_locs = np.arange(-2, 15, .2)
y_locs = np.arange(-2, 15, .2)
rates = np.zeros((y_locs.size, x_locs.size))
pT = time.time()
cT = time.time()
# Cycle through head locations, collecting trans rates
for n,x in enumerate(x_locs):
    for m,y in enumerate(y_locs):
        for i in range(trials):
            rates[m, n] += xb.tran01((x, y))
    # Tell me how much time is left, about
    cT = time.time()
    rT = (cT-pT)*(x_locs.size - (n+1))
    print('On col %(c)04d of %(t)04d, about %(m)02d:%(s)02d left' \
          %{'c':n, 't':x_locs.size, 'm':rT//60, 's':rT%60})
    pT = cT
            
# Save it
dataout = open('trans01_xxCG.pkl', 'wb')
pkl.dump([xb, trials, x_locs, y_locs, rates], dataout, 2)
dataout.close()
# Plot the output
contour.title = ("0->1 transition rate of an xxCG crossbridge as a \n" \
                "function of head location")
contour.xlabel = "Location of XB head (nm)"
contour.ylabel = "Location of XB head (nm)"
contour.levels = np.array([.2, .4, .6, .8, .9])*rates.max()
contour.contour(x_locs, y_locs, rates)
# Display all the graphs produced in this script
contour.show()

# Normalize the hit likelihood
#min = np.min(hits)
#max = np.max(hits)
#hits = (hits - min)/(max-min)
#contour.title = "Probability of an xxCG crossbridge being\n found at a given #head locations"
#contour.xlabel = "Location of XB head (nm)"
#contour.ylabel = "Location of XB head (nm)"
#contour.levels = [.1, .3, .5, .7, .9] 
#contour.contour(x_locs, y_locs, hits)
#print(hits)