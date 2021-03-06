#!/usr/bin/env python
# encoding: utf-8
""" Crossbridge.py
Created by Dave Williams on 2009-07-01.

This file defines the system for the various crossbridge types. 
These crossbridges have four points that may be represented as 
linear or torsional springs. A schematic follows::
           H   - Head of the myosin
           |  
           G   - Globular domain
           |  
           C   - Converter region
          /   
         N     - Neck region
        /     
    ===T====== - Thick filament
"""

import warnings
import numpy.random as random
from scipy.optimize import fmin_bfgs as fmin
from numpy import pi, sqrt, log
import math as m


class Spring(object):
    """A generic spring that handles some accounting"""
    def __init__(self, spring_config):
        """Create the spring with a set of passed values
        
        Takes:
            spring_config: a dictionary containing the keys
                weak: the spring rest value when unbound or weakly bound
                strong: the spring rest value when strongly bound
                spring_konstant: the spring constant
        Returns: 
            None
        """
        self.weak = spring_config['weak']
        self.strong = spring_config['strong']
        self.k = spring_config['spring_konstant']
        # Diffusion related
        temperature = 288 # in K
        boltzman = 1.381 * 10**-23 # Boltzman const (in J/K)
        k_t = boltzman * temperature * 10**21  # kT without pN/nM conversion
        # Normalize is a factor used to normalize the PDF of the segment vals
        self.normalize = sqrt(2*pi*k_t/self.k)
        self.stand_dev = sqrt(k_t/self.k) # of seg vals
    
    def rest(self, state):
        """Return the rest value of the spring at in a given state
        
        Takes:
            state: the spring state of interest [1, 2, or 3]
        Returns:
            None
        """
        if state in [1, 2]:
            return self.weak
        elif state == 3:
            return self.strong
        else:
            warnings.warn("Improper value for spring state")
    
    def energy(self, curr_val, state):
        """Given a current value, return the energy the spring stores"""
        if state in [1, 2]:
            return (0.5 * self.k * m.pow((curr_val-self.weak), 2))
        elif state == 3:
            return (0.5 * self.k * m.pow((curr_val-self.strong), 2))
        else:
            warnings.warn("Improper value for spring state")
    
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
        return (random.normal(self.weak, self.stand_dev))
    

class Crossbridge(object):
    def __init__(self, config = None):
        """A generic cross-bridge, a TNCG one by default
        
        This cross-bridge class provides default functionality, such
        as calculation of energies and transition rate constants, that
        later classes of cross-bridge may build upon.
        Takes:
            config: a dictionary containing the spring definitions, a sample:
                {'T': {'weak': 40*pi/180,
                       'strong': 40*pi/180,       
                       'spring_konstant': 100
                      }, 
                 'N': ...
                 'C': ...
                 'G': ... 
                }
        """
        # Eventually, take out default config, put in GenerateData
        if config == None:
            config = {
                'T': {
                    'weak': 40*pi/180,
                    'strong': 40*pi/180,
                    'spring_konstant': 100
                },
                'N': {
                    'weak':  10.5,
                    'strong':  10.5,
                    'spring_konstant': 10
                },
                'C': {
                    'weak':2*pi - 165*pi/180,
                    'strong': 2*pi - 110*pi/180,
                    'spring_konstant': 40
                },
                'G': {
                    'weak': 9.6,
                    'strong': 9.6,
                    'spring_konstant': 5
                }
            }
        self.config = config
        # Define initial spring values
        self.t = Spring(self.config['T'])
        self.n = Spring(self.config['N'])
        self.c = Spring(self.config['C'])
        self.g = Spring(self.config['G'])
    
    def minimize_energy(self, h_loc, state):
        """The cross-bridge's minimum energy for a given head location, state
        
        Takes:
            h_loc:  [x,y] location of the cross-bridge tip
            state: state of the cross-bridge, 1, 2, or 3
        Returns:
            energy: cross-bridge's minimum energy 
            min_conv: [x,y] converter location yielding the minimum energy 
        """
        rest_conv_loc = (self.n.rest(state) * m.cos(self.t.rest(state)),
                         self.n.rest(state) * m.sin(self.t.rest(state)))
        min_conv = fmin(self.energy, 
            rest_conv_loc, args = (h_loc, state), disp=0)
        return (self.energy(min_conv, h_loc, state), min_conv)
    
    def energy(self, conv_loc, h_loc, state):
        """Return the energy in the xb with the given parameters"""
        (t_ang, n_len, c_ang, g_len) = self.seg_values(conv_loc, h_loc)
        return float(
            self.t.energy(t_ang, state) + 
            self.n.energy(n_len, state) + 
            self.c.energy(c_ang, state) + 
            self.g.energy(g_len, state) 
        )
    
    def free_energy(self, h_loc, state):
        """Return the free energy in the xb with the given parameters"""
        g_0 = 13 #in RT 
        atp_conc = 0.005 # or 5 mM
        adp_conc = 0.00003 # or 30 uM
        phos_conc  = 0.003 # or 3 mM
        g_lib = - g_0 - log(atp_conc / (adp_conc * phos_conc))
        alph = 0.28 #G_lib freed in 0->1 trans, from Bert/Tom/Pate/Cooke
        eta = 0.68 #ditto, for 1->2 trans
        if state is 1:
            return float(0)
        elif state is 2:
            return float(alph * g_lib + self.minimize_energy(h_loc, state)[0])
        elif state is 3:
            return float(eta * g_lib + self.minimize_energy(h_loc, state)[0])
        
    def force(self, h_loc, state):
        """From the head loc, the force vector being exerted by the XB"""
        (energy, conv_loc) = self.minimize_energy(h_loc, state)
        (t_ang, n_len, c_ang, g_len) = self.seg_values(conv_loc, h_loc)
        del(energy, t_ang, n_len) # Not needed
        c_k = self.c.k
        g_k = self.g.k
        c_s = self.c.rest(state)
        g_s = self.g.rest(state)
        f_x = (-g_k * (g_len - g_s) * m.cos(c_ang) + 
                1/g_len * c_k * (c_ang - c_s) * m.sin(c_ang))
        f_y = (-g_k * (g_len - g_s) * m.sin(c_ang) + 
                1/g_len * c_k * (c_ang - c_s) * m.cos(c_ang))
        return [float(f_x), float(f_y)]
    
    def seg_values(self, conv_loc, h_loc):
        """Calculate the values of the segments of the XB"""
        diff = [h_loc[0] - conv_loc[0], h_loc[1] - conv_loc[1]]
        t_ang = m.atan2(conv_loc[1], conv_loc[0])
        n_len = m.hypot(conv_loc[0], conv_loc[1])
        c_ang = m.atan2(diff[1], diff[0]) + m.pi - t_ang
        g_len = m.hypot(diff[0], diff[1])
        return (t_ang, n_len, c_ang, g_len)
    
    def bind_or_not(self, b_site):
        """Given an (x,y) location of an open binding site, bind or not after
        bopping the cross-bridge head to a new location. Return a boolean,
        True for a binding event and False for no binding event.
        """
        ## Bop the springs to get new values
        t_ang = self.t.bop()
        n_len = self.n.bop()
        c_ang = self.c.bop()
        g_len = self.g.bop()
        ## Translate those values to (x,y) postitions
        conv_loc = (n_len * m.cos(t_ang),
                    n_len * m.sin(t_ang))
        h_loc = (conv_loc[0] + g_len * m.cos(c_ang + t_ang - m.pi), 
                 conv_loc[1] + g_len * m.sin(c_ang + t_ang - m.pi))
        ## Find the distance to the binding site
        distance = m.hypot(b_site[0]-h_loc[0], b_site[1]-h_loc[1])
        ## The binding prob is dept on the exp of a dist
        b_prob = 12*m.exp(-distance**2)
        ## Throw a random number to check binding
        return (b_prob > random.rand())
    
    def r12(self, b_site, trials):
        """Give the prob of binding, given a b_site and number of trials """
        # Binds gives us the number of times binding occurs out of all trials
        binds = sum(self.bind_or_not(b_site) for t in range(trials))
        return (float(binds)/float(trials))
    
    def r23(self, b_site):
        """Given a binding site, b_site, to which a myosin head is loosely
        bound, return a probability of transition to a tightly bound state
        """
        state2_energy = self.minimize_energy(b_site, 2)[0]
        state3_energy = self.minimize_energy(b_site, 3)[0]
        rate = .1 * (1 + m.tanh(.4 * (state2_energy - state3_energy)+4))+.001
        # Note that the .001 is just to keep rates above 0.0000 at all times
        return float(rate)
    
    def r31(self, b_site):
        """Given a binding site, b_site, to which a myosin head is tightly
        bound, return a probability of transition to an unbound state
        """
        state3_energy = self.minimize_energy(b_site, 3)[0]
        #rate = m.exp(-1 / (state3_energy + 1e-9)) #1e-9 avoids 1/0 at rest loc
        rate =  m.sqrt(.01 * state3_energy) + 0.02
        return float(rate)
    

class FourSpring(Crossbridge):
    """An instance of the four-spring crossbridge.
        
           H   - Head of the myosin
           |  
           G   - Globular domain,   linear spring
           |  
           C   - Converter region,  torsional spring
          /   
         N     - Neck region,       linear spring
        /     
    ===T====== - Thick filament,    torsional spring
    """
    def __init__(self, config = None):
        """No modification of values needed, just trigger the att calc"""
        Crossbridge.__init__(self, config)
    

class TwoSpring(Crossbridge):
    """An instance of the two-spring crossbridge.
        
           H   - Head of the myosin
           |  
           G   - Globular domain,   linear spring
           |  
           C   - Converter region,  torsional spring
          /   
         N     - Neck region,       fixed length
        /     
    ===T====== - Thick filament,    fixed angle
    """
    def __init__(self, config = None):
        """Modify values for this spring system, trigger the attribute calc"""
        Crossbridge.__init__(self, config)
    
    def minimize_energy(self, h_loc, state):
        """Return the min energy in the XB with the head at the given loc"""
        rest_conv_loc = (self.n.rest(state) * m.cos(self.t.rest(state)),
                         self.n.rest(state) * m.sin(self.t.rest(state)))
        return (self.energy(rest_conv_loc, h_loc, state), rest_conv_loc)
    
    def bind_or_not(self, b_site):
        """Given an (x,y) location of an open binding site, bind or not after
        bopping the cross-bridge head to a new location. Return a boolean,
        True for a binding event and False for no binding event.
        """
        ## Bop the springs to get new values
        t_ang = self.t.rest(1)
        n_len = self.n.rest(1)
        c_ang = self.c.bop()
        g_len = self.g.bop()
        ## Translate those values to (x,y) postitions
        conv_loc = (n_len * m.cos(t_ang),
                    n_len * m.sin(t_ang))
        h_loc = (conv_loc[0] + g_len * m.cos(c_ang + t_ang - m.pi), 
                 conv_loc[1] + g_len * m.sin(c_ang + t_ang - m.pi))
        ## Find the distance to the binding site
        distance = m.hypot(b_site[0]-h_loc[0], b_site[1]-h_loc[1])
        ## The binding prob is dept on the exp of a dist
        b_prob = 72*m.exp(-distance**2) +.00001
        ## Throw a random number to check binding
        return (b_prob > random.rand())
    

class OneSpring(Crossbridge):
    """An instance of the one-spring crossbridge"""
    def __init__(self, config = None):
        Crossbridge.__init__(self, config)
    
    def minimize_energy(self, h_loc, state):
        """Return the min energy of the XB at h_loc, ignore y dimension"""
        # Ignore y dim and only use energy in neck
        return (self.n.energy(h_loc[0], state), (h_loc[0], 0))
    
    def bind_or_not(self, b_site):
        """Given an (x,y) location of an open binding site, bind or not after
        bopping the cross-bridge head to a new location. Return a boolean,
        True for a binding event and False for no binding event.
        """
        ## Bop the spring to get a new value
        n_len = self.n.bop()
        ## Find the distance to the binding site
        distance = abs(b_site[0]-n_len) #Ignore y dim
        ## The binding prob is dept on the exp of a dist
        b_prob = m.exp(-distance)
        ## Throw a random number to check binding
        return (b_prob > random.rand())
    

