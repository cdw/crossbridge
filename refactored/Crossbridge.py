#!/usr/bin/env python
# encoding: utf-8
#
## This file defines the system for the various crossbridge types. 
## These crossbridges have four points that may be represented as 
## linear or torsional springs. A schematic follows.
##       H   - Head of the myosin
##       |  
##       G   - Globular domain
##       |  
##       C   - Converter region
##      /   
##     N     - Neck region
##    /     
## ==T====== - Thick filament


from numpy import pi, sin, cos, arctan2, sqrt, hypot
import numpy as np
from scipy.optimize import fmin_powell as fmin
import time
import cPickle as pkl


class XB:
    """Provide a base class attempting to be agnostic to XB type."""
    def __init__(self):
        self.type = self.__class__.__name__
        """Create the values we'll be referencing for the XB"""
        self.Ts = pi/4 # angle (rad) of the connection to the thick fil
        self.Tk = 100  # spring constant of connection to the thick fil
        self.Ns = 7    # rest length of neck region
        self.Nr = 5    # rigor length of neck region
        self.Nk = 10    # spring constant of neck region
        self.Cs = pi/3+(pi-self.Ts) #rest angle (rad) of the converter domain
        self.Cr = self.Cs + pi/6 #rest angle (rad) in rigor state
        self.Ck = 100  # torsional spring const of converter domain
        self.Gs = 3    # rest length of globular domain
        self.Gr = 3    # rest length of globular domain in rigor state
        self.Gk = 5    # spring constant of globular domain
    
    def calc_attributes(self):
        """Set and maintain a set of an XB's related values.
            
        This is included here as opposed to in the initialization so that
        different xb types can modify the spring values before setting these.
        """
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
        # Free energy related values
        G_0 = 13 #in RT 
        #FIXME: Put above line in Joules or convert other units
        ATP_conc = 0.005 # or 5 mM
        ADP_conc = 0.00003 # or 30 uM
        Pi_conc  = 0.003 # or 3 mM
        self.G_lib = - G_0 - np.log(ATP_conc / (ADP_conc * Pi_conc))
        self.alpha = 0.28 #G_lib freed in 0->1 trans, from Bert/Tom/Pate/Cooke
        self.eta = 0.68 #ditto, for 1->2 trans
        # Current state and identity of XB
        self.rest_conv_and_head_loc() # Set conv_loc and head_loc to rest locs
        self.bound = False
        self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly
    
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
    
    def free_energy(self, state, h_loc=None):
        """What's the free energy in a given state at the specified values?
            
        Note that eta and alpha are the fractions of energy liberted from ATP 
        during this whole sordid affair.
        """
        if state is 0:
            return 0
        elif state is 1:
            return self.alpha * self.G_lib + self.minimize_energy(h_loc)
        elif state is 2:
            return self.eta * self.G_lib + self.minimize_energy(h_loc)
    
    def force(self, h_loc, state=1):
        """From the head loc, get the for vector being exerted by the XB"""
        self.minimize_energy(h_loc, state)
        C = self.conv_ang()
        G = self.glob_len()
        Ck = self.Ck
        Gk = self.Gk
        if state is 2:
            Cs = self.Cr
            Gs = self.Gr
        else:
            Cs = self.Cs
            Gs = self.Gs
        Fx = -Gk * (G - Gs) * cos(C) + 1/G * Ck * (C - Cs) * sin(C)
        Fy = -Gk * (G - Gs) * sin(C) + 1/G * Ck * (C - Cs) * cos(C)
        return (Fx, Fy)
    


class TNCG(XB):
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
    def __init__(self):
        """No modification of values needed, just trigger the att calc"""
        XB.__init__(self)
        self.calc_attributes()
    
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
        T = np.random.normal(self.Ts, self.Tsig)
        N = np.random.normal(self.Ns, self.Nsig)
        C = np.random.normal(self.Cs, self.Csig)
        G = np.random.normal(self.Gs, self.Gsig)
        self.set_conv_and_head_from_segments(T, N, C, G)
        return self.head_loc
    
    def probability(self, h_loc=None):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location
        """
        if h_loc is not None: 
            self.head_loc = h_loc
        U = self.minimize_energy() #sets C to lowest U loc, returns U
        T = self.thic_ang()
        N = self.neck_len()
        C = self.conv_ang()
        G = self.glob_len()
        pT = 1/self.Tz * np.exp(-U / self.kT)
        pN = 1/self.Nz * np.exp(-U / self.kT)
        pC = 1/self.Cz * np.exp(-U / self.kT)
        pG = 1/self.Gz * np.exp(-U / self.kT)
        return pT * pN * pC * pG
    
    def minimize_energy(self, h_loc=None, state=1):
        """Set the conv_loc to minimize the XB's energy 
        and return the newly located minimum energy
        """
        if state is 2:
            Nstore = self.Ns
            Cstore = self.Cs
            Gstore = self.Gs
            self.Ns = self.Nr
            self.Cs = self.Cr
            self.Gs = self.Gr
        if h_loc is not None: 
            self.head_loc = h_loc
        e = lambda (l): self.e_dep_conv(l)
        min_n_len = fmin(e, self.conv_loc, disp=0)
        if state is 2:
            self.Ns = Nstore
            self.Cs = Cstore
            self.Gs = Gstore
        return self.energy()
    
    def e_dep_conv(self, conv):
        """Update the conv_loc and return the xb energy"""
        self.conv_loc = (conv[0], conv[1])
        return self.energy()
    
    def set_state(self, state):
        """Set the state of the XB"""
        if state is 0:
            self.state = 0
            self.Cs = pi/3+(pi-self.Ts)
        elif state is 1:
            self.state = 1
            self.Cs = pi/3+(pi-self.Ts)
        elif state is 2:
            self.state = 2
            self.Cs = pi/2+(pi-self.Ts)
    


class TNCx(XB):
    pass



class xxCG(XB):
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
    def __init__(self):
        """Modify values for this spring system, trigger the attribute calc"""
        XB.__init__(self)
        self.Ns = 1    # rest length of neck region
        self.Gs = 10.5 # rest length of globular domain
        self.calc_attributes()
    
    def tran12(self):
        """Cautiously working on this fellow to get rate of strong binding"""
        print("You haven't written tran12 yet")
        # Bop again
        # FIXME : Talk to Tom about this, I am just concerned as we are bopping from a point that already has some energy stored in it. Should I be using the random normal blah blah or some sort of assymetrical distribution?
        #hLoc = self.bop()
        C = np.random.normal(self.Cs, self.Csig)
        G = np.random.normal(self.Gs, self.Gsig)
        # Calculate energy difference
        EnergyDiff = (self.free_energy(2, self.Ts, self.Ns, C, G) -
            self.free_energy(1, self.Ts, self.Ns, C, G))
        # Transition if energy at new position is less
        if EnergyDiff <= 0:
            self.state = 2
        pass
    
    def tran20(self):
        """Cautiously working on unbinding rates"""
        
    
    def bop(self):
        """Bop the xb to a new location, return the new head location.
            
        See docstring for TNCG.bop() for justification and implementation.
        """
        C = np.random.normal(self.Cs, self.Csig)
        G = np.random.normal(self.Gs, self.Gsig)
        self.set_conv_and_head_from_segments(self.Ts, self.Ns, C, G)
        return self.head_loc
    
    def probability(self, h_loc=None):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location
        """
        if h_loc is not None: 
            self.head_loc = h_loc
        G = self.glob_len()
        C = self.conv_ang()
        U = self.energy()
        pG = 1/self.Gz * np.exp(-U / self.kT)
        pC = 1/self.Cz * np.exp(-U / self.kT)
        return pG * pC
    
    def minimize_energy(self, h_loc=None, state=1):
        """Set the conv_loc to minimize the XB's energy 
        and return the newly located minimum energy
        """
        if h_loc is not None: 
            self.head_loc = h_loc
        # Nothing more is needed as the xxCG XB's converter domain is fixed
        return self.energy()
    
    def set_state(self, state):
        """Set the state of the XB"""
        if state is 0:
            self.state = 0
            self.Cs = pi/3+(pi-self.Ts)
        elif state is 1:
            self.state = 1
            self.Cs = pi/3+(pi-self.Ts)
        elif state is 2:
            self.state = 2
            self.Cs = pi/2+(pi-self.Ts)
        else:
            print("Invalid XB state, must be 0, 1, or 2")
            return
    



class xNxx(XB):
    """An instance of the classic one-spring crossbridge.
        
           H   - Head of the myosin
           |  
           G   - Globular domain,   fixed length
           |  
           C   - Converter region,  fixed angle
          /   
         N     - Neck region,       linear spring
        /     
    ===T====== - Thick filament,    fixed angle
    """
    def __init__(self):
        """Modify values as needed for diff geom, then trigger the att calc"""
        XB.__init__(self)
        self.Ts = 0 # angle (rad) of the connection to the thick fil
        self.Ns = 5    # rest length of neck region
        self.Nr = 3    # rigor length of neck region
        self.Nk = 5    # spring constant of neck region
        self.Cs = pi/2+(pi-self.Ts) #rest angle (rad) of the converter domain
        self.Cr = self.Cs
        self.Gs = 10    # rest length of globular domain
        self.calc_attributes()
    
    def bop(self):
        """Bop the xb to a new location, return the new head location.
            
        See docstring for TNCG.bop() for justification and implementation.
        """
        N = np.random.normal(self.Ns, self.Nsig)
        self.set_conv_and_head_from_segments(self.Ts, N, self.Cs, self.Gs)
        return self.head_loc
    
    def probability(self, h_loc=None):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location
        """
        if h_loc is not None: 
            self.head_loc = h_loc
        N = self.neck_len()
        U = self.energy()
        pN = 1/self.Nz * np.exp(-U / self.kT)
        return pN
    
    def minimize_energy(self, h_loc=None):
        """Set the conv_loc to minimize the XB's energy 
        and return the newly located minimum energy
        """
        if h_loc is not None: 
            self.head_loc = h_loc
        # NOTE Head Y location is ignored, consistant with G being inflexible
        self.set_conv_and_head_from_segments(self.Ts, self.head_loc[0], 
                                             self.Cs, self.Gs)
        return self.energy()
    
        def set_state(self, state):
            """Set the state of the XB"""
            if state is 0:
                self.state = 0
                self.Ns = 5
            elif state is 1:
                self.state = 1
                self.Ns = 5
            elif state is 2:
                self.state = 2
                self.Ns = 0
    



class xNCx(XB):
    """A counterfactual instance of an alternate two-spring crossbridge.
        
           H   - Head of the myosin
           |  
           G   - Globular domain,   fixed length
           |  
           C   - Converter region,  torsional spring
          /   
         N     - Neck region,       linear spring
        /     
    ===T====== - Thick filament,    fixed angle
    """
    ###IMPORTANT: MINEFIELD, DON'T USE TO GENERATE FIGURES, GERDUNKEN ONLY###
    # NOTE This is bogus as of 4-29-09 cdw
    # Spot the problem... give up? Try to get farther than Gs away from the 
    # line that Ts allows the converter location to trave in. That's right, 
    # without G or T being a spring is it not possible.
    # Moral of the story: you need at least one freely moving linear spring
    #                     downstream of one torsional spring in order to have
    #                     freedom of movement that gets you to any (X,Y) loc
    def __init__(self):
        """No modification of values needed, just trigger the att calc"""
        XB.__init__(self)
        self.calc_attributes()
    
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
        N = np.random.normal(self.Ns, self.Nsig)
        C = np.random.normal(self.Cs, self.Csig)
        self.set_conv_and_head_from_segments(self.Ts, N, C, self.Gs)
        return self.head_loc
    
    def probability(self, h_loc=None):
        """Given the location of the XB head,
        return the probability that it is there,
        relative to the probability of being in the rest location
        """
        if h_loc is not None: 
            self.head_loc = h_loc
        self.minimize_energy()
        N = self.neck_len()
        C = self.conv_ang()
        U = self.energy()
        pN = 1/self.Nz * np.exp(-U / self.kT)
        pC = 1/self.Cz * np.exp(-U / self.kT)
        return pG * pC
    
    def minimize_energy(self, h_loc = None):
        """Set the conv_loc to minimize and return the XB's energy 
            
        This is based on taking the equation for energy in the xNCx Xb:
            U(N,C)= Nk (N - Ns) + Ck (C - Cs)
        Writeing C in terms of Ts and the x and y components of G gives:
            C = (pi/2 - T) + pi/2 + arctan[Gy/Gx]
        Where the x and y components of G can be written, using the x and y
        locations of the XB head, as:
            Gx = Hx - N cos[T]
            Gy = Hy - N sin[T]
        Plugging all this back into U we can write U(N) thusly:
            U(N) = Nk (N - Ns) + Ck (C(N) - Cs)
            U(N) = Nk (N - Ns) + Ck ((pi - T + 
                arctan[(Hy - N sin[T]) / (Hx - N cos[T])]) - Cs)
        Which we then take the derivitive of wrt N
            U'(N) = (Nk (Hx^2 + Hy^2 + N^2) + 
                (Hy Ck - 2 Hx Nk N) Cos[T] - (Hx Ck + 2 Hy Nk N) Sin[T]) 
                / (Hx^2 + Hy^2 + N^2 - 2 Hx N Cos[T] - 2 Hy N Sin[T])
        Setting this to zero and solving for N we get:
            N = Hx Cos[T] + Hy Sin[T] +- 
                1/Nk Sqrt[-Nk (Hy Cos[T] - Hx Sin[T]) 
                (Ck + Hy Nk Cos[T] - Hx Nk Sin[T])]
        Which is what we set N to in order to get the minimum energy
        """
        if h_loc is None: 
            h_loc = self.head_loc
        Hx = h_loc[0]
        Hy = h_loc[1]
        T = self.Ts
        N = (Hx * cos(T) + Hy * sin(T) + 1/self.Nk * 
            sqrt(-self.Nk * (Hy * cos(T) - Hx * sin(T)) * 
                 (self.Ck + Hy * self.Nk * cos(T) - Hx * self.Nk * sin(T))))
        # Update conv_loc for ease of use, this temp makes head loc inaccurate
        self.set_conv_and_head_from_segments(self.Ts, N, self.Cs, self.Gs)
        # Find the proper angle for C, then fix the head loc
        C = pi - T + arctan2(Hy - self.conv_loc[1], Hx - self.conv_loc[0])
        self.set_conv_and_head_from_segments(self.Ts, N, C, self.Gs)
        return self.energy()
    





