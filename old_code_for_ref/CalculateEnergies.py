#TODO
# * implement a TNCG class solver           Wednesday
# * test map generation with TNCG class     Wednesday

## HMM Based Spring Energy Profile
# Given a torsional and linear spring constant and rest angle, this spits back a
# matrix of energies at various locations

## Input parameters
# * k_th    - torsional spring constant
# * k_lin   - linear spring constant
# * s_th    - rest angle of torsional spring
# * s_lin   - rest length of linear spring (arm)
# * s_head  - length of head region
# * s_phi   - angle of arm
# * res     - level to resolve energy matrix
# * y_range - 2-vector of maximum and minimum y values for energy matrix
# * x_range - 2-vector of maximum and minimum x values for energy matrix

## Calculate energies

#Crossbridge classes
# Crossbridge classes are determined by a code of which points are springs
# The points are TNCG, so the classic class is xxCG, and a four spring
# class would be TNCG. The classic xxCG class is the default.

##"""A general representation of the unbound myosin crossbridge
##      H   - Head of the myosin
##      |  
##      G   - Globular domain
##      |  
##      C   - Converter region
##     /   
##    N     - Neck region
##   /     
##==T====== - Thick filament attachment
##
##A general representation is needed
##to let a call for comparision be heeded
##we pass instructions
##to define XB construction
##so that bad XB classes can be weeded
##"""

## Change log and misc musings
## 2009.01.08 CDW - I had this written as a class. That didn't work out so hot.
##  Now it is a flat file which allows me to not have to append self. to every
##  dang variable that I want to use.
## 2009.01.13 CDW - But, past Dave, what you didn't realize is that then it
##  is hard to change variables from within functions (hard to reclass an XB)
##  so to mitigate this we take on the xb as a class but have the functions
##  that interact with it remain seperate.


# A little bit of startup organizaiton to provide support for various things
# that I want to treat as built in (putting the batteries in the toy)
from math import pi, sin, cos, tan, atan2, sqrt, hypot
import scipy.optimize.optimize as opt


class xb:
    ## Set the crossbridge values to an initial state (for class xxCG)
    # Define default rest values
    R_T = pi/4 #rest angle (rad) of connection to thick fil
    R_N = 8.5  #rest length (nm) of thick-fil-to-converter rod
    R_C = pi/3+(pi-R_T) #rest angle (rad) of converter domain
    R_G = 6    #rest length (nm) of globular head domain
    # Define the current values of the springs as the rest values
    T = R_T
    N = R_N
    C = R_C
    G = R_G
    # Define default crossbridge as class one
    SPR_T = False
    SPR_N = False
    SPR_C = True
    SPR_G = True
    # Define default spring values
    K_T = -1
    K_N = -1
    K_C = 200 #in pN/(nm*rad)
    K_G = 5   #in pN/nm
    # Define head displacement penalty
    K_H = 10 * reduce(max, [K_T, K_N, K_C, K_G])

    def __repr__(self):
        return ('XB Class: ' + self.give_me_curr_xb_class() +
                '\n Rest head loc: ' + str(self.head_loc(True)) +
                '\n Curr head loc: ' + str(self.head_loc(False)))

    def give_me_curr_xb_class(self):
        """This fellow gives back a little string formated in the
        xxCG format to indicate the current type of xb we have
        """
        x_or_let = lambda boo,abc: boo*abc+(not boo)*"x"
        return (x_or_let(self.SPR_T,"T")+x_or_let(self.SPR_N,"N")+
                x_or_let(self.SPR_C,"C")+x_or_let(self.SPR_G,"G"))
    
    def find_locs(self):
        """I take an XB structure and find the locations of the converter
        and the globular domain, I pass them back as a tuple
        """
        XB = self
        c_loc = (XB.N*cos(XB.T), XB.N*sin(XB.T))
        h_loc = (c_loc[0]+XB.G*cos(XB.C-pi+XB.T),
                 c_loc[1]+XB.G*sin(XB.C-pi+XB.T))
        return (c_loc, h_loc)

    def give_me_badness(self, head):
        """If you give me the desired head location, I will give you
        how far short of it I fall, the badness
        """
        [c_loc, h_loc] = self.find_locs() #suboptimal, should be stored, meh
        badness = pow((pow(head[0]-h_loc[0],2)+pow(head[1]-h_loc[1],2)),.5)
        return badness

    def set_xb_class(self, xb_class='xxCG',
                     rest_vals = (pi/4, 8.5, pi/3, 6),
                     spr_vals=(-1, -1, 200, 5)):
        """Give me a class, some rest values, and spring values
        and I will change the class of crossbridge currently represented
        """
        #Parse the xb_class variable
        let_to_bool = lambda let: xb_class.count(let)==1
        self.SPR_T = let_to_bool('T')
        self.SPR_N = let_to_bool('N')
        self.SPR_C = let_to_bool('C')
        self.SPR_G = let_to_bool('G')
        #Parse rest values
        self.R_T = rest_vals[0] #rest angle (rad) of connection to thick fil
        self.R_N = rest_vals[1] #rest length (nm) of thick-fil-to-converter rod
        self.R_C = rest_vals[2] #rest angle (rad) of converter domain
        self.R_G = rest_vals[3] #rest length (nm) of globular head domain
        #Parse spring values
        self.K_T = spr_vals[0] #in pN/(nm*rad)
        self.K_N = spr_vals[1] #in pN/nm
        self.K_C = spr_vals[2] #in pN/(nm*rad)
        self.K_G = spr_vals[3] #in pN/nm

    def head_loc(self, use_rest_values=False):
        """If you give me the state of the XB component values,
        I'll give you the location of the myosin head.
        """
        #Unpack component values (or their rest positions)
        if use_rest_values == False:
            [T,N,C,G]=[self.T, self.N, self.C, self.G]
        else:
            [T,N,C,G]=[self.R_T, self.R_N, self.R_C, self.R_G]
        #Calculate the x and y location of the head
        #x_T = 0 #this is assumed
        #y_T = 0 #this is assumed
        x_C = cos(T)*N
        y_C = sin(T)*N
        x_H = x_C + cos(C)*G
        y_H = y_C + sin(C)*G
        #Give back a tuple
        return (x_H, y_H)

    def energy_function(self, H=None):
        """Returns a function for the energy in the XB at a given head loc"""
        # If not given an H, go get the current H
        if H is None:
            H = self.head_loc()
        f = lambda T, N, C, G: (self.SPR_T * .5 * (1/N) * pow(T-self.R_T, 2) +
                                self.SPR_N * .5 * pow(N-self.R_N, 2) +
                                self.SPR_C * .5 * (1/G) * pow(C-self.R_C, 2) + 
                                self.SPR_G * .5 * pow(G-self.R_G, 2))
 


def give_me_rest_xb_vals(xb):
    return (xb.R_T, xb.R_N, xb.R_C, xb.R_G)

def head_loc_to_xb_vals(xb, H=None):
    """If you give me the location of H,
    I will give you the lengths and angles of the xb
    In the form of an xb class
    """
    # If not given an H, go get the current H
    if H is None:
        H = xb.head_loc()
    if xb.SPR_N == True and xb.SPR_T == True and xb.SPR_C == True and xb.SPR_G == True:
        #in the TNCG case, do something
        xb_imob = xb
        xb_imob.SPR_T = False
        xb_imob.SPR_N = False
        xb_imob = head_loc_to_xb_vals(xb_imob, H)
        #poop poop poop, the stratagy I was going for didn't allow for boundry conditions
        xb_energy_func(xb)
    elif xb.SPR_C == True and xb.SPR_G == True:
        #Where C is fixed in location, and the system looks like the
        #system from CDW's rotation, solve the old way
        #Unpack some values
        [R_T, R_N] = [xb.R_T, xb.R_N]
        x_C = cos(R_T)*R_N
        y_C = sin(R_T)*R_N
        run_G = H[0]-x_C
        rise_G = H[1]-y_C
        xb.C = atan2(rise_G, run_G)
        xb.G = hypot(run_G, rise_G)
        return xb

def xb_energy_func(xb):
    """If you give me a state of the XB segments,
    I'll give you a function that takes in their values
    and returns the total energy stored in them.
    """
    
    #Currently only worked on for TNCG case
    f = lambda T, N, C, G: (SPR_T*.5*(1/N)*pow(T-R_T,2) +
                            SPR_N*.5*pow(N-R_N,2) +
                            SPR_C*.5*(1/G)*pow(C-R_C,2) + 
                            SPR_G*.5*pow(G-R_G,2))
    return f
        

def xb_energy(xb):
    """If you give me a state of the XB segments,
    I'll give you an energy of the whole XB.
    """
    #Unpack component values
    [T,N,C,G]=[xb.T, xb.N, xb.C, xb.G]
    [R_T,R_N,R_C,R_G]=[xb.R_T, xb.R_N, xb.R_C, xb.R_G]
    [SPR_T, SPR_N, SPR_C, SPR_G] = [xb.SPR_T, xb.SPR_N, xb.SPR_C, xb.SPR_G]
    U_T = SPR_T*.5*(1/N)*pow(T-R_T,2) #energy in T 
    U_N = SPR_N*.5*pow(N-R_N,2)       #energy in N
    U_C = SPR_C*.5*(1/G)*pow(C-R_C,2) #energy in C
    U_G = SPR_G*.5*pow(G-R_G,2)       #energy in G
    return U_T+U_N+U_C+U_G

def xb_forces(xb_vals=None):
    """If you give me a state of the XB segments,
    I'll give you an array of forces on each node?
    """
    pass
