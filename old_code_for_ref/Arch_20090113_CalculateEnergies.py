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
##  that interact with it remain seperate. TODO TODO TODO

# SEE TODO RIGHT ABOVE THIS!



# A little bit of startup organizaiton to provide support for various things
# that I want to treat as built in (putting the batteries in the toy)
from math import pi, sin, cos, tan, atan2, sqrt, hypot
import matplotlib.pylab as plt

class xb:
    ## Set the crossbridge values to an initial state (for class xxCG)
    #Define default rest values
    R_T = pi/4 #rest angle (rad) of connection to thick fil
    R_N = 8.5  #rest length (nm) of thick-fil-to-converter rod
    R_C = pi/3 #rest angle (rad) of converter domain
    R_G = 6    #rest length (nm) of globular head domain
    #Define default crossbridge as class one
    SPR_T = False
    SPR_N = False
    SPR_C = True
    SPR_G = True
    #Define default spring values
    K_T = -1
    K_N = -1
    K_C = 200 #in pN/(nm*rad)
    K_G = 5   #in pN/nm

## # The __init__ def is commented out since CrossBridge is no longer a class.
##def __init__(self):
##    """Starts up the crossbridge as a class one fellow (xxCG)
##    """
##        #Define default rest values
##        self.R_T = pi/4 #rest angle (rad) of connection to thick fil
##        self.R_N = 8.5  #rest length (nm) of thick-fil-to-converter rod
##        self.R_C = pi/3 #rest angle (rad) of converter domain
##        self.R_G = 6    #rest length (nm) of globular head domain
##            ##  #Define the current values
##            ##  self.T = pi/4 
##            ##  self.N = 8.5  
##            ##  self.C = pi/3 
##            ##  self.G = 6    
##        #Define default crossbridge as class one
##        #self.SPR_T = False
##        self.SPR_N = False
##        self.SPR_C = True
##        self.SPR_G = True
##        #Define default spring values
##        self.K_T = -1
##        self.K_N = -1
##        self.K_C = 200 #in pN/(nm*rad)
##        self.K_G = 5   #in pN/nm
##    #Define the default energy calculating equation
##    self.En_Eq = lambda G, C: K_C/(2*G)*pow(C-R_C,2) + .5*K_G*pow(G-R_G,2)

def give_me_rest_xb_vals():
    return (R_T, R_N, R_C, R_G)

def give_me_curr_xb_class():
    """This fellow gives back a little string formated in the
    xxCG format to indicate the current type of xb we have
    """
    x_or_let = lambda boo,abc: boo*abc+(not boo)*"x"
    return (x_or_let(SPR_T,"T")+x_or_let(SPR_N,"N")+
            x_or_let(SPR_C,"C")+x_or_let(SPR_G,"G"))

def set_xb_class(xb_class='xxCG',
                 rest_vals = (pi/4, 8.5, pi/3, 6),
                 spr_vals=(-1, -1, 200, 5)):
    """Give me a class, some rest values, and spring values
    and I will change the class of crossbridge currently represented
    """
    #Parse the xb_class variable
    let_to_bool = lambda let: xb_class.count(let)==1
    SPR_T = let_to_bool('T')
    SPR_N = let_to_bool('N')
    SPR_C = let_to_bool('C')
    SPR_G = let_to_bool('G')
    #Parse rest values
    R_T = rest_vals[0] #rest angle (rad) of connection to thick fil
    R_N = rest_vals[1] #rest length (nm) of thick-fil-to-converter rod
    R_C = rest_vals[2] #rest angle (rad) of converter domain
    R_G = rest_vals[3] #rest length (nm) of globular head domain
    #Parse spring values
    K_T = spr_vals[0] #in pN/(nm*rad)
    K_N = spr_vals[1] #in pN/nm
    K_C = spr_vals[2] #in pN/(nm*rad)
    K_G = spr_vals[3] #in pN/nm    

def xb_vals_to_head_loc(xb_vals=None):
    """If you give me the state of the XB component values,
    I'll give you the location of the myosin head.
    """
    #Unpack component values (or their rest positions)
    if xb_vals is None:
        [T,N,C,G]=self.give_me_rest_xb_vals()
    else:
        [T,N,C,G]=xb_vals
    #Calculate the x and y location of the head
    #x_T = 0 #this is assumed
    #y_T = 0 #this is assumed
    x_C = cos(T)*N
    y_C = sin(T)*N
    x_H = x_C + cos(C)*G
    y_H = y_C + sin(C)*G
    #Give back a tuple
    return (x_H, y_H)

def head_loc_to_xb_vals(H=None):
    """If you give me the location of H,
    I will give you the lengths and angles of the xb
    In the form (T, N, C, G)
    """
    # If not given an H, go get the rest H
    if H is None:
        H = xb_vals_to_head_loc()
    if SPR_N == True or SPR_T == True:
        #in the case that C can move from default location, do something
        pass
    elif SPR_C == True and SPR_G == True:
        #Where C is fixed in location, and the system looks like the
        #system from CDW's rotation, solve the old way
        x_C = cos(R_T)*R_N
        y_C = sin(R_T)*R_N
        run_G = H[0]-x_C
        rise_G = H[1]-y_C
        C = atan2(rise_G, run_G)
        G = hypot(run_G, rise_G)
        return (R_T, R_N, C, G)

def xb_energy(xb_vals=None):
    """If you give me a state of the XB segments,
    I'll give you an energy of the whole XB.
    """
    #Unpack component values (or their rest positions)
    if xb_vals is None:
        [T,N,C,G]=give_me_rest_xb_vals()
    else:
        [T,N,C,G]=xb_vals
    U_T = SPR_T*.5*(1/N)*pow(T-R_T,2) #energy in T 
    U_N = SPR_N*.5*pow(N-R_N,2)       #energy in N
    U_C = SPR_C*.5*(1/G)*pow(C-R_C,2) #energy in C
    U_G = SPR_G*.5*pow(G-R_G,2)       #energy in G
    return U_T+U_N+U_C+U_G

def xb_forces(xb_vals=None):
    """If you give me a state of the XB segments,
    I'll give you an array of forces on each node?
    """
    
#--Define the plotting commands we will use
def plotContour(ContourArray, XVals, YVals, Title):
    im = plt.imshow(ContourArray,
                interpolation='bilinear',
                origin='lower',
                cmap=plt.cm.gray,
                extent=(XVals[0],XVals[-1],YVals[0],YVals[-1]),
                alpha=1,
                vmin=-0.01,
                vmax=1.0)
    plt.colorbar()
##  This commented region creates contours
##    cset = plt.contour(ContourArray, arange(.05,.16,.05),
##                   origin='lower',
##                   linewidths=2,
##                   extent=(XVals[0],XVals[-1],YVals[0],YVals[-1]) )
##  This commented region displays labels on contours
##    plt.clabel(cset,
##	inline=1,
##	fmt='%.2f',
##	fontsize=10)
##    plt.autumn()
    plt.axis('image')
    plt.xlabel('Location along thin filament (nm)')
    plt.ylabel('Distance from thin to thick filament (nm)')
    plt.title(Title)
    #savefig(ParamStr+str(ParamValue)+'.png')
    #show()
