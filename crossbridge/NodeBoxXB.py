import cPickle as pickle
import numpy as np

# Load data
#file_name = "head_and_conv_loc.pkl"
file_name = "NodeBoxData.pkl"
stream = open(file_name, 'r')
data = pickle.load(stream)
x_range = data.get('x_range')
y_range = data.get('y_range')
x_locs = np.arange(x_range[0], x_range[1], x_range[2]) 
y_locs = np.arange(y_range[0], y_range[1], y_range[2])
h_locs = data.get('h_locs')
c_locs = data.get('c_locs')
energy = data.get('energy')
cfg = data.get('config')
rest_loc = (
            cfg['N']['weak'] * np.cos(cfg['T']['weak']) + 
            cfg['G']['weak'] * np.cos(cfg['C']['weak'] - np.pi + cfg['T']['weak']), 
            cfg['N']['weak'] * np.sin(cfg['T']['weak']) + 
            cfg['G']['weak'] * np.sin(cfg['C']['weak'] - np.pi + cfg['T']['weak'])            
            )

# Set canvas size and whatnot
scale_factor = 10
dims = (scale_factor*(x_range[1]-x_range[0]), 
        scale_factor*(y_range[1]-y_range[0]))
size(dims[0], dims[1])
speed(30)

# Calculate offsets, account for
ra = 1*scale_factor # Circle radius
if x_range[0]<0: # or the next line gets wonky
    os = (-x_range[0]*scale_factor,ra) # Offset
else: 
    os = (x_range[0]*scale_factor,ra)

# Function to draw an xb head
def draw_xb(conv, head):
    # Scale up conv and head
    conv = [conv[0]*scale_factor, conv[1]*scale_factor]
    head = [head[0]*scale_factor, head[1]*scale_factor]
    conv = [os[0]+conv[0], os[1]+conv[1]]
    head = [os[0]+head[0], os[1]+head[1]]
    # Thick fil
    strokewidth(2)
    stroke(.8, .2, .2)
    fill(.8, .2, .2)
    rect(0, 0, scale_factor*(x_range[1]-x_range[0]), ra)
    # Attachment site
    stroke(.2, .2, .8)
    fill(.2, .2, .8)
    oval(os[0]-ra, os[1]-ra, 2*ra, 2*ra)
    # Neck
    strokewidth(10)
    line(os[0], os[1], conv[0], conv[1])
    # Converter
    strokewidth(2)
    oval(conv[0]-ra, conv[1]-ra, 2*ra, 2*ra)
    # Globular
    strokewidth(10)
    line(conv[0], conv[1], head[0], head[1])
    # Head
    strokewidth(2)
    oval(head[0]-ra, head[1]-ra, 2*ra, 2*ra)
    
# Function to draw the energy at the current level
def draw_energy(energy):
    font("Helvetica", 18)
    text("Energy:", 10, dims[1]-20)
    text("% 5.1f" % energy, 70, dims[1]-20)

def draw():
    # Get the head location's ideal x and y locs
    x = (MOUSEX - os[0])/scale_factor
    y = MOUSEY/scale_factor

    # Get the closest calculated conv and head locations
    x_ind = x_locs.searchsorted(x)-1
    y_ind = y_locs.searchsorted(y)-1

    draw_xb(c_locs[0, y_ind, x_ind], h_locs[0, y_ind, x_ind])
    draw_energy(energy[0, y_ind, x_ind])
    
    #Draw the rest location
    fill(0.0, 1.0, 0.0)
    stroke(0, 1, 0)
    star(rest_loc[0] * scale_factor + os[0], rest_loc[1] * scale_factor + os[1], 9, 10, 5)


    