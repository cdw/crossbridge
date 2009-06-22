import numpy as np
import matplotlib as matp
import matplotlib.pyplot as plt

title = " "
xlabel = " "
ylabel = " "
levels = [1, 5, 10, 15]
cm = matp.cm.jet

def contour(x,y,z):
    """Take 1-d x and y arrays, and a 2-d z array and plot out a contour"""
    # Set something?
    matp.rcParams['xtick.direction'] = 'out'
    matp.rcParams['ytick.direction'] = 'out'
    # Produce the coordinate matrices from two coordinate vectors
    X, Y = np.meshgrid(x, y)
    # Plot the contours with default labeling
    fig = plt.figure()
    im = plt.imshow(z, interpolation='bilinear', origin='upper',
                cmap=cm, extent=(x.min(),x.max(),y.max(),y.min()), 
                norm=matp.colors.Normalize(vmin=min(levels),
                                           vmax=max(levels)))
    CS = plt.contour(X, Y, z, levels, colors='k', linewidths=2)
    plt.clabel(CS, inline=1, fontsize=10, fmt='%1.3f')
    fig.axes[0].invert_yaxis()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # Display our contour

def hexbin(x,y, extremes=None):
    """Feed me an x and a y vector, each of the same length such that 
    (x[i], y[i]) is a point that we would want to count as a 'hit'
    """
    if extremes is None:
        extremes = [x.min(), x.max(), y.max(), y.min()]    
    
    plt.hexbin(x,y, cmap=cm)
    plt.axis(extremes)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cb = plt.colorbar()


def quiver(x,y,Vec):
    plt.quiver(x, y, Vec[:,:,0], Vec[:,:,1])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    
def show():
    plt.show()