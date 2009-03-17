import numpy as np
import matplotlib as matp
import matplotlib.pyplot as plt


def contour(x,y,z):
    """Take 1-d x and y arrays, and a 2-d z array and plot out a contour"""
    # Set something?
    matp.rcParams['xtick.direction'] = 'out'
    matp.rcParams['ytick.direction'] = 'out'
    # Produce the coordinate matrices from two coordinate vectors
    X, Y = np.meshgrid(x, y)
    # Plot the contours with default labeling
    plt.figure()
    im = plt.imshow(z, interpolation='bilinear', origin='lower',
                cmap=matp.cm.jet, extent=(x.min(),x.max(),y.min(),y.max()), norm=matp.colors.Normalize(vmin=0, vmax=30))
    CS = plt.contour(X, Y, z, [1, 3, 5, 7, 9, 11], colors='k',origin='lower',
                 linewidths=2)
    plt.clabel(CS, inline=1, fontsize=9)
    plt.title('xxCG Energy vs Location')
    # Display our contour
    plt.show()