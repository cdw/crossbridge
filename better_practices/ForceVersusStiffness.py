#!/usr/bin/env python
# encoding: utf-8
"""
ForceVersusStiffness.py

Created by Dave Williams on 2009-11-08.
Copyright (c) 2009 Dave Williams. All rights reserved.
"""
import Crossbridge
import numpy as np
from numpy import pi, radians
import matplotlib.pyplot as plt

def main():
    """Plot the axial and radial forces produced by crossbridges with 
    various levels of globular domain stiffness.
    """
    print "Starting config and whatnot... ",
    ## Configure the variables that we are investigating
    xbs = 3                       # Crossbridge state
    ls = 15.77                    # Lattice spacing
    svals = np.arange(1, 10, 0.5) # Stiffness values
    ## Set parameters that we will be changing less often
    x_locs = np.arange(0, 20, 0.1) 
    c4sXB = {
        'T': {
            'weak': radians(40),
            'strong': radians(40),
            'spring_konstant': 100
        },
        'N': {
            'weak': 10.5,
            'strong': 10.5,
            'spring_konstant': 10
        },
        'C': {
            'weak': 2*pi - radians(165), #pi/3 + (pi - radians(40)), #40deg from T weak
            'strong': 2*pi - radians(110), # pi/2 + (pi - pi/4),
            'spring_konstant': 40
        },
        'G': {
            'weak': 9.6,
            'strong': 9.6,
            'spring_konstant': 5
        }
    }
    c2sXB = {
         'T': {
             'weak': 0,
             'strong': 0,
             'spring_konstant': 1
         },
         'N': {
             'weak': 0,
             'strong': 0,
             'spring_konstant': 1
         },
         'C': {
             'weak': radians(47.16) + pi,
             'strong': radians(73.20) + pi,
             'spring_konstant': 40
         },
         'G': {
             'weak': 19.93,
             'strong': 16.47,
             'spring_konstant': 2
         }
    }
    ## Process the setup variables
    xb4 = Crossbridge.FourSpring(c4sXB)
    xb2 = Crossbridge.TwoSpring(c2sXB)
    d10 = 1.5 * (ls + 6.90) # See fil_sep_to_d10 in CreateData
    f4sXB = np.zeros((svals.size, x_locs.size, 2))
    f2sXB = np.zeros((svals.size, x_locs.size, 2))
    print "done."
    print "Generating data...",
    ## Generate data
    i = 0
    for k in svals:
        print str(i),
        xb4.g.k = k
        f4sXB[i, :] = [xb4.force((x, ls), xbs) for x in x_locs]
        xb2.g.k = k
        f2sXB[i, :] = [xb2.force((x, ls), xbs) for x in x_locs]
        i += 1
    print " done."
    ## Set up the plot for the data
    fig = plt.figure(1, figsize=(8, 4))
    axe = ([fig.add_subplot(2, 2, g+1) for g in range(2*2)])
    grey_steps = ('.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8')
    force_lvls = (-15, -10, -5, 0, 5, 10, 15, 20, 25)
    x_grid, y_grid = np.meshgrid(x_locs, svals)
    ## Plot the contours
    axe[0].set_title("4sXB axial force")
    axe[0].contourf(x_grid, y_grid, f4sXB[:, :, 0], 
        force_lvls, colors = grey_steps)
    axe[2].set_title("4sXB radial force")
    c = axe[2].contourf(x_grid, y_grid, f4sXB[:, :, 1], 
        force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[3], shrink = 0.6)
    cb.ax.set_position([.94, .11, .2, .29])
    axe[1].set_title("2sXB axial force")
    axe[1].contourf(x_grid, y_grid, f2sXB[:, :, 0], 
        force_lvls, colors = grey_steps)
    axe[3].set_title("2sXB radial force")
    c = axe[3].contourf(x_grid, y_grid, f2sXB[:, :, 1], 
        force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[3], shrink = 0.6)
    cb.ax.set_position([.94, .11, .2, .29])
    ## Fix the limits and annotate axes
    for a in axe:
        a.set_xlim([x_locs[0], x_locs[-1]])
        a.set_ylim([svals[0], svals[-1]])
        a.set_xlabel("Binding site offset (nm)")
        a.set_ylabel("Linear Spring\n Stiffness (pN)")
    ## Display the results
    fig.subplots_adjust(wspace=0.30, hspace=0.70,
                  left=0.10, right=0.92,
                  top=0.92, bottom=0.11) 
    plt.show()


if __name__ == '__main__':
    main()

