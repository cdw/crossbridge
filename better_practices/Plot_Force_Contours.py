#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Force_Contours.py
Plot contours of the force magnitude and direction (orientation) for the 
4sXB and the 2sXB.
Created by Dave Williams on 2009-09-11.
"""

import sys
import Storage
from numpy import sqrt, exp, pi, tanh, log
import numpy as np
import matplotlib.pyplot as plt
import pylab

def main():
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)]
    force3 = np.array([s.get("force3") for s in store])
    mag = np.array([[[np.hypot(entry[0], entry[1]) for entry in row] 
        for row in xb] for xb in force3])
    ore = np.array([[[np.arctan2(entry[1], entry[0]) for entry in row] 
        for row in xb] for xb in force3])
    # Load and process x/y related values
    assert(store[0].get("x_range")==store[1].get("x_range"))
    assert(store[0].get("y_range")==store[1].get("y_range"))
    x_range = store[0].get("x_range")
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    x_grid, y_grid = np.meshgrid(x_locs, y_locs)
    # Set up to the plot
    fig = plt.figure(1, figsize=(8, 4))
    axe = ([fig.add_subplot(2, 2, g+1) for g in range(2*2)])
    colors = ('#1F1E24', '#76D753', '#FF466F', '#F6D246', '#32298F')
    grey_steps = ('.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8')
    #ore_colors = ('#e00000', '#ff5000', '#ffd000', '#97ff39', '#00ffad', 
    #'#00c4ff', '#003dff', '#0000c9',
    #'#003dff', '#00c4ff', '#00ffad', '#97ff39', '#ffd000', '#ff5000')
    ore_colors = ('#e00000', '#ffd000', '#00ffad', 
    '#003dff', '#0000c9', '#003dff', '#00ffad', '#ffd000')
    ore_steps = np.arange(-np.pi, np.pi+np.pi*2/len(ore_colors),
                                  np.pi*2/len(ore_colors))
    # # Plot contours
    # force_lvls = (0, 3, 6, 9, 12, 15)
    # # 4sXB Magnitude Contour
    # axe[0].contourf(x_grid, y_grid, mag[1], 
    #     force_lvls, colors = grey_steps)
    # axe[0].contour(x_grid, y_grid, mag[1], 
    #     force_lvls, colors = 'white')
    # axe[0].set_title("4sXB force magnitude")
    # # 2sXB Magnitude Contour
    # axe[1].contourf(x_grid, y_grid, mag[0], 
    #     force_lvls, colors = grey_steps)
    # axe[1].contour(x_grid, y_grid, mag[0], 
    #     force_lvls, colors = 'white')
    # axe[1].set_title("2sXB force magnitude")
    # # 4sXB Directional Contour
    # axe[2].contourf(x_grid, y_grid, ore[0], ore_steps, colors=ore_colors)
    # axe[2].set_title("4sXB force direction")
    # # 2sXB Directional Contour
    # contour = axe[3].contourf(x_grid, y_grid, ore[1], ore_steps, 
    #                           colors=ore_colors)
    # colorbar = plt.colorbar(contour, ax=axe[3], shrink = 0.6, 
    #     ticks=[-3.2, -1.6, 0, 1.6, 3.2])
    # colorbar.ax.set_position([.93, .11, .2, .29])
    # axe[3].set_title("2sXB force direction")
    
    #Plot other contours
    force_lvls = (-15, -10, -5, 0, 5, 10, 15, 20, 25)
    axe[0].set_title("4sXB axial force")
    axe[0].contourf(x_grid, y_grid, force3[1,:,:,0], 
        force_lvls, colors = grey_steps)
    axe[1].set_title("4sXB radial force")
    c = axe[1].contourf(x_grid, y_grid, force3[1,:,:,1], 
        force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[3], shrink = 0.6)
    cb.ax.set_position([.94, .11, .2, .29])
    axe[2].set_title("2sXB axial force")
    axe[2].contourf(x_grid, y_grid, force3[0,:,:,0], 
        force_lvls, colors = grey_steps)
    axe[3].set_title("2sXB radial force")
    c = axe[3].contourf(x_grid, y_grid, force3[0,:,:,1], 
        force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[3], shrink = 0.6)
    cb.ax.set_position([.94, .11, .2, .29])
    
    # Fix the limits and annotate
    for a in axe:
        a.set_xlim([x_locs[0], x_locs[-1]])
        a.set_ylim([y_locs[0], y_locs[-1]])
        a.set_xlabel("Binding site offset (nm)")
        a.set_ylabel("$d_{1,0}$ (nm)")
    # Display the results
    fig.subplots_adjust(wspace=0.30, hspace=0.70,
                        left=0.10, right=0.92,
                        top=0.92, bottom=0.11)
    plt.show()
    


if __name__ == '__main__':
    main()

