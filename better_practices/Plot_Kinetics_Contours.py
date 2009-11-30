#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Kinetics_Contours.py
Created by Dave Williams on 2009-09-09.
"""

import sys
import Storage
from numpy import sqrt, exp, pi, tanh, log
import numpy as np
import matplotlib.pyplot as plt
import pylab

def main():
    # Load properties that will be needed
    print "Loading stored data...",
    store = [Storage.Storage(2), Storage.Storage(4)] 
    free_energy = [s.get("free_energy") for s in store]
    r12 = [s.get("r12") for s in store]
    r23 = [s.get("r23") for s in store]
    r31 = [s.get("r31") for s in store]
    print "done."
    # Load and process x/y related values
    assert(store[0].get("x_range")==store[1].get("x_range"))
    assert(store[0].get("y_range")==store[1].get("y_range"))
    x_range = store[0].get("x_range")
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    x_grid, y_grid = np.meshgrid(x_locs, y_locs)
    y_steps = [26, 30, 34, 38]
    y_steps = [y_steps, [y_locs.searchsorted(y) for y in y_steps]]
    # Set up to the plot
    fig = plt.figure(1, figsize=(8, 10))
    axe = ([fig.add_subplot(4, 2, g+1) for g in range(4*2)])
    colors = ('#1F1E24', '#76D753', '#FF466F', '#F6D246', '#32298F')
    grey_steps = ('.2', '.25', '.3', '.35', '.4', '.45', '.5', '.55', '.6', '.65', '.7', '.75', '.8')
    RT = 3.97
    ener_lvls = (-100*RT, -1*RT, 1*RT, 3*RT, 5*RT, 100*RT)
    prob_lvls = (-0.001, .1, .2, .4, .6, .8, 1.001)
    prob_ticks = (0, .5, 1)
    # Plot the contours
    # Energy of the 4sXB
    axe[0].annotate("4sXB free energy", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    axe[0].contourf(x_grid, y_grid, free_energy[1], 
        ener_lvls, colors = grey_steps)
    # Energy of the 2sXB
    axe[1].annotate("2sXB free energy", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    contour = axe[1].contourf(x_grid, y_grid, free_energy[0], 
        ener_lvls, colors = grey_steps)
    colorbar = plt.colorbar(contour, ax=axe[1], ticks = ener_lvls[1:-1], 
                            format='% 4.1f')
    colorbar.ax.set_position([.92, .80, .07, .13])    
    colorbar.set_label("Energy")
    # Binding rate of the 4sXB
    axe[2].annotate("4sXB attachment rate", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    axe[2].contourf(x_grid, y_grid, r12[1], 
        prob_lvls, colors = grey_steps)
    print "============================"
    print "4sXB max attachment rate loc"
    x_steps = [x_locs[r12[1][y].index(max(r12[1][y]))] for y in y_steps[1]]
    for i in range(len(x_steps)):
        print "@ ", y_steps[0][i], "nm, max x loc is ", x_steps[i]
    # Binding rate of the 2sXB
    axe[3].annotate("2sXB attachment rate", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    contour = axe[3].contourf(x_grid, y_grid, r12[0], 
        prob_lvls, colors = grey_steps)
    colorbar = plt.colorbar(contour, ax=axe[3], format='%.1f', 
                            ticks=prob_ticks, spacing='proportional')
    colorbar.ax.set_position([.92, .57, .07, .13])
    colorbar.set_label("Binding Rate")
    print "============================"
    print "2sXB max attachment rate loc"
    x_steps = [x_locs[r12[0][y].index(max(r12[0][y]))] for y in y_steps[1]]
    for i in range(len(x_steps)):
        print "@ ", y_steps[0][i], "nm, max x loc is ", x_steps[i]
    # Power stroke rate of the 4sXB
    axe[4].annotate("4sXB powerstroke rate", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    axe[4].contourf(x_grid, y_grid, r23[1], 
        prob_lvls, colors = grey_steps)
    print "============================"
    print "4sXB max powerstroke rate loc"
    x_steps = [x_locs[r23[1][y].index(max(r23[1][y]))] for y in y_steps[1]]
    for i in range(len(x_steps)):
        print "@ ", y_steps[0][i], "nm, max x loc is ", x_steps[i]    
    # Power stroke rate of the 2sXB
    axe[5].annotate("2sXB powerstroke rate", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    contour = axe[5].contourf(x_grid, y_grid, r23[0], 
        prob_lvls, colors = grey_steps)
    colorbar = plt.colorbar(contour, ax=axe[5], format='%.1f', 
                            ticks=prob_ticks, spacing='proportional')
    colorbar.ax.set_position([.92, .33, .07, .13])
    colorbar.set_label("Powerstroke Rate")
    print "============================"
    print "2sXB max powerstroke rate loc"
    x_steps = [x_locs[r23[0][y].index(max(r23[0][y]))] for y in y_steps[1]]
    for i in range(len(x_steps)):
        print "@ ", y_steps[0][i], "nm, max x loc is ", x_steps[i]    
    # Detachment rate of the 4sXB
    axe[6].annotate("4sXB detachment rate", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    axe[6].contourf(x_grid, y_grid, r31[1], 
        prob_lvls, colors = grey_steps)
    # Detachment rate of the 2sXB
    axe[7].annotate("2sXB detachment rate", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    contour = axe[7].contourf(x_grid, y_grid, r31[0], 
        prob_lvls, colors = grey_steps)
    colorbar = plt.colorbar(contour, ax=axe[7], format='%.1f', 
                            ticks=prob_ticks, spacing='proportional')
    colorbar.ax.set_position([.92, .09, .07, .13])
    colorbar.set_label("Detachment Rate")
    # Fix and annotate
    titles = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"]
    for a in axe:
        #a.set_yticks((32, 34, 36, 38, 40, 42))
        #a.set_xlim([x_locs[0], x_locs[-1]])
        a.set_ylim([y_locs[0], y_locs[-1]])
        a.set_yticks([26, 30, 34, 38])
        a.set_xlabel("Axial offset (nm)")
        a.set_ylabel("$d_{1,0}$ (nm)")
        a.set_title(titles.pop(0), x=-0.18, y=0.95, weight="demi")
        a.xaxis.set_ticks_position('bottom')
        a.yaxis.set_ticks_position('left')
    # Display the results
    fig.subplots_adjust(wspace=0.35, hspace=0.60,
                        left=0.10, right=0.91,
                        top=0.94, bottom=0.08)
    plt.show()
    

if __name__ == '__main__':
    main()

