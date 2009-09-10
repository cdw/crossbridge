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
    # Set up to the plot
    fig = plt.figure(1, figsize=(7, 7))
    axe = ([fig.add_subplot(4, 2, g+1) for g in range(4*2)])
    colors = ('#1F1E24', '#76D753', '#FF466F', '#F6D246', '#32298F')
    grey_steps = ('.2', '.3', '.4', '.5', '.6', '.7', '.8')
    RT = 3.97
    ener_lvls = (-1e1000*RT, -1*RT, 0*RT, 3*RT, 6*RT, 1e1000*RT)
    prob_lvls = (-1, .05, .1, .5, .9, .999, 2)
    # Plot the contours
    # Energy of the 4sXB
    contour = axe[0].contourf(x_grid, y_grid, free_energy[1], 
        ener_lvls, colors = grey_steps)
    axe[0].contour(x_grid, y_grid, free_energy[1], 
        ener_lvls, colors = 'white')
    # Energy of the 2sXB
    axe[1].contourf(x_grid, y_grid, free_energy[0], 
        ener_lvls, colors = grey_steps)
    axe[1].contour(x_grid, y_grid, free_energy[0], 
        ener_lvls, colors = 'white')
    # Binding rate of the 4sXB
    axe[2].contourf(x_grid, y_grid, r12[1], 
        prob_lvls, colors = grey_steps)
    axe[2].contour(x_grid, y_grid, r12[1], 
        prob_lvls, colors = 'white')
    # Binding rate of the 2sXB
    axe[3].contourf(x_grid, y_grid, r12[0], 
        prob_lvls, colors = grey_steps)
    axe[3].contour(x_grid, y_grid, r12[0], 
        prob_lvls, colors = 'white')
    # Power stroke rate of the 4sXB
    axe[4].contourf(x_grid, y_grid, r23[1], 
        prob_lvls, colors = grey_steps)
    axe[4].contour(x_grid, y_grid, r23[1], 
        prob_lvls, colors = 'white')
    # Power stroke rate of the 2sXB
    axe[5].contourf(x_grid, y_grid, r23[0], 
        prob_lvls, colors = grey_steps)
    axe[5].contour(x_grid, y_grid, r23[0], 
        prob_lvls, colors = 'white')
    # Detachment rate of the 4sXB
    axe[6].contourf(x_grid, y_grid, r31[1], 
        prob_lvls, colors = grey_steps)
    axe[6].contour(x_grid, y_grid, r31[1], 
        prob_lvls, colors = 'white')
    # Detachment rate of the 2sXB
    axe[7].contourf(x_grid, y_grid, r31[0], 
        prob_lvls, colors = grey_steps)
    axe[7].contour(x_grid, y_grid, r31[0], 
        prob_lvls, colors = 'white')
    # Fix the and annotate
    titles = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"]
    for a in axe:
        a.set_yticks((32, 34, 36, 38, 40, 42))
        a.set_xlim([x_locs[0], x_locs[-1]])
        a.set_ylim([y_locs[0], y_locs[-1]])
        a.set_xlabel("Binding site offset (nm)")
        a.set_ylabel("$d_{1,0}$ (nm)")
        a.set_title(titles.pop(0), x=-0.18, y=0.95, weight="demi")
    # Display the results
    fig.subplots_adjust(wspace=0.30, hspace=0.40,
                        left=0.10, right=0.95,
                        top=0.94, bottom=0.08)
    plt.show()
    

if __name__ == '__main__':
    main()

