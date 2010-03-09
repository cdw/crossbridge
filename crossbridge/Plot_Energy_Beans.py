#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Energy_Beans.py
Created by Dave Williams on 2010-02-07.
"""

# TODO: Add in calculation of 1sXB energies
# TODO: Figure out how to plot without frame if easy
# TODO: Plot cuts through energies at rest lattice spacing.

import sys
import Storage
from numpy import sqrt, exp, pi, tanh, log
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
import pylab

def main():
    # Load properties that will be needed
    print "Loading stored data...",
    store = [Storage.Storage(2), Storage.Storage(4)] 
    energy = [s.get("energy") for s in store]
    print "done."
    # Process limit values
    x_range = store[0].get("x_range")
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    x_grid, y_grid = np.meshgrid(x_locs, y_locs)
    # Create data for the single spring cross-bridge
    e1sxb = single_spring_energy(x_locs, y_locs)
    # Set up the plot
    fig = plt.figure(1, figsize=(32, 24))
    axe = ([fig.add_subplot(3, 2, g+1) for g in range(3*2)])
    lnwdth = 2 # width of plot lines in points
    lnclrs = ('#1F1E24', '#76D753', '#FF466F')
    colors = ('#4444BB','#2288FF','#77DD55','#FFDD44','#FF7722','#FF4466')
    RT = 3.97
    ener_lvls = (-1*RT, .05*RT, .1*RT, .2*RT, .4*RT, .6*RT, .8*RT)
    cut_loc = np.searchsorted(y_locs, 34)
    
    # Plot the energy beans
    ## 1sXB
    axe[0].contourf(x_grid, y_grid, e1sxb, ener_lvls, colors=colors)
    axe[0].axis('off')
    ## 4sXB
    axe[2].contourf(energy[1], ener_lvls, colors=colors, origin='upper')
    axe[2].axis('off')
    ## 2sXB
    contour = axe[4].contourf(energy[0], ener_lvls, colors=colors, origin='upper')
    axe[4].axis('off')
    colorbar = plt.colorbar(contour, ax=axe[4], orientation='horizontal')
    
    # Plot a cut through their energies at rest lattice spacing
    axe[1].plot(x_locs, e1sxb[cut_loc], color=lnclrs[0], lw=lnwdth)
    axe[1].set_ylabel('1sXB energy') 
    axe[3].plot(x_locs, energy[1][cut_loc], color=lnclrs[1], lw=lnwdth)
    axe[3].set_ylabel('2sXB energy') 
    axe[5].plot(x_locs, energy[0][cut_loc], color=lnclrs[2], lw=lnwdth)
    axe[5].set_ylabel('4sXB energy') 
    [a.set_ylim((-.1, 35)) for a in axe[1::2]]
    [a.set_yticks((0, 10, 20, 30)) for a in axe[1::2]]
    [a.set_xlabel('Binding Site Offset (nm)') for a in axe[1::2]]
    # Display the results
    fig.subplots_adjust(wspace=0.35, hspace=0.50,
                        left=0.10, right=0.88,
                        top=0.94, bottom=0.08)
    plt.show()

def single_spring_energy(x_locs, y_locs=[0]):
    """Generate energy in 1sXB as in the days of yore"""
    ## Define parameters
    k_xb = 5 / 3.976 # From Mathematica
    # 3.976 is provided by entering the following into Mathematica:
    #   << PhysicalConstants`
    #   << Units` 
    #   Convert[(MolarGasConstant * 288 Kelvin * AvogadroConstant^-1)
    #       /(Nano Meter)^2, (Pico Newton)/(Nano Meter)]
    xb_0 =13.55 # adjusted from sqrt(eta * DeltaG / k_xb) to match moderns
    # Return energy at x locations 
    x_energies = .5 * k_xb * np.power(x_locs-xb_0, 2) 
    return np.tile(x_energies, (len(y_locs), 1))

if __name__ == '__main__':
    main()
