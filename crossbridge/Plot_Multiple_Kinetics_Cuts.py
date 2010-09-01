#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Multiple_Kinetics_Cuts.py
Plots the free energy and forward transition rates at multiple lattice 
spacings for the 4sXB and the 2sXB models.
Created by Dave Williams on 2010-10-26.
"""

import sys
import Storage
from numpy import sqrt, exp, pi, tanh, log
import numpy as np
import matplotlib.pyplot as plt

class Usage(Exception):
    """Passes mesages back to the command line"""
    def __init__(self, msg):
        self.msg = msg


def main():
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)] 
    free_energy = [s.get("free_energy") for s in store]
    r12 = [np.multiply(s.get("r12"), 1000) for s in store]
    r23 = [np.multiply(s.get("r23"), 1000) for s in store]
    r31 = [np.multiply(s.get("r31"), 1000) for s in store]
    force3 = np.array([s.get("force3") for s in store])
    # Parse out rest lattice spacing related values
    # This is a todo
    # Load and process x/y related values
    x_range = store[0].get("x_range") #better be same in the 4 and 2 cases
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    ## Recreate Bert's Data
    # This should be integrated as an external file at some point?
    # Set contour levels and cut positions
    #RT = 3.97
    #cut_locs = np.searchsorted(y_locs, (36, 37, 38))
    cut_vals = (30, 32, 34, 36, 38)
    cut_locs = np.searchsorted(y_locs, cut_vals)
    ## Set up
    fig = plt.figure(1, figsize=(4.86, 6.06))
    axe = ([fig.add_subplot(4, 2, g+1) for g in range(4*2)])
    # colors = ('#1F1E24', '#76D753', '#FF466F', '#F6D246', '#32298F')
    colors = ('#000000', '#0077BB', '#00B840', '#FF5500', '#F6E010')
    lnw=2 # width of plot lines in points
    # Set up a plotting shortcut
    def plot_cuts(axis, prop): 
        for cut,color in zip(cut_locs, colors):
            axis.plot(x_locs, prop[cut], color=color, lw=lnw)
        
    
    ## Plot free energy and transition rates
    # Basic plots 
    plot_cuts(axe[0], free_energy[0]) # Free energy
    plot_cuts(axe[1], free_energy[1])
    plot_cuts(axe[2], r12[0]) # Attachment rate
    plot_cuts(axe[3], r12[1])
    plot_cuts(axe[4], r23[0]) # Power stroke rate
    plot_cuts(axe[5], r23[1])
    plot_cuts(axe[6], r31[0]) # Detachment rate
    plot_cuts(axe[7], r31[1])
    # Axis Labels
    set_x = lambda a,t: a.set_xlabel(t, size=10)
    set_x(axe[6], "Binding site offset (nm)")
    set_x(axe[7], "Binding site offset (nm)")
    set_y = lambda a,t: a.set_ylabel(t, labelpad=14, size=10, 
                                     ha="center", linespacing=1)
    set_y(axe[0], "Free Energy (RT) \n")
    set_y(axe[2], "Binding Rate \n Constant (s$^{-1}$)")
    set_y(axe[4], "Powerstroke Rate \n Constant (s$^{-1}$)")
    set_y(axe[6], "Detachment Rate \n Constant (s$^{-1}$)")
    # Y Axis Limits and Ticks
    for axis in axe[0:2]:
        tics = (-10, -5, 0, 5, 10, 15)
        axis.set_ybound(tics[0], tics[-1])
        axis.set_yticks(tics)
        axis.set_yticklabels(tics, size=9)
        axis.yaxis.set_ticks_position('left')
    for axis in axe[2:]:
        tics = range(0,1001,200)
        axis.set_ybound(tics[0], tics[-1])
        axis.set_yticks(tics)
        axis.set_yticklabels(tics, size=9)
        axis.yaxis.set_ticks_position('left')
    # X Axis Ticks
    for axis in axe:
        tics = range(0,21,5)
        axis.set_xbound(tics[0], tics[-1])
        axis.set_xticks(tics)
        axis.set_xticklabels(tics, size=9)
        axis.xaxis.set_ticks_position('bottom')
    # Titles
    mr_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H"]
    for axis, letter in zip(axe, mr_alphabet):
        axis.set_title(letter, x=-0.20, y=1.04, size=12, weight="demi")
    fig.text(.34,.96, "2sXB Cuts", ha='center', size=12)
    fig.text(.80,.96, "4sXB Cuts", ha='center', size=12)
    # Legend
    leg_vals = [str(val)+" nm" for val in cut_vals]
    axe[5].legend(leg_vals, loc=9, borderpad=0.3,
                  handlelength=1.4, handletextpad=0.1, labelspacing=0.1,
                  fancybox=True, ncol=2, columnspacing=.8, 
                  prop={"size":9})
    ## Display
    fig.subplots_adjust(wspace=0.34, hspace=0.44,
                        left=0.17, right=0.97,
                        top=0.92, bottom=0.08)
    plt.show()


if __name__ == '__main__':
    sys.exit(main())
