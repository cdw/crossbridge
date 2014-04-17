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
    store = [Storage.Storage(2), Storage.Storage(4)] 
    free_energy = [np.flipud(s.get("free_energy")) for s in store]
    r12 = [np.flipud(np.multiply(s.get("r12"), 1000)) for s in store]
    r23 = [np.flipud(np.multiply(s.get("r23"), 1000)) for s in store]
    r31 = [np.flipud(np.multiply(s.get("r31"), 1000)) for s in store]
    # Load and process x/y related values
    assert(store[0].get("x_range")==store[1].get("x_range"))
    assert(store[0].get("y_range")==store[1].get("y_range"))
    x_range = store[0].get("x_range")
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[1], y_range[0], -y_range[2])
    x_grid, y_grid = np.meshgrid(x_locs, y_locs)
    # Set up to the plot
    fig = plt.figure(1, figsize=(4.68, 5.85))
    axe = ([fig.add_subplot(4, 2, g+1) for g in range(4*2)])
    colors = ('#4444BB','#2288FF','#77DD55','#FFDD44','#FF7722','#FF4466')
    grey_steps = ('.2', '.25', '.3', '.35', '.4', '.45', '.5', '.55', 
                  '.6', '.65', '.7', '.75', '.8')
    RT = 3.97
    ener_lvls = (-100*RT, -1*RT, 1*RT, 3*RT, 5*RT, 100*RT)
    ener_ticks = ("-1RT", "RT", "3RT", "5RT")
    prob_lvls = (-1, 100, 200, 400, 600, 800, 1001)
    prob_ticks = (0, 500, 1000)
    # Plotting shortcut
    def plot_contour(axis, prop, levels):
        cont = axis.contourf(x_grid, y_grid, prop, levels, colors=colors,
                            origin='upper')
        return cont 
    
    # Basic Plots
    plot_contour(axe[0], free_energy[1], ener_lvls)
    energy_contour = plot_contour(axe[1], free_energy[0], ener_lvls)
    plot_contour(axe[2], r12[1], prob_lvls)
    r12_contour = plot_contour(axe[3], r12[0], prob_lvls)
    plot_contour(axe[4], r23[1], prob_lvls)
    r23_contour = plot_contour(axe[5], r23[0], prob_lvls)
    plot_contour(axe[6], r31[1], prob_lvls)
    r31_contour = plot_contour(axe[7], r31[0], prob_lvls)
    # Colorbars
    def format_colorbar(colorbar, title, tic_labels, pos):
        colorbar.set_label(title, size=9, ha='center',
                          linespacing=0.9, labelpad=13)
        colorbar.ax.set_yticklabels(tic_labels, size=8)
        colorbar.ax.set_position(pos) 
    colorbar = plt.colorbar(energy_contour, ax=axe[1], 
                            ticks = ener_lvls[1:-1], format='% 4.1f')
    format_colorbar(colorbar, "Energy", ener_ticks,
                    [.85, .805, .07, .13]) 
    colorbar = plt.colorbar(r12_contour, ax=axe[3], format='%d', 
                            ticks=prob_ticks, spacing='proportional')
    format_colorbar(colorbar, "Binding \n Rate Constant (s$^{-1}$)",
                 prob_ticks, [.85, .565, .07, .13]) 
    colorbar = plt.colorbar(r23_contour, ax=axe[5], format='%d', 
                            ticks=prob_ticks, spacing='proportional')
    format_colorbar(colorbar, "Power stroke \n Rate Constant (s$^{-1}$)",
                 prob_ticks, [.85, .325, .07, .13]) 
    colorbar = plt.colorbar(r31_contour, ax=axe[7], format='%d', 
                            ticks=prob_ticks, spacing='proportional')
    format_colorbar(colorbar, "Detachment \n Rate Constant (s$^{-1}$)",
                 prob_ticks, [.85, .085, .07, .13]) 
    #colorbar.ax.set_position([.89, .80, .07, .13]) 
    #colorbar.ax.set_yticklabels(ener_ticks, size=8)
    #colorbar.set_label("Energy")
    
    #colorbar.ax.set_position([.89, .57, .07, .13])
    #colorbar.set_label("Binding \n Rate Constant (s$^-1)")
    
    #colorbar.ax.set_position([.89, .33, .07, .13])
    #colorbar.set_label("Power stroke  \n Rate Constant (s$^-1)")
    #colorbar = plt.colorbar(r31_contour, ax=axe[7], format='%d', 
    #                        ticks=prob_ticks, spacing='proportional')
    #colorbar.ax.set_position([.89, .09, .07, .13])
    #colorbar.set_label("Detachment  \n Rate Constant (s$^-1)")
    # Annotate the axes
    annotate = lambda axis, title: axis.annotate(title, (0.5, 1.05),
                                 xycoords='axes fraction', ha='center',
                                 size=10)
    annotate(axe[0], "4sXB free energy")
    annotate(axe[1], "2sXB free energy")
    annotate(axe[2], "4sXB attachment \n rate constant")
    annotate(axe[3], "2sXB attachment \n rate constant")
    annotate(axe[4], "4sXB power stroke \n rate constant")
    annotate(axe[5], "2sXB power stroke \n rate constant")
    annotate(axe[6], "4sXB detachment \n rate constant")
    annotate(axe[7], "2sXB detachment \n rate constant")
    # Axis Labels
    set_x = lambda a,t: a.set_xlabel(t, size=11, linespacing=1.0)
    [set_x(axis, "Axial offset (nm)") for axis in axe[6:]]
    set_y = lambda a,t: a.set_ylabel(t, labelpad=9, size=9.5, 
                                     ha="center", linespacing=1)
    set_y(axe[0], "$d_{1,0}$ (nm)")
    set_y(axe[2], "$d_{1,0}$ (nm)")
    set_y(axe[4], "$d_{1,0}$ (nm)")
    set_y(axe[6], "$d_{1,0}$ (nm)")
    # Y Axis Limits and Ticks
    for axis in axe:
        limits = [y_range[0], y_range[1]]
        tics = [26, 30, 34, 38] 
        axis.set_yticks(tics)
        axis.set_yticklabels(tics, size=8)
        axis.yaxis.set_ticks_position('left')
        axis.set_ylim(limits[1], limits[0])
    # X Axis Ticks
    for axis in axe:
        tics = range(0,21,5)
        axis.set_xbound(tics[0], tics[-1])
        axis.set_xticks(tics)
        axis.set_xticklabels(tics, size=8)
        axis.xaxis.set_ticks_position('bottom')
    # Sub-figure Titles
    mr_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H"]
    for axis, letter in zip(axe, mr_alphabet):
        axis.set_title(letter, x=-0.20, y=1.04, size=12, weight="demi")
    # Display the results
    fig.subplots_adjust(wspace=0.32, hspace=0.70,
                        left=0.10, right=0.84,
                        top=0.94, bottom=0.08)
    plt.show()
    

if __name__ == '__main__':
    main()

