#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Kinetics_Cuts.py
Plots the free energy and forward transition rates at a resting lattice 
spacing for the 4sXB, the 2sXB, and the historical single spring crossbridge.
Created by Dave Williams on 2009-09-09.
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


def bert_data(kinetic_prop, x_locs):
    """Generate kinetic properties as in the days of yore"""
    ## Initial constants
    #R = 8.314
    #T = 288
    Gatp = 13 # In units of RT
    atp = 5  * 10**-3
    adp = 30 * 10**-6
    phos  = 3  * 10**-3
    DeltaG = abs(-Gatp - log(atp / (adp * phos)))
    alpha = 0.28
    eta = 0.68
    ## Define parameters
    A = 2000  # From Tanner, 2008 Pg 1209
    B = 100   # Ditto for C through P
    C = 1
    D = 1
    M = 3600
    N = 40
    P = 20
    ## Calculate crossbridge spring values
    k_xb = 5 / 3.976 # From Mathematica
    # 3.976 is provided by entereing the following into Mathematica:
    #   << PhysicalConstants`
    #   << Units` 
    #   Convert[(MolarGasConstant * 288 Kelvin * AvogadroConstant^-1)
    #       /(Nano Meter)^2, (Pico Newton)/(Nano Meter)]
    xb_0 =13.55 # adjusted from sqrt(eta * DeltaG / k_xb) to match moderns
    xb_1 =5 # adjusted from sqrt(eta * DeltaG / k_xb) to match moderns
    ## Create functions for free energy in a given location and state:
    g_1 = lambda x: 0 * x
    g_2 = lambda x: alpha * -DeltaG + k_xb * (x - xb_0)**2
    g_3 = lambda x: eta * -DeltaG + k_xb * x**2
    ## Create functions to yield transition rate at a given location
    r_12 = lambda x: (
        A * sqrt(k_xb / (2 * pi)) * exp(-.5 * k_xb * (x - xb_0)**2))
    r_23 = lambda x: (
        B / sqrt(k_xb) * (1 - tanh(C * sqrt(k_xb) * (x - xb_0))) + D)
    r_31 = lambda x: (
        sqrt(k_xb) * (sqrt(M * (x-4.76)**2) - N * (x-4.76)) + P)
    ## Create functions to derive reverse transition rates using the
    ## thermo equation:
    ## r_forward(x)/r_backward(x) = exp(G_next(x) - G_prev(x))
    r_21 = lambda x: r_12(x) / exp(g_1(x) - g_2(x))
    r_32 = lambda x: r_23(x) / exp(g_2(x) - g_3(x))
    r_13 = lambda x: 0 * x # See Tanner, 2007 Pg 1209 for justification]
    ## Create functions to return the force generated
    af = lambda x: - k_xb * (xb_1 - x)
    rf = lambda x: x * 0
    ## Create the data
    #x_locs = np.arange(-5, 15, .1)
    if kinetic_prop == "energy_1":
        return g_1(x_locs)  # Free energy state 1
    elif kinetic_prop == "energy_2":
        return g_2(x_locs)  # Free energy state 2
    elif kinetic_prop == "energy_3":
        return g_3(x_locs)  # Free energy state 3
    elif kinetic_prop == "rates_12":
        return r_12(x_locs) # Forward rate 12
    elif kinetic_prop == "rates_21":
        return r_21(x_locs) # Reverse rate 21
    elif kinetic_prop == "rates_23":
        return r_23(x_locs) # Forward rate 23
    elif kinetic_prop == "rates_32":
        return r_32(x_locs) # Reverse rate 32
    elif kinetic_prop == "rates_31":
        return r_31(x_locs) # Forward rate 31
    elif kinetic_prop == "rates_13":
        return r_13(x_locs) # Reverse rate 13
    elif kinetic_prop == "axial_force":
        return af(x_locs)    # Axial force
    elif kinetic_prop == "radial_force":
        return rf(x_locs)    # Radial force
    else:
        raise Usage("Unhandled crossbridge property")


def main():
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)] 
    free_energy = [s.get("free_energy") for s in store]
    r12 = [np.multiply(s.get("r12"), 1000) for s in store]
    r23 = [np.multiply(s.get("r23"), 1000) for s in store]
    r31 = [np.multiply(s.get("r31"), 1000) for s in store]
    force3 = np.array([s.get("force3") for s in store])
    axial_force = list(force3[:,:,:,0])
    radial_force = list(force3[:,:,:,1])
    # Parse out rest lattice spacing related values
    # This is a todo
    # Load and process x/y related values
    x_range = store[0].get("x_range") #better be same in the 4 and 2 cases
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    ## Recreate Bert's Data
    # This should be integrated as an external file at some point?
    free_energy.append(bert_data("energy_2", x_locs))
    r12.append(bert_data("rates_12", x_locs))
    r23.append(bert_data("rates_23", x_locs))
    r31.append(bert_data("rates_31", x_locs))
    axial_force.append(bert_data("axial_force", x_locs))  
    radial_force.append(bert_data("radial_force", x_locs)) 
    # Set contour levels and cut positions
    #RT = 3.97
    cut_locs = np.searchsorted(y_locs, [31, 34])
    cut = np.searchsorted(y_locs, 34)
    ## Set up
    fig = plt.figure(1, figsize=(6.83, 3.30))
    axe = ([fig.add_subplot(2, 3, g+1) for g in range(3*2)])
    colors = ('#1F1E24', '#76D753', '#FF466F', '#F6D246', '#32298F')
    lnw=2 # width of plot lines in points
    # Set up a plotting shortcut
    def plot_cuts(axis, prop): 
        axis.plot(x_locs, prop[2], color=colors[0], lw=lnw)
        axis.plot(x_locs, prop[0][cut], color=colors[1], lw=lnw)
        axis.plot(x_locs, prop[1][cut], color=colors[2], lw=lnw)
        
    
    ## Plot free energy and transition rates
    # Basic plots
    plot_cuts(axe[0], free_energy)
    plot_cuts(axe[1], r12)
    plot_cuts(axe[2], r23)
    plot_cuts(axe[3], r31)
    plot_cuts(axe[4], axial_force)
    plot_cuts(axe[5], radial_force)
    # Axis Labels
    set_x = lambda a,t: a.set_xlabel(t, size=10)
    [set_x(axis, "Binding site offset (nm)") for axis in axe[3:]]
    set_y = lambda a,t: a.set_ylabel(t, labelpad=14, size=9.5, 
                                     ha="center", linespacing=1)
    set_y(axe[0], "Free Energy (RT) \n")
    set_y(axe[1], "Binding Rate \n Constant (s$^{-1}$)")
    set_y(axe[2], "Powerstroke Rate \n Constant (s$^{-1}$)")
    set_y(axe[3], "Detachment Rate \n Constant (s$^{-1}$)")
    set_y(axe[4], "Axial force (pN) \n")
    set_y(axe[5], "Radial Force (pN) \n")
    # Y Axis Limits and Ticks
    def set_y_axis(axis, tics):
        axis.set_ybound(tics[0], tics[-1])
        axis.set_yticks(tics)
        axis.set_yticklabels(tics, size=8)
        axis.yaxis.set_ticks_position('left')
    
    set_y_axis(axe[0], np.arange(-10, 16, 5))
    for axis in axe[1:4]:
        set_y_axis(axis, np.arange(0, 1001, 200))
    set_y_axis(axe[4], np.arange(-5, 21, 5))
    set_y_axis(axe[5], np.arange(-5, 16, 5))
    # X Axis Ticks
    for axis in axe:
        tics = range(0,21,5)
        axis.set_xbound(tics[0], tics[-1])
        axis.set_xticks(tics)
        axis.set_xticklabels(tics, size=8)
        axis.xaxis.set_ticks_position('bottom')
    # Titles
    mr_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H"]
    for axis, letter in zip(axe, mr_alphabet):
        axis.set_title(letter, x=-0.20, y=1.04, size=12, weight="demi")
    # Legend 
    legend_vals = ("1sXB", "2sXB", "4sXB")
    axe[2].legend(legend_vals , loc=0, borderpad=0.3,
                  handlelength=1.8, handletextpad=0.3, labelspacing=0.2,
                  fancybox=True, prop={'size':9})
    ## Display
    fig.subplots_adjust(left=0.12,   bottom=0.13, 
                        right=0.96,  top=0.88, 
                        wspace=0.68, hspace=0.50)
    plt.show()


if __name__ == '__main__':
    sys.exit(main())
