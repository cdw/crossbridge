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
    r_13 = lambda x: 0 * x # See Tanner, 2007 Pg 1209 for justification
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
    else:
        raise Usage("Unhandled crossbridge property")


def main():
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)] 
    free_energy = [s.get("free_energy") for s in store]
    r12 = [np.multiply(s.get("r12"), 1000) for s in store]
    r23 = [np.multiply(s.get("r23"), 1000) for s in store]
    r31 = [np.multiply(s.get("r31"), 1000) for s in store]
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
    cut_locs = np.searchsorted(y_locs, (31, 34, 37))
    ## Set up
    fig = plt.figure(1, figsize=(8, 6))
    axe = ([fig.add_subplot(2, 2, g+1) for g in range(2*2)])
    colors = ('#1F1E24', '#76D753', '#FF466F', '#F6D246', '#32298F')
    lnwdth=2
    ## Plot free energy and transition rates
    # Energy of the loosely bound state
    axe[0].plot(x_locs, bert_data("energy_2", x_locs), color=colors[0],
                lw=lnwdth)
    axe[0].plot(x_locs, free_energy[0][cut_locs[1]], color=colors[1],
                lw=lnwdth)
    axe[0].plot(x_locs, free_energy[1][cut_locs[1]], color=colors[2],
                lw=lnwdth)
    axe[0].set_xlabel("Binding site offset (nm)")
    axe[0].set_ylabel("Free Energy (RT)")
    axe[0].set_ylim((-10, 15))
    axe[0].set_yticks((-10, -5, 0, 5, 10, 15))
    axe[0].set_title("a)", x=-0.20, y=0.97, weight="demi")
    axe[0].legend(("1sXB", "2sXB", "4sXB"), loc=0, borderpad=0.3,
                  handlelength=1.8, handletextpad=0.3, labelspacing=0.2,
                  fancybox=True)
    # Binding rates
    axe[1].plot(x_locs, bert_data("rates_12", x_locs), color=colors[0],
                lw=lnwdth)
    axe[1].plot(x_locs, r12[0][cut_locs[1]], color=colors[1], lw=lnwdth)
    axe[1].plot(x_locs, r12[1][cut_locs[1]], color=colors[2], lw=lnwdth)
    axe[1].set_title("b)", x=-0.20, y=0.97, weight="demi")
    axe[1].set_ylabel("Binding Rate (s$^{-1}$)")
    # Powerstroke rates
    axe[2].plot(x_locs, bert_data("rates_23", x_locs), color=colors[0],
                lw=lnwdth)
    axe[2].plot(x_locs, r23[0][cut_locs[1]], color=colors[1], lw=lnwdth)
    axe[2].plot(x_locs, r23[1][cut_locs[1]], color=colors[2], lw=lnwdth)
    axe[2].set_title("c)", x=-0.20, y=0.97, weight="demi") 
    axe[2].set_ylabel("Powerstroke Rate (s$^{-1}$)")
    # Detachment rates
    axe[3].plot(x_locs, bert_data("rates_31", x_locs), color=colors[0],
                lw=lnwdth)
    axe[3].plot(x_locs, r31[0][cut_locs[1]], color=colors[1], lw=lnwdth)
    axe[3].plot(x_locs, r31[1][cut_locs[1]], color=colors[2], lw=lnwdth)
    axe[3].set_title("d)", x=-0.20, y=0.97, weight="demi")
    axe[3].set_ylabel("Detachment Rate (s$^{-1}$)")
    # Add lables and limits
    for a in axe[1:]:
        a.set_xlabel("Binding site offset (nm)")
        a.set_ylim((-.1, 1000.1))
        a.set_yticks((0, 200, 400, 600, 800, 1000))
    ## Display
    fig.subplots_adjust(wspace=0.35, hspace=0.40,
                        left=0.10, right=0.95,
                        top=0.94, bottom=0.08)
    plt.show()


if __name__ == '__main__':
    sys.exit(main())
