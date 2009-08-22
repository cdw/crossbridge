#!/usr/bin/env python
# encoding: utf-8
"""
FigureS1.py
Created by Dave Williams on 2009-08-20.
"""

import sys
import getopt
import Storage
import FigureConstructor
from numpy import sqrt, exp, pi, tanh, log
import numpy as np



__help_message__ = '''
Rolling stones gather no moss, also here are some options, pass some along:
Option                Results
-h, --help            You see this message
-r, --reverse         Include reverse rates (not implemented)
-e, --energy          Include all free energies (not implemented)
-s, --save            Write output to disk instead of display
'''


class Usage(Exception):
    """Passes mesages back to the command line"""
    def __init__(self, msg):
        self.msg = msg

def bert_data(kinetic_prop, x_locs):
    """Generate kinetic properties as in the days of yore"""
    ## Initial constants
    R = 8.314
    T = 288
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
    xb_0 = sqrt(eta * DeltaG / k_xb)
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
        sqrt(k_xb) * (sqrt(M * x**2) - N * x) + P)
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



def main(argv=None):
    """Parse options and generate figures for the SingleXB paper"""
    if argv is None:
        argv = sys.argv
    try:
        try:
            short_opts = "hres"
            long_opts = ["help" , "reverse", "energy", "save"]
            opts, args = getopt.getopt(argv[1:], short_opts, long_opts)
        except getopt.error, msg:
            raise Usage(msg)
        ## Option processing
        save = False # Save to default file name
        show_reverse_rates = False
        show_all_energies = False
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(__help_message__)
            elif option in ("-r", "--reverse"):
                show_reverse_rates = True
            elif option in ("-e", "--energy"):
                show_all_energies = True
            elif option in ("-s", "--save"):
                print("Using default values")
            else:
                raise Usage("Unhandled option")
        ## Now for the meat of the issue
        # Load properties that will be needed
        store = Storage.Storage(2) # Read data for xb type 2
        free_energy = store.get("free_energy")
        r12 = store.get("r12")
        r23 = store.get("r23")
        r31 = store.get("r31")
        # Parse out rest lattice spacing related values
        # This is a todo
        # Load and process x/y related values
        x_range = store.get("x_range")
        x_locs = np.arange(x_range[0], x_range[1], x_range[2])
        y_range = store.get("y_range")
        y_locs = np.arange(y_range[0], y_range[1], y_range[2])
        ## Recreate Bert's Data
        # This is done but needs integration
        ## Configure matplotlib
        f_c = FigureConstructor.FigureConstructor((2, 2))
        ## Plot 
        # Set contour levels and cut positions
        RT = 3.97
        energy_levels = (-1e1000*RT, -1*RT, 0*RT, 3*RT, 6*RT, 1e1000*RT)
        cuts = (36, 37, 38)
        # Plot energy
        f_c.quick_plot(0, (x_locs, bert_data("energy_2", x_locs)))
        f_c.cut_plot(0, (x_locs, y_locs), free_energy, cuts, {
            "title": "Energy of the loosely bound state",
            "y_lab": "Free Energy (RT)",
            "y_ticks": (-10, -5, 0, 5, 10, 15),
            "y_limits": (-10,15)
        })
        # Plot r12 contour
        f_c.quick_plot(1, (x_locs, bert_data("rates_12", x_locs)))
        f_c.cut_plot(1, (x_locs, y_locs), np.multiply(r12, 1000), cuts, {
            "title": "Binding Rates",
            "y_lab": "Transition Rate (s$^{-1}$)",
            "y_ticks": (0, 200, 400, 600, 800, 1000),
            "y_limits": (-.1, 1000.1)
        })
        # Plot r23 contour
        f_c.quick_plot(2, (x_locs, bert_data("rates_23", x_locs)))
        f_c.cut_plot(2, (x_locs, y_locs), np.multiply(r23, 1000), cuts, {
            "title": "Strong binding rates",
            "y_lab": "Transition Rate (s$^{-1}$)",
            "y_ticks": (0, 200, 400, 600, 800, 1000),
            "y_limits": (-.1, 1000.1)
        })
        # Plot r31 contour
        f_c.quick_plot(3, (x_locs, bert_data("rates_31", x_locs)))
        f_c.cut_plot(3, (x_locs, y_locs), np.multiply(r31, 1000), cuts, {
            "title": "Detachment Rates",
            "x_lab": "Binding site offset (nm)",
            "y_lab": "Transition Rate (s$^{-1}$)",
            "y_ticks": (0, 200, 400, 600, 800, 1000),
            "y_limits": (-.1, 1000.1)
        })
        if save is False:
            # Show plot
            f_c.show_plot()
        else:
            saveloc = "./FigureS1.pdf"
            f_c.save_plot(location=saveloc, trans=False)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2

if __name__ == '__main__':
    sys.exit(main())

