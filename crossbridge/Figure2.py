#!/usr/bin/env python
# encoding: utf-8
"""
Figure2.py
Created by Dave Williams on 2009-07-06.
"""

import sys
import getopt
import Storage
import FigureConstructor
import numpy as np


__help_message__ = '''
Vertical stripes are slimming, also here are some options, pass some along:
Option                Values         Results
-h, --help                           You see this message
-x, --crossbridge     1,2,4          Picks crossbridge to plot
-s, --save            Exists or not  Write output to disk instead of display

Defaults: x=2 and s=Nope
'''


class Usage(Exception):
    """Passes mesages back to the command line"""
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    """Parse options and generate figures for the SingleXB paper"""
    if argv is None:
        argv = sys.argv
    try:
        try:
            short_opts = "hx:s"
            long_opts = ["help" , "crossbridge=", "save"]
            opts, args = getopt.getopt(argv[1:], short_opts, long_opts)
        except getopt.error, msg:
            raise Usage(msg)
        ## Option processing
        save = False # Save to default file name
        xbtype = 2   # Read data for xb type 2
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(__help_message__)
            elif option in ("-x", "--crossbridge"):
                if value in ("1", "2", "4"):
                    xbtype = int(value)
                else:
                    raise Usage("Allowed xb types are 1, 2 and 4 (spring)")
            elif option in ("-s", "--save"):
                print("Using default values")
            else:
                raise Usage("Unhandled option")
        print("Graphing crossbridge with "+str(xbtype)+" springs")
        ## Now for the meat of the issue
        # Load properties that will be needed
        store = Storage.Storage(xbtype)
        free_energy = store.get("free_energy")
        r12 = store.get("r12")
        r23 = store.get("r23")
        r31 = store.get("r31")
        # Load and process x/y related values
        x_range = store.get("x_range")
        x_locs = np.arange(x_range[0], x_range[1], x_range[2])
        y_range = store.get("y_range")
        y_locs = np.arange(y_range[0], y_range[1], y_range[2])
        ## Configure matplotlib
        f_c = FigureConstructor.FigureConstructor((4, 2))
        ## Plot 
        # Set contour levels and cut positions
        RT = 3.97
        energy_levels = (-1e1000*RT, -1*RT, 0*RT, 3*RT, 6*RT, 1e1000*RT)
        trans_prob_levels = (-1, .05, .1, .5, .9, .999, 2)
        cuts = (36, 37, 38)
        # Plot energy two contour
        f_c.two_contour(0, (x_locs, y_locs), free_energy, energy_levels, {
            "title": "Energy of the loosely bound state",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
        })
        # Plot energy cuts
        f_c.cut_plot(1, (x_locs, y_locs), free_energy, cuts, {
            "title": "Energy Slice",
            "y_limits": (-6, 3),
            "y_ticks": (-6, -3, 0, 3)
        })
        # Plot r12 contour
        f_c.two_contour(2, (x_locs, y_locs), r12, trans_prob_levels, {
            "title": "r$_{12}$: unbound to loosely bound",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
        })
        # Plot r12 cuts
        f_c.cut_plot(3, (x_locs, y_locs), r12, cuts, {
            "title": "r$_{12}$ Slice",
            "y_limits": (-.1, 1.1),
            "y_ticks": (0, .5, 1)
        })
        # Plot r23 contour
        f_c.two_contour(4, (x_locs, y_locs), r23, trans_prob_levels, {
            "title": "r$_{23}$: loosely to strongly bound",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
        })
        # Plot r23 cuts
        f_c.cut_plot(5, (x_locs, y_locs), r23, cuts, {
            "title": "r$_{23}$ Slice",
            "y_limits": (-.1, 1.1),
            "y_ticks": (0, .5, 1)
        })
        # Plot r31 contour
        f_c.two_contour(6, (x_locs, y_locs), r31, trans_prob_levels, {
            "title": "r$_{31}$: strongly to unbound",
            "x_lab": "Binding site offset (nm)",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
        })
        # Plot r31 cuts
        f_c.cut_plot(7, (x_locs, y_locs), r31, cuts, {
            "title": "r$_{31}$ Slice",
            "x_lab": "Binding site offset (nm)",
            "y_limits": (-.1, 1.1),
            "y_ticks": (0, .5, 1)
        })
        if save is False:
            # Show plot
            f_c.show_plot()
        else:
            saveloc = "./"+str(xbtype)+"springFig.pdf"
            f_c.save_plot(location=saveloc, trans=False)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2

if __name__ == '__main__':
    sys.exit(main())

