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
        xbtype = 2   # Read data for xb type 2
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
        store = Storage.Storage(xbtype)
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
        f_c.two_contour(0, (x_locs, y_locs), free_energy, energy_levels, {
            "title": "Energy of the loosely bound state",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
        })
        # Plot r12 contour
        f_c.two_contour(1, (x_locs, y_locs), r12, trans_prob_levels, {
            "title": "r$_{12}$: unbound to loosely bound",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
        })
        # Plot r23 contour
        f_c.two_contour(2, (x_locs, y_locs), r23, trans_prob_levels, {
            "title": "r$_{23}$: loosely to strongly bound",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
        })
        # Plot r31 contour
        f_c.two_contour(3, (x_locs, y_locs), r31, trans_prob_levels, {
            "title": "r$_{31}$: strongly to unbound",
            "x_lab": "Binding site offset (nm)",
            "y_lab": "d$_{1,0}$ (nm)",
            "y_ticks": (34, 36, 38, 40, 42),
            "y_limits": (34,42),
            "cuts": {"cut_locs":cuts}
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

