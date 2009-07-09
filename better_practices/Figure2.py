#!/usr/bin/env python
# encoding: utf-8
"""
Figure1.py
Created by Dave Williams on 2009-07-06.
"""

#import sys
#import os
import Storage
import FigureConstructor
import numpy as np


def main():
    """Generate the second figure of the SingleXB paper"""
    # Load properties that will be needed
    store = Storage.Storage(2)
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
    cuts = (8, 10, 12)
    # Plot energy two contour
    f_c.two_contour(0, (x_locs, y_locs), free_energy, energy_levels, {
        "title": "Energy of the loosely bound state",
        "y_lab": "Spacing"
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
        "y_lab": "Spacing"
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
        "y_lab": "Spacing"
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
        "y_lab": "Spacing"
    })
    # Plot r31 cuts
    f_c.cut_plot(7, (x_locs, y_locs), r31, cuts, {
        "title": "r$_{31}$ Slice",
        "x_lab": "Binding site offset (nm)",
        "y_limits": (-.1, 1.1),
        "y_ticks": (0, .5, 1)
    })
    # Show plot
    f_c.show_plot()

if __name__ == '__main__':
    main()

