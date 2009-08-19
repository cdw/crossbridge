#!/usr/bin/env python
# encoding: utf-8
"""
FigureS2.py
Separated axial and radial components of force exerted by the 
two spring crossbridge at several lattice spacings.
Created by Dave Williams on 2009-07-09.
"""

import Storage
import FigureConstructor
import numpy as np

def main():
    """Generate figure S2 of the SingleXB paper"""
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)]
    force3 = np.array([s.get("force2") for s in store])
    # Note: force3[0,:,:] is all 2 spring data, 0->1 for 4 spring data
    # Select the cuts to take 
    cuts = (36, 37, 38)
    # Load and process x/y related values
    x_range = [s.get("x_range") for s in store]
    x_locs = np.array([np.arange(x[0], x[1], x[2]) for x in x_range])
    y_range = [s.get("y_range") for s in store]
    y_locs = np.array([np.arange(y[0], y[1], y[2]) for y in y_range])
    y_ticks = (-50,0,50)
    # If we plot bits outside our y window, the external arrows can stick into 
    # the window of interest, so trim the data
    #y_ind = (np.searchsorted(y_locs[0, ::2], y_limits[0])-1, np.searchsorted(y_locs[0, ::2], y_limits[1]))
    
    ## Configure matplotlib
    f_c = FigureConstructor.FigureConstructor((2, 2), (5,10))
    ## Plot
    f = 5 #horrible hack to cut off bad values FIXME at source, 4 spring force
    # Need this fella quiver_cut(self, sub, x_y_values, j_k_grid, cut_component, cut_locs, labels_n_limits={}):
    f_c.quiver_cut(0, (x_locs[0, ::2], y_locs[0, f::2]), 
                    force3[0, f::2, ::2], 0, cuts, {
                    "title": "X component of force for 2 spring crossbridge",
                    "x_lab": "Binding site offset (nm)",
                    "y_lab": "x/axial component of force",
                    "y_limits": (y_ticks[0], y_ticks[-1]),
                    "y_ticks": y_ticks
                    })
    f_c.quiver_cut(1, (x_locs[0, ::2], y_locs[0, f::2]), 
                    force3[0, f::2, ::2], 1, cuts, {
                    "title": "Y component of force for 2 spring crossbridge",
                    "x_lab": "Binding site offset (nm)",
                    "y_lab": "y/radial component of force",
                    "y_limits": (y_ticks[0], y_ticks[-1]),
                    "y_ticks": y_ticks
                    })
    f_c.quiver_cut(2, (x_locs[0, ::2], y_locs[0, f::2]), 
                    force3[1, f::2, ::2], 0, cuts, {
                    "title": "X component of force for 4 spring crossbridge",
                    "x_lab": "Binding site offset (nm)",
                    "y_lab": "x/axial component of force",
                    "y_limits": (y_ticks[0], y_ticks[-1]),
                    "y_ticks": y_ticks
                    })
    f_c.quiver_cut(3, (x_locs[0, ::2], y_locs[0, f::2]), 
                    force3[1, f::2, ::2], 1, cuts, {
                    "title": "Y component of force for 4 spring crossbridge",
                    "x_lab": "Binding site offset (nm)",
                    "y_lab": "y/radial component of force",
                    "y_limits": (y_ticks[0], y_ticks[-1]),
                    "y_ticks": y_ticks
                    })
    # Show plot
    f_c.show_plot()



if __name__ == '__main__':
    main()

