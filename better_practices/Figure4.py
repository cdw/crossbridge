#!/usr/bin/env python
# encoding: utf-8
"""
Figure4.py
Created by Dave Williams on 2009-07-09.
"""

#import sys
#import os
import Storage
import FigureConstructor
import numpy as np

def main():
    """Generate the fourth figure of the SingleXB paper"""
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)]
    force3 = np.array([s.get("force2") for s in store])
    # Load and process x/y related values
    x_range = [s.get("x_range") for s in store]
    x_locs = np.array([np.arange(x[0], x[1], x[2]) for x in x_range])
    y_range = [s.get("y_range") for s in store]
    y_locs = np.array([np.arange(y[0], y[1], y[2]) for y in y_range])
    y_ticks = (34, 34, 36, 38, 40, 42)
    y_limits = (34,42)
    # If we plot bits outside our y window, the external arrows can stick into 
    # the window of interest, so trim the data
    #y_ind = (np.searchsorted(y_locs[0, ::2], y_limits[0])-1, np.searchsorted(y_locs[0, ::2], y_limits[1]))
    
    ## Configure matplotlib
    f_c = FigureConstructor.FigureConstructor((3, 1), (5,10))
    ## Plot
    force_scale = 300
    f = 5 #horrible hack to cut off bad values FIXME at source, 4 spring force
    f_c.quiver_plot(0, (x_locs[0, ::2], y_locs[0, f::2]), 
                    force3[0, f::2, ::2], {
                    "title": "Force of a loosely bound 2 spring crossbridge",
                    "y_lab": "d$_{1,0}$ spacing (nm)",
                    "y_ticks": y_ticks,
                    "y_limits": y_limits,
                    "scale": force_scale
                    })
    f_c.quiver_plot(1, (x_locs[1, ::2], y_locs[1, f::2]), 
                    force3[1, f::2, ::2], {
                    "title": "Force of a loosely bound 4 spring crossbridge",
                    "y_lab": "d$_{1,0}$ spacing (nm)",
                    "y_ticks": y_ticks,
                    "y_limits": y_limits,
                    "scale": force_scale
                    })
    f_c.quiver_plot(2, (x_locs[0, ::2], y_locs[0, f::2]), 
                    force3[0, f::2, ::2] - force3[1, f::2, ::2], {
                    "title": "Force difference between a 2 spring crossbridge \n and a 4 spring crossbridge",
                    "x_lab": "Binding site offset (nm)",
                    "y_lab": "d$_{1,0}$ spacing (nm)",
                    "y_ticks": y_ticks,
                    "y_limits": y_limits,
                    "scale": force_scale,
                    "quiv_key": {"x":0.85, "y":0.03, 
                                 "scale": 15, "label": "3 (pn)"}
                    })
    # Show plot
    f_c.show_plot()



if __name__ == '__main__':
    main()

