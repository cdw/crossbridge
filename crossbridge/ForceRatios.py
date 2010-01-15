#!/usr/bin/env python
# encoding: utf-8
"""
ForceRatios.py

Created by Dave Williams on 2009-11-10.
Copyright (c) 2009 Dave Williams.
"""

import Storage
import numpy as np
import matplotlib.pyplot as plt

def main():
    ## Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)] 
    force3 = np.array([s.get("force3") for s in store])
    r23 = np.array([s.get("r23") for s in store])
    r31 = np.array([s.get("r31") for s in store])
    ## Load and process x/y related values
    assert(store[0].get("x_range")==store[1].get("x_range")) # XBs should have
    assert(store[0].get("y_range")==store[1].get("y_range")) # same x/y ranges
    x_range = store[0].get("x_range")
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    x_grid, y_grid = np.meshgrid(x_locs, y_locs)
    ## Clip low force regions for non-exploding ratios
    force_limit = 0.5
    cliped_f = np.array([[[[v if abs(v)>force_limit else 1 for v in xb] \
        for xb in row] for row in col] for col in force3[:, :, :, :]])
    ## Calculate the radial/axial ratios
    axrad = np.divide(cliped_f[:, :, :, 1], cliped_f[:, :, :, 0])
    ## Set up the plotting 'paper'
    fig = plt.figure(1, figsize=(8, 2.5))
    axe = ([fig.add_subplot(1, 2, g+1) for g in range(1*2)])
    ## Plot the 4sXB radial/axial ratio
    axe[0].annotate("4sXB |radial/axial| force ratio", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    c = axe[0].contourf(x_grid, y_grid, np.abs(axrad[0, :, :]), np.arange(0,2.25,.25))
    cb = plt.colorbar(c, ax=axe[0], shrink = 0.7, fraction=.16, format='%.1f')
    cb.ax.set_position([.42, .23, .13, .55])
    cb.set_label("radial/axial ratio")
    ## Plot the 2sXB radial/axial ratio
    axe[1].annotate("2sXB |radial/axial| force ratio", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    c = axe[1].contourf(x_grid, y_grid, np.abs(axrad[1, :, :]), np.arange(0,2.25,.25))
    cb = plt.colorbar(c, ax=axe[1], shrink = 0.7, fraction=.16, format='%.1f')
    cb.ax.set_position([.92, .23, .13, .55])
    cb.set_label("radial/axial ratio")
    ## Set the titles and axis labels
    titles = ["a)", "b)"]
    for a in axe:
        a.set_ylim([y_locs[0], y_locs[-1]])
        a.set_xlabel("Axial offset (nm)")
        a.set_ylabel("$d_{1,0}$ (nm)")
        a.set_title(titles.pop(0), x=-0.15, y=0.99, weight="demi")
        a.xaxis.set_ticks_position('bottom')
        a.yaxis.set_ticks_position('left')
    ## Display the results
    fig.subplots_adjust(wspace=0.60, hspace=0.60,
                        left=0.10, right=0.91,
                        top=0.80, bottom=0.22)
    plt.show()


if __name__ == '__main__':
    main()

