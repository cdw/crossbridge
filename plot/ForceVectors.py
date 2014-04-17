#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Force_Vectors.py

Created by Dave Williams on 2009-09-19.
Copyright (c) 2009 Dave Williams.
"""

import Storage
import numpy as np
import matplotlib.pyplot as plt

def main():
    """Plot vectors of the forces the 4sXB and 2sXB generate, followed by 
    contour plots of the axial and radial forces for the same.
    """
    ## Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)] 
    force3 = np.array([s.get("force3") for s in store])
    r23 = np.array([s.get("r23") for s in store])
    r31 = np.array([s.get("r31") for s in store])
    ## Print force properties
    print_force_properties(force3)
    print_clipped_force_properties(force3)
    chance = get_chance(r31, r23)
    ## Make chance all or nothing
    chance = np.array([[[1 if xb > -0.3 else 0 for xb in r] for r in c] for c in chance])
    ## Clip unlikely locations from forces for ratios
    # cliped_f = np.array([[[[v if abs(v)>force_limit else 1 for v in xb]
    #       for xb in row] for row in col] for col in force3[:, :, :, :]])
    ## Load and process x/y related values
    # Cross-bridges should have the same x and y ranges.
    assert(store[0].get("x_range")==store[1].get("x_range")) 
    assert(store[0].get("y_range")==store[1].get("y_range")) 
    x_range = store[0].get("x_range")
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    x_grid, y_grid = np.meshgrid(x_locs, y_locs)
    ## Set x/y locs of the big/sparse and detailed/magnified views
    sparse_locs = ((2, 4, 6, 8, 10, 12, 14, 16, 18), 
                   (26, 28, 30, 32, 34, 36, 38))
    magnified_locs = ((3.5, 4, 4.5, 5, 5.5, 6, 6.5), (32, 33, 34, 35, 36))
    ## Set x/y indices of x/y locs for sparse and magnified plots
    sparse_inds = [np.searchsorted(x_locs, sparse_locs[0]), 
                   np.searchsorted(y_locs, sparse_locs[1])]
    magnified_inds = [np.searchsorted(x_locs, magnified_locs[0]), 
                      np.searchsorted(y_locs, magnified_locs[1])]
    ## Create coordinate grids for both sparse and magnified plots
    sparse_grid = np.meshgrid(x_locs[sparse_inds[0]], y_locs[sparse_inds[1]])
    magnified_grid = np.meshgrid(x_locs[magnified_inds[0]], 
                                 y_locs[magnified_inds[1]]) 
    ## Set up the plotting 'paper'
    fig = plt.figure(1, figsize=(8, 10))
    axe = ([fig.add_subplot(4, 2, g+1) for g in range(4*2)])
    ## Force of the 4sXB
    axe[0].annotate("4sXB force vectors", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    axe[0].quiver(sparse_grid[0], sparse_grid[1], 
        force3[1, :, :, 0][sparse_inds[1], :][:, sparse_inds[0]], 
        force3[1, :, :, 1][sparse_inds[1], :][:, sparse_inds[0]],
        chance[1, :, :][sparse_inds[1], :][:, sparse_inds[0]],
        cmap=plt.get_cmap("Greys"))
    axe[0].set_xticks([0, 5, 10, 15, 20])
    axe[0].set_yticks([26, 30, 34, 38])
    ## Magnified force of the 4sXB
    axe[2].annotate("Small offset 4sXB force vectors", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center') 
    axe[2].quiver(magnified_grid[0], magnified_grid[1], 
            force3[1, :, :, 0][magnified_inds[1], :][:, magnified_inds[0]], 
            force3[1, :, :, 1][magnified_inds[1], :][:, magnified_inds[0]],
                  #chance[1, :, :][magnified_inds[1], :][:, magnified_inds[0]],
            cmap=plt.get_cmap("Greys"))
    axe[2].set_xticks(magnified_locs[0][1::2])
    axe[2].set_yticks(magnified_locs[1])
    ## Force of the 2sXB
    axe[1].annotate("2sXB force vectors", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    Sparse2sXB = axe[1].quiver(sparse_grid[0], sparse_grid[1], 
        force3[0, :, :, 0][sparse_inds[1], :][:, sparse_inds[0]], 
        force3[0, :, :, 1][sparse_inds[1], :][:, sparse_inds[0]],
        chance[0, :, :][sparse_inds[1], :][:, sparse_inds[0]],
        cmap=plt.get_cmap("Greys"))
    axe[1].set_xticks([0, 5, 10, 15, 20])
    axe[1].set_yticks([26, 30, 34, 38])
    axe[1].quiverkey(Sparse2sXB, 
                    1.12,           # x coord 
                    0.56,           # y coord
                    10,             # scale
                    "10 pn",        # label
                    labelpos='S',   # place text below the arrow
                    coordinates='axes')
    ## Magnified force of the 2sXB
    axe[3].annotate("Small offset 2sXB force vectors", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    Magnified2sXB = axe[3].quiver(magnified_grid[0], magnified_grid[1], 
            force3[0, :, :, 0][magnified_inds[1], :][:, magnified_inds[0]], 
            force3[0, :, :, 1][magnified_inds[1], :][:, magnified_inds[0]],
                                  #chance[0, :, :][magnified_inds[1], :][:, magnified_inds[0]],
            cmap=plt.get_cmap("Greys"))
    axe[3].set_xticks(magnified_locs[0][1::2])
    axe[3].set_yticks(magnified_locs[1])
    axe[3].quiverkey(Magnified2sXB, 
                    1.12,           # x coord 
                    0.56,           # y coord
                    2,              # scale
                    "2 pn",         # label
                    labelpos='S',   # place text below the arrow
                    coordinates='axes') 
    ## Plot other contours 
    grey_steps = ('.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8')
    colors = ('#4444BB','#333388','#2288FF','#77DD55',
              '#FFDD44','#FF7722','#AA4422', '#FF4466','#FF0000' )
    force_lvls = (-15, -10, -5, 0, 5, 10, 15, 20, 35)
    ## 4sXB axial force
    axe[4].annotate("4sXB axial force", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    c = axe[4].contourf(x_grid, y_grid, force3[1, :, :, 0], 
        force_lvls, colors = colors)
    ## 4sXB radial force
    axe[6].annotate("4sXB radial force", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    c = axe[6].contourf(x_grid, y_grid, force3[1,:,:,1], 
        force_lvls, colors = colors)
    ## 2sXB axial force
    axe[5].annotate("2sXB axial force", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    c = axe[5].contourf(x_grid, y_grid, force3[0, :, :, 0], 
       force_lvls, colors = colors)
    cb = plt.colorbar(c, ax=axe[5], shrink = 0.7, fraction=.16, format='% d',
                           ticks=force_lvls[1:-1:2])
    cb.ax.set_position([.92, .33, .07, .13])
    cb.set_label("axial force (pn)")
    ## 2sXB radial force
    axe[7].annotate("2sXB radial force", (0.5, 1.05), 
                    xycoords='axes fraction', horizontalalignment='center')
    c = axe[7].contourf(x_grid, y_grid, force3[0, :, :, 1], 
       force_lvls, colors = colors)
    cb = plt.colorbar(c, ax=axe[7], shrink = 0.7, fraction=.16, format='% d',
                        ticks=force_lvls[1:-1:2])
    cb.ax.set_position([.92, .09, .07, .13])
    cb.set_label("radial force (pn)")
    ## Fix the limits
    for i in range(4, 8):
        axe[i].set_xlim([x_locs[0], x_locs[-1]])
        axe[i].set_ylim([y_locs[0], y_locs[-1]])
        axe[i].set_yticks([26, 30, 34, 38])
        axe[i].set_xticks([0, 5, 10, 15, 20])
    ## Set the titles and axis labels
    titles = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"]
    for a in axe:
        a.set_xlabel("Axial offset (nm)")
        a.set_ylabel("$d_{1,0}$ (nm)")
        a.set_title(titles.pop(0), x=-0.15, y=0.99, weight="demi")
        a.xaxis.set_ticks_position('bottom')
        a.yaxis.set_ticks_position('left')
    ## Display the results
    fig.subplots_adjust(wspace=0.35, hspace=0.60,
                        left=0.10, right=0.91,
                        top=0.94, bottom=0.08)
    plt.show()

def print_force_properties(force3):
    """Print force max, min, etc for 4sXB and 2sXB"""
    print "## Force properties"
    print "4sXB max/min axial force: ", np.max(force3[1, :, :, 0]), "/", \
            np.min(force3[1, :, :, 0])
    print "4sXB max/min radial force: ", np.max(force3[1, :, :, 1]), "/", \
            np.min(force3[1, :, :, 1])
    print "2sXB max/min axial force: ", np.max(force3[0, :, :, 0]), "/", \
            np.min(force3[0, :, :, 0])
    print "2sXB max/min radial force: ", np.max(force3[0, :, :, 1]), "/", \
            np.min(force3[0, :, :, 1])

def print_clipped_force_properties(force3, force_limit = 0.5):
    """Clip low force regions for non-exploding ratios"""
    cliped_f = np.array([[[[v if abs(v)>force_limit else 1 for v in xb] 
                        for xb in row] for row in col] 
                        for col in force3[:, :, :, :]])
    print "4sXB max radial/axial ratio: ", \
        np.max(np.divide(cliped_f[1, :, :, 1], cliped_f[1, :, :, 0]))
    print "4sXB min radial/axial ratio: ", \
        np.min(np.divide(cliped_f[1, :, :, 1], cliped_f[1, :, :, 0]))    
    print "2sXB max radial/axial ratio: ", \
        np.max(np.divide(cliped_f[0, :, :, 1], cliped_f[0, :, :, 0]))
    print "2sXB min radial/axial ratio: ", \
        np.min(np.divide(cliped_f[0, :, :, 1], cliped_f[0, :, :, 0]))

def get_chance(r31, r23):
    """Calculate the likelihood of force being exerted at all locations"""
    #chance = np.clip( \
            #    [np.add(np.subtract(1, r31[i]), r23[i]) for i in [0,1]], 0, 1)
    chance = np.array([np.subtract(r23[i], r31[i]) for i in [0,1]])
    print "## Chance scales"
    print "4sXB chance max is ", np.max(chance[0]), 
    print " and min is ", np.min(chance[0])
    print "2sXB chance max is ", np.max(chance[1]), 
    print " and min is ", np.min(chance[1])
    return chance 

if __name__ == '__main__':
    main()

