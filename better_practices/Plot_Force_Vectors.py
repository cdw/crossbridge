#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Force_Vectors.py

Created by Dave Williams on 2009-09-19.
Copyright (c) 2009 Dave Williams. All rights reserved.
"""

#import sys
#import os
import Storage
import numpy as np
import matplotlib.pyplot as plt
import pylab


def main():
    ## Load properties that will be needed
    print "Loading stored data...",
    store = [Storage.Storage(2), Storage.Storage(4)] 
    force3 = np.array([s.get("force3") for s in store])
    r23 = np.array([s.get("r23") for s in store])
    r31 = np.array([s.get("r31") for s in store])
    print "4sXB r31 max: ", np.max(r31[0])
    print "4sXB r31 min: ", np.min(r31[0])
    print "2sXB r31 max: ", np.max(r31[1])
    print "2sXB r31 min: ", np.min(r31[1])
    chance = np.clip([np.add(np.subtract(1, r31[i]), r23[i]) for i in [0,1]], 0, 1)
    print "4sXB chance max is ", np.max(chance[0]), " and min is ", np.min(chance[0])
    print "2sXB chance max is ", np.max(chance[1]), " and min is ", np.min(chance[1])
    #chance = np.array([np.divide(3, r31[i]) for i in [0,1]])
    #chance = np.array([np.add(s.get("r23"), np.divide(1, s.get("r31"))) for s in store])
    print "done."
    
    ## Load and process x/y related values
    assert(store[0].get("x_range")==store[1].get("x_range"))
    assert(store[0].get("y_range")==store[1].get("y_range"))
    x_range = store[0].get("x_range")
    x_locs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    x_grid, y_grid = np.meshgrid(x_locs, y_locs)
    # Set x/y locs of the large scale sparse and the detailed magnified plots
    sparse_locs = [[2, 4, 6, 8, 10, 12, 14, 16, 18], [28, 30, 32, 34, 36, 38]]
    magnified_locs = [[3.5, 4, 4.5, 5, 5.5, 6, 6.5], [32, 33, 34, 35, 36]]
    # Set x/y indices of the x/y locs for the sparse and magnified plots
    sparse_inds = [np.searchsorted(x_locs, sparse_locs[0]), 
                   np.searchsorted(y_locs, sparse_locs[1])]
    magnified_inds = [np.searchsorted(x_locs, magnified_locs[0]), 
                      np.searchsorted(y_locs, magnified_locs[1])]
    # Create coordinate grids for both sparse and magnified plots
    sparse_grid = np.meshgrid(x_locs[sparse_inds[0]], y_locs[sparse_inds[1]])
    magnified_grid = np.meshgrid(x_locs[magnified_inds[0]], 
                                 y_locs[magnified_inds[1]]) 
    
    ## Set up the plotting 'paper'
    fig = plt.figure(1, figsize=(8, 10))
    axe = ([fig.add_subplot(4, 2, g+1) for g in range(4*2)])
    
    ## Force of the 4sXB
    Sparse4sXB = axe[0].quiver(sparse_grid[0], sparse_grid[1], 
        force3[1,:,:,0][sparse_inds[1], :][:, sparse_inds[0]], 
        force3[1,:,:,1][sparse_inds[1], :][:, sparse_inds[0]],
        chance[1,:,:][sparse_inds[1], :][:, sparse_inds[0]],
        cmap=plt.get_cmap("Greys"))
    axe[0].set_xticks([0, 5, 10, 15, 20])
    axe[0].quiverkey(Sparse4sXB, 
                    0.86,           # x coord 
                    1.05,           # y coord
                    10,             # scale
                    "10pn",         # label
                    labelpos='E',   # place text right of arrow
                    coordinates='axes') 
     
    ## Magnified force of the 4sXB
    Magnified4sXB = axe[2].quiver(magnified_grid[0], magnified_grid[1], 
            force3[1,:,:,0][magnified_inds[1], :][:, magnified_inds[0]], 
            force3[1,:,:,1][magnified_inds[1], :][:, magnified_inds[0]],
            chance[1,:,:][magnified_inds[1], :][:, magnified_inds[0]],
            cmap=plt.get_cmap("Greys"))
    axe[2].set_xticks(magnified_locs[0][1::2])
    axe[2].set_yticks(magnified_locs[1])
    axe[2].quiverkey(Magnified4sXB, 
                    0.89,           # x coord 
                    1.05,           # y coord
                    2,              # scale
                    "2pn",          # label
                    labelpos='E',   # place text right of arrow
                    coordinates='axes') 
    
    ## Force of the 2sXB
    Sparse2sXB = axe[1].quiver(sparse_grid[0], sparse_grid[1], 
        force3[0,:,:,0][sparse_inds[1], :][:, sparse_inds[0]], 
        force3[0,:,:,1][sparse_inds[1], :][:, sparse_inds[0]],
        chance[0,:,:][sparse_inds[1], :][:, sparse_inds[0]],
        cmap=plt.get_cmap("Greys"))
    #axe[1].set_xticks(sparse_locs[0][1::2])
    axe[1].set_xticks([0, 5, 10, 15, 20])
       
    ## Magnified force of the 2sXB
    Magnified2sXB = axe[3].quiver(magnified_grid[0], magnified_grid[1], 
            force3[0,:,:,0][magnified_inds[1], :][:, magnified_inds[0]], 
            force3[0,:,:,1][magnified_inds[1], :][:, magnified_inds[0]],
            chance[0,:,:][magnified_inds[1], :][:, magnified_inds[0]],
            cmap=plt.get_cmap("Greys"))
    axe[3].set_xticks(magnified_locs[0][1::2])
    axe[3].set_yticks(magnified_locs[1])
    
    ## Plot other contours 
    grey_steps = ('.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8')
    force_lvls = (-15, -10, -5, 0, 5, 10, 15, 20, 25)
    
    # 4sXB axial force
    plt.subplot(4, 2, 5)
    c = axe[4].contourf(x_grid, y_grid, force3[1,:,:,0], 
       force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[4], shrink = 0.7, fraction=.16, format='% d')
    cb.ax.set_position([.42, .33, .07, .13])
    cb.set_label("4sXB axial force")
    
    # 4sXB radial force
    plt.subplot(4, 2, 6)
    c = axe[6].contourf(x_grid, y_grid, force3[1,:,:,1], 
       force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[6], shrink = 0.7, fraction=.16, format='% d')
    cb.ax.set_position([.42, .09, .07, .13])
    cb.set_label("4sXB radial force")
    
    # 2sXB axial force
    plt.subplot(4, 2, 7)
    c = axe[5].contourf(x_grid, y_grid, force3[0,:,:,0], 
       force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[5], shrink = 0.7, fraction=.16, format='% d')
    cb.ax.set_position([.92, .33, .07, .13])
    cb.set_label("2sXB axial force")
    
    # 2sXB radial force
    plt.subplot(4, 2, 8)
    c = axe[7].contourf(x_grid, y_grid, force3[0,:,:,1], 
       force_lvls, colors = grey_steps)
    cb = plt.colorbar(c, ax=axe[7], shrink = 0.7, fraction=.16, format='% d')
    cb.ax.set_position([.92, .09, .07, .13])
    cb.set_label("2sXB radial force")
    for i in range(4,8):
        axe[i].set_xlim([x_locs[0], x_locs[-1]])
        axe[i].set_ylim([y_locs[0], y_locs[-1]])
    
    ## Fix and annotate
    titles = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"]
    #titles = ["4sXB", "Section from 4sXB", "2sXB", "Section from 2sXB"]
    for a in axe:
        a.set_xlabel("Binding site offset (nm)")
        a.set_ylabel("$d_{1,0}$ (nm)")
        a.set_title(titles.pop(0), x=-0.15, y=0.99, weight="demi")
    
    ## Display the results
    fig.subplots_adjust(wspace=0.58, hspace=0.60,
                        left=0.10, right=0.91,
                        top=0.94, bottom=0.08)
    plt.show()


if __name__ == '__main__':
    main()

