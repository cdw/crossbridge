#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Energy_Minima.py
Plot the step size of the multi-spring cross-bridges.
Created by Dave Williams on 2010-09-23.
"""

import Storage
import numpy as np
from scipy import signal as sig
import matplotlib.pyplot as plt

def main():
    """Will calculate and plot step size versus lattice spacing"""
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)] 
    pre_energy = [s.get("free_energy") for s in store]
    post_energy = [s.get("post_energy") for s in store]
    x_range = store[0].get("x_range")
    xlocs = np.arange(x_range[0], x_range[1], x_range[2])
    y_range = store[0].get("y_range")
    ylocs = np.arange(y_range[0], y_range[1], y_range[2])
    # Calculate energy minimums
    xb2mins = min_energies(pre_energy[0], post_energy[0], xlocs) 
    xb4mins = min_energies(pre_energy[1], post_energy[1], xlocs) 
    # Set up the figure
    fig = plt.figure(1, figsize=(7.5,2.5)) 
    axe = (fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2))
    # Plot the results
    axe[0].plot(ylocs, xb4mins[0], color='#FF466F', lw=4)
    axe[0].plot(ylocs, xb4mins[1], color='#FF466F', lw=4)
    axe[1].plot(ylocs, xb2mins[0], color='#76D753', lw=4)
    axe[1].plot(ylocs, xb2mins[1], color='#76D753', lw=4)
    # Annotate the plots
    axe[0].set_title("4sXB Energy Minima")
    axe[0].set_xlabel("Lattice spacing (nm)") 
    axe[0].set_ylabel("Minimum location (nm)")
    axe[0].set_xlim((25.5, 39))
    axe[1].set_title("2sXB Energy Minima")
    axe[1].set_xlabel("Lattice spacing (nm)") 
    axe[1].set_ylabel("Minimum location (nm)")
    axe[1].set_xlim((25.5, 39))
    # Y Axis Limists and Ticks
    def set_y_axis(axis, tics):
        axis.set_ylim(tics[0], tics[-1])
        axis.set_yticks(tics)
        axis.set_yticklabels(tics, size=8)
        axis.yaxis.set_ticks_position('left')
    set_y_axis(axe[0], range(1,9))
    set_y_axis(axe[1], range(1,9))
    # X Axis Ticks
    # TODO: Left off here!
    for axis in axe:
        tics = range(0,21,5)
        axis.set_xlim(tics[0], tics[-1])
        axis.set_xticks(tics)
        axis.set_xticklabels(tics, size=8)
        axis.xaxis.set_ticks_position('bottom')
    set
    # Display the plots
    fig.subplots_adjust(wspace=0.25, hspace=0.48,
                        left=0.08, right=0.98,
                        top=0.85, bottom=0.21)
    plt.show()

def min_energies(pre, post, xlocs):
    """Give back the locations of all the minimum energies"""
    pre_inds = [np.argmin(row) for row in pre]
    post_inds = [lmin(row, start) for row, start in zip(post,pre_inds)]
    #post_inds = [np.argmin(row) for row in post]
    pre_loc = [xlocs[ind] for ind in pre_inds]
    post_loc =[xlocs[ind] for ind in post_inds] 
    return (smooth(pre_loc), smooth(post_loc))

def smooth(x, sam=10, type='gaussian'):
    """Smooth the landscape, x, with a window of given length and type.
    x : 1d array or list in need of smoothing
    sam : number of samples to smooth over, default is 10
    type : type of smoothing window, currently gaussian or flat
    """
    if type is 'gaussian':
        win = sig.gaussian(sam, sam*.3) #the std is a shot in the dark 
    elif type is 'flat':
        win = np.ones(sam)
    else:
        raise ValueError, "type name not recognized"
    scape = np.r_[2*x[0]-x[sam:1:-1],x,2*x[-1]-x[-1:-sam:-1]] 
    smoothed = sig.convolve(np.divide(win, win.sum()), scape, mode='same')
    return smoothed[sam-1:-(sam-1)]

def lmin(scape, start):
    """Return the local minimum's index given a stating point.
    Note that this is not robust for noisy data and will always 
    move down the indices if placed at a peak. The extra numerical
    factors require a minimum slope.
    """
    i = start
    while scape[i - 1] < scape[i] - 0.06:
        i -= 1
    while scape[i + 1] < scape[i] - 0.06:
        i += 1
    return i

if __name__ == '__main__':
    main()

