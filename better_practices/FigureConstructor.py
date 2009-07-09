#!/usr/bin/env python
# encoding: utf-8
"""
FigureSupport.py
Created by Dave Williams on 2009-07-06.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#from pylab import figure, show, colorbar

class FigureConstructor:
    """docstring for FigureConstructor"""
    def __init__(self, rowcol):
        # Use Helvetica as default font
        matplotlib.rcParams['font.family'] = 'sans-serif'
        matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
        # Allow the use of latex
        matplotlib.rcParams['text.usetex'] = True
        # Make negative lines in contours not dashed
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        # Institate figure and subplots
        self.fig = plt.figure(1)
        self.axe = ([self.fig.add_subplot(
                        rowcol[0], 
                        rowcol[1],
                        row * rowcol[1] + col + 1) 
                    for row in range(rowcol[0])
                    for col in range(rowcol[1])])
        # Configure colors
        self.colors = {
            "grey_steps": ('.2', '.3', '.4', '.5', '.6', '.7', '.8'),
            "white": ('white')
            }
        
    
    def two_contour(self, sub, x_y_values, z_grid, levels, labels={}):
        """Plot a grayscale contour of the passed z_vals with lines at the 
        given levels. Let the graph be bounded by the passed x_y_values. 
        Mark up the plot with the passed labels, if any. Labels may include:
        "title", "x_lab", and "y_lab".
        """ 
        
        # Create the grid of x and y vals
        x_grid, y_grid = np.meshgrid(x_y_values[0], x_y_values[1])
        # Plot the solid color filled contours
        self.axe[sub].contourf(x_grid, y_grid, z_grid, levels, 
            colors = self.colors["grey_steps"])
        # Plot the highlighting white lines
        self.axe[sub].contour(x_grid, y_grid, z_grid, levels, 
            colors = self.colors["white"])
        # Fix the limits
        self.axe[sub].set_xlim([x_y_values[0][0], x_y_values[0][-1]])
        self.axe[sub].set_ylim([x_y_values[1][0], x_y_values[1][-1]])
        # Apply labels
        if labels.has_key("title"):
            self.axe[sub].set_title(labels["title"])
        if labels.has_key("x_lab"):
            self.axe[sub].set_xlabel(labels["x_lab"]) 
        if labels.has_key("y_lab"):
            self.axe[sub].set_ylabel(labels["y_lab"]) 
    
    def cut_plot(self, sub, x_y_values, z_grid, cut_locs, labels_n_limits={}):
        """Plot a single, or series of, cuts"""
        # Find the indices at which to take the cuts
        cut_indices = np.searchsorted(x_y_values[1], cut_locs)
        # Plot the cuts
        for cut in cut_indices:
            self.axe[sub].plot(x_y_values[0], z_grid[cut], color='0.7')
        # Fix the limits
        #self.axe[sub].set_yticks([x_y_values[1][0],
        #    .5*(x_y_values[1][-1]-x_y_values[1][0])
        #    ,x_y_values[1][-1]])
        # Apply labels n' limits
        if labels_n_limits.has_key("title"):
            self.axe[sub].set_title(labels_n_limits["title"])
        if labels_n_limits.has_key("x_lab"):
            self.axe[sub].set_xlabel(labels_n_limits["x_lab"]) 
        if labels_n_limits.has_key("y_lab"):
            self.axe[sub].set_ylabel(labels_n_limits["y_lab"]) 
        if labels_n_limits.has_key("y_limits"):
            self.axe[sub].set_ylim(labels_n_limits["y_limits"])
        if labels_n_limits.has_key("y_ticks"):
            self.axe[sub].set_yticks(labels_n_limits["y_ticks"])
    
    def quiver_plot(self, sub, x_y_values, j_k_grid, labels_n_limits={}):
        """Plot an arrows graphy"""
        self.axe[sub].quiver(x_y_values[0], x_y_values[1], 
                             j_k_grid[:,:,0], j_k_grid[:,:,1])
        if labels_n_limits.has_key("title"):
            self.axe[sub].set_title(labels_n_limits["title"])
        if labels_n_limits.has_key("x_lab"):
            self.axe[sub].set_xlabel(labels_n_limits["x_lab"]) 
        if labels_n_limits.has_key("y_lab"):
            self.axe[sub].set_ylabel(labels_n_limits["y_lab"]) 
        if labels_n_limits.has_key("y_limits"):
            self.axe[sub].set_ylim(labels_n_limits["y_limits"])
        if labels_n_limits.has_key("y_ticks"):
            self.axe[sub].set_yticks(labels_n_limits["y_ticks"])
    
    def save_plot(self, location="./image.pdf", trans=True):
        """Save the plot to the specified location and type"""
        self.fig.subplots_adjust(wspace=0.25, hspace=0.50,
                            left=0.10, right=0.95,
                            top=0.94, bottom=0.08)
        self.fig.savefig(location, transparent=trans)
    
    def show_plot(self):
        """Show the plot in an interactive window"""
        self.fig.subplots_adjust(wspace=0.25, hspace=0.50,
                            left=0.10, right=0.95,
                            top=0.94, bottom=0.08)
        plt.show()


# setup_figure(rows, columns)
# two_contour(subplotloc, vals=[x,y,z], levels, cut_loc, labels)
# cut_plot(subplotloc, vals=[x,y,z], cut_loc, labels)
# quiver_plot(subplotloc, vals=[x,y,[fx,fy]], labels)
# finish(outputtype, outputfile)