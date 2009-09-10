#!/usr/bin/env python
# encoding: utf-8
"""
FigureSupport.py
Created by Dave Williams on 2009-07-06.
"""

import numpy as np
import matplotlib
# Choose backend
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
#from pylab import figure, show, colorbar
#import pdb

class FigureConstructor:
    """docstring for FigureConstructor"""
    def __init__(self, rowcol, fig_size=(8, 8)):
        # Use Helvetica as default font
        matplotlib.rcParams['font.family'] = 'sans-serif'
        matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
        # Allow the use of latex
        matplotlib.rcParams['text.usetex'] = True
        # Make negative lines in contours not dashed
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        # Institate figure and subplots
        self.fig = plt.figure(1)
        self.fig.set_size_inches(fig_size) #, forward=True) only on mac?
        self.axe = ([self.fig.add_subplot(
                        rowcol[0], 
                        rowcol[1],
                        subgraph + 1) 
                    for subgraph in range(rowcol[0]*rowcol[1])])
        # Configure colors
        self.colors = {
            "grey_steps": ('.2', '.3', '.4', '.5', '.6', '.7', '.8'),
            "white": ('white'),
            "cuts": ('#1E90FF', '#22FF22', '#FF4444', 
                     '#13E0E0', '#FFFF66', '#DB73FF')
            }
        # self.line_styles = ["--", "-.", ":", "."]
        self.line_styles = ["solid", "solid", "solid", "solid", "solid", 
                            "solid", "solid", "solid"] #try out solids
        
    
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
        # Deal with any passed cuts
        if labels.has_key("cuts"):
            curr_cut = 0
            for cut in labels["cuts"]["cut_locs"]:
                self.axe[sub].plot([x_y_values[0][0], x_y_values[0][-1]], 
                    [cut, cut], ls = self.line_styles[curr_cut], 
                    color=self.colors["cuts"][curr_cut], linewidth=1.5)
                curr_cut += 1
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
        # Create line differentiations
        curr_line = 0
        # Plot the cuts
        for cut in cut_indices:
            self.axe[sub].plot(x_y_values[0], z_grid[cut], 
                color=self.colors["cuts"][curr_line], 
                ls = self.line_styles[curr_line], lw=2.0)
            curr_line += 1
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
    
    def quick_plot(self, sub, x_y_values):
        """Add a simple plot to the specified subfigure"""
        self.axe[sub].plot(x_y_values[0], x_y_values[1], color='0.0', lw=2.0)
    
    def quiver_plot(self, sub, x_y_values, j_k_grid, labels_n_limits={}):
        """Plot an arrows graphy"""
        # Get arrow scale if it is passed
        if labels_n_limits.has_key("scale"):
            k_words = {"scale": labels_n_limits["scale"]}
        # Manually mesh x and y into grid, see http://tinyurl.com/km8ut4
        x_mesh, y_mesh = np.meshgrid(x_y_values[0], x_y_values[1])
        # Plot the passed quivers
        quiv = self.axe[sub].quiver(x_mesh, y_mesh, 
            j_k_grid[:,:,0], j_k_grid[:,:,1], **k_words)
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
        if labels_n_limits.has_key("quiv_key"):
                self.axe[sub].quiverkey(quiv, 
                    labels_n_limits["quiv_key"]["x"], 
                    labels_n_limits["quiv_key"]["y"],  
                    labels_n_limits["quiv_key"]["scale"], 
                    labels_n_limits["quiv_key"]["label"],
                        labelpos='E',
                        coordinates='figure',
                        fontproperties={'weight': 'bold'})
    
    def quiver_cut(self, sub, x_y_values, j_k_grid, cut_component, cut_locs, labels_n_limits={}):
        '''Hack out the component of the j_k grid that you want to have 
        a cut of (0 for j, 1 for k) and then pass on the remainder to 
        cut_plot
        '''
        # j_k_grid -> z_grid
        # can this work by just indexing as done two lines below?
        # Pass that buck
        self.cut_plot(sub, x_y_values, j_k_grid[:, :, cut_component],
            cut_locs, labels_n_limits)
    
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

## Example Useage
# setup_figure(rows, columns)
# two_contour(subplotloc, vals=[x,y,z], levels, cut_loc, labels)
# cut_plot(subplotloc, vals=[x,y,z], cut_loc, labels)
# quiver_plot(subplotloc, vals=[x,y,[fx,fy]], labels)
# finish(outputtype, outputfile)
    
