#!/usr/bin/env python
# encoding: utf-8
"""
CalculateStepSize.py
Created by Dave Williams on 2009-11-23.
"""

import Storage 
import numpy as np

store = [Storage.Storage(i) for i in (2,4)]
en2 = [s.get("free_energy") for s in store]
r23 = np.array([s.get("r23") for s in store]) 
en3 = [s.get("post_energy") for s in store]

# Load and process the x and y values
assert(store[0].get("x_range")==store[1].get("x_range"))
assert(store[0].get("y_range")==store[1].get("y_range"))
x_range = store[0].get("x_range")
x_locs = np.arange(x_range[0], x_range[1], x_range[2])
y_range = store[0].get("y_range")
y_locs = np.arange(y_range[0], y_range[1], y_range[2])
y_steps = [26, 30, 34, 38]
y_steps = [y_steps, [y_locs.searchsorted(y) for y in y_steps]]

# Calculate the minimum energy positions and step sizes
print "# Crossbridge Step Size" 
print " "
print "## Raw values"
print "### 4sXB min loosely bound energy loc"
x4e2 = [x_locs[en2[1][y].index(min(en2[1][y]))] for y in y_steps[1]]
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x4e2[i]
print "### 4sXB min strongly bound energy loc"
x4e3 = [x_locs[en3[1][y].index(min(en3[1][y]))] for y in y_steps[1]]
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x4e3[i]
print "### 4sXB powerstroke inflection loc"
i_lvl = (max(r23[1][0])-min(r23[1][0]))/2 + min(r23[1][0]) 
x4ps = [x_locs[(-r23[1][y]).searchsorted(-i_lvl)] for y in y_steps[1]]
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x4ps[i]
print "### 2sXB min loosely bound energy loc"
x2e2 = [x_locs[en2[0][y].index(min(en2[0][y]))] for y in y_steps[1]]
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x2e2[i]
print "### 2sXB min strongly bound energy loc"
x2e3 = [x_locs[en3[0][y].index(min(en3[0][y]))] for y in y_steps[1]]
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x2e3[i]
print "### 2sXB powerstroke inflection loc"
i_lvl = (max(r23[0][0])-min(r23[0][0]))/2 + min(r23[0][0]) 
x2ps = [x_locs[(-r23[0][y]).searchsorted(-i_lvl)] for y in y_steps[1]]
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x2ps[i]
print " "
print "## Calculated values"
print "### 4sXB energy step size"
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x4e2[i]-x4e3[i]
print "### 4sXB inflection step size" 
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x4ps[i]-x4e3[i]
print "### 2sXB energy step size"
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x2e2[i]-x2e3[i] 
print "### 2sXB inflection step size" 
for i in range(len(y_steps[0])):
    print "* @", y_steps[0][i], "nm, x is", x2ps[i]-x2e3[i]

