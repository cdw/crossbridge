#!/usr/bin/env python
# encoding: utf-8
"""
ParameterEstimation.py
Created by Dave Williams on 2009-09-09.
Copyright (c) 2009 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import numpy as np
from numpy import pi, cos, sin

def head_and_conv_locs(xb, state):
    xb = xb[state]
    conv_loc = (xb['neck_length'] * cos(xb['thic_angle']),
                xb['neck_length'] * sin(xb['thic_angle']))
    # need to calc the head loc, return both
    

def main():
    four = [{'state':1}, {'state':2}, {'state':3}]
    two = [{'state':1}, {'state':2}, {'state':3}]
    # From Davis and Epstein, 2009, pg6140
    four[1]['conv_angle'] = np.radians(125)
    four[2]['conv_angle'] = np.radians(125)
    four[3]['conv_angle'] = np.radians(70)
    
    # From Morel:1997
    # Rest lattice spacing with respect to radial force is 34nm d10
    # From Millman:1998
    # Rest lattice spacing with respect to compression is 37nm d10
    
    
    
    
    # Needed
    # four[1]['thic_angle']
    # four[1]['neck_length']
    # four[1]['glob_length']
    # four[1]['thic_spring'] # May be unattainable
    # four[1]['neck_spring'] # May be unattainable
    # four[1]['glob_spring'] # May be unattainable
    
    
    
    
    

if __name__ == '__main__':
    main()

