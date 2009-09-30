#!/usr/bin/env python
# encoding: utf-8
"""
NodeBoxDataGen.py

Created by Dave Williams on 2009-09-29.
Copyright (c) 2009 Dave Williams. All rights reserved.
"""

import cPickle as pickle
import numpy as np
from numpy import pi, radians
import Crossbridge

def main():
    config = {
        'T': {
            'weak': radians(40),
            'strong': radians(40),
            'spring_konstant': 100
        },
        'N': {
            'weak': 10.5,
            'strong': 10.5,
            'spring_konstant': 10
        },
        'C': {
            'weak': 2*pi - radians(165), # 40deg from T weak
            'strong':2*pi - radians(110),
            'spring_konstant': 100
        },
        'G': {
            'weak': 9.6,
            'strong': 9.6,
            'spring_konstant': 5
        }
    }
    # config = {
    #     'T': {
    #         'weak': pi/4,
    #         'strong': pi/4,
    #         'spring_konstant': 100
    #     },
    #     'N': {
    #         'weak': 5,
    #         'strong': 5,
    #         'spring_konstant': 10
    #     },
    #     'C': {
    #         'weak': pi/3 + (pi - pi/4), #40deg from T weak
    #         'strong': pi/2 + (pi - pi/4),
    #         'spring_konstant': 100
    #     },
    #     'G': {
    #         'weak': 3,
    #         'strong': 3,
    #         'spring_konstant': 5
    #     }
    # }
    xb = Crossbridge.FourSpring(config)
    x_range = [0, 30, .2] 
    y_range = [0, 30, .2]
    x_locs = np.arange(x_range[0], x_range[1], x_range[2]) 
    y_locs = np.arange(y_range[0], y_range[1], y_range[2])
    # c_loc and h_loc assignment is confusing. The indices go like this: 
    # (a, b, c, d)
    # Where a is state, 0 for first two, 1 for strongly bound
    #       b is the row, corresponding to the y location in y_locs
    #       c is the col, corresponding to the x location in x_locs
    #       d is the x/y output value selector, 0 for x, 1 for y
    c_locs = np.zeros((2, y_locs.size, x_locs.size, 2))
    h_locs = np.zeros((2, y_locs.size, x_locs.size, 2))
    energy = np.zeros((2, y_locs.size, x_locs.size))
    springs = 4
    
    xit, yit, sit = 0, 0, 0
    for state in [2]:
        for y in y_locs:
            for x in x_locs:
                h_locs[sit, yit, xit] = [x, y]
                energy[sit, yit, xit], c_locs[sit, yit, xit] = \
                                xb.minimize_energy((x, y), state)
                xit += 1
            yit += 1
            xit = 0
        sit += 1
        yit = 0
        xit = 0
    
    
    file_name = "NodeBoxData.pkl"
    stream = open(file_name, 'w')
    #p_out = pickle.Pickler(stream) #can add ,1 after stream for smaller files
    #p_out.dump(x_range)
    #p_out.dump(y_range)
    #p_out.dump(h_locs)
    # p_out.dump(c_locs)
    # p_out.dump(energy)
    # p_out.dump(config)
    # p_out.dump(springs)
    pickle.dump({
        "x_range": x_range, 
        "y_range": y_range, 
        "h_locs": h_locs, 
        "c_locs": c_locs, 
        "energy": energy, 
        "config": config, 
        "springs": springs
        }, stream)
    stream.close()
    return


if __name__ == '__main__':
    main()

