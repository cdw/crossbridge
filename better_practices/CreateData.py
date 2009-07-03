#!/usr/bin/env python
# encoding: utf-8
"""
CreateData.py

Created by Dave Williams on 2009-07-01.
Copyright (c) 2009 __MyCompanyName__. All rights reserved.
"""

import Storage
import Crossbridge
import warnings
import sys
import getopt
import numpy as np
from numpy import pi


help_message = '''
Don't eat crumbs in bed
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            elif option in ("-h", "--help"):
                raise Usage(help_message)
            elif option in ("-o", "--output"):
                output = value
            elif option in ("-x", "--crossbridge"):
                if value in (2, 4):
                    xbtype = value
                else:
                    raise Usage("Allowed xb types are 2 and 4 (spring)")
            elif option in ("-t", "--trials"):
                trials = value
            elif option in ("-d", "--defaults"):
                print("Using default values")
                xbtype = 4
                trials = 10
            else:
                raise Usage("Unhandled option")
        # If something wasn't passed, fill in a default
        xbtype = xbtype if vars().has_key["xbtype"] else 4
        trials = trials if vars().has_key["trials"] else 10
        # Set ranges used
        x_range = (-5, 15, 10) 
        y_range = (5, 15, 2.5)
        # Choose config based on xbtype
        if xbtype == 4:
            config = {
                'T': {
                    'weak': pi/4,
                    'strong': pi/4,
                    'spring_konstant': 100
                },
                'N': {
                    'weak': 7,
                    'strong': 5,
                    'spring_konstant': 10
                },
                'C': {
                    'weak': pi/3 + (pi - pi/4), #pi/4 from T weak
                    'strong': pi/2 + (pi - pi/4),
                    'spring_konstant': 100
                },
                'G': {
                    'weak': 3,
                    'strong': 3,
                    'spring_konstant': 5
                }
            }
        elif xbtype == 2:
            config = {
                 'T': {
                     'weak': 0,
                     'strong': 0,
                     'spring_konstant': 1
                 },
                 'N': {
                     'weak': 0,
                     'strong': 0,
                     'spring_konstant': 1
                 },
                 'C': {
                     'weak': pi/3 + pi,
                     'strong': pi/2 + pi,
                     'spring_konstant': 100
                 },
                 'G': {
                     'weak': 11,
                     'strong': 11,
                     'spring_konstant': 5
                 }
             }
        
        store = Storage.Storage(xbtype, config, x_range, y_range)
        
        xb = Crossbridge.FourSpring(config)
        energies = calc_values(store, xb, x_range, y_range, 'energy', state=1)
        store.write('energy', energies)
        store.save()
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    

def calc_values(s_inst, xb_inst, x_r, y_r, val_type, state =1, trials =1):
    """Calculate values and write parameters to storage instance"""
    x_locs = np.arange(x_r[0], x_r[1], x_r[2]) 
    y_locs = np.arange(y_r[0], y_r[1], y_r[2])
    # Select the value to generate with a dict and some lambdas
    try:
        value_gen_func = {
            'energy': lambda x, y: xb_inst.minimize_energy((x,y), state)[0], 
            'r12': lambda x, y: xb_inst.r12((x,y), trials),
            'r23': lambda x, y: xb_inst.r23((x,y)),
            'r31': lambda x, y: xb_inst.r31((x,y)),
            'force': lambda x, y: xb_inst.force((x,y), state)
        }[val_type]
    except KeyError:
        warnings.warn("Invalid value type, can't calculate request values")
        return
    # Create the new values in format [[row1], [row2], ...]
    new_vals = [[float(value_gen_func(x, y)) for x in x_locs] for y in y_locs]
    return new_vals

if __name__ == "__main__":
    sys.exit(main())
