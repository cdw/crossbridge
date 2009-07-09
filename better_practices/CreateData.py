#!/usr/bin/env python
# encoding: utf-8
"""
CreateData.py
Created by Dave Williams on 2009-07-01.
"""

import Storage
import Crossbridge
import warnings
import sys
import getopt
import numpy as np
from numpy import pi


__help_message__ = '''
Don't eat crumbs in bed, also here are some options:
Option                Values         Results
-h, --help                           You see this message
-x, --crossbridge     2,4            Picks number of springs to emulate
-t, --trials          Interger       How many trials for each r12 point
-p, --property        Some strs      No value chooses all props, some value
                                       selects a given property
-d, --defaults        Exists or no   Chooses all default values       
'''


class Usage(Exception):
    """Passes mesages back to the command line"""
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    """Parse options and manage the creation of desired data"""
    if argv is None:
        argv = sys.argv
    try:
        try:
            short_opts = "ho:vx:t:p:d"
            long_opts = ["help" , "output=", "crossbridge=", 
                         "trials=", "property=", "defaults"]
            opts, args = getopt.getopt(argv[1:], short_opts, long_opts)
        except getopt.error, msg:
            raise Usage(msg)
        # Default values, retained for non-passed options
        xbtype = 4
        trials = 10
        prop_to_gen = None # Triggers generation of all properties
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            elif option in ("-h", "--help"):
                raise Usage(__help_message__)
            elif option in ("-o", "--output"):
                output = value
            elif option in ("-x", "--crossbridge"):
                if value in ("2", "4"):
                    xbtype = int(value)
                else:
                    raise Usage("Allowed xb types are 2 and 4 (spring)")
            elif option in ("-t", "--trials"):
                trials = int(value)
            elif option in ("-p", "--property"):
                prop_to_gen = value
            elif option in ("-d", "--defaults"):
                print("Using default values")
            else:
                raise Usage("Unhandled option")
        # Set ranges used
        x_range = [-5, 15, 1] 
        y_range = [5, 15, .5]
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
            xb = Crossbridge.FourSpring(config)
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
            xb = Crossbridge.TwoSpring(config)
        # Make or load a place to store results
        store = Storage.Storage(xbtype, config, x_range, y_range)
        # Generate some properties, or all of them
        if prop_to_gen is None:
            energy = calc_values(xb, x_range, y_range, 'energy', state=1)
            free_e = calc_values(xb, x_range, y_range, 'free_energy', state=2)
            r12 = calc_values(xb, x_range, y_range, 'r12', trials = trials)
            r23 = calc_values(xb, x_range, y_range, 'r23')
            r31 = calc_values(xb, x_range, y_range, 'r31')
            force1 = calc_values(xb, x_range, y_range, 'force', state=1)
            force2 = calc_values(xb, x_range, y_range, 'force', state=2)
            force3 = calc_values(xb, x_range, y_range, 'force', state=3)
            store.write('energy', energy)
            store.write('free_energy', free_e)
            store.write('r12', r12)
            store.write('r23', r23)
            store.write('r31', r31)
            store.write('force1', force1)
            store.write('force2', force2)
            store.write('force3', force3)
        else:
            prop_gen_func = {
                'energy': lambda:
                calc_values(xb, x_range, y_range, 'energy', state=1), 
                'free_energy': lambda:
                calc_values(xb, x_range, y_range, 'free_energy', state=2), 
                'r12': lambda:
                calc_values(xb, x_range, y_range, 'r12', trials = trials),
                'r23': lambda:
                calc_values(xb, x_range, y_range, 'r23'),
                'r31': lambda:
                calc_values(xb, x_range, y_range, 'r31'),
                'force1': lambda:
                calc_values(xb, x_range, y_range, 'force', state=1),
                'force2': lambda:
                calc_values(xb, x_range, y_range, 'force', state=2),
                'force3': lambda:
                calc_values(xb, x_range, y_range, 'force', state=3)
            }[prop_to_gen]
            prop_vals = prop_gen_func()
            store.write(prop_to_gen, prop_vals)
        # Save results to disk
        store.save()
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    

def calc_values(xb_inst, x_r, y_r, val_type, state =1, trials =1):
    """Calculate values and return them"""
    x_locs = np.arange(x_r[0], x_r[1], x_r[2]) 
    y_locs = np.arange(y_r[0], y_r[1], y_r[2])
    # Select the value to generate with a dict and some lambdas
    try:
        value_gen_func = {
            'energy': lambda x, y: xb_inst.minimize_energy((x,y), state)[0], 
            'free_energy': lambda x, y: xb_inst.free_energy((x,y), state), 
            'r12': lambda x, y: xb_inst.r12((x,y), trials),
            'r23': lambda x, y: xb_inst.r23((x,y)),
            'r31': lambda x, y: xb_inst.r31((x,y)),
            'force': lambda x, y: xb_inst.force((x,y), state)
        }[val_type]
    except KeyError:
        warnings.warn("Invalid value type, can't calculate request values")
        return
    # Create the new values in format [[row1], [row2], ...]
    new_vals = [[value_gen_func(x, y) for x in x_locs] for y in y_locs]
    return new_vals

if __name__ == "__main__":
    sys.exit(main())
