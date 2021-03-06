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
import time
import getopt
import numpy as np
from numpy import pi, radians


__help_message__ = '''
Don't eat crumbs in bed, also here are some options, pass some along:
Option                Values         Results
-h, --help                           You see this message
-x, --crossbridge     1,2,4          Picks number of springs to emulate
-t, --trials          Interger       How many trials for each r12 point
-p, --property        Some strs      No value chooses all props, some value
                                       selects a given property
-d, --defaults        Exists or not  Chooses all default values       
'''


class Usage(Exception):
    """Passes mesages back to the command line"""
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    """Parse options and manage the creation of desired data"""
    tic = time.time()
    if argv is None:
        argv = sys.argv
    try:
        try:
            short_opts = "hx:t:p:d"
            long_opts = ["help" , "crossbridge=", 
                         "trials=", "property=", "defaults"]
            opts, args = getopt.getopt(argv[1:], short_opts, long_opts)
        except getopt.error, msg:
            raise Usage(msg)
        # Default values, retained for non-passed options
        xbtype = 4
        trials = 10
        prop_to_gen = None # Triggers generation of all properties
        # option processing
        if len(opts) == 0:
            raise Usage(__help_message__)
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(__help_message__)
            elif option in ("-x", "--crossbridge"):
                if value in ("1", "2", "4"):
                    xbtype = int(value)
                else:
                    raise Usage("Allowed xb types are 1, 2 and 4 (spring)")
            elif option in ("-t", "--trials"):
                trials = int(value)
            elif option in ("-p", "--property"):
                prop_to_gen = value
            elif option in ("-d", "--defaults"):
                print("Using default values")
            else:
                raise Usage("Unhandled option")
        # Set ranges used
        x_range = [0, 20, 0.1] 
        y_range = [10, 20, 0.1]
        d10_range = [fil_sep_to_d10(y_range[0]), 
                    fil_sep_to_d10(y_range[1]), 1.5 * y_range[2]]
            # The y range is chosen based on values of SL length range for 
            # vertabrate cardiac muscle, actin thickness, myosin thickness, 
            # and an assumption of constant volume over those ranges the basic 
            # equation goes like this:
            #   filcenter_dist = face_dist + .5 * dia_actin + .5 * dia_myosin
            #   d10 = 1.5 * filcenter_dist
            #   vol = (2/np.sqrt(3)) * s_l * d10**2
            #   d10 = np.sqrt(np.sqrt(3)/2) * (vol/sl)
            # Millman 1998, pg375 says cardiac d10 at 2.2um is 37nm, and that 
            # cardiac SL range is 1.9-2.5um
            # So, plugging this in (and using the initial d10 at 2.2um to get 
            # the volume which remains constant thereafter) we get extreme 
            # values for the filcenter_dist of 8.8 nm and 5.4 nm. (These are 
            # extreme values of ~40nm and ~34.7nm if given in d10). Hence the 
            # choice of 5 and 10 for the limits of the y_range
        # Choose config based on xbtype
        if xbtype == 4:
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
                    'weak': 2*pi - radians(165), #pi/3 + (pi - radians(40)), #40deg from T weak
                    'strong': 2*pi - radians(110), # pi/2 + (pi - pi/4),
                    'spring_konstant': 40
                },
                'G': {
                    'weak': 9.6,
                    'strong': 9.6,
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
                     'weak': radians(47.16) + pi,
                     'strong': radians(73.20) + pi,
                     'spring_konstant': 40
                 },
                 'G': {
                     'weak': 19.93,
                     'strong': 16.47,
                     'spring_konstant': 2
                 }
            }
            xb = Crossbridge.TwoSpring(config)
        elif xbtype == 1:
            config = {
                 'T': {
                     'weak': 0,
                     'strong': 0,
                     'spring_konstant': 1
                 },
                 'N': {
                     'weak': 5,
                     'strong': 0,
                     'spring_konstant': 5
                 },
                 'C': {
                     'weak': pi,
                     'strong': pi,
                     'spring_konstant': 1
                 },
                 'G': {
                     'weak': 0,
                     'strong': 0,
                     'spring_konstant': 1
                 }
            }
            xb = Crossbridge.OneSpring(config)
        # Make or load a place to store results
        store = Storage.Storage(xbtype, config, x_range, d10_range)
        # Generate some properties, or all of them
        if prop_to_gen is None:
            runing_tic = time.time()
            print "Calculating energies... "
            energy = calc_values(xb, x_range, y_range, 'energy', state=1)
            free_e = calc_values(xb, x_range, y_range, 'free_energy', state=2)
            post_e = calc_values(xb, x_range, y_range, 'free_energy', state=3)
            print "done. Took " + str(time.time()-runing_tic) + " seconds."
            runing_tic = time.time()
            print "Calculating binding rate... "
            r12 = calc_values(xb, x_range, y_range, 'r12', trials = trials)
            print "done. Took " + str(time.time()-runing_tic) + " seconds."
            runing_tic = time.time()
            print "Calculating all other values... "
            r23 = calc_values(xb, x_range, y_range, 'r23')
            r31 = calc_values(xb, x_range, y_range, 'r31')
            force1 = calc_values(xb, x_range, y_range, 'force', state=1)
            force2 = calc_values(xb, x_range, y_range, 'force', state=2)
            force3 = calc_values(xb, x_range, y_range, 'force', state=3)
            print "done. Took " + str(time.time()-runing_tic) + " seconds."
            runing_tic = time.time()
            print "Storing output, will exit when done."
            store.write('energy', energy)
            store.write('free_energy', free_e)
            store.write('post_energy', post_e)
            store.write('r12', r12)
            store.write('trials', trials)
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
                'post_energy': lambda:
                calc_values(xb, x_range, y_range, 'free_energy', state=3),
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
        print "The whole thing took " + str(time.time()-tic) + " seconds."
        
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

def fil_sep_to_d10(face_to_face):
    """Convert filament seperation values from filament-face-to-filament face
    to d10 values that folks are used to seeing in x-ray diffraction studies
    """
    # # Filament paramter setup
    #    thick_dia = 26 # in nm, from Millman, 1998, pg378 
    #    # Note that 31 nm value is in Kensler & Harris, 2008, PMCID: PMC2242758
    #    thin_dia = 9.5 # in nm, from Millman, 1998 on pg378
    #    # The dist from fil center to fil center must consider these diameters
    #    cent_to_cent = face_to_face + 0.5 * thick_dia + 0.5 * thin_dia
    #    # From Millman 1998, pg362, we derive
    #    # d10 = 2*(np.cos(30*(np.pi/180)) * cent_to_cent) * np.cos(30*(np.pi/180))
    #    # or, since np.cos(30*(np.pi/180)) = np.sqrt(3/4)
    #    # d10 = 2*(np.sqrt(3/4) * cent_to_cent) * np.sqrt(3/4) or
    #    # d10 = 2* np.sqrt(3/4) * np.sqrt(3/4) * cent_to_cent  or
    #    # d10 = 3/2 * cent_to_cent
    
    # Based on values from ParamterEstimation.py
    cent_to_cent = face_to_face + 6.90
    return (1.5 * cent_to_cent)

if __name__ == "__main__":
    sys.exit(main())
