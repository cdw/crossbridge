#!/usr/bin/env python
# encoding: utf-8
"""
ControlXB.py

Created by Dave Williams
"""

import numpy as np
import time
import Crossbridge
import Graph

def create_forces(xb_type, xb_state=0, limits=(3,.5,15,3,.5,15)):
    x_locs = np.arange(limits[0], limits[2], limits[1]) 
    y_locs = np.arange(limits[3], limits[5], limits[4])
    forces = np.zeros((y_locs.size, x_locs.size, 2))
    # Instantiate the xb, dict provides a switch between XB types
    try:
        xb = {'TNCG':Crossbridge.TNCG(), 'xxCG':Crossbridge.xxCG(),
              'xNxx':Crossbridge.xNxx()}[xb_type]
    except KeyError:
        print("Invalid XB type, available types are TNCG, xxCG and xNxx")
        return
    #xb = Crossbridge.xNCx()
    xb.set_state(xb_state)
    print(xb)
    pT = time.time()
    cT = time.time()
    # Cycle through and collect all the probabilities
    for n,x in enumerate(x_locs):
        for m,y in enumerate(y_locs):
            forces[m, n] = xb.force((x,y),xb_state)
        # Tell me how much time is left, about
        cT = time.time()
        rT = (cT-pT)*(x_locs.size - (n+1))
        print('On col %(c)04d of %(t)04d, about %(m)02d:%(s)02d left' \
              %{'c':n, 't':x_locs.size, 'm':rT//60, 's':rT%60})
        pT = cT
    return forces

def create_energies(xb_type, xb_state=0, limits=(-2,.5,15,-2,.5,15)):
    x_locs = np.arange(limits[0], limits[2], limits[1]) 
    y_locs = np.arange(limits[3], limits[5], limits[4])
    energs = np.zeros((y_locs.size, x_locs.size))
    # Instantiate the xb, dict provides a switch between XB types
    try:
        xb = {'TNCG':Crossbridge.TNCG(), 'xxCG':Crossbridge.xxCG(),
              'xNxx':Crossbridge.xNxx()}[xb_type]
    except KeyError:
        print("Invalid XB type, available types are TNCG, xxCG and xNxx")
        return
    #xb = Crossbridge.xNCx()
    xb.set_state(xb_state)
    print(xb)
    pT = time.time()
    cT = time.time()
    # Cycle through and collect all the probabilities
    for n,x in enumerate(x_locs):
        for m,y in enumerate(y_locs):
            energs[m, n] = xb.minimize_energy(h_loc=(x,y))
        # Tell me how much time is left, about
        cT = time.time()
        rT = (cT-pT)*(x_locs.size - (n+1))
        print('On col %(c)04d of %(t)04d, about %(m)02d:%(s)02d left' \
              %{'c':n, 't':x_locs.size, 'm':rT//60, 's':rT%60})
        pT = cT
    return energs


def create_free_energies(xb_type, xb_state=0, limits=(-2,.5,15,-2,.5,15)):
    x_locs = np.arange(limits[0], limits[2], limits[1]) 
    y_locs = np.arange(limits[3], limits[5], limits[4])
    energs = np.zeros((y_locs.size, x_locs.size))
    # Instantiate the xb, dict provides a switch between XB types
    try:
        xb = {'TNCG':Crossbridge.TNCG(), 'xxCG':Crossbridge.xxCG(),
              'xNxx':Crossbridge.xNxx()}[xb_type]
    except KeyError:
        print("Invalid XB type, available types are TNCG, xxCG and xNxx")
        return
    #xb = Crossbridge.xNCx()
    xb.set_state(xb_state)
    print(xb)
    pT = time.time()
    cT = time.time()
    # Cycle through and collect all the probabilities
    for n,x in enumerate(x_locs):
        for m,y in enumerate(y_locs):
            energs[m, n] = xb.free_energy(xb_state,h_loc=(x,y))
        # Tell me how much time is left, about
        cT = time.time()
        rT = (cT-pT)*(x_locs.size - (n+1))
        print('On col %(c)04d of %(t)04d, about %(m)02d:%(s)02d left' \
              %{'c':n, 't':x_locs.size, 'm':rT//60, 's':rT%60})
        pT = cT
    return energs


def tran01(xb_type, limits=(-2,.5,15,-2,.5,15), trials=1000):
    x_locs = np.arange(limits[0], limits[2], limits[1]) 
    y_locs = np.arange(limits[3], limits[5], limits[4])
    trans01 = np.zeros((y_locs.size, x_locs.size))
    # Instantiate the xb, dict provides a switch between XB types
    try:
        xb = {'TNCG':Crossbridge.TNCG(), 'xxCG':Crossbridge.xxCG(),
              'xNxx':Crossbridge.xNxx()}[xb_type]
    except KeyError:
        print("Invalid XB type, available types are TNCG, xxCG and xNxx")
        return
    #xb = Crossbridge.xNCx()
    print(xb)
    pT = time.time()
    cT = time.time()
    # Cycle through and collect all the probabilities
    for n,x in enumerate(x_locs):
        for m,y in enumerate(y_locs):
            for t in range(trials):
                trans01[m, n] += xb.tran01((x,y))
        # Tell me how much time is left, about
        cT = time.time()
        rT = (cT-pT)*(x_locs.size - (n+1))
        print('On col %(c)04d of %(t)04d, about %(m)02d:%(s)02d left' \
              %{'c':n, 't':x_locs.size, 'm':rT//60, 's':rT%60})
        pT = cT
    # Normalize and return
    trans01 = trans01/trials
    return trans01

def graph_forces(forces, limits=(3,.5,15,3,.5,15)):
    x_locs = np.arange(limits[0], limits[2], limits[1]) 
    y_locs = np.arange(limits[3], limits[5], limits[4])
    Graph.title = "Forces"
    Graph.xlabel = "Offset (nm)"
    Graph.ylabel = "LS (nm)"
    Graph.quiver(x_locs, y_locs, forces)
    Graph.show()

def graph_energies(energs, limits=(-2,.5,15,-2,.5,15)):
    x_locs = np.arange(limits[0], limits[2], limits[1]) 
    y_locs = np.arange(limits[3], limits[5], limits[4])
    Graph.title = "Energy level of an TNCG crossbridge at different head locations"
    Graph.xlabel = "Location of XB head (nm)"
    Graph.ylabel = "Location of XB head (nm)"
    Graph.levels = [1, 4, 7, 10, 13]
    Graph.contour(x_locs, y_locs, energs)
    Graph.show()
