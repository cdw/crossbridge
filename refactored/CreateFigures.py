#!/usr/bin/env python
# encoding: utf-8
"""
CreateFigures.py

Created by Dave Williams on 2009-05-25.
"""

import matplotlib
from pylab import figure, show, colorbar
import numpy as np
import scipy.io
import ControlXB
import Crossbridge

## Configure matplotlib
# Use Helvetica as default font
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text', usetex=True)

MakexxCGfig = True
MakeTNCGfig = False
UseStoredDiffusion = True
Cuts = False
SeriesOfCuts = False

Range = (-5, .2, 15, 5, .2, 15)
T01Trials = 6000


if MakexxCGfig is True:
    ## Calculations
    # Energies
    xxCGe0 = ControlXB.create_energies("xxCG", xb_state=0, limits=Range)
    xxCGe1 = ControlXB.create_energies("xxCG", xb_state=1, limits=Range)
    xxCGe2 = ControlXB.create_energies("xxCG", xb_state=2, limits=Range)
    xxCGfe0 = ControlXB.create_free_energies("xxCG", xb_state=0, limits=Range)
    xxCGfe1 = ControlXB.create_free_energies("xxCG", xb_state=1, limits=Range)
    xxCGfe2 = ControlXB.create_free_energies("xxCG", xb_state=2, limits=Range)
    # Forward rates
    if UseStoredDiffusion is True:
        StoredMat = scipy.io.loadmat('./xxCGr01.mat')
        if np.all(StoredMat['Range']==np.array(Range)) is False:
            print('Error: Range update since last r01 calc')
        else:
            xxCGr01 = StoredMat['xxCGr01']
    else:
        # FIXME : Scaling factor present
        xxCGr01 = 0.001+ 3 * ControlXB.tran01("xxCG", limits=Range, trials=T01Trials)
        ToSave = {"Range":Range, "xxCGr01":xxCGr01}
        f = open('./xxCGr01.mat', 'w')
        scipy.io.savemat(f, ToSave)
        f.close()
    xxCGr12 = .5 * (1 + np.tanh(.6 * (xxCGe1 - xxCGe2)))+.001
    xxCGr20 = np.exp(-1 / xxCGe2)
    # Reverse rates
    xxCGr10 = xxCGr01 / np.exp(xxCGfe0 - xxCGfe1)
    xxCGr21 = xxCGr12 / np.exp(xxCGfe1 - xxCGfe2)
    xxCGr02 = 0 * xxCGr20 # See Tanner, 2007 Pg 1209 for justification 
    
    ## Plotting
    # For levels, axes, etc
    RT = 3.97
    EnergyLevels = (-3*RT, -1*RT, 1*RT, 3*RT, 5*RT, 7*RT, RT*10**10)
    TransLevels = (-0.001, .1 ,.3, .5, .7, .9, 1,  10**100)
    clrs = ('blue','cyan','lightgreen','yellow','orange','red','darkred')
    X = np.arange(Range[0], Range[2], Range[1])
    Y = np.arange(Range[3], Range[5], Range[4])
    # Make the figure and subplots
    fig0 = figure()
    ax00 = fig0.add_subplot(422)
    ax01 = fig0.add_subplot(423)
    ax02 = fig0.add_subplot(425)
    ax03 = fig0.add_subplot(427)
    ax04 = fig0.add_subplot(424)
    ax05 = fig0.add_subplot(426)
    ax06 = fig0.add_subplot(428)
    # Plot and annotate
    if Cuts is False and SeriesOfCuts is False:
        C00 = ax00.contourf(X, Y, xxCGfe1, EnergyLevels, 
                colors=clrs)
        ax00.set_title('Energy Profile')
        ax00.set_ylabel('Lattice Spacing (nm)')    
        C01 = ax01.contourf(X, Y, xxCGr01, TransLevels, 
                colors=clrs)
        #ax01.clabel(C01, inline=1, fontsize=10)
        ax01.set_title('r$_{\sffamily 01}$: unbound to loosely bound')
        ax01.set_ylabel('Lattice Spacing (nm)')
        C02 = ax02.contourf(X, Y, xxCGr12, TransLevels, 
                colors=clrs)       
        ax02.set_title('r$_{12}$ Transition rates')
        ax02.set_ylabel('Lattice Spacing (nm)')
        #fig0.colorbar(C02, shrink=0.8, extend='both')
        C03 = ax03.contourf(X, Y, xxCGr20, TransLevels, 
                colors=clrs)
        ax03.set_title('r$_{20}$ Transition rates')
        ax03.set_ylabel('Lattice Spacing (nm)')
        ax03.set_xlabel('Binding site offset (nm)')
        C04 = ax04.contourf(X, Y, xxCGr10, TransLevels, 
                colors=clrs)
        ax04.set_title('r$_{10}$ Transition rates')
        C05 = ax05.contourf(X, Y, xxCGr21, TransLevels, 
                colors=clrs)
        ax05.set_title('r$_{21}$ Transition rates')
        C06 = ax06.contourf(X, Y, xxCGr02, TransLevels, 
                colors=clrs)
        ax06.set_title('r$_{02}$ Transition rates')
        ax06.set_xlabel('Binding site offset (nm)')
    elif SeriesOfCuts is False:
        C00 = ax00.plot(X, xxCGfe1[xxCGfe1.shape[0]/2,:])
        ax00.set_title('Energy Profile at 10 nm')
        ax00.set_ylabel('Energy (RT)')    
        C01 = ax01.plot(X, xxCGr01[xxCGr01.shape[0]/2,:])
        #ax01.clabel(C01, inline=1, fontsize=10)
        ax01.set_ylim([-0.1,1.1])
        ax01.set_title('r$_{01}$ Transition rates at 10 nm')
        ax01.set_ylabel('Transition probability')
        C02 = ax02.plot(X, xxCGr12[xxCGr12.shape[0]/2,:])
        ax02.set_ylim([-0.1,1.1])
        ax02.set_title('r$_{12}$ Transition rates at 10 nm')
        ax02.set_ylabel('Transition probability')
        #fig0.colorbar(C02, shrink=0.8, extend='both')
        C03 = ax03.plot(X, xxCGr20[xxCGr20.shape[0]/2,:])
        ax03.set_ylim([-0.1,1.1])
        ax03.set_title('r$_{20}$ Transition rates at 10 nm')
        ax03.set_ylabel('Transition probability')
        ax03.set_xlabel('Binding site offset (nm)')
        C04 = ax04.plot(X, xxCGr10[xxCGr10.shape[0]/2,:])
        ax04.set_ylim([-0.1,1.1])
        ax04.set_title('r$_{10}$ Transition rates at 10 nm')
        C05 = ax05.plot(X, xxCGr21[xxCGr21.shape[0]/2,:])
        ax05.set_ylim([-0.1,1.1])
        ax05.set_title('r$_{21}$ Transition rates at 10 nm')
        C06 = ax06.plot(X, xxCGr02[xxCGr02.shape[0]/2+20,:])
        ax06.set_ylim([-0.1,1.1])
        ax06.set_title('r$_{02}$ Transition rates at 10 nm')
        ax06.set_xlabel('Binding site offset (nm)')
    elif SeriesOfCuts is True:
        #C00 = ax00.plot(X, xxCGfe1[xxCGfe1.shape[0]/2,:])
        #ax00.set_title('Energy Profile at 10 nm')
        #ax00.set_ylabel('Energy (RT)')    
        C01 = ax01.plot(X, xxCGr01[3*xxCGr01.shape[0]/10,:],
                        X, xxCGr01[4*xxCGr01.shape[0]/10,:],
                        X, xxCGr01[5*xxCGr01.shape[0]/10,:])
        ax01.legend(('8 nm', '9 nm', '10 nm'), 'upper right')
        ax01.set_ylim([-0.1,1.1])
        ax01.set_title('r$_{01}$ Transition rates at 8, 9, and 10 nm')
        ax01.set_ylabel('Transition probability')
        C02 = ax02.plot(X, xxCGr12[3*xxCGr12.shape[0]/10,:],
                        X, xxCGr12[5*xxCGr12.shape[0]/10,:],
                        X, xxCGr12[7*xxCGr12.shape[0]/10,:])
        ax02.legend(('8 nm', '9 nm', '10 nm'), 'upper right')
        ax02.set_ylim([-0.1,1.1])
        ax02.set_title('r$_{12}$ Transition rates at 8, 9, and 10 nm')
        # FIXME This is really at 8, 10, and 12 nm
        ax02.set_ylabel('Transition probability')
        #fig0.colorbar(C02, shrink=0.8, extend='both')
        C03 = ax03.plot(X, xxCGr20[3*xxCGr20.shape[0]/10,:],
                        X, xxCGr20[5*xxCGr20.shape[0]/10,:],
                        X, xxCGr20[7*xxCGr20.shape[0]/10,:])
        ax03.legend(('8 nm', '9 nm', '10 nm'), 'lower right')
        ax03.set_ylim([-0.1,1.1])
        ax03.set_title('r$_{20}$ Transition rates at 8, 9, and 10 nm')
        ax03.set_ylabel('Transition probability')
        ax03.set_xlabel('Binding site offset (nm)')
        #C04 = ax04.plot(X, xxCGr10[xxCGr10.shape[0]/2,:])
        #ax04.set_ylim([-0.1,1.1])
        #ax04.set_title('r$_{10}$ Transition rates at 10 nm')
        #C05 = ax05.plot(X, xxCGr21[xxCGr21.shape[0]/2,:])
        #ax05.set_ylim([-0.1,1.1])
        #ax05.set_title('r$_{21}$ Transition rates at 10 nm')
        #C06 = ax06.plot(X, xxCGr02[xxCGr02.shape[0]/2+20,:])
        #ax06.set_ylim([-0.1,1.1])
        #ax06.set_title('r$_{02}$ Transition rates at 10 nm')
        #ax06.set_xlabel('Binding site offset (nm)')
            
        

    

if MakeTNCGfig is True:
    ## Calculations
    # Energies
    TNCGe0 = ControlXB.create_energies("TNCG", xb_state=0, limits=Range)
    TNCGe1 = ControlXB.create_energies("TNCG", xb_state=1, limits=Range)
    TNCGe2 = ControlXB.create_energies("TNCG", xb_state=2, limits=Range)
    TNCGfe0 = ControlXB.create_free_energies("TNCG", xb_state=0, limits=Range)
    TNCGfe1 = ControlXB.create_free_energies("TNCG", xb_state=1, limits=Range)
    TNCGfe2 = ControlXB.create_free_energies("TNCG", xb_state=2, limits=Range)
    # Forward rates
    if UseStoredDiffusion is True:
        StoredMat = scipy.io.loadmat('./TNCGr01.mat')
        if np.all(StoredMat['Range']==np.array(Range)) is False:
            print('Error: Range update since last r01 calc')
        else:
            TNCGr01 = StoredMat['TNCGr01']
    else:
        # FIXME : Scaling factor present
        TNCGr01 = 0.001+ 3 * ControlXB.tran01("TNCG", limits=Range, trials=T01Trials)
        ToSave = {"Range":Range, "TNCGr01":TNCGr01}
        f = open('./TNCGr01.mat', 'w')
        scipy.io.savemat(f, ToSave)
        f.close()
    TNCGr12 = .5 * (1 + np.tanh(.6 * (TNCGe1 - TNCGe2)))+.001
    TNCGr20 = np.exp(-1 / TNCGe2)
    # Reverse rates
    TNCGr10 = TNCGr01 / np.exp(TNCGfe0 - TNCGfe1)
    TNCGr21 = TNCGr12 / np.exp(TNCGfe1 - TNCGfe2)
    TNCGr02 = 0 * TNCGr20 # See Tanner, 2007 Pg 1209 for justification 

    ## Plotting
    # For levels, axes, etc
    RT = 3.97
    EnergyLevels = (-3*RT, -1*RT, 1*RT, 3*RT, 5*RT, 7*RT, RT*10**10)
    TransLevels = (-0.001, .1 ,.3, .5, .7, .9, 1,  10**100)
    clrs = ('blue','cyan','lightgreen','yellow','orange','red','darkred')
    X = np.arange(Range[0], Range[2], Range[1])
    Y = np.arange(Range[3], Range[5], Range[4])
    # Make the figure and subplots
    fig1 = figure()
    ax10 = fig1.add_subplot(422)
    ax11 = fig1.add_subplot(423)
    ax12 = fig1.add_subplot(425)
    ax13 = fig1.add_subplot(427)
    ax14 = fig1.add_subplot(424)
    ax15 = fig1.add_subplot(426)
    ax16 = fig1.add_subplot(428)
    # Plot and annotate
    if Cuts is False:
        C10 = ax10.contourf(X, Y, TNCGfe1, EnergyLevels, 
                colors=clrs)
        ax10.set_title('Energy Profile')
        ax10.set_ylabel('Lattice Spacing (nm)')    
        C11 = ax11.contourf(X, Y, TNCGr01, TransLevels, 
                colors=clrs)
        #ax11.clabel(C01, inline=1, fontsize=10)
        ax11.set_title('r$_{01}$ Transition rates')
        ax11.set_ylabel('Lattice Spacing (nm)')
        C12 = ax12.contourf(X, Y, TNCGr12, TransLevels, 
                colors=clrs)
        ax12.set_title('r$_{12}$ Transition rates')
        ax12.set_ylabel('Lattice Spacing (nm)')
        #fig1.colorbar(C02, shrink=0.8, extend='both')
        C13 = ax13.contourf(X, Y, TNCGr20, TransLevels, 
                colors=clrs)
        ax13.set_title('r$_{20}$ Transition rates')
        ax13.set_ylabel('Lattice Spacing (nm)')
        ax13.set_xlabel('Binding site offset (nm)')
        C14 = ax14.contourf(X, Y, TNCGr10, TransLevels, 
                colors=clrs)
        ax14.set_title('r$_{10}$ Transition rates')
        C15 = ax15.contourf(X, Y, TNCGr21, TransLevels, 
                colors=clrs)
        ax15.set_title('r$_{21}$ Transition rates')
        C16 = ax16.contourf(X, Y, TNCGr02, TransLevels, 
                colors=clrs)
        ax16.set_title('r$_{02}$ Transition rates')
        ax16.set_xlabel('Binding site offset (nm)')
    else:
        C10 = ax10.plot(X, TNCGfe1[TNCGfe1.shape[0]/2-20,:])
        ax10.set_title('Energy Profile at 10 nm')
        ax10.set_ylabel('Energy (RT)')    
        C11 = ax11.plot(X, TNCGr01[TNCGr01.shape[0]/2-20,:])
        #ax11.clabel(C11, inline=1, fontsize=10)
        ax11.set_ylim([-0.1,1.1])
        ax11.set_title('r$_{01}$ Transition rates at 10 nm')
        ax11.set_ylabel('Transition probability')
        C12 = ax12.plot(X, TNCGr12[TNCGr12.shape[0]/2-20,:])
        ax12.set_ylim([-0.1,1.1])
        ax12.set_title('r$_{12}$ Transition rates at 10 nm')
        ax12.set_ylabel('Transition probability')
        #fig1.colorbar(C12, shrink=0.8, extend='both')
        C13 = ax13.plot(X, TNCGr20[TNCGr20.shape[0]/2-20,:])
        ax13.set_ylim([-0.1,1.1])
        ax13.set_title('r$_{20}$ Transition rates at 10 nm')
        ax13.set_ylabel('Transition probability')
        ax13.set_xlabel('Binding site offset (nm)')
        C14 = ax14.plot(X, TNCGr10[TNCGr10.shape[0]/2-20,:])
        ax14.set_ylim([-0.1,1.1])
        ax14.set_title('r$_{10}$ Transition rates at 10 nm')
        C15 = ax15.plot(X, TNCGr21[TNCGr21.shape[0]/2-20,:])
        ax15.set_ylim([-0.1,1.1])
        ax15.set_title('r$_{21}$ Transition rates at 10 nm')
        C16 = ax16.plot(X, TNCGr02[TNCGr02.shape[0]/2-20,:])
        ax16.set_ylim([-0.1,1.1])
        ax16.set_title('r$_{02}$ Transition rates at 10 nm')
        ax16.set_xlabel('Binding site offset (nm)')
        


# Draw what we've written and show it on the screen
if MakexxCGfig is True and MakeTNCGfig is True:    
    fig0.canvas.draw()
    fig1.canvas.draw()
    show()
elif MakexxCGfig is True:
    fig0.canvas.draw()
    show()
elif MakeTNCGfig is True:
    fig1.canvas.draw()
    show()
    
