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
Cuts = True
SeriesOfCuts = False

Range = (-5, .2, 15, 5, .2, 15)
T01Trials = 600


if MakexxCGfig is True:
    ## Calculations
    # Energies
    e0 = ControlXB.create_energies("xxCG", xb_state=0, limits=Range)
    e1 = ControlXB.create_energies("xxCG", xb_state=1, limits=Range)
    e2 = ControlXB.create_energies("xxCG", xb_state=2, limits=Range)
    fe0 = ControlXB.create_free_energies("xxCG", xb_state=0, limits=Range)
    fe1 = ControlXB.create_free_energies("xxCG", xb_state=1, limits=Range)
    fe2 = ControlXB.create_free_energies("xxCG", xb_state=2, limits=Range)
    # Forward rates
    if UseStoredDiffusion is True:
        StoredMat = scipy.io.loadmat('./xxCGr01.mat')
        if np.all(StoredMat['Range']==np.array(Range)) is False:
            print('Error: Range update since last r01 calc')
        else:
            r01 = StoredMat['xxCGr01']
    else:
        # FIXME : Scaling factor present
        r01 = 0.001+ 3 * ControlXB.tran01("xxCG", 
                            limits=Range, trials=T01Trials)
        ToSave = {"Range":Range, "xxCGr01":r01}
        f = open('./xxCGr01.mat', 'w')
        scipy.io.savemat(f, ToSave)
        f.close()
    r12 = .5 * (1 + np.tanh(.6 * (e1 - e2)))+.001
    r20 = np.exp(-1 / e2)
    # Reverse rates
    r10 = r01 / np.exp(fe0 - fe1)
    r21 = r12 / np.exp(fe1 - fe2)
    r02 = 0 * r20 # See Tanner, 2007 Pg 1209 for justification 
elif MakeTNCGfig is True:
    ## Calculations
    # Energies
    e0 = ControlXB.create_energies("TNCG", xb_state=0, limits=Range)
    e1 = ControlXB.create_energies("TNCG", xb_state=1, limits=Range)
    e2 = ControlXB.create_energies("TNCG", xb_state=2, limits=Range)
    fe0 = ControlXB.create_free_energies("TNCG", xb_state=0, limits=Range)
    fe1 = ControlXB.create_free_energies("TNCG", xb_state=1, limits=Range)
    fe2 = ControlXB.create_free_energies("TNCG", xb_state=2, limits=Range)
    # Forward rates
    if UseStoredDiffusion is True:
        StoredMat = scipy.io.loadmat('./TNCGr01.mat')
        if np.all(StoredMat['Range']==np.array(Range)) is False:
            print('Error: Range update since last r01 calc')
        else:
            r01 = StoredMat['TNCGr01']
    else:
        # FIXME : Scaling factor present
        r01 = 0.001+ 3 * ControlXB.tran01("TNCG", 
                            limits=Range, trials=T01Trials)
        ToSave = {"Range":Range, "TNCGr01":r01}
        f = open('./TNCGr01.mat', 'w')
        scipy.io.savemat(f, ToSave)
        f.close()
    r12 = .5 * (1 + np.tanh(.6 * (e1 - e2)))+.001
    r20 = np.exp(-1 / e2)
    # Reverse rates
    r10 = r01 / np.exp(fe0 - fe1)
    r21 = r12 / np.exp(fe1 - fe2)
    r02 = 0 * r20 # See Tanner, 2007 Pg 1209 for justification 

## Plotting Functions
def twoContour(figure, loc, X, Y, Z , levels, 
               xlab=None, ylab=None, title=None):
    '''Take care of contour plotting'''
    axis = figure.add_subplot(loc)
    axis.contourf(X, Y, Z, levels, colors=fclrs)
    axis.contour(X, Y, Z, levels, colors=clrs)
    if title is not None:
        axis.set_title(title)
    if ylab is not None:
        axis.set_ylabel(ylab)
    if xlab is not None:
        axis.set_xlabel(xlab)

def cutPlot(figure, loc, X, Y, Z, CutLoc=10,
            xlab=None, ylab=None, title=None):
    '''Take care of cut plotting'''
    axis = figure.add_subplot(loc)
    CutInd = np.flatnonzero(np.logical_and(Y>CutLoc-0.001, Y<CutLoc+0.001))[0]
    axis.plot(X, Z[CutInd, :],'k')
    axis.set_ylim([-0.1,1.1])
    if title is not None:
        axis.set_title(title)
    if ylab is not None:
        axis.set_ylabel(ylab)
    if xlab is not None:
        axis.set_xlabel(xlab)
    
## Plotting
# For levels, axes, etc
RT = 3.97
EnergyLevels = (-100*RT, -1*RT, 0*RT, 3*RT, 6*RT, RT*10**10)
TransLevels = (-0.001, .2 ,.4, .6, .8, 1,  10**100)
#clrs = ('blue','cyan','lightgreen','yellow','orange','red','darkred')
fclrs = ('.2', '.3', '.4', '.5', '.6', '.7', '.8')
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
clrs = ('white')
X = np.arange(Range[0], Range[2], Range[1])
Y = np.arange(Range[3], Range[5], Range[4])
Extent = (Range[0], Range[2], Range[3], Range[5])
Cmap = matplotlib.cm.gray
# Make the figure
fig0 = figure()
fig0.subplots_adjust(wspace=0.25, hspace=0.50,
                    left=0.10, right=0.95,
                    top=0.94, bottom=0.08)
# Plot and annotate
if Cuts is False and SeriesOfCuts is False:
    twoContour(fig0, 422, X, Y, fe1 , EnergyLevels, 
                    title='Loosely bound energy', 
                    ylab='Lattice Dist (nm)', 
                    xlab=None)
    twoContour(fig0, 423, X, Y, r01 , TransLevels, 
                    title='r$_{01}$: unbound to loosely bound', 
                    ylab='Lattice Dist (nm)', 
                    xlab=None)
    twoContour(fig0, 425, X, Y, r12 , TransLevels, 
                    title='r$_{12}$: loosely to strongly bound', 
                    ylab='Lattice Dist (nm)', 
                    xlab=None)
    twoContour(fig0, 427, X, Y, r20 , TransLevels, 
                    title='r$_{20}$: strongly bound to unbound', 
                    ylab='Lattice Dist (nm)', 
                    xlab='Binding site offset (nm)')
    twoContour(fig0, 424, X, Y, r10 , TransLevels, 
                    title='r$_{10}$: loosely bound to unbound', 
                    ylab=None, 
                    xlab=None)
    twoContour(fig0, 426, X, Y, r21 , TransLevels, 
                    title='r$_{21}$: strongly to loosely bound', 
                    ylab=None, 
                    xlab=None)
    twoContour(fig0, 428, X, Y, r02 , TransLevels, 
                    title='r$_{02}$: unbound to strongly bound', 
                    ylab=None, 
                    xlab='Binding site offset (nm)')
elif SeriesOfCuts is False:
    CutHere = 10.0
    cutPlot(fig0, 422, X, Y, fe1, 
                CutLoc=CutHere,
                title='Energy Profile at 10 nm',
                ylab='Energy (RT)',
                xlab=None)
    cutPlot(fig0, 423, X, Y, r01, 
                CutLoc=CutHere,
                title='r$_{01}$ Transition rates at 10 nm', 
                ylab='Transition probability', 
                xlab=None)        
    cutPlot(fig0, 425, X, Y, r12, 
                CutLoc=CutHere,
                title='r$_{12}$ Transition rates at 10 nm', 
                ylab='Transition probability', 
                xlab=None)
    cutPlot(fig0, 427, X, Y, r20, 
                CutLoc=CutHere,
                title='r$_{20}$ Transition rates at 10 nm', 
                ylab='Transition probability', 
                xlab='Binding site offset (nm)')
    cutPlot(fig0, 424, X, Y, r10, 
                CutLoc=CutHere,
                title='r$_{10}$ Transition rates at 10 nm', 
                ylab=None, 
                xlab=None)
    cutPlot(fig0, 426, X, Y, r21, 
                CutLoc=CutHere,
                title='r$_{21}$ Transition rates at 10 nm', 
                ylab=None, 
                xlab=None)
    cutPlot(fig0, 428, X, Y, r02, 
                CutLoc=CutHere,
                title='r$_{02}$ Transition rates at 10 nm', 
                ylab=None, 
                xlab='Binding site offset (nm)')
elif SeriesOfCuts is True:
    #C00 = ax00.plot(X, xxCGfe1[xxCGfe1.shape[0]/2,:])
    #ax00.set_title('Energy Profile at 10 nm')
    #ax00.set_ylabel('Energy (RT)')    
    C01 = ax01.plot(X, r01[3*r01.shape[0]/10,:],
                    X, r01[4*r01.shape[0]/10,:],
                    X, r01[5*r01.shape[0]/10,:])
    ax01.legend(('8 nm', '9 nm', '10 nm'), 'upper right')
    ax01.set_ylim([-0.1,1.1])
    ax01.set_title('r$_{01}$ Transition rates at 8, 9, and 10 nm')
    ax01.set_ylabel('Transition probability')
    C02 = ax02.plot(X, r12[3*r12.shape[0]/10,:],
                    X, r12[5*r12.shape[0]/10,:],
                    X, r12[7*r12.shape[0]/10,:])
    ax02.legend(('8 nm', '9 nm', '10 nm'), 'upper right')
    ax02.set_ylim([-0.1,1.1])
    ax02.set_title('r$_{12}$ Transition rates at 8, 9, and 10 nm')
    # FIXME This is really at 8, 10, and 12 nm
    ax02.set_ylabel('Transition probability')
    #fig0.colorbar(C02, shrink=0.8, extend='both')
    C03 = ax03.plot(X, r20[3*r20.shape[0]/10,:],
                    X, r20[5*r20.shape[0]/10,:],
                    X, r20[7*r20.shape[0]/10,:])
    ax03.legend(('8 nm', '9 nm', '10 nm'), 'lower right')
    ax03.set_ylim([-0.1,1.1])
    ax03.set_title('r$_{20}$ Transition rates at 8, 9, and 10 nm')
    ax03.set_ylabel('Transition probability')
    ax03.set_xlabel('Binding site offset (nm)')
    


# Draw what we've written and show it on the screen
fig0.canvas.draw()
show()



    