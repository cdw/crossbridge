#!/usr/bin/env python
# encoding: utf-8
#
# Replicate the transition rate figure from Bert's 2007 PLoS paper
# CDW 20090504


from pylab import figure, show
import numpy as np
from numpy import sqrt, exp, pi, tanh, log

## Calculate energies
R = 8.314
T = 288
Gatp = 13 # In units of RT
atp = 5  * 10**-3
adp = 30 * 10**-6
phos  = 3  * 10**-3
DeltaG = abs(-Gatp - log(atp / (adp * phos)))
alpha = 0.28
eta = 0.68


## Define parameters
A = 2000  # From Tanner, 2008 Pg 1209
B = 100   # Ditto for C through P
C = 1
D = 1
M = 3600
N = 40
P = 20


## Calculate crossbridge spring values
k_xb = 5 / 3.976 # From Mathematica
# 3.976 is provided by entereing the following into Mathematica:
#   << PhysicalConstants`
#   << Units` 
#   Convert[(MolarGasConstant * 288 Kelvin * AvogadroConstant^-1)
#       /(Nano Meter)^2, (Pico Newton)/(Nano Meter)]
xb_0 = sqrt(eta * DeltaG / k_xb)


## Create functions to provide free energy in a given location and state:
g_1 = lambda x: 0 * x
g_2 = lambda x: alpha * -DeltaG + k_xb * (x - xb_0)**2
g_3 = lambda x: eta * -DeltaG + k_xb * x**2


## Create functions to yield transition rate at a given location
r_12 = lambda x: A * sqrt(k_xb / (2 * pi)) * exp(-.5 * k_xb * (x - xb_0)**2)
r_23 = lambda x: B / sqrt(k_xb) * (1 - tanh(C * sqrt(k_xb) * (x - xb_0))) + D
r_31 = lambda x: sqrt(k_xb) * (sqrt(M * x**2) - N * x) + P


## Create functions to derive reverse transition rates using the thermo eqn:
## r_forward(x)/r_backward(x) = exp(G_next(x) - G_prev(x))
r_21 = lambda x: r_12(x) / exp(g_1(x) - g_2(x))
r_32 = lambda x: r_23(x) / exp(g_2(x) - g_3(x))
r_13 = lambda x: 0 * x # See Tanner, 2007 Pg 1209 for justification


## Create the data for plotting
locs = np.arange(-5, 15, .1)
energy_1 = g_1(locs)  # Free energy state 1
energy_2 = g_2(locs)  # Free energy state 2
energy_3 = g_3(locs)  # Free energy state 3
rates_12 = r_12(locs) # Forward rate 12
rates_21 = r_21(locs) # Reverse rate 21
rates_23 = r_23(locs) # Forward rate 23
rates_32 = r_32(locs) # Reverse rate 32
rates_31 = r_31(locs) # Forward rate 31
rates_13 = r_13(locs) # Reverse rate 13


## Truncate some values for proper vector file generation
energy_axes = [-5, 15, -30, 5]
rate_axes = [-5, 15, -20, 1000]
eamax = energy_axes[3]*2
ramax = rate_axes[3]*2
energy_2 = [v*(v<=eamax)+eamax*(v>eamax) for v in energy_2]
energy_3 = [v*(v<=eamax)+eamax*(v>eamax) for v in energy_3]
rates_21 = [v*(v<=ramax)+ramax*(v>ramax) for v in rates_21]
rates_32 = [v*(v<=ramax)+ramax*(v>ramax) for v in rates_32]


## Produce the actual figure
# Make the figure and subplots
fig = figure()
ax0 = fig.add_subplot(221)
ax1 = fig.add_subplot(222)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(224)
# Plot/annotate the free energies
ax0.plot(locs, energy_1, 'k-', label='$G_1(x)$')
ax0.plot(locs, energy_2, 'k--')
ax0.plot(locs, energy_3, 'k:')
ax0.axis(energy_axes)
ax0.annotate('$G_1(x)$', [10, -4])
ax0.annotate('$G_2(x)$', [5.5, -8])
ax0.annotate('$G_3(x)$', [1.5, -18])
ax0.set_title('Energy states')
ax0.set_ylabel('Free Energy (RT)')
# Plot/annotate the r_12 and r_21 transitions
ax1.plot(locs, rates_12, 'k-')
ax1.plot(locs, rates_21, 'k--')
ax1.axis(rate_axes)
ax1.annotate('$r_{x,12}(x)$', [5, 600], alpha=1.0, backgroundcolor='w')
ax1.annotate('$r_{x,21}(x)$', [7, 200])
ax1.set_title('Binding rates')
ax1.set_ylabel('Transition Rate (s$^{-1}$)')
# Plot/annotate the r_23 and r_32 transitions
ax2.plot(locs, rates_23, 'k-')
ax2.plot(locs, rates_32, 'k--')
ax2.axis(rate_axes)
ax2.annotate('$r_{x,23}(x)$', [5, 50])
ax2.annotate('$r_{x,32}(x)$', [3.5, 700])
ax2.set_title('Strong binding rates')
ax2.set_xlabel('x (nm)')
ax2.set_ylabel('Transition Rate (s$^{-1}$)')
# Plot/annotate the r_31 and r_13 transitions
ax3.plot(locs, rates_31, 'k-')
ax3.plot(locs, rates_13, 'k--')
ax3.axis(rate_axes)
ax3.annotate('$r_{x,31}(x)$', [-2.5, 300])
ax3.annotate('$r_{x,13}(x)$', [9, 10])
ax3.set_title('Detachment rates')
ax3.set_xlabel('x (nm)')
ax3.set_ylabel('Transition Rate (s$^{-1}$)')
# Draw what we've written and show it on the screen
fig.canvas.draw()
show()