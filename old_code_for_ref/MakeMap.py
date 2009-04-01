## Create a map of the energy density for the xxCG class crossbridge as a
## test of the CalculateEnergies functions

#Give us basic math
from math import pi, sin, cos, tan, atan2, sqrt, hypot

import CalculateEnergies as ce
print("CalculateEnergies is imported as ce")

import numpy as np
print("NumPy is imported as np")


x_vals = range(5,16)
y_vals = range(5,16)
#change the crossbridge type
xb = ce.xb()
xb.set_xb_class('TNCG', (pi/4, 8.5, pi/3, 6), (200, 10, 200, 5))

#energy_map = [[ce.xb_energy(ce.head_loc_to_xb_vals(xb, (x,y))) for x in x_vals] for y in y_vals]

#ce.plotContour(energy_map, x_vals, y_vals, "Energy of xxCG")
#ce.plt.show()

import XBWindow as win
can = win.new_root()
win.draw_XB(can, xb, (10,9))
