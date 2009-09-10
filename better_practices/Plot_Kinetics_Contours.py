#!/usr/bin/env python
# encoding: utf-8
"""
Plot_Kinetics_Contours.py
Created by Dave Williams on 2009-09-09.
"""

import sys
import Storage
from numpy import sqrt, exp, pi, tanh, log
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Load properties that will be needed
    store = [Storage.Storage(2), Storage.Storage(4)] 
    free_energy = [s.get("free_energy") for s in store]
    r12 = [np.multiply(s.get("r12"), 1000) for s in store]
    r23 = [np.multiply(s.get("r23"), 1000) for s in store]
    r31 = [np.multiply(s.get("r31"), 1000) for s in store]
    


if __name__ == '__main__':
    main()

