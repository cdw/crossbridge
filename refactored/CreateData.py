#!/usr/bin/env python
# encoding: utf-8
"""
CreateData.py
Created by Dave Williams on 2009-06-28.

Creates sets of data for a given crossbridge type

Crossbridge types may be specified by their spring number: 1, 2, or 4
Data produced may be specified by their range, precision and type
"""
#
#import sys
#import getopt
#
#
#help_message = '''
#Useage: CreateData [OPTIONS] CrossBridgeType
#Where CrossBridgeType may be specified by a spring number: 1, 2, or 4
#
#Options include:
#Data produced may be specified by their range, precision and type
#'''
#
#
#class Usage(Exception):
#    def __init__(self, msg):
#        self.msg = msg
#
#
#def main(argv=None):
#    if argv is None:
#        argv = sys.argv
#    try:
#        try:
#            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
#        except getopt.error, msg:
#            raise Usage(msg)
#    
#        # option processing
#        for option, value in opts:
#            if option == "-v":
#                verbose = True
#            if option in ("-h", "--help"):
#                raise Usage(help_message)
#            if option in ("-o", "--output"):
#                output = value
#    
#    except Usage, err:
#        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
#        print >> sys.stderr, "\t for help use --help"
#        return 2
#
#
#if __name__ == "__main__":
#    sys.exit(main())

CrossBrigeType = 1 #1, 2, or 4
DataRange = (-5, 1, 15, 5, .5, 15) #(xmin, xstep, xmax, ymin, ystep, ymax)
ProduceEnergies = True
ProduceForces = True
ProduceTransitions = True
T01Trials = 1000
#        config={'Ts': pi/4, 'Ns': 7,  'Cs': pi/3+(pi-pi/4),      'Gs': 3, 
#                'Tr': pi/4, 'Nr': 7,  'Cr': pi/3+(pi-pi/4)+pi/6, 'Gr': 3, 
#                'Tk': 100,  'Nk': 10, 'Ck': 100,                 'Gk': 10}


def fname():
    """docstring for fname"""
    pass