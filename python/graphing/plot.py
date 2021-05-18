#!/bin/env/python

# Sanjeev Kolli and Jonas Kaufman
# May 17, 2018
# Script to label sites in Monte Carlo supercells

from __future__ import print_function


# TincArgumentParser to allow command line args or json config as only command line flag
import sys, os

from tinc import *

parser = TincArgumentParser(description='Graph generator')

parser.add_argument('temp_interest', type=float, help='value to highlight', default = 300.0,nargs='?')
parser.add_argument('xlabel', type=str, help='x-axis label', default = 'x',nargs='?')
parser.add_argument('ylabel', type=str, help='y-axis label', default = 'y',nargs='?')
parser.add_argument('miny', type=float, help='minimum y value to display', default = 0.0,nargs='?')
parser.add_argument('maxy', type=float, help='maximum y value to display',default = 0.0, nargs='?')
parser.add_argument('inxFile', type=str, help='input filename', default = 'inx.bin', nargs='?')
parser.add_argument('inyFile', type=str, help='input filename', default = 'iny.bin', nargs='?')
parser.add_argument('__output_dir', type=str, help='Output directory', default = 'cached_output', nargs='?')
parser.add_argument('__output_names', type=str, help='output filename', default = 'out.png', nargs='?')

args = parser.get_args()

#print(args)

# Draw graph
import numpy as np
import matplotlib.pyplot as plt
import sys

datax = np.fromfile(args['inxFile'], dtype=np.float64)
datay = np.fromfile(args['inyFile'], dtype=np.float64)

#print(datay, file=sys.stderr)
#print(datax, file=sys.stderr)
fig, ax1 = plt.subplots(figsize=(4,3))

ax1.set_xlabel(args['xLabel'])
ax1.set_ylabel(args['yLabel'])
ax1.plot(datax, datay)

plt.ylim((float(args['miny'])* 0.95,float(args['maxy'])* 1.05))

ax1.vlines(float(args['temp_interest']), float(args['miny'])* 0.95, float(args['maxy'])* 1.05, 'r')

fig.tight_layout()  # otherwise the right y-label is slightly clipped

print("Saving:" + args['__output_dir'] + args['__output_names'][0])
plt.savefig(args['__output_dir'] + args['__output_names'][0], dpi=300)


