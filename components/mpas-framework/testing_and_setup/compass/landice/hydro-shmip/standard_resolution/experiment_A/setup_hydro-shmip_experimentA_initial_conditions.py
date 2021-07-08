#!/usr/bin/env python
'''
Generate initial conditions for hydro-SHMIP land ice test case A
Details here: http://shmip.bitbucket.org/
'''

from __future__ import absolute_import, division, print_function, unicode_literals

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import sys

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file to setup", metavar="FILE")
parser.add_option("-n", "--number", dest="number", type='int', help="test variant to set up, 1-6", metavar="NUMBER")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print('No file specified.  Attempting to use landice_grid.nc')

# Setup dictionaries of parameter values for each experiment
a_params = {1:7.93e-11,  2:1.59e-09, 3:5.79e-09, 4:2.5e-8, 5:4.5e-8, 6:5.79e-7}

# Open the file, get needed dimensions
gridfile = NetCDFFile(options.filename,'r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
# Get variables
xCell = gridfile.variables['xCell']
yCell = gridfile.variables['yCell']
xEdge = gridfile.variables['xEdge']
yEdge = gridfile.variables['yEdge']
xVertex = gridfile.variables['xVertex']
yVertex = gridfile.variables['yVertex']
thickness = gridfile.variables['thickness']
bedTopography = gridfile.variables['bedTopography']
layerThicknessFractions = gridfile.variables['layerThicknessFractions']
SMB = gridfile.variables['sfcMassBal']
waterFluxMask = gridfile.variables['waterFluxMask']


# adjust mesh origin if it hasn't been done yet
if yVertex[:].min() != 0.0:
   # x=100km should be value on right vertices
   # this would be the average of the last two vertices
   # which is the second to last edge
   xShift = np.unique(xEdge[:])[-2] - 100.0e3
   xCell[:] = xCell[:] - xShift
   xEdge[:] = xEdge[:] - xShift
   xVertex[:] = xVertex[:] - xShift

   # y origin should have yVertex of bottom row at y=0
   yShift = yVertex[:].min()
   yCell[:] = yCell[:] - yShift
   yEdge[:] = yEdge[:] - yShift
   yVertex[:] = yVertex[:] - yShift
print("Min/Max yVertex: {}  {}".format(yVertex[:].min(), yVertex[:].max()))
if np.absolute(np.unique(yVertex[:])[-2] - np.unique(yVertex[:])[1] - 20000.0) > gridfile.variables['dcEdge'][:].min() * 0.866:
   sys.exit('Error: Domain width is too far from 20000.0! Adjust ny in periodic_hex namelist and/or --remove_extra_y argument to mark_periodic_boundaries_for_culling.py.  Target is the difference between the second to last vertices in north and south to be within 0.866*dcEdge of 20000.0.')

# thickness
x = 100.0e3 - np.absolute(xCell[:] - 100.0e3)
thickness[0,:] = 6.0 *( np.sqrt( x[:] + 5.0e3) - np.sqrt(5.0e3) ) + 1.0
thickness[0, x<0.0]=0.0  # make a margin

# flat bed
bedTopography[:] = 0.0

# Set up no flux mask along north and south and east
waterFluxMask = gridfile.variables['waterFluxMask'][0,:]
waterFluxMask[np.nonzero(yEdge[:]==np.unique(yEdge[:])[0])] = 2
waterFluxMask[np.nonzero(yEdge[:]==np.unique(yEdge[:])[-1])] = 2
waterFluxMask[np.nonzero(xEdge[:]==np.unique(xEdge[:])[-3])] = 2  # last 3 edge values
waterFluxMask[np.nonzero(xEdge[:]==np.unique(xEdge[:])[-2])] = 2  # last 3 edge values
waterFluxMask[np.nonzero(xEdge[:]==np.unique(xEdge[:])[-1])] = 2  # last 3 edge values
gridfile.variables['waterFluxMask'][0,:] = waterFluxMask


print("Setting up test number: {}".format(options.number))

# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

print("Using water input value (m/s) of: ".format(a_params[options.number])
waterInput = gridfile.variables['externalWaterInput'][0,:]
waterInput[:] = a_params[options.number] * 1000.0  # Convert from m/s to kg/m2/s
#scale down the source term in the 'tents' in the N and S rows - there is no x-dir throughflow here so excess water is introduced
scaleFactor = 2.0/3.0 # for perfect hexagon
waterInput[np.nonzero(yCell[:]==np.unique(yCell[:])[0])] = a_params[options.number] * 1000.0 * scaleFactor
waterInput[np.nonzero(yCell[:]==np.unique(yCell[:])[-1])] = a_params[options.number] * 1000.0 * scaleFactor
gridfile.variables['externalWaterInput'][0,:] = waterInput


# melt
gridfile.variables['basalMeltInput'][:] = 0.0 * 1000.0 # no basal melting

# velocity
gridfile.variables['uReconstructX'][:] = 1.0e-6 # all tests use this velocity

# initial thickness
gridfile.variables['waterThickness'][0,:] = 0.05 * np.absolute(xCell[:] - 100.0e3)/100.0e3 # start with something to avoid weird adapative time steps initially

# initial waterPressure
#gridfile.variables['waterPressure'][:] = 0.5 * 9.81 * 910.0 * thickness[:] * (0.97 + 0.03 * np.random.rand(thickness[:].size)) # start with half of Pice.  Last factor adds some random noise to disrupt symmetry.
gridfile.variables['waterPressure'][:] = 0.5 * 9.81 * 910.0 * thickness[:]

gridfile.close()

print('Successfully added initial conditions to: ' + options.filename)

