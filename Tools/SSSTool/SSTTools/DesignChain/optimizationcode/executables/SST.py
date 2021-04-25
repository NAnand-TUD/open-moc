#!/usr/bin/env python3
# 
################## FILE NAME: RST.py ##########################
#==========================================================================
# author: Jozef Stuijt				                                      |
# 	: Master Student,                                                     |
#	: Process and Energy Departmemt,                                      |
#	: TU Delft,                                                           |
#	: The Netherlands                                                     |
#                                                                         |
# email : joe.stuijt@gmail.com                                            | 
#																		  |
# 								                               			  |
# Description: MAIN CODE TO RUN RST: Basic Job-> Data Handeling           |
#																		  |
#=========================================================================|
#
#
# MAIN CODE
#
#---------------------------------------------------------------------------------------------#
## START: Initializing Packages 
import sys				            # for sys.exit(-1)
import time				            # for timing the code
#import matplotlib.pyplot as plt		# for plotting routines
import pdb				            # Debugging Module
import numpy as np
import math
from scipy.optimize import minimize
from srcSSS.ClassRST import *
from srcSSS.IO import *
try:                                       # check if import possible related to cluster
   import matplotlib.pyplot as plt         # for plotting routines
   imprtplt = True
except ImportError: 
   print('\nWarining: Plotting module import unsuccessful\n') 
   imprtplt = False
## END: Initializing Packages

# Read User Input file 'RST_Config.in' <Sample Input file>
#

#---------------------------------------------------------------------------------------------#
## START: Reading Input File
DIR					  = os.getcwd()+'/'

t0              	  = time.time()
try:    INFile		  = DIR+sys.argv[1]
except: INFile        = DIR+'RST_Config.in'

IN					  = ReadUserInput(INFile)
RADIAL			  = False
if IN['STATOR_KIND']== "RADIAL":
	RADIAL			  = True
	nBlades           = int(IN['nBlades'])
	pitch             = 360.0/nBlades
elif IN['STATOR_KIND']== "AXIAL":
	pitch 			  = float(IN['AxialPitch'])

flowAngle             = float(IN['flowAngle'])
radiusOut             = float(IN['radiusOut'])
outFileMoc            = IN['outFileMoc']
kernelRadiusRatio     = float(IN['kernelRadiusRatio'])
ActualRout            = float(IN['ActualRout'])
ActualRin             = float(IN['ActualRin'])
RealRout              = float(IN['RealRout'])
TEminT                = float(IN['TEminT'])
Mode 				  = IN['Mode']

try: ScaleMax             = int(IN['ScaleMax'])
except KeyError: ScaleMax = 50

try: 
	AreaRatio = float(IN['AreaRatio'])
	AR = True
except KeyError: AR = False#pass#Throat = 1e-3

try: ScaleNoz              = float(IN['ScaleNoz'])
except KeyError: ScaleNoz  = 1

try:
	ScaleMesh 			  = float(IN['ScaleMesh'])
	nPS					  = int(IN['nPointsScale'])
	UMG2Name 			  = IN['UMG2Name']
	SpecsName			  = IN['SpecsName']
	CoordsName			  = IN['CoordsName']
except KeyError: pass

try:
	plotName 			     = IN['plotName']
	if plotName != 'show': pltFrmt                  = plotName.split('.')[1]
	try: nBldPlt             = int(IN['nBldPlt'])
	except KeyError: nBldPlt = nBlades
except KeyError: pass


Stator = SuperSonicStator(IN)
Stator.NewBladeBuild()

print('\nStart: Building Blade\n')
if Mode=='Sim':
	writeCoords(Stator.Vane,CoordsName)
	print('\nStart: Creating Meshing Boundaries\n')
	Stator.MeshBoundariesBlade()
	if imprtplt:
		plotBoundaries(Stator)	
		Stator.plotStatorStage()
		plt.show()
	print('\nStart: Writting Output Files\n')
	writeUMG2out(Stator, UMG2Name, nPS, ScaleMesh)
	writeBladeSpecs(Stator, ScaleNoz, SpecsName)
	print('\nCompleted SST procedure in {0:.2f} seconds\n'.format(time.time() - t0))

if Mode=='Blade':
	print('\nStart: Plotting Blade and Writting Files \n')
	if imprtplt:
		#Stator.MeshBoundariesBlade()
		#plotBoundaries(Stator)
		#plotStatorStage(Stator.Vane,nBldPlt,pitch)
		print('\nCompleted SST procedure in {0:.2f} seconds\n'.format(time.time() - t0))
		Stator.plotStatorStage()
		if RADIAL:
			if (nBlades==nBldPlt): RadiusLimitPlot([RealRout,ActualRin],['--k','--k'],Stator,pitch)
			else: RadiusLimitPlot([RealRout,ActualRin],['--k','--k'],Stator,-pitch,nBldPlt)
		else:
			axes = plt.gca()
			ll,ul =axes.get_ylim()
			plt.plot([-RealRout,-RealRout],[ll,ul],'--k')
			plt.plot([-ActualRin,-ActualRin],[ll,ul],'--k')

if imprtplt:
	plt.axis('equal')
	if plotName == 'show': plt.show()
	else: plt.savefig(plotName,format=pltFrmt,dpi=2000)
#if not warn: print('\nCompleted RST procedure in {0:.2f} seconds\n'.format(time.time() - t0))
#else: print('\nCompleted RST procedure in {0:.2f} seconds WITHOUT CONVERGENCE\n'.format(time.time() - t0))
	
