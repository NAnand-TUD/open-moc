#!/usr/bin/env python3
################## FILE NAME: MAIN.py ##########################
#================================================================
# author: Nitish Anand    & Jozef Stuijt                         |
#     :Master Student,                                           |
#    :Process and Energy Departmemt,                             |
#    :TU Delft,                                                  |
#    :The Netherlands                                            |
#                                                                |
# email : nitish.ug2@gmail.com                                   |
#																 |
# Description: MAIN CODE TO RUN MOC: Basic Job-> Data Handeling	 |
#																 |
#================================================================

#********************************************************************
#*							*** READ ME***							*
#********************************************************************
#																	*
# Limitation: (*) Always run the code in Python3					*
#	    : (*) Specify correct gas properties [Else M<1 error]		*
#	    : (*) if n is very large, dtau should be small				*
#																	*
# MODULES   : Name are in CAPS		eg: CLASSES.py					*
#																	*
# FUNCTION  : 1st char CAPS		eg:SauerAnalysis()					*
#								    								*
#********************************************************************
#
#
# MAIN CODE
#
#---------------------------------------------------------------------------------------------#
## START: Initializing Packages 
import sys				                    # for sys.exit(-1)
import time				                    # for timing the code
import os				                    # Check Existance of File
import pdb				                    # Debugging Module
import numpy as np
from srcMOC.CLASSES import *
from srcMOC.INTERIOR import *
from srcMOC.MASSCALC import *
from srcMOC.IVP2KERNAL import *				# NITISH:" Check if this is still REQUIRED. I THink it is not Required
from srcMOC.IO import *
from copy import deepcopy
## END: Initializing Packages

DIR				= os.getcwd()+'/'
t0              = time.time()
try:    INFile  = DIR+sys.argv[1]
except: INFile  = DIR+'MOC_Config.in'

IN              = ReadUserInput(INFile)
GasEqu          = IN['GAS_EQU']
rho_d           = float(IN['rho_d'])
Noz_Mach        = float(IN['Noz_Design_Mach'])
n_reflex        = float(IN['n_ref'])
y_t             = float(IN['y_t'])
n               = int(IN['n'])
y               = ()
y1              = ()
Sauer           = ()
IVP             = ()
IVP_s           = ()
IVP_loc         = ()
IVP1            = ()
IVP_New         = ()
Nozzle_Wall     = []
IVP_MAIN        = []						# 0.a Main LIST (LIST containing list of Objects of each IVP)
FileName_NCoods = DIR+IN['File_Name_NCods']
FileName_NProp  = DIR+IN['File_Name_NProp']

try: 
	File_Name_Plot = IN['File_Name_Plot']
	PlotResults = True
	pltFrmt = File_Name_Plot.split('.')[1]
except KeyError: PlotResults = False


#---------------------------------------------------------------------------------------------#
## 1. START: Sauer Analysis
print('Start: Sauer Analysis\n')
y = np.linspace(-y_t,y_t,(2*(n)-1))			# 1.a Initializing "y" on Sauer Line
Real_Sauer = Points()
if (GasEqu=='RefProp'or GasEqu=='RefProp_TAB' or GasEqu=="CoolProp"): Real_Sauer = Sauer_RealGas(Real_Sauer)

for i in y:
	P=deepcopy(Real_Sauer)
	P.y=i
	Sauer = Sauer + tuple([P])				# 1.b Creating Objects for Sauer Analysis
	
Pcount=0
for P in Sauer:								# 1.c Computing properties using Sauer Analysis and Saving as IVP
	P.SauerMain()
	IVP_s = IVP_s + tuple([P])
	Pcount=Pcount+1
	printProgress(Pcount,len(Sauer),'1 of 4')
IVP = IVP + tuple([IVP_s[-1]])				# 1.d IVP point for IVP2KERNAL
print('Start: Sauer 2 IVP\n')

#---------------------------------------------------------------------------------------------#
## 2. START: Sauer 2 IVP
for j in range(0,(n-1)*2):					# 2.a Exending the Sauer to IVP
	printProgress(j,(n-1)*2-1,'2 of 4')
	k=len(IVP_s)
	IVP_New = ()
	for i in range(0,int(k)-1):			    # 2.b Interior Point from Sauer and Line after that
		pa=InteriorPoint(IVP_s[i],IVP_s[i+1],IN)
		IVP_New = IVP_New+tuple([pa])
	del(IVP_s)								# 2.c del() is used as destructor
	IVP_s = deepcopy(IVP_New)
	IVP_MAIN.append(IVP_s)
	del(IVP_New)
	IVP = IVP + tuple([IVP_s[-1]])			# 2.d Generating IVP line from Sauer

print('\n')
IVP_MAIN.append(IVP)
Nozzle_Wall.append(IVP[0])

#---------------------------------------------------------------------------------------------#
# 3. START: IVP 2 KERNAL
tau = np.linspace(0,90,90/float(IN['dtau'])+1)# 3.a Wall angles at the expansion regime (max = 90)
tau = tau [1:-1:1]
wall = Points()
print('\n\n\t***** Developing Kernal Region *****\n\n')
print('\n\tProgress\t\tAngle\t\tMa_D\t\tMa_C\t\tx\t\tMass\t\tClockTime\n')
c=0
for ang in tau:
	c=c+1
	IVP_New = ()
	pt3=deepcopy(IVP[0])
	pt1 = deepcopy(IVP[1])
	(wall.x,wall.y) = Rot_abt_center(y_t,ang,rho_d)	
											# 3.b Finding points on the Nozzle Wall
	flag = 0				
	k=0
	while flag!=1:
		try:								# 3.c Try to find Nozzle Wall point property
			[wall.u,wall.v,kk] = Inv_Wall_pnt(wall,pt3,pt1,ang,IN)
			flag = 1
		except NameError:					# 3.c Else move to the next point on IVP
			k=k+1
			pt3 = deepcopy(IVP[0+k])
			pt1 = deepcopy(IVP[1+k])
	IVP_New = IVP_New + tuple([wall])		# 3.c Add Nozzle Point to List
	for j in range(0,int(len(IVP)-1)-k):	# 3.d Interior Wall point for rest of the points
		if j==0:pt3=deepcopy(wall)
		else: pt3 = deepcopy(kk)
		pt2 = deepcopy(IVP[j+1+k])
		IVP_New = IVP_New + tuple([InteriorPoint(pt2,pt3,IN)])
		kk = IVP_New[j]
	pt2 = deepcopy(IVP_New[-1])				# 3.e Calculating center Point
	pt3 = deepcopy(IVP_New[-1])
	pt3.v = -pt3.v
	pt3.y = -pt3.y
	IVP_New = IVP_New + tuple([InteriorPoint(pt3,pt2,IN)])
	IVP_MAIN.append(IVP_New)
	del(IVP)
	IVP = deepcopy(IVP_New)
	Nozzle_Wall.append(IVP[0])
	t2 = time.time()
	if c%10 == 0:print('\n\tProgress\t\tAngle\t\tMa_D\t\tMa_C\t\tx\t\tMass\t\tClockTime')
	printProgress(IVP_New[-1].M-1,Noz_Mach-1,'3 of 4','',2,15)
	print("\t\t\t\t%2.2f\t\t%1.2f\t\t%1.2f\t\t%1.4f\t\t%.4f\t\t%.2f\t\t" %(ang,Noz_Mach,IVP_New[-1].M,IVP_New[-1].x,MassFlowCalc(IVP_New,IN)/2,t2-t0,))
	
	if IVP_New[-1].M > Noz_Mach:
		print('\n\n\t***** Kernal Region Done *****\n\n')
		break
	del(IVP_New)

#---------------------------------------------------------------------------------------------#
## 4. START: Reflex Zone
print('\nStart: Reflex Zone\n')
# A: Finding Boundary
length = []
for i in range(1,len(IVP)):
	Ref_L = IVP[0:i+1]						# 4.a Mass in Line
	Ref_P = tuple([IVP[i]])					# 4.b Mass out point
	MassIn = MassFlowCalc(Ref_L,IN)
	MassOut = MassFlowCalc( Ref_P,IN)
	l = MassIn/MassOut						# 4.c Length of Nozzle
	LocWall =  Ref_P[0].ReflexLine(1,l)
	length.append(l)
	Nozzle_Wall.append(LocWall)
	printProgress(i,len(IVP)-1,'4A of 4')
	# Added by Jozef to change plot settings
	if (PlotResults and imprtplt):
		plt.figure(0)
		plt.plot(LocWall.x,LocWall.y,'.k')

	
Nozzle_Wall = [Points(0,0)] + Nozzle_Wall + [Points(Nozzle_Wall[-1].x,0)]
WriteNozzleDim(Nozzle_Wall,FileName_NCoods)	# 4.d Open Write and close Nozzle Wall
# B: Finding Properties
l = max(length)/n_reflex
n=1

while n!=n_reflex+1:
	Ref_Line = ()
	LN_p = IVP[-1].ReflexLine(n,l)			# 4.e Calculate 1st Point at center
	Ref_Line = Ref_Line + tuple([InteriorPoint(IVP[-1],LN_p,IN,1)])
	for Pt in IVP[-2:0:-1]:					# 4.f Calculate rest of the line in a backward fashion
		Ref_Line = Ref_Line + tuple([InteriorPoint(Pt,Ref_Line[-1],IN,1)])
		if not Ref_Line[-1].PointInPoly(Nozzle_Wall):
			IVP_MAIN.append(Ref_Line)
			break
	n=n+1
	printProgress(n,n_reflex+1,'4B of 4')

#---------------------------------------------------------------------------------------------#
## 5. START: Writing data file				# 5.a Write Data File
try:os.remove(os.getcwd()+'/'+FileName_NProp)
except: print ("Continue: No File to Delete <Prop File> \n Location:",os.getcwd())
WriteDataFile('w',IVP_MAIN,FileName_NProp)
WriteDataFile('a',IVP_MAIN,FileName_NProp)
print('\nCompleted in %f seconds\n'%(time.time() - t0))
# Added by Jozef to toggle plotting option
if (PlotResults and imprtplt): plt.savefig(DIR+File_Name_Plot,format=pltFrmt,dpi=300) #plt.show()
##
## END
##---------------------------------------------------------------------------------------------#
