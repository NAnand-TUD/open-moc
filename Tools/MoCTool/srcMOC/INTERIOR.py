##################FILE NAME: INTERIOR.py########################
#================================================================
# author: Nitish Anand    & Jozef Stuijt                         |
#     :Master Student,                                           |
#    :Process and Energy Departmemt,                             |
#    :TU Delft,                                                  |
#    :The Netherlands                                            |
#                                                                |
# email : nitish.ug2@gmail.com                                   |
# 								                                 |
# Description: MOC Routine: position and velocity calculator. 	 |
#	     : Refer Thesis or Zucrow & Hoffman for the 	         |
#	     : derivation of equation				                 |
# MODULE Function:						                         |
#	a. InteriorPoint()					                         |
#	b. InteriorPoint_Prop()					                     |
#	c. AxialPoint_Prop()					                     |
#	d. Inv_Wall_pnt()				                             |
#	e. FlowProp()						                         |
#	f. EulerPredictorCorrector()				                 |
#	g. Direct_Wall_pnt()					                     |
#================================================================
#
# ShortHands
#	EPC: Euler Predictor Corrector
#	PisS: PLUS is simple
#	MisS: MINUS is simple

from srcMOC.CLASSES import *
import time
import pdb
import math as calc
from srcMOC.MASSCALC import *
try:                                       # check if import possible related to cluster
   import matplotlib.pyplot as plt         # for plotting routines
   imprtplt = True
except ImportError: 
   print('\nWarining: Plotting module import unsuccesful\n')
   imprtplt = False
#import matplotlib.pyplot as plt
import sys

#---------------------------------------------------------------------------------------------#
#
# Main code to compute the flow property of the third point from two points(P:Plus,M:Minus)
#	includes euler corrector predictor step
#

def InteriorPoint(paramP1,paramM1,IN,PisS=0,MisS=0,EPC1=1):  
    # WORKS COMPLETELY FINE: COMPLIMETS THE SOLUTION OF BOOK Ex: 12.5 ZH-Vol2

    # Added to read display
    #try:INFile = sys.argv[1]
    #except:INFile = 'MOC_Config.in'
    #IN     = ReadUserInput(INFile)
    #disp   = IN['PlotResults']
    try:  
        IN['File_Name_Plot']
        disp = True
    except KeyError: disp = False
    #disp = IN['PlotResults']
    #if (disp!='YES' and disp!='NO'):
    #    raise IOError('Invalid Input: PlotResults')
    #    sys.exit()

    flag = 0;
    #Check if Point on AXIS
    if abs(abs(paramP1.y)-abs(paramM1.y))<1e-5:
        param1=AxisPoint_Prop(paramP1,paramM1,IN)
    else:
        param1=InteriorPoint_Prop(paramP1,paramM1,IN,'p',PisS,MisS)
    # Check: Euler Predictor Corrector Yes/No
    if EPC1==1:
        param2 = EulerPredCorr(paramP1,paramM1,param1,IN)
    else:
        param2 = param1
    # PLOTTING ROUTINE: Red: Positive Characteristics, Green: Negative Characterisitcs
    # Jozef: Added plot toggle
    #if (disp == 'YES'):
    if (disp and imprtplt):
        plt.figure(0)
        if (paramP1.y>0 or paramM1.y>0 or param2.y>0):
            plt.plot([paramP1.x,param2.x],[paramP1.y,param2.y],'r')
            plt.plot([paramM1.x,param2.x],[paramM1.y,param2.y],'g')
    
    return param2

#---------------------------------------------------------------------------------------------#
#
# calcualtes the flow properties of the 3rd point (intersection of 2 characterisitcs)
#

def InteriorPoint_Prop(paramP,paramM,IN,flag='c',PisS=0,MisS=0):
    NozTyp = str(IN['NozzleType'])
    param = Points()
    paramP.ThermoProp(IN)
    # For now always PLANAR
    if NozTyp == 'PLANAR':
        delta = 0
    else:
        delta = 1
    paramP.theta = calc.atan(paramP.v/paramP.u)
    try:paramP.alpha = calc.asin(paramP.a/paramP.V)
    except:
        #print('*** Error: Mach less than 1 ***\nCHECK INPUT VARIABLE')
        #pdb.set_trace()
        raise ValueError('Mach less than 1 ***\nCHECK INPUT VARIABLE')
    #Check if PLUS Characterisitcs is Simple Wave
    if PisS == 0:
        paramP.lam = calc.tan(paramP.theta+paramP.alpha)
    else:
        paramP.lam = calc.tan(0+paramP.alpha)
    paramP.Q = paramP.u**2 - paramP.a**2
    paramP.R1 = (2*paramP.u*paramP.v) - (paramP.Q*paramP.lam)
    #If Axial Point
    if abs(paramP.y)<=1e-9:
        paramP.S = 0
    else:
        paramP.S = delta*(paramP.a**2 * paramP.v)/paramP.y
    paramM.ThermoProp(IN)
    paramM.theta = calc.atan(paramM.v/paramM.u)
    paramM.alpha = calc.asin(paramM.a/paramM.V)
    #Check if MINUS Characterisitcs is Simple Wave
    if MisS == 0:
        paramM.lam = calc.tan(paramM.theta-paramM.alpha)
    else:
        paramM.lam = calc.tan(0-paramM.alpha)
    paramM.Q = paramM.u**2 - paramM.a**2
    paramM.R1 = (2*paramM.u*paramM.v) - (paramM.Q*paramM.lam)
    #If Axial Point
    if abs(paramM.y)<=1e-9:
        paramM.S=0
    else:
        paramM.S = delta*(paramM.a**2 * paramM.v)/paramM.y
    #Calculate Interior Point using P and M Point
    if abs(paramM.y)<=1e-9:
        param.x=paramM.x
        param.y = (paramP.lam*(param.x-paramP.x))+paramP.y
    elif abs(paramP.y)<=1e-9:
        param.x=paramP.x
        param.y=(paramM.lam*(param.x-paramM.x))+paramM.y
    else:
        [param.x,param.y] = Intersection([paramM.x,paramM.y],[paramP.x,paramP.y],[paramM.lam,paramP.lam])
    #Calculating U and V for solution point
    paramP.T = paramP.S*(param.x-paramP.x)+(paramP.Q*paramP.u)+(paramP.R1*paramP.v)
    paramM.T = paramM.S*(param.x-paramM.x)+(paramM.Q*paramM.u)+(paramM.R1*paramM.v)
    if abs(paramM.Q)<=1e-9: 
        param.u = paramP.u
        param.v = 0    
    elif abs(paramP.Q)<=1e-9:
        param.u = paramM.u
        param.v = 0
    else:
        if abs(paramP.y)<=1e-9:
            param.v=0
        elif abs(paramM.y)<=1e-9:
            param.v=0
        else:
            if abs(((paramP.T*paramM.Q)-(paramM.T*paramP.Q)))<=1e-9:
                param.v = 0        
            elif abs(((paramM.T*paramP.Q)-(paramP.T*paramM.Q)))<=1e-9:
                param.v=0
            else:
                param.v = ((paramP.T*paramM.Q)-(paramM.T*paramP.Q))/((paramM.Q*paramP.R1)-(paramP.Q*paramM.R1))
        param.u = (paramM.T - (paramM.R1*param.v))/paramM.Q
    #Return x,y,u,v for solution point
    return param

#---------------------------------------------------------------------------------------------#
#
# Calculates the flow properties of the Axis Point
#

def AxisPoint_Prop(paramP,paramM,IN):
    NozTyp = str(IN['NozzleType'])
    if NozTyp == 'PLANAR':
        delta = 0
    else:
        delta = 1
    param = Points()
    # PLUS Characteristics
    paramP.ThermoProp(IN)
    paramP.theta = calc.atan(paramP.v/paramP.u)
    paramP.alpha = calc.asin(paramP.a/paramP.V)
    paramP.lam = calc.tan(paramP.theta+paramP.alpha)
    paramP.Q = paramP.u**2 - paramP.a**2
    paramP.R1 = (2*paramP.u*paramP.v) - (paramP.Q*paramP.lam)
    paramP.S = delta*(paramP.a**2 * paramP.v)/paramP.y

    # MINUS Characteristics
    paramM.ThermoProp(IN)
    paramM.theta = calc.atan(paramM.v/paramM.u)
    paramM.alpha = calc.asin(paramM.a/paramM.V)
    paramM.lam = calc.tan(paramM.theta-paramM.alpha)
    paramM.Q = paramM.u**2 - paramM.a**2
    paramM.R1 = (2*paramM.u*paramM.v) - (paramM.Q*paramM.lam)
    paramM.S = delta*(paramM.a**2 * paramM.v)/paramM.y

    # Calculating x , y , u , v for AXIAL Solution Point
    [param.x,param.y] =Intersection([paramM.x,paramM.y],[paramP.x,paramP.y],[paramM.lam,paramP.lam])
    param.y=0
    paramM.T = paramM.S*(param.x-paramM.x)+(paramM.Q*paramM.u)+(paramM.R1*paramM.v)
    param.v=0
    param.u = (paramM.T)/paramM.Q
    return param

#---------------------------------------------------------------------------------------------#
#
# Inverse Wall Point calcualtor for the inital Nozzle expansion area [Kernel Region]
#

def Inv_Wall_pnt(wall,pt33,pt11,wall_ang,IN):
    Corr_v=float(IN['tolerance_v'])
    flagv=0
    flagu=0
    pt2 =Points()
    pt4 = Points()
    pt44 = Points()
    pt44.x = pt33.x
    #Determining PLUS/MINUS Characterisitcs
    if (pt33.y)<0:
        C = 'n'
    else:
        C = 'p'
    # Ref: Zucrow and Hoffman
    #Finding point on line pt1 & pt3
    #Convergence Criteria : du $ dv < Corr_v
    while (abs(pt2.v-flagv)<Corr_v) & (abs(pt2.u-flagu)<Corr_v):     
        lam31 = (pt33.y-pt11.y)/(pt33.x-pt11.x)
        pt2.u = pt33.u
        pt2.v = pt33.v  
        lam42 =FlowProp(pt2,1,C,IN)
        [pt2.x,pt2.y]=Intersection([pt33.x,pt33.y],[wall.x,wall.y],[lam31,lam42])
	#Rasing Error: Important to Move On to next Point in the MAIN.py(Ln:129)
        if pt2.x > pt11.x:
             raise NameError('BeyondPointp33xIWP')
        flagv=pt2.v
        flagu=pt2.u
        #Interpolate the U and V value
        pt2.u=Inter(pt33.x,pt11.x,pt33.u,pt11.u,pt2.x)
        pt2.v=Inter(pt33.x,pt11.x,pt33.v,pt11.v,pt2.x)
    #Now Finding properties at the WALL(4)
    pt4.m = lam42
    pt4.x = wall.x
    pt4.y = wall.y
    # Special FlowProp Calculator (since, Not all value has to be calcualted)
    pt2 = FlowProp(pt2,2,C,IN)
    if C=='p':
        factor = 1
    else:
        factor = -1
    T = pt2.S*(pt4.x-pt2.x)+(pt2.Q*pt2.u)+(pt2.R1*pt2.v)
    pt4.u= T/(pt2.Q+(pt2.R1*calc.tan(calc.radians(factor*wall_ang))))
    pt4.v=calc.tan(calc.radians(factor*wall_ang))*pt4.u
    flagu=0
    flagv=0
    #RETURN: Wall velocities
    return pt4.u,pt4.v,pt2

#---------------------------------------------------------------------------------------------#
#
# Calcualtes flow properties (limited)
#	- Used in Inv Wall Points where calc of all prop is not required

def FlowProp(paramP,flag,Ctype,IN):
    paramP.ThermoProp(IN)
    NozTyp = str(IN['NozzleType'])
    if NozTyp == 'PLANAR':
        delta = 0
    else:
        delta = 1
    paramP.theta = calc.atan(paramP.v/paramP.u)
    paramP.alpha = calc.asin(paramP.a/paramP.V)
    if Ctype == 'p':
        paramP.lam = calc.tan(paramP.theta+paramP.alpha)
    else:
        paramP.lam = calc.tan(paramP.theta-paramP.alpha)
    if flag == 1:
        return paramP.lam
    elif flag ==2:
        paramP.Q = paramP.u**2 - paramP.a**2
        paramP.R1 = (2*paramP.u*paramP.v) - (paramP.Q*paramP.lam)
        if abs(paramP.v)<1e-9:
            paramP.S = 0
        else:
            paramP.S = delta*(paramP.a**2 * paramP.v)/paramP.y      
    return paramP

#---------------------------------------------------------------------------------------------#
#
# UNDER INVESTIGATION: Scrapped Code!!
# NOT USED IN MAIN CODE

def Direct_Wall_pnt(wall1,pt2d,wall_ang,IN):
    import math as calc
    import EulerPredCorr as EPC
    import matplotlib.pyplot as plt
    Corr_v=float(IN['tolerance_v'])
    pt4 =Points()
    if (pt2d.y)<0:
        C = 'n'
    else:
        C = 'p' 
    lam42 = FlowProp(pt2d,2,C,IN)
    pt4.m = lam42
    pt4.x = wall1.x
    pt4.y = wall1.y
    if C=='p':
        factor = 1
    else:
        factor = -1
    T = pt2d.S*(pt4.x-pt2d.x)+(pt2d.Q*pt2d.u)+(pt2d.R1*pt2d.v)
    pt4.u= T/(pt2d.Q+(pt2d.R1*calc.tan(calc.radians(factor*wall_ang))))
    pt4.v=calc.tan(calc.radians(factor*wall_ang))*pt4.u
    Corr = int(IN['Corr_n_inv'])
    plt.figure(0)
    plt.title("Figure 0")
    plt.plot(pt4.x,pt4.y,'>r')
    return pt4.u,pt4.v,pt2d
#---------------------------------------------------------------------------------------------#
#
# EulerPredCorr()
#	-Takes care of the corrector step in EPC

def EulerPredCorr(paramP,paramM,param6,IN):
    #Property Calculator for InteriorPoint[Corrector]
    #-----------------------------------------------------------------------------------------#
    def LocalProp(paramP,paramM,P,M):
        NozTyp = str(IN['NozzleType'])
        if NozTyp == 'PLANAR':
            delta = 0
        else:
            delta = 1
        param =Points()

	#PLUS Characterisitcs
        paramP.ThermoProp(IN)
        paramM.ThermoProp(IN)
        paramP.theta = calc.atan(paramP.v/paramP.u)
        paramP.alpha = calc.asin(paramP.a/paramP.V)
        paramP.lam = calc.tan(paramP.theta+paramP.alpha)
        paramP.Q = paramP.u**2 - paramP.a**2
        paramP.R1 = (2*paramP.u*paramP.v) - (paramP.Q*paramP.lam)
        if abs(paramP.y)<=1e-9:
            paramP.S = 0
        else:
            paramP.S = delta*(paramP.a**2 * paramP.v)/paramP.y

	#MINUS Characterisitcs
        paramM.theta = calc.atan(paramM.v/paramM.u)
        paramM.alpha = calc.asin(paramM.a/paramM.V)
        paramM.lam = calc.tan(paramM.theta-paramM.alpha)
        paramM.Q = paramM.u**2 - paramM.a**2
        paramM.R1 = (2*paramM.u*paramM.v) - (paramM.Q*paramM.lam)
        if abs(paramM.y)<=1e-9:
            paramM.S=0
        else:
            paramM.S = delta*(paramM.a**2 * paramM.v)/paramM.y

	#Intersection of Points
        if abs(paramM.y)<=1e-9:
            param.x=paramM.x
            param.y = (paramP.lam*(param.x-paramP.x))+paramP.y
        elif abs(paramP.y)<=1e-9:
            param.x=paramP.x
            param.y=(paramM.lam*(param.x-paramM.x))+paramM.y
        else:
            [param.x,param.y] =Intersection([M.x,M.y],[P.x,P.y],[paramM.lam,paramP.lam])
            paramP.T = paramP.S*(param.x-paramP.x)+(paramP.Q*P.u)+(paramP.R1*P.v)
            paramM.T = paramM.S*(param.x-paramM.x)+(paramM.Q*M.u)+(paramM.R1*M.v)
            if abs(paramM.Q)<=1e-9:
                param.u = P.u
                param.v = 0    
            elif abs(paramP.Q)<=1e-9:
                param.u = M.u
                param.v = 0
            else:
                if abs(paramP.y)<=1e-9:
                    param.v=0
                elif abs(paramM.y)<=1e-9:
                    param.v=0
                else:
                    if abs(((paramP.T*paramM.Q)-(paramM.T*paramP.Q)))<=1e-9:
                        param.v = 0        
                    elif abs(((paramM.T*paramP.Q)-(paramP.T*paramM.Q)))<=1e-9:
                        param.v=0
                    else:
                        param.v = ((paramP.T*paramM.Q)-(paramM.T*paramP.Q))/((paramM.Q*paramP.R1)-(paramP.Q*paramM.R1))
                param.u = (paramM.T - (paramM.R1*param.v))/paramM.Q
        return param
    #Property Calculator for AxisPoint[Corrector]
    #-----------------------------------------------------------------------------------------#
    def LocalAxisProp(paramP,paramM,P,M):
        NozTyp = str(IN['NozzleType'])
        if NozTyp == 'PLANAR':
            delta = 0
        else:
            delta = 1
        param = Points()
	
	#PLUS Characterisitcs
        paramP.ThermoProp(IN)
        paramM.ThermoProp(IN)
        paramP.theta = calc.atan(paramP.v/paramP.u)
        paramP.alpha = calc.asin(paramP.a/paramP.V)
        paramP.lam = calc.tan(paramP.theta+paramP.alpha)
        paramP.Q = paramP.u**2 - paramP.a**2
        paramP.R1 = (2*paramP.u*paramP.v) - (paramP.Q*paramP.lam)
        paramP.S = delta*(paramP.a**2 * paramP.v)/paramP.y

        #MINUS Characteristics
        paramM.theta = calc.atan(paramM.v/paramM.u)
        paramM.alpha = calc.asin(paramM.a/paramM.V)
        paramM.lam = calc.tan(paramM.theta-paramM.alpha)
        paramM.Q = paramM.u**2 - paramM.a**2
        paramM.R1 = (2*paramM.u*paramM.v) - (paramM.Q*paramM.lam)
        paramM.S = delta*(paramM.a**2 * paramM.v)/paramM.y

	#Intersection of Points
        [param.x,param.y] =Intersection([M.x,M.y],[P.x,P.y],[paramM.lam,paramP.lam])
        paramP.T = paramP.S*(param.x-paramP.x)+(paramP.Q*P.u)+(paramP.R1*P.v)
        paramM.T = paramM.S*(param.x-paramM.x)+(paramM.Q*M.u)+(paramM.R1*M.v)
        param.v = 0
        param.u = (paramM.T - (paramM.R1*param.v))/paramM.Q
        return param
    #MAIN CODE START HERE
    Corr_n = 1000
    tol_x = float(IN['tolerance_x'])
    tol_v = float(IN['tolerance_v'])
    paramP2 =Points()
    paramM2 =Points()
    param5 = Points()
    for i in range(0,int(Corr_n)):
	#[CHECK]: velocity Convergence
        if (abs(param5.y-param6.y)<tol_x)&(abs(param5.u-param6.u)<tol_v)&(abs(param5.v-param6.v)<tol_v):
            break
        else:
            param5 = param6
            # Finding average properties [PLUS]
            paramP2.x = paramP.x
            paramP2.u = (paramP.u+param5.u)/2
            paramP2.v = (paramP.v+param5.v)/2
            paramP2.y = (paramP.y+param5.y)/2
            # Finding average properties [MINUS]
            paramM2.x = paramM.x
            paramM2.u = (paramM.u+param5.u)/2
            paramM2.v = (paramM.v+param5.v)/2
            paramM2.y = (paramM.y+param5.y)/2
	    #[CHECK]: Spacial Convergence
            if abs(abs(paramP.y)-abs(paramM.y))<1e-5:
                param6 = LocalAxisProp(paramP2,paramM2,paramP,paramM)
            else:
                param6 = LocalProp(paramP2,paramM2,paramP,paramM)
    param6.ThermoProp(IN)
    return param6

##
## END
##
