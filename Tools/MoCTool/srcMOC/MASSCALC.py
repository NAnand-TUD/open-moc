#####################FILE NAME: MASSCALC.py#####################
#================================================================
#================================================================
# author: Nitish Anand    & Jozef Stuijt                         |
#     :Master Student,                                           |
#    :Process and Energy Departmemt,                             |
#    :TU Delft,                                                  |
#    :The Netherlands                                            |
#                                                                |
# email : nitish.ug2@gmail.com                                   |
#                                                                |
# Description: This Module is responsible for calculating the    |
#	Mass in the Nozzle. Mass conservation is checked using       |
#	these functions. Other supporting function are also          |
#	included in this section.				                     |
# MASSCALC MODULE:						                         |
#	: a. MassFlowCalc()					                         |
#	: b. DistCalc()						                         |
#	: c. UnitVec()						                         |
#	: d. vecDotVec()					                         |
#	: e. MakeVec()						                         |
#	: f. RotVec()						                         |
#	: g: Intersection()					                         |
#	: h: Inter()						                         |
#================================================================

import pdb
import math as calc
import sys

#---------------------------------------------------------------------------------------------#
#
# Mass Flow Calculator 
#
def MassFlowCalc(IO,IN):
    flag = 1
    m=len(IO)
    mass = 0
    NozTyp = str(IN['NozzleType'])
    if NozTyp == 'PLANAR':
        delta = 0
    else:
        delta = 1
    if m==1:
        IO[0].ThermoProp(IN)
        M = IO[0].M
        ang = (calc.asin(1/M)+calc.atan(IO[0].v/IO[0].u)-(calc.pi/2))
        V =  vecDotVec(IO[0].u,IO[0].v,calc.cos(ang),calc.sin(ang))
        mass = mass + IO[0].rho*V
    else:
        for i in range(1,m):
            l=DistCalc(IO[i-flag].x,IO[i-flag].y,IO[i].x,IO[i].y)
            try:vec =MakeVec(IO[i-flag].x,IO[i-flag].y,IO[i].x,IO[i].y)
            except: pdb.set_trace()
            rotVec =  RotVec(vec,-90)
            IO[i-flag].ThermoProp(IN)
            IO[i].ThermoProp(IN)
            rho = (IO[i-flag].rho + IO[i].rho)/2
            u = (IO[i-flag].u + IO[i].u)/2
            v = (IO[i-flag].v + IO[i].v)/2
            x = (IO[i-flag].x + IO[i].x)/2
            y = (IO[i-flag].y + IO[i].y)/2
            V =  vecDotVec(u,v,rotVec[0],rotVec[1])
            if delta:
                R = calc.sqrt(x**2 + y**2)
                A = calc.pi*2*R
            else:
                A = 1
            mass = mass + (rho * l * V * A)
    return mass
#
#FUNCTION HERE AFTER ARE SUPPORTING FUNCTION TO CALCULATE MASS
# Calculates, distance betweek points
#

def DistCalc(x1,y1,x2,y2):
    import math as calc
    dist = calc.sqrt((x2-x1)**2+(y2-y1)**2)
    return dist

#---------------------------------------------------------------------------------------------#
#
# Finding Unit Vector
#

def UnitVec(u,v):
    mag = calc.sqrt(u**2+v**2)
    u1 = u/mag
    v1 = v/mag
    return [u1,v1]

#---------------------------------------------------------------------------------------------#
#
# Finding dot Product
#

def vecDotVec(u0,u1,x0,x1):
    ret = u0*x0+u1*x1
    return ret

#---------------------------------------------------------------------------------------------#
#
# Converting cordinates into vector
#

def MakeVec(x1,y1,x2,y2):
    a = [None]*2
    a[0]=x1-x2
    a[1]=y1-y2
    mag = calc.sqrt(a[0]**2+a[1]**2)
    a[0]=a[0]/mag
    a[1]=a[1]/mag
    return a

#---------------------------------------------------------------------------------------------#
#
# Rotate Vector
#

def RotVec(x,angle):
    a = [None]*2
    a[0] = x[0]*calc.cos(calc.radians(angle))-x[1]*calc.sin(calc.radians(angle))
    a[1] = x[0]*calc.sin(calc.radians(angle))+x[1]*calc.cos(calc.radians(angle))
    return a

#---------------------------------------------------------------------------------------------#
#
# Intersection(): Finds position of the intersection of two lines defined by x,y and slope
#

def Intersection(x1,x2,m):
    x3 = (m[0]*x1[0] - m[1]*x2[0] + x2[1] - x1[1])/(m[0]-m[1])
    y3 = m[0]*(x3-x1[0])+x1[1]
    return x3,y3

#---------------------------------------------------------------------------------------------#
#
# Inter(): Liner Interpolation function
#

def Inter(x1,x2,y1,y2,x):
    k=y1+(x-x1)*((y2-y1)/(x2-x1));
    return(k);

#---------------------------------------------------------------------------------------------#
##
## END
##---------------------------------------------------------------------------------------------#
