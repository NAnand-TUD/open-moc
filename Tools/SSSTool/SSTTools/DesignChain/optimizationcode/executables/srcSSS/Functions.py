#! /usr/bin/python

import numpy as np
import math
import sys

#---------------------------------------------------------------------------------------------#
# HELP FUNCTIONS #

def UnitVec(u,v):
	mag = math.sqrt(u**2+v**2)
	u1 = u/mag
	v1 = v/mag
	return [u1,v1]
#-----------------------------------------------------------------------------------------------#
#
#
#
# Liner Interpolation function for Table
def Inter(X,x1,x2,y1,y2):
    return(y1+(X-x1)*((y2-y1)/(x2-x1)))
#-----------------------------------------------------------------------------------------------#
#
#
#
def Rotate_about_center(y_t,ang,rho_d):
	x = [y_t[0],y_t[1]]
	center = [rho_d[0],rho_d[1]]
	x=np.array(x)-np.array(center)
	ang = math.radians(ang)
	M = np.matrix([[math.cos(ang), -math.sin(ang)],[math.sin(ang),math.cos(ang)]])
	x = np.transpose(np.matrix(x))
	rot_x = np.matrix(M) * np.matrix(x)
	rot_x = np.transpose(rot_x)
	rot_x = np.mat(rot_x) + center    
	return rot_x[0,0],rot_x[0,1], 0.0

# Intersection(): Finds position of the intersection of two lines defined by x,y and slope
#-----------------------------------------------------------------------------------------------#
#
#
#
def Intersection(x1,x2,m):
	x3 = (m[0]*x1[0] - m[1]*x2[0] + x2[1] - x1[1])/(m[0]-m[1])
	y3 = (m[0]*(x3-x1[0]))+x1[1]
	return x3,y3
#-----------------------------------------------------------------------------------------------#
#
#
#
def LinePtsInter(p1,p2,p3,p4):
	m1,k1 = Line2Pts(p1,p2)
	m2,k2 = Line2Pts(p3,p4)
	x = (k2-k1)/(m1-m2)
	y = m1*x+k1
	return [x,y]
#-----------------------------------------------------------------------------------------------#
#
#
#
# Find the gradient and y-intercept of line between two points
def Line2Pts(p1,p2):
	if (p1[0]>p2[0]): m = (p1[1]-p2[1])/(p1[0]-p2[0])
	else: m = (p2[1]-p1[1])/(p2[0]-p1[0])
	k = p1[1]-m*p1[0]
	return m,k
#-----------------------------------------------------------------------------------------------#
#
#
#
# Find the intersection of line between two points and circle of radius r
def InterLineCirc(p1,p2,r):
	[m,k] = Line2Pts(p1,p2)
	A = m**2+1
	B = 2*m*k
	C = k**2-r**2
	x1 = (-B+(B**2-4*A*C)**0.5)/(2*A)
	y1 = x1*m+k
	x2 = (-B-(B**2-4*A*C)**0.5)/(2*A)
	y2 = x2*m+k
	p3 = [x1,y1]
	p4 = [x2,y2]
	return p3,p4

def FindMidCurve(array):
	length = []
	length.append(0)
	for i in range(len(array)):
		if (i==0): pass
		else:
			d = ((array[i,0]-array[i-1,0])**2+(array[i,1]-array[i-1,1])**2)**0.5
			length.append(d+length[i-1])
	return np.argmin(np.abs([length-(length[-1]/2)]))
#-----------------------------------------------------------------------------------------------#
#
#
#
# find point at x/L (0<x<1) length of curve
def FindXCurve(array,x):
	length = []
	length.append(0)
	for i in range(len(array)):
		if (i==0): pass
		else:
			d = ((array[i,0]-array[i-1,0])**2+(array[i,1]-array[i-1,1])**2)**0.5
			length.append(d+length[i-1])
	return np.argmin(np.abs([length-(length[-1]*x)]))

def closest_point(pt, array):
	deltas = array - pt
	dist_2 = np.einsum('ij,ij->i', deltas, deltas)
	return np.argmin(dist_2)
#-----------------------------------------------------------------------------------------------#
#
#
#
# https://stackoverflow.com/questions/8954326/how-to-calculate-the-mirror-point-along-a-line#8954454
def mirror_point(pt,P1,P2):
	A = P2[1]-P1[1]
	B = -(P2[0]-P1[0])
	C = -A*P1[0] - B*P1[1]

	M = math.sqrt(A*A+B*B)

	Ap = A/M
	Bp = B/M
	Cp = C/M

	D = Ap*pt[0]+Bp*pt[1]+Cp
	return [pt[0]-2*Ap*D,pt[1]-2*Bp*D]

### --- Optimization functions
#-----------------------------------------------------------------------------------------------#
#
#
#
def OptimizeScaling(args, Blade):
	Blade.SetRadiusOut(args[0])
	Blade.SetScale(args[1])
	Blade.Build_Diverging()
	sys.stdout.write("Tip Radius: %.10f Scaling Factor: %.10f RESIDUAL: %.10f \r"%(args[0],args[1],Blade.GetDist()))
	return(Blade.GetDist())
	
