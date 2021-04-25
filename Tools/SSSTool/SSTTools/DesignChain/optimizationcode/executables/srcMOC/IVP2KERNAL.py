###############FILE NAME:IVP2KERNEL.py###########################
#================================================================
#================================================================
# author: Nitish Anand    & Jozef Stuijt                         |
#     :Master Student,                                           |
#    :Process and Energy Departmemt,                             |
#    :TU Delft,                                                  |
#    :The Netherlands                                            |
#                                                                |
# email : nitish.ug2@gmail.com                                   |
# 								                                 |
# Description: Initial Expansion is defined with a circle:       |
#	     : Function in this module calculates the point on       |
#	     : the cicle                                             |
# IVP2KERNAL MODULE:                                             |
#	:Rot_abt_center()                                            |
#================================================================

import math as calc
import pdb
import numpy as np

#---------------------------------------------------------------------------------------------#
#
# Rotates point about a center (Kernel Regions)
#
def Rot_abt_center(y_t,ang,rho_d):
    x = [0,y_t]
    center = [0,rho_d+y_t]
    x=np.array(x)-np.array(center)
    ang = calc.radians(ang)
    M = np.matrix([[calc.cos(ang), -calc.sin(ang)],[calc.sin(ang),calc.cos(ang)]])
    x = np.transpose(np.matrix(x))
    rot_x = np.matrix(M) * np.matrix(x)
    rot_x = np.transpose(rot_x)
    rot_x = np.mat(rot_x) + center    
    return rot_x[0,0],rot_x[0,1]

##
## END




#---------------------------------------------------------------------------------------------#
