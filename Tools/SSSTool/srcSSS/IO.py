######################FILE NAME: IO.py###########################
# =============================================================================================
# author: Jozef Stuijt & Nitish Anand                                                        |
#       : Master Student,                                                                    |
#           : Process and Energy Department,                                                 |
#           : TU Delft,                                                                      |
#           : The Netherlands                                                                |
#                                                                                            |
# email : joe.stuijt@gmail.com                                                               |
#                                                                                            |
# Description: Plots, reads and writes data files                                            |
#                                                                                            |
# MODULE IO: Input/Output                                                                    |
#       :I. Read User Input                                                                  |
#       :II. File Writting                                                                   |
#       :a.     writeCoords                                                                  |
#           :b.	writeUMG2out                                                                 |
#           :b.	writeBladeSpecs	                                                             |
#       :III. Plotting                                                                       |
#               :a. RadiusLimitPlot                                                          |
#               :b. plotStatorStage                                                          |
#               :c. plotBoundaries                                                           |
# =============================================================================================

import re
import sys
import os
import math
import pdb
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from srcSSS.Functions import *


# ---------------------------------------------------------------------------------------------#
#
# User Input file 'MOC_Config.in' <Sample Input file>
#
# ---------------------------------------------------------------------------------------------#


def ReadUserInput(name):
    IN = {}
	# infile = open(name,'r')
    with open(name, 'r') as infile:
        for line in infile:
            words = re.split(r'=| |%|\n|#', line)
            if not any(words[0] in s for s in ['\n', '%', ' ', '#']):
                words = list(filter(None, words))
                IN[words[0]] = words[1]
    return IN


# ---------------------------------------------------------------------------------------------#
#
# Writes the vane coordinates in the text file <file_name.out>
#
# ---------------------------------------------------------------------------------------------#

def writeCoords(Vane, name):
    with open(name, 'w') as file:
        for i in Vane: file.write("%f\t%f\t%f\n" % (i[0], i[1], 0.0))


# ---------------------------------------------------------------------------------------------#
#
# Writes the geometry output file in UMG2 format
#
# ---------------------------------------------------------------------------------------------#

def writeUMG2out(stator, name, factor=1, scale=1, angle=0, ncurv=8):
    "Deep copy all the parts"
    PRESSURE = deepcopy(stator.PRESSURE)
    SUCTION = deepcopy(stator.SUCTION)
    INLET = deepcopy(stator.INLET)
    OUTLET = deepcopy(stator.OUTLET)
    PERIO1 = deepcopy(stator.PERIO1)
    PERIO2 = deepcopy(stator.PERIO2)

    "Rotate all the points to the first quadrant"
    for i in range(len(PRESSURE)): PRESSURE[i] = Rotate_about_center(PRESSURE[i], angle, [0, 0])
    for i in range(len(SUCTION)): SUCTION[i] = Rotate_about_center(SUCTION[i], angle, [0, 0])
    for i in range(len(INLET)): INLET[i] = Rotate_about_center(INLET[i], angle, [0, 0])
    for i in range(len(OUTLET)): OUTLET[i] = Rotate_about_center(OUTLET[i], angle, [0, 0])
    for i in range(len(PERIO1)): PERIO1[i] = Rotate_about_center(PERIO1[i], angle, [0, 0])
    for i in range(len(PERIO2)): PERIO2[i] = Rotate_about_center(PERIO2[i], angle, [0, 0])

    #Extend Outlet further downstream
    if stator._AXIAL:
        for i in range(len(OUTLET)):
            OUTLET[i,0] += stator._pitch
            OUTLET[i,1] += stator._pitch

        PERIO1[-1,0] +=stator._pitch
        PERIO1[-1, 1] += stator._pitch

        PERIO2[-1,0] +=stator._pitch
        PERIO2[-1, 1] += stator._pitch

    stator.TEpt = Rotate_about_center(stator.TEpt, angle, [0, 0])
    stator.TEpt = [scale * stator.TEpt[0], scale * stator.TEpt[1], 0]

    stator.thrtStrtUp = Rotate_about_center(stator.thrtStrtUp, angle, [0, 0])
    stator.thrtEndUp = Rotate_about_center(stator.thrtEndUp, angle, [0, 0])
    stator.thrtStrtD = Rotate_about_center(stator.thrtStrtD, angle, [0, 0])
    stator.thrtEndD = Rotate_about_center(stator.thrtEndD, angle, [0, 0])

    stator.nozEndUp = Rotate_about_center(stator.nozEndUp, angle, [0, 0])
    stator.nozStrtUp = Rotate_about_center(stator.nozStrtUp, angle, [0, 0])
    stator.nozStrtD = Rotate_about_center(stator.nozStrtD, angle, [0, 0])
    stator.nozEndD = Rotate_about_center(stator.nozEndD, angle, [0, 0])

    "Split the periodic boundaries for more control"
    PERIO3 = PERIO1[0:2]
    PERIO4 = PERIO2[0:2]

    if INLET[0,0] == INLET[-1,0]:
        PERIO2_n = np.delete(PERIO1, [0,2,5], axis=0)
        PERIO1_n = np.delete(PERIO2, [0,2,5], axis=0)
    else:
        PERIO2_n = np.delete(PERIO1, [0], axis=0)
        PERIO1_n = np.delete(PERIO2, [0], axis=0)
    if (len(PRESSURE) % factor == 0):
        PRESSURE_C = np.zeros(((len(PRESSURE) // factor, 3)))
    else:
        PRESSURE_C = np.zeros(((len(PRESSURE) // factor + 1, 3)))
    # pdb.set_trace()
    for i in range(len(PRESSURE)):
        if (i % factor == 0): PRESSURE_C[i // factor] = PRESSURE[i]
    PRESSURE_C[-1] = PRESSURE[-1]

    if (len(SUCTION) % factor == 0):
        SUCTION_C = np.zeros(((len(SUCTION) // factor, 3)))
    else:
        SUCTION_C = np.zeros(((len(SUCTION) // factor + 1, 3)))
    for i in range(len(SUCTION)):
        if (i % factor == 0): SUCTION_C[i // factor] = SUCTION[i]
    SUCTION_C[-1] = SUCTION[-1]

    plt.plot(INLET[:, 0], INLET[:, 1], 'g')
    plt.plot(OUTLET[:, 0], OUTLET[:, 1], 'r')
    plt.plot(PRESSURE[:, 0], PRESSURE[:, 1], 'k')
    plt.plot(SUCTION[:, 0], SUCTION[:, 1], 'k')
    plt.plot(PERIO1[:, 0], PERIO1[:, 1], 'm')
    plt.plot(PERIO2[:, 0], PERIO2[:, 1], 'm')

    # ncurv = 6
    if os.path.isdir('Db'):
        filename = 'Db/geometry.' + name
    else:
        filename = 'geometry.' + name
    curvetype = '               \' S \'\n'
    dim_np = '                 dim                 np\n'
    val_dim_np = '                   2                  '
    x_y = '\n                   x                   y\n'

    with open(filename, 'w') as f:
        # f = open(filename,'w')

        # f.write('  Number of surfaces\n                  '+str(ncurv)+'\n')
        f.write('  Number of surfaces\n                  {}\n'.format(ncurv))

        # Pressure side
        # f.write('PRESSURE SIDE\n'+curvetype+dim_np+val_dim_np+str(len(PRESSURE_C))+x_y+'\n')
        # for i in PRESSURE_C: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('PRESSURE SIDE\n' + curvetype + dim_np + val_dim_np + '{}'.format(len(PRESSURE_C)) + x_y + '\n')
        for i in PRESSURE_C: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))

        # Suction side
        # f.write('SUCTION SIDE\n'+curvetype+dim_np+val_dim_np+str(len(SUCTION_C))+x_y+'\n')
        # for i in SUCTION_C: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('SUCTION SIDE\n' + curvetype + dim_np + val_dim_np + str(len(SUCTION_C)) + x_y + '\n')
        for i in SUCTION_C: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))

        # Inflow
        # f.write('INFLOW\n'+curvetype+dim_np+val_dim_np+str(len(INLET))+x_y+'\n')
        # for i in INLET: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('INFLOW\n' + curvetype + dim_np + val_dim_np + '{}'.format(len(INLET)) + x_y + '\n')
        for i in INLET: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))

        # Periodic Boundary 2
        # f.write('PERIO2\n'+curvetype+dim_np+val_dim_np+str(len(PERIO2_n))+x_y+'\n')
        # for i in PERIO2_n: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('PERIO2\n' + curvetype + dim_np + val_dim_np + '{}'.format(len(PERIO2_n)) + x_y + '\n')
        for i in PERIO2_n: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))

        ##Periodic Boundary 3
        # f.write('PERIO3\n'+curvetype+dim_np+val_dim_np+str(len(PERIO3))+x_y+'\n')
        # for i in PERIO3: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('PERIO3\n' + curvetype + dim_np + val_dim_np + '{}'.format(len(PERIO3)) + x_y + '\n')
        for i in PERIO3: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))

        # Outflow
        # f.write('OUTFLOW\n'+curvetype+dim_np+val_dim_np+str(len(OUTLET))+x_y+'\n')
        # for i in OUTLET: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('OUTFLOW\n' + curvetype + dim_np + val_dim_np + '{}'.format(len(OUTLET)) + x_y + '\n')
        for i in OUTLET: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))

        # Periodic Boundary 4
        # f.write('PERIO4\n'+curvetype+dim_np+val_dim_np+str(len(PERIO4))+x_y+'\n')
        # for i in PERIO4: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('PERIO4\n' + curvetype + dim_np + val_dim_np + '{}'.format(len(PERIO4)) + x_y + '\n')
        for i in PERIO4: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))

        # Periodic Boundary 1
        # f.write('PERIO1\n'+curvetype+dim_np+val_dim_np+str(len(PERIO1_n))+x_y+'\n')
        # for i in PERIO1_n: f.write("%e\t%e\n" %(i[0]*scale,i[1]*scale))
        f.write('PERIO1\n' + curvetype + dim_np + val_dim_np + '{}'.format(len(PERIO1_n)) + x_y + '\n')
        for i in PERIO1_n: f.write('{0:e}\t{1:e}\n'.format(i[0] * scale, i[1] * scale))


# f.close()

# ---------------------------------------------------------------------------------------------#
#
# Writes file with important blade specifications
#
# ---------------------------------------------------------------------------------------------#

def writeBladeSpecs(RSS, ScaleNoz, name):
    with open(name, 'w') as f:
        # f = open(name,'w')
        f.write('### Stator Blade Specifications ###\n\n\n')
        # f.write('# Warnings during blade optimization procedure\nWarning = '+str(warn)+'\n\n')
        # f.write('# Trailing edge coordinates\nTE_coords = ['+str(RSS.TEpt[0])+','+str(RSS.TEpt[1])+']\n\n')
        f.write('# Trailing edge coordinates\nTE_coords = [{0},{1}]\n\n'.format(RSS.TEpt[0], RSS.TEpt[1]))
        # f.write('# Throat coordinate points\nThrt_coords = ['+str(RSS.thrtStrtUp[0])+','+str(RSS.thrtStrtUp[1])+','+str(RSS.thrtEndUp[0])+','+str(RSS.thrtEndUp[1]))
        # f.write(','+str(RSS.thrtStrtUp[0])+','+str(RSS.thrtStrtUp[1])+','+str(RSS.thrtEndUp[0])+','+str(RSS.thrtEndUp[1])+']\n\n')
        f.write(
            '# Throat coordinate points\nThrt_coords_up = [{},{},{},{}]\n'.format(RSS.thrtStrtUp[0], RSS.thrtStrtUp[1],
                                                                                  RSS.thrtEndUp[0], RSS.thrtEndUp[1]))
        f.write('Thrt_coords_down = [{},{},{},{}]\n\n'.format(RSS.thrtStrtD[0], RSS.thrtStrtD[1], RSS.thrtEndD[0],
                                                              RSS.thrtEndD[1]))
        # f.write('# Nozzle coordinate points\nNoz_coords = ['+str(RSS.nozStrtUp[0])+','+str(RSS.nozStrtUp[1])+','+str(RSS.nozEndUp[0])+','+str(RSS.nozEndUp[1]))
        # f.write(','+str(RSS.nozStrtUp[0])+','+str(RSS.nozStrtUp[1])+','+str(RSS.nozEndUp[0])+','+str(RSS.nozEndUp[1])+']\n\n')
        # f.write('# Nozzle coordinate points\nNoz_coords = [{0},{1},{2},{3},{4},{5},{6},{7}]\n\n'.format(RSS.nozStrtUp[0],RSS.nozStrtUp[1],
        #	RSS.nozEndUp[0],RSS.nozEndUp[1],RSS.nozStrtD[0],RSS.nozStrtD[1],RSS.nozEndD[0],RSS.nozEndD[1]))
        f.write('# Nozzle coordinate points\nNoz_coords_up = [{},{},{},{}]\n'.format(RSS.nozStrtUp[0], RSS.nozStrtUp[1],
                                                                                     RSS.nozEndUp[0], RSS.nozEndUp[1]))
        f.write('Noz_coords_down = [{},{},{},{}]\n\n'.format(RSS.nozStrtD[0], RSS.nozStrtD[1], RSS.nozEndD[0],
                                                             RSS.nozEndD[1]))
        # f.write('# Nozzle cross sectional exit length:\nnozArea = '+str(RSS.nozArea)+'\n\n')
        f.write('# Nozzle cross sectional exit length:\nnozArea = {}\n\n'.format(RSS.nozArea))
        # f.write('# Nozzle half throat length:\nHalfThroat = '+str(RSS.scaledHalfThroat)+'\n\n')
        f.write('# Nozzle throat length:\nThroat = {}\n\n'.format(2 * RSS.scaledHalfThroat))
        # f.write('# Nozzle scaling value:\nScaleNoz = '+str(ScaleNoz)+'\n\n')
        f.write('# Nozzle scaling value:\nScaleNoz = ' + str(ScaleNoz) + '\n\n')
        # f.write('# Blade outlet radius:\nRbld = '+str(RSS.BladeRout)+'\n\n')
        f.write('# Blade outlet radius:\nRbld = {}\n\n'.format(RSS.BladeRout))
        f.write('# Blade TE thickness:\nTE_thk = {}\n\n'.format(RSS.TE_thk))
        # f.write('# Blade pitch length:\np = '+str(RSS.pitchLen)+'\n\n')
        f.write('# Blade pitch length:\np = {}\n\n'.format(RSS.pitchLen))
        # f.write('# Blade chord length:\nC = '+str(RSS.chordLen)+'\n\n\n')
        f.write('# Blade chord length:\nC = {}\n\n\n'.format(RSS.chordLen))
        f.write('#Stator inlet boundaries:\nB_IN_1 = {},{}\nB_IN_2 = {},{}\n\n\n'.format(RSS.INLET[0,0],RSS.INLET[0,1],RSS.INLET[-1,0],RSS.INLET[-1,1]))
        f.write(
            '#Stator outlet boundaries:\nB_OUT_1 = {},{}\nB_OUT_2 = {},{}\n\n\n'.format(RSS.OUTLET[0, 0], RSS.OUTLET[0, 1], RSS.OUTLET[-1, 0],
                                                                   RSS.OUTLET[-1, 1]))
        f.write('### END OF FILE ###')


# f.close()

# ---------------------------------------------------------------------------------------------#
#
# Plots the dimentional constraints
#
# ---------------------------------------------------------------------------------------------#

def RadiusLimitPlot360(radii, colors):
    rad = np.linspace(0, 2 * math.pi, 181)
    for i in range(len(radii)):
        x = []
        y = []
        for j in rad:
            x = x + [radii[i] * math.cos(j)]
            y = y + [radii[i] * math.sin(j)]
        plt.plot(x, y, colors[i])


def RadiusLimitPlot(radii, colors, RSS, pitch, nBlades=-1):
    if (nBlades == -1):
        RadiusLimitPlot360(radii, colors)
    else:
        if (nBlades % 2 == 0):
            strtpt = Rotate_about_center(RSS.LEcoords, nBlades / 2 * pitch, [0, 0])
            ang1 = math.atan2(strtpt[1], strtpt[0])
            endpt = Rotate_about_center(RSS.TEpt, (-nBlades / 2 + 1) * pitch, [0, 0])
            ang2 = math.atan2(endpt[1], endpt[0])
        else:
            strtpt = Rotate_about_center(RSS.LEcoords, math.floor(nBlades / 2) * pitch, [0, 0])
            # plt.plot(strtpt[0],strtpt[1],'*r')
            ang1 = math.atan2(strtpt[1], strtpt[0])
            endpt = Rotate_about_center(RSS.TEpt, math.ceil(-nBlades / 2) * pitch, [0, 0])
            # plt.plot(endpt[0],endpt[1],'*r')
            ang2 = math.atan2(endpt[1], endpt[0])

        if ang1 < 0: ang1 = 2 * math.pi - abs(ang1)
        if ang2 < 0: ang2 = 2 * math.pi - abs(ang2)

        if ang2 < ang1:
            temp = ang1
            ang1 = ang2
            ang2 = temp

        # pdb.set_trace()

        rad = np.linspace(ang1, ang2, int(math.degrees(abs(ang2 - ang1))))

        for i in range(len(radii)):
            x = []
            y = []
            for j in rad:
                x = x + [radii[i] * math.cos(j)]
                y = y + [radii[i] * math.sin(j)]
            plt.plot(x, y, colors[i])


# ---------------------------------------------------------------------------------------------#
#
# Plots the complete stator stage including all the vanes
#
# ---------------------------------------------------------------------------------------------#

def plotStatorStage(Vane, nBlades, pitch):
    # nBlades *= nBlades
    pltVane = np.empty((len(Vane), 3))
    # nBlades = 5
    if (nBlades % 2 == 0):
        rng = np.linspace(-nBlades / 2 + 1, nBlades / 2, nBlades)
    else:
        rng = np.linspace(math.ceil(-nBlades / 2), math.floor(nBlades / 2), nBlades)
    # pdb.set_trace()
    for i in rng:
        for j in range(len(Vane)): pltVane[j, :] = Rotate_about_center(Vane[j, :], i * pitch, [0, 0])
        plt.plot(pltVane[:, 0], pltVane[:, 1], 'k')


# def plotStatorStage(Vane,nBlades,pitch):
#    for i in range(nBlades*2):
#        for j in range(len(Vane)): Vane[j,:] = Rotate_about_center(Vane[j,:],pitch,[0,0])
#        plt.plot(Vane[:,0],Vane[:,1],'k')

# ---------------------------------------------------------------------------------------------#
#
# Plots the boundaries for the meshing
#
# ---------------------------------------------------------------------------------------------#

def plotBoundaries(RSS):
    plt.plot(RSS.PERIO2[:, 0], RSS.PERIO2[:, 1], '--c')
    plt.plot(RSS.PERIO1[:, 0], RSS.PERIO1[:, 1], '--c')
    plt.plot(RSS.INLET[:, 0], RSS.INLET[:, 1], '--m')
    plt.plot(RSS.OUTLET[:, 0], RSS.OUTLET[:, 1], '--g')
    plt.plot(RSS.PRESSURE[:, 0], RSS.PRESSURE[:, 1], '--r')
    plt.plot(RSS.SUCTION[:, 0], RSS.SUCTION[:, 1], '--y')

##---------------------------------------------------------------------------------------------#
## END
##---------------------------------------------------------------------------------------------#

