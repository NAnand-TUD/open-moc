######################FILE NAME: IO.py###########################
#================================================================
# author: Nitish Anand                                           |
#     :Master Student,                                           |
#    :Process and Energy Departmemt,                             |
#    :TU Delft,                                                  |
#    :The Netherlands                                            |
#                                                                |
# email : nitish.ug2@gmail.com                                   |
# 								                                 |
# Description: Reads and Writes data file.			             |
#								                                 |
# MODULE IO: Input/Output					                     |
#	:a. Read User Input					                         |
#	:b. WriteDataFile					                         |
#	:c. WriteNozzleDim					                         |
#	:d. printProgress					                         |
#================================================================

import re
import sys
import os

#---------------------------------------------------------------------------------------------#
#
# Read User Input file 'MOC_Config.in' <Sample Input file>
#

def ReadUserInput(name):
    IN = {}
    infile = open(name,'r')
    for line in infile:
      words = re.split(r'=| |%|\n|#',line)
      if not any(words[0] in s for s in ['\n','%',' ','#']):
        words = list(filter(None,words))
        IN[words[0]] = words[1]
    return IN

#---------------------------------------------------------------------------------------------#
#
# Writes the Nozzle properties in the text file <Nozzle_prop.out>
#

def WriteDataFile(Operation,VAR1,File_Name,flag='Nrm'):

    f=open(File_Name,Operation)
    if Operation == 'w':
        f.write('it\tx\ty\tu\tv\tV\tM\trho\ta\tT\tP\n')
    elif Operation == 'a':
        if flag == 'Nrm':
            for i in VAR1:
                for j in i:
                    f.write("%f %f %f %f %f %f %f %f %f %f\n" %(j.x,j.y,j.u,j.v,j.V,j.M,j.rho,j.a,j.Temp,j.Press))
        elif flag == 'TAB':
            f.write("%d %f %f %f %f %f %f %f\n" %(VAR1.V,VAR1.M,VAR1.rho,VAR1.a,VAR1.Temp,VAR1.Press,VAR1.gamma))
        else:
            print('*** Unknow Flag ***')
            sys.exit(-1)
    f.close()

#---------------------------------------------------------------------------------------------#
#
# Writes the Nozzle co-ordinates in the text file <Nozzle_coords.out>
#

def WriteNozzleDim(Coords,File_Name):
    try:os.remove(os.getcwd()+'/'+File_Name)
    except: print ("Continue: No File to Delete <Co-ordinate file>\n Location:",os.getcwd())
    f=open(File_Name,'w')
    for i in Coords:
        f.write("%f %f %f\n" %(i.x,i.y,0.0))
    f.close()

#---------------------------------------------------------------------------------------------#
#
# Prints the progress bar on the screen [shows the percentage calcualtion left]
#

def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    if percents>100.0: percents=100.0
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")

##
## END
##---------------------------------------------------------------------------------------------#

