#!/usr/bin/python3

#==========================================================================================/
# author: Nitish Anand				                                           |
# 	: PhD Candidate,                                                                   |
#	: Power and Propulsion,                                                            |
#	: TU Delft,                                                                        |
#	: The Netherlands                                                                  |
#                                                                                          |
# email : nitish.ug2@gmail.com                                                             |
#                                                                                          |
#                                                                                          |
# Description: CODE Sending Runs to different Features                                     |
#                                                                                          |
#==========================================================================================|

#---------------------------------------------------------------------------------------------#
## START: Load Modules
#---------------------------------------------------------------------------------------------#
import os
import sys
import errno, subprocess
import pdb
import re

#---------------------------------------------------------------------------------------------#
## TEST PACKAGES: Test if Packages are available
#---------------------------------------------------------------------------------------------#
try:
    import matplotlib.pyplot as plt
except:
    print("MatplotLib not availabe do: \n\t sudo pip3 install matplotlib\n\t sudo apt-get install python-tk (I guess) \n\t See README.md")
try:
    import numpy as np
except:
    print("Numpy not availabe do: \n\t sudo pip3 install numpy\n\t See README.md")
try:
    import scipy as sc
except:
    print("scipy not availabe do: \n\t sudo pip3 install scipy\n\t See README.md")
try:
    import sympy as sc
except:
    print("sympy not availabe do: \n\t sudo pip3 install sympy\n\t See README.md")
try:
    import CoolProp as sc
except:
    print("CoolProp not availabe do: \n\t sudo pip3 install coolprop\n\t See README.md")

#---------------------------------------------------------------------------------------------#
## IMPORT SPL PACKAGES: Function to Read Input Files <Should be used from IO file>
#---------------------------------------------------------------------------------------------#
#TODO Make IO common
MOC_HOME = os.environ['MOC_HOME']
sys.path.append(MOC_HOME+'/Tools/MoCTool/srcMOC/')

from IO import *

#---------------------------------------------------------------------------------------------#
## PROCESS: Load the file and Save the Directory
#---------------------------------------------------------------------------------------------#
DIR                       = os.getcwd()+'/'
try:    INFile            = sys.argv[1]
except: 
    print("\n\n *** CONFIGURATION FILE NOT PROVIDED *** \n\n\tCheckout TestCases.. \n\tExiting... \n")
    sys.exit()

IN                        = ReadUserInput(INFile)

#---------------------------------------------------------------------------------------------#
## PROCESS: Save the location of the tools
#---------------------------------------------------------------------------------------------#	
SSS_Exec = MOC_HOME+'/Tools/SSSTool/SST.py'
MOC_Exec = MOC_HOME+'/Tools/MoCTool/moctool.py'
SSR_Exec = MOC_HOME+'/Tools/SSRTool/moctool.py'
SIM_Exec = MOC_HOME+'/Tools/SimTool/moctool.py'


#---------------------------------------------------------------------------------------------#
## PROCESS: Create the command depending on the information in the config file
#---------------------------------------------------------------------------------------------#
if IN['TOOL_TYPE'] == 'MOC':
        Command  = MOC_Exec + " " + INFile
elif IN['TOOL_TYPE'] == "SSS":
        Command  = SSS_Exec + " " + INFile
elif IN['TOOL_TYPE'] == "SSR":
        Command  = SSR_Exec + " " + INFile
elif IN['TOOL_TYPE'] == "SimulationLoop":
        Command  = SIM_Exec + " " + INFile
else:        
        print('*** OPTION MISSING ??  ***\n\n ')
        print('\tTOOL_TYPE not specified in the configuration file.. \n\tPlease check test cases or Documentation\n\t')
        print('Exiting .... ')
        sys.exit(0)

#---------------------------------------------------------------------------------------------#
## PROCESS: Run the command
#---------------------------------------------------------------------------------------------#

try:
    proc = subprocess.Popen( Command, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE )
    return_code = proc.wait()
    message = (proc.stderr.read())
    print(message.decode('ascii'))
except:
    print("Did subprocess crash ??? ")

#---------------------------------------------------------------------------------------------#
## END: END OF PYTHON FILE
#---------------------------------------------------------------------------------------------#        
