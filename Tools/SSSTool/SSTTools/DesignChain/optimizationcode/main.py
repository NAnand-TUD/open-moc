# 
################## FILE NAME: main.py ##########################
#================================================================
# author: Jozef Stuijt                                          |
#       : Master Student,                                       |
#	    : Process and Energy Department,                    |
#	    : TU Delft,						                        |
#	    : The Netherlands					                    |
#								                                |
# email : joe.stuijt@gmail.com					                |
# 								                                |
# Description: Main script needed to run design & simulation    |
#              chain, based on the original work of PhD         |
#              candidate Stephan Smit (3mE, TU Delft)           |
#                                                               |
#================================================================
#
#
# MAIN CODE
#
import logging
import sys
import re
from os import path
import os
import numpy as np
import pdb
from executionplan.jozefexecutionplan import *
import multiprocessing

from jobs.slurmjob import *
from jobs.localjob import *
from CoolProp.CoolProp import PropsSI

def set_log_setting():
    logging.basicConfig(filename='main.log', level=logging.DEBUG)

def remove_completed_sims(var_rng, simdir):
    vals = []
    for dirname in os.listdir(simdir):
        try:
            vals.append(float(dirname.split('_')[-1].replace('-','.')))
        except:
            pass
    var_not_dup = []
    for v in var_rng:
        if v not in vals:
            var_not_dup.append(v)

    return np.array(var_not_dup)

def mesh_convergence(var_rng, var_type, nupr=7, cluster=True, numit=None, blfac=None, tefac=None, nemel=None, Ma=None):
    mocconfig, sstconfig, umg2config, su2per_config, su2cfd_config, su2sol_config = config_init()

    if Ma != None:
        mocconfig.Noz_Design_Mach = Ma

    if numit != None:
        su2per_config.EXT_ITER = numit + 1
        su2cfd_config.EXT_ITER = numit + 1
        su2sol_config.EXT_ITER = numit + 1


    if blfac != None:
        umg2config.BL_FAC = blfac

    if tefac != None:
        umg2config.RAD_FAC = tefac

    if nemel != None:
        umg2config.MIN_NELEM = nemel

    if var_type=='numit':
        sim_id = 'Mesh_Convergence_NumIt'
        job_name = 'Iter_{}'.format(var_rng[0]).replace('.', '-')

        su2per_config.EXT_ITER = var_rng[0] + 1
        su2cfd_config.EXT_ITER = var_rng[0] + 1
        su2sol_config.EXT_ITER = var_rng[0] + 1

    if var_type=='BL':
        sim_id = 'Mesh_Convergence_BL'

    if var_type=='TE':
        sim_id = 'Mesh_Convergence_TE'

    if var_type=='numel':
        sim_id = 'Mesh_Convergence_NemEl'    

    for v in var_rng:

        #if var_type=='iter':
        if var_type=='BL':
            umg2config.BL_FAC = v
            job_name = 'BL_FAC_{:.2f}'.format(v).replace('.', '-')
        if var_type=='TE':
            umg2config.RAD_FAC = v
            job_name = 'TE_FAC_{:.2f}'.format(v).replace('.', '-')
        if var_type=='numel':
            umg2config.MIN_NELEM = v
            job_name = 'MIN_NELEM_{:.2e}'.format(v).replace('.', '-')

        ep = JozefExecutionPlan_Full(mocconfig, sstconfig, umg2config, su2per_config, su2cfd_config, su2sol_config, nupr)
        if cluster:
            slurmjob = SlurmJob(nupr, 1, 2, executionplan=ep, simulation_id=sim_id)
            slurmjob.name = job_name
            slurmjob.initialize_jobfolder(ep)
            slurmjob.execute()
        else:
            localjob = LocalJob(executionplan=ep, simulation_id=sim_id)
            localjob.name = job_name
            localjob.initialize_jobfolder(ep)
            localjob.execute()

def sim_range_sweep(var_rng, var_nam = 'Mach', nupr = 7, cluster=True, Pback = None, AreaRatio = None, NozMach = None, RealRout = None, phi= None, Rout=None, mocconfig = None, 
    ThroatMoc = None, BC = 'Giles', Build = 'full', filter_dups = True):

    sstconfig=None
    umg2config=None
    su2per_config=None
    su2cfd_config=None
    su2sol_config=None

    add_id = []

    if BC == 'Riemann':
        mocconfig=MOCConfig('')
        su2per_config=SU2Config('')
        su2cfd_config=SU2Config('')
        su2sol_config=SU2Config('')
        Pout = su2per_config.MARKER_GILES.split(', ')[11]
        MARKER_RIEMANN = '(inflow, TOTAL_CONDITIONS_PT, {}, {}, -1.0, 0.0, 0.0, outflow, STATIC_PRESSURE, {}, 0.0, 0.0, 0.0, 0.0)'.format(mocconfig.Po, mocconfig.To, Pout)
        su2per_config.MARKER_RIEMANN = MARKER_RIEMANN
        su2cfd_config.MARKER_RIEMANN = MARKER_RIEMANN
        su2sol_config.MARKER_RIEMANN = MARKER_RIEMANN
        su2per_config.RAMP_OUTLET_PRESSURE = 'NO'
        su2cfd_config.RAMP_OUTLET_PRESSURE = 'NO'
        su2sol_config.RAMP_OUTLET_PRESSURE = 'NO'
        su2per_config.MARKER_GILES = '( NONE )'
        su2cfd_config.MARKER_GILES = '( NONE )'
        su2sol_config.MARKER_GILES = '( NONE )'

        su2per_config.excludeconfig.append('RAMP_OUTLET_PRESSURE_COEFF')
        su2cfd_config.excludeconfig.append('RAMP_OUTLET_PRESSURE_COEFF')
        su2sol_config.excludeconfig.append('RAMP_OUTLET_PRESSURE_COEFF')

        add_id.append('RiemannBC')


    if Pback != None:
        if BC == 'Giles':
            su2per_config=SU2Config('')
            su2cfd_config=SU2Config('')
            su2sol_config=SU2Config('')
            su2per_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[11],str(Pback))
            su2cfd_config.MARKER_GILES = su2cfd_config.MARKER_GILES.replace(su2cfd_config.MARKER_GILES.split(', ')[11],str(Pback))
            su2sol_config.MARKER_GILES = su2sol_config.MARKER_GILES.replace(su2sol_config.MARKER_GILES.split(', ')[11],str(Pback))
        if BC == 'Riemann':
            su2per_config.MARKER_RIEMANN = su2per_config.MARKER_RIEMANN.replace(su2per_config.MARKER_RIEMANN.split(', ')[11],str(Pback))
            su2cfd_config.MARKER_RIEMANN = su2cfd_config.MARKER_RIEMANN.replace(su2cfd_config.MARKER_RIEMANN.split(', ')[11],str(Pback))
            su2sol_config.MARKER_RIEMANN = su2sol_config.MARKER_RIEMANN.replace(su2sol_config.MARKER_RIEMANN.split(', ')[11],str(Pback))

        mocconfig = MOCConfig('')
        #pdb.set_trace()
        add_id.append('PR_{:.1f}'.format(mocconfig.Po/Pback).replace('.','-'))

    mocconfig, sstconfig, umg2config, su2per_config, su2cfd_config, su2sol_config = config_init(mocconfig, sstconfig, umg2config, su2per_config, su2cfd_config, su2sol_config)

    if NozMach != None:
        mocconfig.Noz_Design_Mach = NozMach
        add_id.append('NozMach_{:.2f}'.format(mocconfig.Noz_Design_Mach).replace('.','-'))

    if phi != None:
        sstconfig.flowAngle = phi
        add_id.append('FlowAngle_{:.1f}'.format(sstconfig.flowAngle).replace('.','-'))

    if RealRout != None:
        sstconfig.RealRout = RealRout
        add_id.append('RoutMesh_{:.3f}'.format(sstconfig.RealRout).replace('.','-'))

    if Rout != None:
        sstconfig.ActualRout = Rout
        sstconfig.radiusOut = Rout
        add_id.append('Rout_{:.3f}'.format(sstconfig.ActualRout).replace('.','-'))

    if AreaRatio != None:
       del sstconfig.excludeconfig[-1]
       sstconfig.AreaRatio = AreaRatio
       add_id.append('CnstAR')

    if ThroatMoc !=None:
        mocconfig.y_t = ThroatMoc
        add_id.append('ThroatMoC_{:.0f}'.format(mocconfig.y_t))

    if var_nam == 'Mach':
        sim_id = 'Mach_Range'

    elif var_nam == 'FlowAngle':
        i=0
        RadIn = np.linspace(0.25,sstconfig.ActualRin,len(var_rng))
        sim_id = 'FlowAngle_Range'

    elif var_nam == 'Rout':
        sim_id = 'Rout_Range'
        dR = sstconfig.ActualRin-sstconfig.ActualRout

    elif var_nam == 'BackPressure':
        sim_id = 'BackPressure_Range'

    elif var_nam == 'RealRout':
        sim_id = 'RealRout_Range'

    elif var_nam == 'Visc':
        sim_id = 'Viscosity_Range'


    if len(add_id)==0:
        sim_id = sim_id+'_BaseCase'
    else:
        for j, add in enumerate(add_id):
            sim_id = sim_id+'_{}'.format(add)

    try:
        oldmoc.Po
    except NameError:
        pass        
    else:
        sim_id = 'NewInit_'+sim_id

    if filter_dups:
        try:
            var_rng = remove_completed_sims(var_rng, os.path.join('simulations',sim_id))
        except OSError:
            pass
 
    for v in var_rng:

        if var_nam == 'Mach':
            mocconfig.Noz_Design_Mach = v
            job_name = 'Mach_{:.2f}'.format(v).replace('.', '-')

        if var_nam == 'FlowAngle':
            sstconfig.ActualRin = RadIn[i]
            sstconfig.flowAngle = v
            job_name = 'FlowAngle_{:.1f}'.format(v).replace('.', '-')
            i+=1

        if var_nam == 'Rout':
            sstconfig.ActualRout = v
            sstconfig.radiusOut = v
            job_name = 'Rout_{:.3f}'.format(v).replace('.', '-')

        if var_nam == 'BackPressure':
            Pback = mocconfig.Po/v 
            su2per_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[11],str(Pback))
            su2cfd_config.MARKER_GILES = su2cfd_config.MARKER_GILES.replace(su2cfd_config.MARKER_GILES.split(', ')[11],str(Pback))
            su2sol_config.MARKER_GILES = su2sol_config.MARKER_GILES.replace(su2sol_config.MARKER_GILES.split(', ')[11],str(Pback))
            
            job_name = 'PR_{:.2f}'.format(v).replace('.', '-')

        if var_nam == 'RealRout':
            sstconfig.RealRout = v
            job_name = 'RealRout_{1:.3f}'.format(v).replace('.', '-')

        if var_nam == 'Visc':
            job_name = 'Mu_{:1.1e}'.format(v).replace('.', '-')

            su2per_config.MU_CONSTANT = v
            su2cfd_config.MU_CONSTANT = v
            su2sol_config.MU_CONSTANT = v
            
        if Build == 'full':
            ep = JozefExecutionPlan_Full(mocconfig, sstconfig, umg2config, su2per_config, su2cfd_config, su2sol_config, nupr)
        if Build == 'blade':
            ep = JozefExecutionPlan_Blade(mocconfig, sstconfig) 
        
        if cluster:
            slurmjob = SlurmJob(nupr, 1, 1.5, executionplan=ep, simulation_id=sim_id)
            slurmjob.name = job_name
            slurmjob.initialize_jobfolder(ep)
            slurmjob.execute()
        else:
            localjob = LocalJob(executionplan=ep, simulation_id=sim_id)
            localjob.name = job_name
            localjob.initialize_jobfolder(ep)
            localjob.execute()

def config_init(mocconfig=None, sstconfig=None, umg2config=None, su2per_config=None, su2cfd_config=None, su2sol_config=None):
    
    if mocconfig==None:
        mocconfig = MOCConfig('')
    if sstconfig==None:
        sstconfig = SSTConfig('')
    if umg2config==None:
        umg2config = UMG2Config('')
    if su2per_config==None:
        su2per_config = SU2Config('')
    if su2cfd_config==None:
        su2cfd_config = SU2Config('')
    if su2sol_config==None:
        su2sol_config = SU2Config('')

    Fluid = mocconfig.FLDNAME
    Po = mocconfig.Po
    To = mocconfig.To

    #pdb.set_trace()

    try:
        Pback = float(su2per_config.MARKER_GILES.split(', ')[11])
    except IndexError:
        Pback = float(su2per_config.MARKER_RIEMANN.split(', ')[9])

    so = PropsSI('S','P',Po,'T',To,Fluid)
    Ho = PropsSI('H','P',Po,'T',To,Fluid)

    GAMMA_IN= PropsSI('C','P',Po,'T',To,Fluid)/PropsSI('CVMASS','P',Po,'T',To,Fluid)
    GAMMA_OUT= PropsSI('C','P',Pback,'S',so,Fluid)/PropsSI('CVMASS','P',Pback,'S',so,Fluid)
    GAMMA = 0.5*(GAMMA_IN+GAMMA_OUT)
    GAS_CONSTANT= PropsSI('GAS_CONSTANT',Fluid)/PropsSI('M',Fluid)
    CRITICAL_TEMPERATURE= PropsSI('TCRIT',Fluid)
    CRITICAL_PRESSURE= PropsSI('PCRIT',Fluid)
    Vc = 1/PropsSI('DMOLAR','P',CRITICAL_PRESSURE,'T',CRITICAL_TEMPERATURE,Fluid)
    MU_IN = PropsSI('V','P',Po,'T',To,Fluid)
    MU_OUT = PropsSI('V','P',Pback,'S',so, Fluid)
    MU_CONSTANT = 0.5*(MU_IN+MU_OUT)
    KT_IN = PropsSI('L','P',Po,'T',To, Fluid)
    KT_OUT = PropsSI('L','P',Pback,'S',so, Fluid)
    KT_CONSTANT = 0.5*(KT_IN+KT_OUT)

    mocconfig.Tc = CRITICAL_TEMPERATURE
    mocconfig.Pc = CRITICAL_PRESSURE
    mocconfig.Vc = Vc
    mocconfig.M = GAS_CONSTANT/1000
    mocconfig.Ho = Ho/1000
    mocconfig.so = so/1000
    mocconfig.gamma = GAMMA

    sstconfig.outFileMoc = mocconfig.File_Name_NCods
    sstconfig.kernelRadiusRatio = mocconfig.rho_t
    sstconfig.UMG2Name = umg2config.FILE_NAME

    umg2config.SPECS_FILE = sstconfig.SpecsName

    try:
        su2per_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[2],str(Po))
        su2per_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[3],str(To))
    except IndexError:
        su2per_config.MARKER_RIEMANN = su2per_config.MARKER_RIEMANN.replace(su2per_config.MARKER_RIEMANN.split(', ')[2],str(Po))
        su2per_config.MARKER_RIEMANN = su2per_config.MARKER_RIEMANN.replace(su2per_config.MARKER_RIEMANN.split(', ')[3],str(To))
    su2per_config.FREESTREAM_PRESSURE = Po
    su2per_config.FREESTREAM_TEMPERATURE = To
    su2per_config.FREESTREAM_DENSITY = PropsSI('D','P',Po,'T',To,Fluid)
    su2per_config.GAMMA_VALUE = GAMMA
    su2per_config.GAS_CONSTANT = GAS_CONSTANT
    su2per_config.CRITICAL_TEMPERATURE = CRITICAL_TEMPERATURE
    su2per_config.CRITICAL_PRESSURE = CRITICAL_PRESSURE
    su2per_config.ACENTRIC_FACTOR = PropsSI('ACENTRIC',Fluid)
    su2per_config.MU_CONSTANT = MU_CONSTANT
    su2per_config.KT_CONSTANT = KT_CONSTANT

    try:
        su2cfd_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[2],str(Po))
        su2cfd_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[3],str(To))
    except IndexError:
        su2cfd_config.MARKER_RIEMANN = su2cfd_config.MARKER_RIEMANN.replace(su2cfd_config.MARKER_RIEMANN.split(', ')[2],str(Po))
        su2cfd_config.MARKER_RIEMANN = su2cfd_config.MARKER_RIEMANN.replace(su2cfd_config.MARKER_RIEMANN.split(', ')[3],str(To))
    su2cfd_config.FREESTREAM_PRESSURE = Po
    su2cfd_config.FREESTREAM_TEMPERATURE = To
    su2cfd_config.FREESTREAM_DENSITY = PropsSI('D','P',Po,'T',To,Fluid)
    su2cfd_config.GAMMA_VALUE = GAMMA
    su2cfd_config.GAS_CONSTANT = GAS_CONSTANT
    su2cfd_config.CRITICAL_TEMPERATURE = CRITICAL_TEMPERATURE
    su2cfd_config.CRITICAL_PRESSURE = CRITICAL_PRESSURE
    su2cfd_config.ACENTRIC_FACTOR = PropsSI('ACENTRIC',Fluid)
    su2cfd_config.MU_CONSTANT = MU_CONSTANT
    su2cfd_config.KT_CONSTANT = KT_CONSTANT

    try:
        su2sol_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[2],str(Po))
        su2sol_config.MARKER_GILES = su2per_config.MARKER_GILES.replace(su2per_config.MARKER_GILES.split(', ')[3],str(To))
    except IndexError:
        su2sol_config.MARKER_RIEMANN = su2sol_config.MARKER_RIEMANN.replace(su2sol_config.MARKER_RIEMANN.split(', ')[2],str(Po))
        su2sol_config.MARKER_RIEMANN = su2sol_config.MARKER_RIEMANN.replace(su2sol_config.MARKER_RIEMANN.split(', ')[3],str(To))
    su2sol_config.FREESTREAM_PRESSURE = Po
    su2sol_config.FREESTREAM_TEMPERATURE = To
    su2sol_config.FREESTREAM_DENSITY = PropsSI('D','P',Po,'T',To,Fluid)
    su2sol_config.GAMMA_VALUE = GAMMA
    su2sol_config.GAS_CONSTANT = GAS_CONSTANT
    su2sol_config.CRITICAL_TEMPERATURE = CRITICAL_TEMPERATURE
    su2sol_config.CRITICAL_PRESSURE = CRITICAL_PRESSURE
    su2sol_config.ACENTRIC_FACTOR = PropsSI('ACENTRIC',Fluid)
    su2sol_config.MU_CONSTANT = MU_CONSTANT
    su2sol_config.KT_CONSTANT = KT_CONSTANT

    su2cfd_config.config_filename = 'su2_cfd' + su2cfd_config.name + '.cfg'
    su2cfd_config.MESH_FILENAME = su2per_config.MESH_OUT_FILENAME
    su2cfd_config.log_filename = 'su2_cfd' + su2cfd_config.name + '.log'
    su2cfd_config.config_filename = 'su2_cfd' + su2cfd_config.name + '.cfg'

    su2sol_config.config_filename = 'su2_sol' + su2cfd_config.name + '.cfg'
    su2sol_config.MESH_FILENAME = su2per_config.MESH_OUT_FILENAME
    su2sol_config.log_filename = 'su2_sol' + su2cfd_config.name + '.log'
    su2sol_config.config_filename = 'su2_sol' + su2cfd_config.name + '.cfg'

    return mocconfig, sstconfig, umg2config, su2per_config, su2cfd_config, su2sol_config


if __name__ == "__main__":

    "Final Runs"

    ### Geometry sensitivity ###
    #sim_range_sweep(np.linspace(45,85,41),var_nam='FlowAngle',Build='blade',NozMach=2.2,nupr=1)
    #sim_range_sweep(np.linspace(45,85,41),var_nam='FlowAngle',Build='blade',NozMach=1.4,nupr=1)
    #sim_range_sweepnp.linspace(45,85,41),var_nam='FlowAngle',Build='blade',NozMach=3.0,nupr=1)
    #sim_range_sweep(np.linspace(1.4,3.0,33),Build='blade',phi=70,nupr=1)
    #sim_range_sweep(np.linspace(1.4,3.0,33),Build='blade',phi=75,nupr=1)
    #sim_range_sweep(np.linspace(0.09,0.14,11),var_nam='Rout', Build='blade',nupr=1, NozMach=1.4)
    #sim_range_sweep(np.linspace(0.09,0.14,11),var_nam='Rout', Build='blade',nupr=1, NozMach=2.2)
    #sim_range_sweep(np.linspace(0.09,0.14,11),var_nam='Rout', Build='blade',nupr=1, NozMach=3.0)

    ### Mesh Convergence ###
    #mesh_convergence([2e4], 'numit', nupr=14, blfac=0.2, tefac=0, nemel=5e4, Ma=3.0)
    #mesh_convergence(np.linspace(0,0.3,31), 'BL', numit=10e3, tefac=5, nemel=1e5, Ma=3.0)
    #mesh_convergence(np.linspace(0,20,21), 'TE', numit=10e3, blfac=0.15, nemel=1e5, Ma=3.0)
    #mesh_convergence(np.hstack((np.linspace(4e4,2e5,17),np.linspace(2.5e5,3.5e5,3))), 'numel', numit=10e3, blfac=0.2, tefac=0, Ma=3.0)

    ### Reynolds Throat ###
    #sim_range_sweep(np.logspace(-12, -1, 23), var_nam='Visc')

    ### Main Run for different BCs ###
    #sim_range_sweep(np.linspace(1.4,3.0,33))
    #sim_range_sweep(np.linspace(1.4,3.0,33), BC='Riemann')

    ### Main run different mesh outlet radius ###
    #sim_range_sweep(np.linspace(1.4,3.0,33), RealRout=0.104)
    #sim_range_sweep(np.linspace(1.4,3.0,33), RealRout=0.116)

    ### Different pressure ratios ###
    #Pin = 32e5
    #sim_range_sweep(np.linspace(1.4,3.0,33), Pback=Pin/10)
    #sim_range_sweep(np.linspace(1.4,3.0,33), Pback=Pin/20)
    #sim_range_sweep(np.linspace(1.4,3.0,33), Pback=Pin/80)
    #sim_range_sweep(np.linspace(1.4,3.0,33), Pback=Pin/120)
    #sim_range_sweep(np.linspace(1.4,3.0,33), Pback=Pin/160)

    ### Fixed blade different back pressures
    #sim_range_sweep(np.linspace(2.5,20,36), var_nam='BackPressure', NozMach=1.8)
    #sim_range_sweep(np.linspace(5,40,36), var_nam='BackPressure', NozMach=2.2)
    #sim_range_sweep(np.linspace(41,80,40), var_nam='BackPressure', NozMach=2.2)
    #sim_range_sweep(np.linspace(10,60,26), var_nam='BackPressure', NozMach=2.4)
    #sim_range_sweep(np.linspace(10,110,21), var_nam='BackPressure', NozMach=2.8)
    #sim_range_sweep(np.linspace(115,160,10), var_nam='BackPressure', NozMach=2.8)

    "Tests"

    #mesh_convergence([51], 'numit', nupr=2, cluster=False, blfac=0.2, tefac=0, nemel=5e4)
    #mesh_convergence([501], 'numit', nupr=14, blfac=0.2, tefac=0, nemel=5e4, Ma=3.0)
    sim_range_sweep([2.3], cluster=False)
    sim_range_sweep(np.linspace(2.5,3.0,2), cluster=False)


    logging.info('running the simulation cases')
logging.info('finished all the cases')
