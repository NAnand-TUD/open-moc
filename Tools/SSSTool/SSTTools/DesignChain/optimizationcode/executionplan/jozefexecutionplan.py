#from os import system
from ..executionunit.executionunit import Python3WithConfExecutionUnit
from ..executionunit.executionunit import SU2_PERExecutionUnit
from ..executionunit.su2executionunit import *
from executionplan import ExecutionPlan
from ..config.mocconfig import MOCConfig
from ..config.sstconfig import SSTConfig
from ..config.umg2config import UMG2Config
from ..config.su2_config_jozef import SU2Config


class JozefExecutionPlan_Blade(ExecutionPlan):
    def __init__(self, mocconfig, sstconfig):

        # ---------------------------- EXECUTION UNITS -----------------------#
        ExecutionPlan.__init__(self)

        # generate the supersonic conv/div nozzle coords
        # mocconfig.output_filename = "test"
        moctool = Python3WithConfExecutionUnit(mocconfig, 'moctool.py')

        # Build radial stator blade
        ssttool = Python3WithConfExecutionUnit(sstconfig, 'SST.py')
        ssttool.req_data = ['Db']

        # ---------------------------- EXECUTION PLAN -----------------------#
        self.executionpath =        [
                                     moctool,
                                     ssttool,
                                    ]

class JozefExecutionPlan_Sim(ExecutionPlan):
    def __init__(self, umg2config, su2per_config, su2cfd_config, su2sol_config, np, mesh_dat = False):

        # ---------------------------- EXECUTION UNITS -----------------------#
        ExecutionPlan.__init__(self)

        meshgen = Python3WithConfExecutionUnit(umg2config, 'UMG2.py')
        if mesh_dat:
            meshgen.req_data = ['Db']

        #su2 per
        meshper = SU2_PERExecutionUnit(su2per_config, 'SU2_PER', True)
        #if su2_dat:
        #    meshper.req_data = ['su2mesh.su2']

        #su2 cfd
        su2cfd = SU2MultiExecutionUnit(np, su2cfd_config, "SU2_CFD")

        #su2solve
        su2sol = SU2SingleExecutionUnit(su2sol_config, "SU2_SOL")

        # ---------------------------- EXECUTION PLAN -----------------------#
        self.executionpath =        [
                                     meshgen,
                                     meshper,
                                     su2cfd,
                                     su2sol
                                    ]

class JozefExecutionPlan_SU2(ExecutionPlan):
    def __init__(self, su2per_config, su2cfd_config, su2sol_config, np):

        # ---------------------------- EXECUTION UNITS -----------------------#
        ExecutionPlan.__init__(self)

        #su2 per
        meshper = SU2_PERExecutionUnit(su2per_config, 'SU2_PER', True)
        meshper.req_data = ['su2mesh.su2']

        #su2 cfd
        su2cfd = SU2MultiExecutionUnit(np, su2cfd_config, "SU2_CFD")

        #su2solve
        su2sol = SU2SingleExecutionUnit(su2sol_config, "SU2_SOL")

        # ---------------------------- EXECUTION PLAN -----------------------#
        self.executionpath =        [
                                     meshper,
                                     su2cfd,
                                     su2sol
                                    ]

class JozefExecutionPlan_Full(ExecutionPlan):
    def __init__(self, mocconfig, sstconfig, umg2config, su2per_config, su2cfd_config, su2sol_config, np):

        # ---------------------------- EXECUTION UNITS -----------------------#
        ExecutionPlan.__init__(self)

        # generate the supersonic conv/div nozzle coords
        # mocconfig.output_filename = "test"
        moctool = Python3WithConfExecutionUnit(mocconfig, 'moctool.py')

        # Build radial stator blade
        ssttool = Python3WithConfExecutionUnit(sstconfig, 'SST.py')
        ssttool.req_data = ['Db']

        # Mesh generation with SU2
        meshgen = Python3WithConfExecutionUnit(umg2config, 'UMG2.py')

        #su2 per
        meshper = SU2_PERExecutionUnit(su2per_config, 'SU2_PER', True)

        #su2 cfd
        su2cfd = SU2MultiExecutionUnit(np, su2cfd_config, "SU2_CFD")

        #su2solve
        su2sol = SU2SingleExecutionUnit(su2sol_config, "SU2_SOL")

        # ---------------------------- EXECUTION PLAN -----------------------#
        self.executionpath =        [
                                     moctool,
                                     ssttool,
                                     meshgen,
                                     meshper,
                                     su2cfd,
                                     su2sol
                                    ]