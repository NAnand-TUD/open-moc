from executionunit import *

class SU2ExecutionUnit(ExecutionUnit):
    def __init__(self, su2config, exec_name):
        ExecutionUnit.__init__(self, su2config, exec_name)
        #self.installdir="~/su2-turbomachinery-installation/bin/"
        #self.installdir="~/su2-feature-turbomachinery-test/bin/"
        self.installdir="~/Programs/SU2/bin"


class SU2SingleExecutionUnit(SU2ExecutionUnit):
    def __init__(self, su2config, exec_name):
        SU2ExecutionUnit.__init__(self, su2config, exec_name)
        self.set_execution_commands()

    def set_execution_commands(self):
        self.execution_commands.append(
            "PATH=$PATH:" +
            self.installdir + " "+
            self.exec_name+ " " +
            self.config_filename + " > " +
            self.log_filename
            )

class SU2MultiExecutionUnit(SU2ExecutionUnit):
    def __init__(self, ncores, su2cfdconfig, exec_name):
        SU2ExecutionUnit.__init__(self, su2cfdconfig, exec_name)
        self.ncores = ncores
        self.set_execution_commands()


    def set_execution_commands(self):
        self.execution_commands.append(
            "PATH=$PATH:" +self.installdir + " mpirun.openmpi -np " + str(self.ncores) +
            " "+self.exec_name+" " + self.config_filename + " > " + self.log_filename
            )