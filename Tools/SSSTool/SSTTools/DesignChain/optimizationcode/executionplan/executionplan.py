
from ..executionunit.su2executionunit import SU2MultiExecutionUnit
class ExecutionPlan(object):
    """The ExecutionPlan Class is a collection of different ExecutionUnit stringed together to form an executionpath"""
    def __init__(self):
        """
        :ivar executionpath: initial value: []
        :vartype executionpath: list[:class:`optimizationcode.executionunit.executionunit.ExecutionUnit`]
        """
        self.executionpath = []

    def add_execution_units(self, executionunit):
        """
        Add the execution units to the executionpath

        :param executionunit: ExecutionUnit object
        :type executionunit: :class:`code.executionunit.executionunit.ExecutionUnit`
        """
        self.executionpath.extend(executionunit)

    def get_last_executionunit_file(self):
        """
        Get the output file of the last ExecutionUnit in the the executionpath

        :return: output_filename
        :rtype: str
        """
        executionunit = self.executionpath[-1]
        return executionunit.output_filename

    def get_su2logfile(self):
        """
        Get the name of the log file of the SU2-Simulation
        
        :return: log_filename
        :rtype: str
        """
        for eu in self.executionpath:
            if isinstance(eu, SU2MultiExecutionUnit):
                return eu.log_filename
