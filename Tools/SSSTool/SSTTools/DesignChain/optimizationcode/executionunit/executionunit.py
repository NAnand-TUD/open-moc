from os import getcwd


class ExecutionUnit(object):
    """
    The ExecutionUnit class is created using a configuration object and will implement the correct bash-commands to execute the different executables
    """
    def __init__(self, config, exec_name):
        """
        Initialisation

        :param config: Configuration object
        :type config: :class:`code.config.config.Config`
        :param exec_name: Name of the executable
        :type exec_name: str

        :ivar str execdir: initial value: current working directory + optimizationcode/executables
        :ivar str exec_name: initial value: exec_name
        :ivar str log_dir: initial value: 'log/'
        :ivar str cfg_dir: initial value: 'cfg/'
        :ivar list[str] req_exec: []
        :ivar str config_filename: cfg_dir + config.config_filename
        :ivar str configfilecontent: config.get_config_filecontent()
        :ivar str log_filename: log_dir + config.log_filename
        :ivar str output_filename: config.output_filename


        """
        self.execdir =  "/".join([getcwd(), 'optimizationcode/executables'])
        self.exec_name = exec_name
        self.log_dir = "log/"
        self.cfg_dir = "cfg/"
        self.result_dir = 'results/'
        self.execution_commands=[]
        self.req_exec=[]
        self.req_data=[]
        self.req_mesh=[]
        self.config_filename = self.cfg_dir +config.config_filename
        self.configfilecontent = config.get_config_filecontent()
        self.log_filename = self.log_dir + config.log_filename
        self.output_filename = config.output_filename

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :raise
        """
        raise NotImplementedError

class waitExecutionUnit():
    def __init__(self):
        self.config_filename='1'
        self.configfilecontent='1'
        self.execution_commands=[]
        self.req_exec=[]
        self.req_mesh=[]
        self.req_data=[]


        self.set_execution_commands()

    def set_execution_commands(self):
        self.execution_commands.append('wait')


class NormalExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'exec_name.exe config_filename > log_filename'
    """
    def __init__(self, config, exec_name, wait):
        ExecutionUnit.__init__(self, config, exec_name)
        self.wait = wait
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :return:
        """
        if not self.wait:
            self.execution_commands.append(self.execdir + "/" +
                                           self.exec_name + " " +
                                           self.config_filename + " > " +
                                           self.log_filename + ' &')
        else:
            self.execution_commands.append(self.execdir + "/" +
                                           self.exec_name + " " +
                                           self.config_filename + " > " +
                                           self.log_filename)

class SU2_PERExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'exec_name.exe < config_filename > log_filename'
    """
    def __init__(self, config, exec_name, wait):
        ExecutionUnit.__init__(self, config, exec_name)
        self.wait = wait
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commandfile_existss to execution commands variable
        :return:
        """
        if not self.wait:
            self.execution_commands.append(self.execdir + "/"+
                                           self.exec_name + " < " +
                                           self.config_filename + " > " +
                                           self.log_filename+ ' &')
        else:
            self.execution_commands.append(self.execdir + "/"+
                                           self.exec_name + " < " +
                                           self.config_filename + " > " +
                                           self.log_filename)

class PythonWithConfExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'python exec_name.py config_filename > log_filename'
    """
    def __init__(self, config, exec_name):
        ExecutionUnit.__init__(self, config, exec_name)
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :return:
        """
        self.execution_commands.append("python " + self.execdir + "/" +
                                       self.exec_name + " " +
                                       self.config_filename + " > " +
                                       self.log_filename)

class Python3WithConfExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'python exec_name.py config_filename > log_filename'
    """
    def __init__(self, config, exec_name):
        ExecutionUnit.__init__(self, config, exec_name)
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :return:
        """
        self.execution_commands.append("python3 " + self.execdir + "/" +
                                       self.exec_name + " " +
                                       self.config_filename + " > " +
                                       self.log_filename)



class PythonWithoutConfExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'python exec_name.py > log_filename'
    """
    def __init__(self, config, exec_name):
        ExecutionUnit.__init__(self, config, exec_name)
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :return:
        """
        self.execution_commands.append("python " + self.execdir + "/" +
                                       self.exec_name + " > " +
                                       self.log_filename)

class PythonWithoutConfWithOutputExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'python exec_name.py > log_filename'
    """
    def __init__(self, config, exec_name):
        ExecutionUnit.__init__(self, config, exec_name)
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :return:
        """
        self.execution_commands.append("python " + self.execdir + "/" +
                                       self.exec_name + " "+ self.output_filename + " > " +
                                       self.log_filename)

class PVBatchWithoutConfExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'pvbatch exec_name.py argv[1] > log_filename'
    """
    def __init__(self, config, exec_name, wait):
        ExecutionUnit.__init__(self, config, exec_name)
        self.wait = wait
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :return:
        """
        if not self.wait:
            self.execution_commands.append("pvbatch " + self.execdir + "/" +
                                           self.exec_name + " "+ self.output_filename + " > " +
                                           self.log_filename + ' &')
        else:
            self.execution_commands.append("pvbatch " + self.execdir + "/" +
                                           self.exec_name + " "+ self.output_filename + " > " +
                                           self.log_filename )

class PVBatchWithConfExecutionUnit(ExecutionUnit):
    """Execution unit in the shape of: 'pvbatch exec_name.py argv[1] > log_filename'
    """
    def __init__(self, config, exec_name, wait):
        ExecutionUnit.__init__(self, config, exec_name)
        self.wait = wait
        self.set_execution_commands()

    def set_execution_commands(self):
        """
        Appends execution commands to execution commands variable
        :return:
        """
        if not self.wait:
            self.execution_commands.append("pvbatch " + self.execdir + "/" +
                                           self.exec_name + " "+ self.config_filename + " > " +
                                           self.log_filename + ' &')
        else:
            self.execution_commands.append("pvbatch " + self.execdir + "/" +
                                           self.exec_name + " "+ self.config_filename + " > " +
                                           self.log_filename )

