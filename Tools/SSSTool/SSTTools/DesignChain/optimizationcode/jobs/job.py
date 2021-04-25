from uuid import uuid4
from os import system, getcwd, path
from time import sleep
import logging
import multiprocessing
from os import path

class Job(object):
    """This is the abstract class
    """
 
    def __init__(self, executionplan, simulation_id):
        """

        :param executionplan: executionplan
        :type executionplan: :class: `optimizationcode.executionplan.executionplan.ExecutionPlan`
        :param simulation_id: simulation id
        :type simulation_id: str

        :ivar str simulation_id: initial value simulation_id
        :ivar str workingdirectory: current working directory
        :ivar str executiondirectory: initial value: cwd/optimizationcode/executables
        :ivar str gen_sim_dir: simulations directory name
        :ivar str simulation_dir: initial value: cwd/simulations
        :ivar list[str] job_header: initial value: []
        :ivar list[str] job_commands: initial value: []
        :ivar list[str] req_exec: initial value: []
        """
        self.simulation_id = simulation_id
        self.workingdirectory = getcwd()
        self.execdir="/".join([self.workingdirectory, 'optimizationcode/executables'])
        self.meshdir="/".join([self.workingdirectory, 'optimizationcode/meshes'])
        self.datadir="/".join([self.workingdirectory, 'optimizationcode/data'])

        self.executionplan = executionplan
        self.gen_sim_dir = 'simulations'
        self.simulation_dir="/".join([self.gen_sim_dir,self.simulation_id])
        self.job_header =[]
        self.job_commands=[]
        self.req_exec=[]
        self.req_data=[]
        self.req_mesh=[]




    def set_jobdir(self):
        """Sets the job directory
        """
        self.jobdir="/".join([self.simulation_dir,self.name])


    def initialize_jobfolder(self, executionplan):
        """
        Initializes the jobfolder where the simulation is going to be run

        :param executionplan: execution plan to be executed
        :type executionplan: :class: `optimizationcode.executionplan.executionplan.ExecutionPlan`
        """
        self.jobdir="/".join([self.simulation_dir,self.name])
        self.log_directory = "/".join([self.jobdir,"log"])
        self.config_directory = "/".join([self.jobdir, "cfg"])
        self.results_directory = "/".join([self.jobdir, "results"])
        self.create_directory(self.gen_sim_dir)
        self.create_directory(self.simulation_dir)
        self.create_directory(self.jobdir)
        self.create_directory(self.log_directory)
        self.create_directory(self.config_directory)
        self.create_directory(self.results_directory)
        self.add_dependent_executables(executionplan)
        self.add_dependent_data(executionplan)
        self.add_dependent_mesh(executionplan)
        self.copy_exec_in_jobdir()
        self.copy_data_in_jobdir()
        self.copy_mesh_in_jobdir()

    def convert_content_in_dict(self, content):
        """
        Converts the content of a file in the shape of "variable=value" in a dictionary

        :param list[str] content: content of the file
        :return: dictionary of content
        :rtype: dict
        """
        dictt = dict()
        for line in content:
            elements = line.strip().split("=")
            dictt[elements[0].strip()] = float(elements[1].strip())
        return dictt

    def write_configurationfiles(self,executionplan):
        """
        Writes all the configuration files of the individual ExecutionUnits in the executionpath of the ExecutionPlan into the correct directory

        :param executionplan: execution plan
        :type executionplan: :class:`code.executionplan.executionplan.ExecutionPlan`
        """
        for executionunit in executionplan.executionpath:
            self.write_file_jobdir(executionunit.config_filename, executionunit.configfilecontent)
        
    def add_dependent_executables(self,executionplan):
        """
        Make one list of the the required executables of the individual ExecutionUnits in the executionpath of the ExecutionPlan

        :param executionplan: execution plan
        :type executionplan: :class:`code.executionplan.executionplan.ExecutionPlan`
        """

        for executionunit in executionplan.executionpath:
            self.req_exec.extend(executionunit.req_exec)

    def add_dependent_data(self,executionplan):
        """
        Make one list of the the required executables of the individual ExecutionUnits in the executionpath of the ExecutionPlan

        :param executionplan: execution plan
        :type executionplan: :class:`code.executionplan.executionplan.ExecutionPlan`
        """

        for executionunit in executionplan.executionpath:
            self.req_data.extend(executionunit.req_data)

    def add_dependent_mesh(self,executionplan):
        """
        Make one list of the the required executables of the individual ExecutionUnits in the executionpath of the ExecutionPlan

        :param executionplan: execution plan
        :type executionplan: :class:`code.executionplan.executionplan.ExecutionPlan`
        """

        for executionunit in executionplan.executionpath:
            self.req_mesh.extend(executionunit.req_mesh)


    def create_uid(self):
        """
        Generate a unique user id and set the name of the job to this userid
        """
        self.uid=uuid4()
        self.name = str(self.uid)

    
    def set_job_commands(self, executionplan):
        """
        Make one list of the job_commands of all individual ExecutionUnits in the executionpath of the ExecutionPlan

        :param executionplan:  execution plan object
        :type executionplan: :class:`code.executionplan.executionplan.ExecutionPlan`
        """
        for executionunit in executionplan.executionpath:
            self.job_commands.extend(executionunit.execution_commands)

    def create_directory(self,directory_name):
        """
        create the directory with a specific name

        :param directory_name: name of the directory
        :type directory_name: str
        :return:
        """
        system('mkdir -p '+ directory_name)

    def set_execution_directory(self, executionplan):
        for executionunit in executionplan.executionpath:
            executionunit.execdir=self.execdir
            executionunit.execution_commands=[]
            executionunit.set_execution_commands()

    def copy_exec_in_jobdir(self):
        """
        Copy all the executables in the req_exec variable into the job directory
        """
        for item in self.req_exec:
            system('cp ' + '/'.join([self.execdir,item]) + ' ' + self.jobdir + '/')   

    def copy_data_in_jobdir(self):

        """
        Copy all the data in the req_data variable into the job directory
        """
        for item in self.req_data:
            if path.isdir('/'.join([self.datadir,item])):
                system('cp -r ' + '/'.join([self.datadir, item]) + ' ' + self.jobdir + '/')
            else:
                system('cp ' + '/'.join([self.datadir, item]) + ' ' + self.jobdir + '/')

            #print('cp -r ' + '/'.join([self.datadir,item]) + ' ' + self.jobdir + '/')



    def copy_mesh_in_jobdir(self):
        """
        Copy all the data in the req_data variable into the job directory
        """
        for item in self.req_mesh:
            system('cp ' + '/'.join([self.meshdir,item]) + ' ' + self.jobdir + '/')

    def create_job_file(self,executionplan):
        """
        Create job file in the job directory

        :param executionplan: execution plan object
        :type executionplan: :class:`code.executionplan.executionplan.ExecutionPlan`
        """
        logging.info(multiprocessing.current_process().name+ ": writing configuration files for: "+ self.name)
        self.write_configurationfiles(executionplan)

        logging.info(multiprocessing.current_process().name + ": creating jobfile for: " + self.name)
        self.set_job_header()
        self.set_job_commands(executionplan)
        jobheader = '\n'.join(self.job_header)
        jobcommands = '\n'.join(self.job_commands)
        jobcontent = '\n'.join([jobheader, jobcommands])
        filename = self.name + ".job"
        self.write_file_jobdir(filename,jobcontent)

    def write_file_jobdir(self, filename, content):
        """
        Write file with specific filename in job directory

        :param filename: name of the file
        :type filename: str
        :param content: content of the file
        :type content: str
        """

        self.job_filename ='/'.join([self.jobdir, filename])
        with open(self.job_filename, 'w') as f:
            f.write(content)

    def set_job_header(self):
        raise NotImplementedError

    def execute(self,executionplan):
        raise NotImplementedError

    def file_exists(self,filename):
        """
        Check if a file exists

        :param filename: name of the file
        :type filename: str
        :return: True if the file exist otherwise False
        :rtype: bool
        """
        return path.exists(filename)

