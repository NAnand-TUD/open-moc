from uuid import uuid4
from os import system, WEXITSTATUS, popen
from job import *
import logging
import multiprocessing

class LocalJob(Job):
    """This is the local job
    """
    def __init__(self, executionplan, simulation_id):
        Job.__init__(self,executionplan, simulation_id)
    def set_job_header(self):
        self.job_header.append('#!/bin/bash')

    def execute(self):
        logging.info(multiprocessing.current_process().name + ": creating local_job "+ self.name)
        self.create_job_file(self.executionplan)
        filename = self.job_filename.split('/')[-1]
        directory = self.jobdir
        logging.info(multiprocessing.current_process().name + ": executing local_job "+ self.name)
        call = system('cd '+ directory +' && bash '+ filename + " & ")
        if (WEXITSTATUS(call) == 0):
            logging.info(multiprocessing.current_process().name +": succesfully executed the local_job: " + self.name)

    def execute_and_wait(self):
        logging.info(multiprocessing.current_process().name + ": creating local_job "+ self.name)
        self.create_job_file(self.executionplan)
        filename = self.job_filename.split('/')[-1]
        directory = self.jobdir
        logging.info(multiprocessing.current_process().name + ": executing local_job "+ self.name)
        call = system('cd '+ directory +' && bash '+ filename)
        if (WEXITSTATUS(call) == 0):
            logging.info(multiprocessing.current_process().name +": succesfully executed the local_job: " + self.name)

    def process_is_running(self, processname):
        return not popen("ps -aef | grep -i '"+processname+"' | grep -v 'grep' | awk '{ print $3 }'").read().strip().split('\n') == ['']

    def return_job(self):
        result_file = self.executionplan.get_last_executionunit_file()
        log_file = self.executionplan.get_su2logfile()

        logging.info(multiprocessing.current_process().name +
                     ": waiting for local_job: "+ self.name)
        result_filename = "/".join([self.jobdir, result_file])
        log_filename = "/".join([self.jobdir, log_file])
        try:
            content = self.wait_for_file_and_return_content(result_filename, log_filename, "SU2_CFD")
            logging.info(multiprocessing.current_process().name +
                         ": received file for local_job: " + self.name)
            dict = self.convert_content_in_dict(content)
        except:
            dict={"totaltostatic_efficiency": 0, "outletpressure_error": 0}


        return dict

   # def return_blade_specs(self):


    #def return_blade_specs(filename):
    #    def ReadFile(name):
    #        IN = {}
    #        infile = open(name, 'r')
    #        for line in infile:
    #            words = re.split(r'=| |%|\n|#', line)
    #            if not any(words[0] in s for s in ['\n', '%', ' ', '#']):
    #                words = list(filter(None, words))
    #                IN[words[0]] = words[1]
    #        return IN

    #   blade_specs = ReadFile(filename)

    #    return blade_specs['Throat'], blade_specs['TE_coords']

    # def return_blade_specs(self):

    def su2cfd_solution_converged(self, su2cfd_log_filename):
        fob = open(su2cfd_log_filename, 'r')
        content = fob.readlines()
        return "Succes" in content[-2]


    def wait_for_file_and_return_content(self, result_filename, log_filename, processname):

        logging.info(multiprocessing.current_process().name +
                     ": is waiting for: "+ log_filename)
        while not self.file_exists(log_filename):
            sleep(10)


        logging.info(multiprocessing.current_process().name +
                     ": is waiting for the process: " + processname)
        while self.process_is_running(processname):
            sleep(10)

        logging.info(multiprocessing.current_process().name +
                     ": checks if the solution is converged: " + processname)
        if not self.su2cfd_solution_converged(log_filename):
            raise ValueError("%s the process has converged" % result_filename)

        logging.info(multiprocessing.current_process().name +
                     ": waits for the postprocessing results: " + processname)

        while not self.file_exists(result_filename):
            sleep(10)

        logging.info(multiprocessing.current_process().name +
                     ": reads the postprocessing results: " + processname)
        if path.isfile(result_filename):
            fob = open(result_filename, 'r')
            content = fob.readlines()
            return content
        else:
            raise ValueError("%s the postprocess has failed" % result_filename)

    #def return_stator_specs(self,specs_filename):

