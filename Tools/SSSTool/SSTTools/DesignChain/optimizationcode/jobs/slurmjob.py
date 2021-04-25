from os import system, popen, path
from job import *
from secrets import *
import re

class SlurmJob(Job):
    def __init__(self, ncores, nnodes, esttime, executionplan,simulation_id ):
        Job.__init__(self, executionplan,simulation_id)
        self.nnodes = nnodes
        self.ncores = ncores
        self.esttime = esttime
        self.username = username
        self.create_uid()

              
    def set_job_header(self):
        self.job_header.append('#!/bin/bash')
        self.job_header.append('#SBATCH -n ' + str(self.ncores))
        self.job_header.append('#SBATCH -N ' + str(self.nnodes))
        self.job_header.append('#SBATCH -t ' + self.get_slurmjob_timeformat())


    def get_slurmjob_timeformat(self):
        hours = int(self.esttime)
        minutes = (self.esttime*60) % 60
        seconds = (self.esttime*3600) % 60
        return "%d:%02d:%02d" % (hours, minutes, seconds)
       
    def execute(self):
        self.create_job_file(self.executionplan)
        filename = self.job_filename.split('/')[-1]

        directory = self.jobdir
        #print('cd '+ directory + ' && sbatch '+ filename)
        system('cd '+ directory + ' && sbatch '+ filename)

    def parse_slurmjob_status(self,returnstring):
        listofstring = returnstring.split("\n")
        try:
            return listofstring[-2].split()[5]
        except:
            return "job not found"

    def get_slurmjob_status(self):
        call = popen("sacct -u "+self.username+" --name="+self.name+".job").read()
        logging.debug(multiprocessing.current_process().name +
                     ": gets the following status back from slurm for job " + self.name +": " + call)
        return self.parse_slurmjob_status(call)

    def slurmjob_running(self):
        status = self.get_slurmjob_status()
        logging.debug(multiprocessing.current_process().name +
                     ": gets the following status for job " + self.name +": " + status)
        return (status == "RUNNING" or status == "PENDING")

    def slurmjob_succesful(self):
        status = self.get_slurmjob_status()
        logging.debug(multiprocessing.current_process().name +
                     ": gets the following status for job " + self.name +": " + status)
        return (status == "COMPLETED")

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
            dict={"totaltostatic_efficiency": -1000, "outletpressure_error": 1000, "theta_velocity": 1000}


        return dict

    def su2cfd_solution_converged(self, su2cfd_log_filename):
        fob = open(su2cfd_log_filename, 'r')
        content = fob.readlines()
        return "Succes" in content[-2]

    def wait_for_file_and_return_content(self, result_filename, log_filename, processname):


        logging.info(multiprocessing.current_process().name +
                     ": is waiting for the slurmjob to finish: " + self.name)
        sleep(5)
        while self.slurmjob_running():
            sleep(10)

        logging.info(multiprocessing.current_process().name +
                     ": checks if the job was succesful: "+  self.name)
        if not self.slurmjob_succesfull():
            raise ValueError("%s the job has failed" % self.name)

        logging.info(multiprocessing.current_process().name +
                     ": waits for the log file of the cfd simulation: " + self.name)

        while not self.file_exists(log_filename):
            sleep(10)

        logging.info(multiprocessing.current_process().name +
                     ": checks if the simulation has converged: " + self.name)
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

    def return_blade_specs(self, filename):
        logging.info(multiprocessing.current_process().name +
                     ": is waiting for the slurmjob to finish: " + self.name)
        sleep(5)
        while self.slurmjob_running():
            sleep(10)

        specs = {}
        with open(path.join(self.jobdir,filename), 'r') as specfile:
            for line in specfile:
                words = re.split(r'=| |%|\n|#', line)
                if not any(words[0] in s for s in ['\n', '%', ' ', '#']):
                    words = list(filter(None, words))
                    specs[words[0]] = words[1]

        return specs['Throat'], specs['TE_coords'], specs['ScaleNoz']
