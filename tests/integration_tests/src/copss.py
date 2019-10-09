import os
import sys
import subprocess
import multiprocessing
import time


class COPSS:
    def __init__(self, copss_dir, compile_tool):
        self.copss_dir = copss_dir
        self.compile_tool = compile_tool
        self.validate_cpu()

    def validate_cpu(self):
        """validate the number of CPU available"""
        if multiprocessing.cpu_count() < 4:
            print("Error: Available number of CPUs ({0}) is less than required.".format(multiprocessing.cpu_count()))
            sys.exit(-1)

    def compile(self, package='POINTPARTICLE', decompile=False):
        """ compile COPSS

        :param package: package name
        :param decompile: if decompile
        :return: None
        """
        if decompile:
            subprocess.call("bash " + self.compile_tool + " -p " + package + " "
                            "-a " + "clean_only", shell=True, stdout =
                            subprocess.PIPE, stderr = subprocess.PIPE)
            return None

        returnvalue = subprocess.call("bash " + self.compile_tool + " -p "+
                                      package, shell=True, stdout =
                                      subprocess.PIPE, stderr = subprocess.PIPE)
        if returnvalue != 0:
            print("--" + "Error: Compilation exited with non-zero exit code "
                         "{0}".format(returnvalue))
            subprocess.call("bash " + self.compile_tool + " -p " + package + " "
                            "-a " + "clean_only", shell=True, stdout =
                            subprocess.PIPE, stderr = subprocess.PIPE)
            sys.exit(returnvalue)
        print("Succeed!")

    def run(self, test_dir):
        # run simulation
        os.chdir(test_dir)

        start = time.time()
        returnvalue = subprocess.call("bash run.sh", shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        end = time.time()
        # if COPSS end with nonzero exit code, it means there is a bug in the code
        if returnvalue != 0:
            print("\nError: COPSS exited with non-zero exit code {}".format(
                returnvalue))
            subprocess.call('bash zclean.sh', shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sys.exit(returnvalue)
        print("Succeed! execution time = {0:.0f} s".format(end-start))

    @staticmethod
    def clean(test_dir):
        os.chdir(test_dir)
        subprocess.call('bash zclean.sh', shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)




