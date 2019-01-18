#!/usr/bin/python
###############################################
# Integration test for COPSS-Hydrodynamics
# Author: Jiyuan Li, Xikai Jiang
#
# Required dependencies:
#    1. COPSS executables in $COPSS_SRC
#    2. Python 2.7
###############################################
import sys
import subprocess
import os
import json
import time
import filecmp
import multiprocessing

##################################################
# Helper functions
##################################################

def compile(package, decompile=False):
    """
    :param package: package name to compile, i.e., "POINTPARTICLE", "RIGIDPARTICLE", etc.
    :param decompile: if True, clean all compiled output
    """
    if decompile:
        subprocess.call("bash " + compile_tool + " -p " + package + " -a " + "clean_only",
                        shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        return None
    
    print("--"*2 + "Recompiling package '{0}' if changes were made in src folder".format(package))
    returnvalue = subprocess.call("bash " + compile_tool + " -p "+ package, shell=True,
                    stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if returnvalue != 0:
            print("--"*2 + "Error: Compilation exited with non-zero exit code {}".format(returnvalue))
            subprocess.call("bash " + compile_tool + " -p " + package + " -a " + "clean_only",
                            shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            sys.exit(returnvalue)
    print("--"*2 + "Compilation done !!!")
    
def compare_output(validation_outputs):
    cwd = os.getcwd()
    for validation_output in validation_outputs:
        benchmark_file = cwd + '/output/' + validation_output
        file = cwd + "/" + validation_output
        if not filecmp.cmp(benchmark_file, file):
            return False, file
    
    return True, None

##################################################
# define global parameter
##################################################
print("Notice: Results generated using different number of CPUs are slightly different due to\
 numerical causes. For simplicity, we only test the results using 4 cpus (number of cpus is\
 defined in 'run.sh' in each test folder).")
if multiprocessing.cpu_count() < 4:
    print("Error: Available number of CPUs ({0}) is less than required.".format(multiprocessing.cpu_count()))
    sys.exit(0)
    
###################################################
# Load tests from test.json
###################################################
print("")
print("Loading tests from 'test.json' ...")
if not os.path.isfile("test.json"):
    print("--"*2 + "Error: 'test.json' does not exist in current folder")
    sys.exit(1)

with open("test.json", "r") as read_file:
    tests = json.load(read_file)

##################################################
# Prepare test environment
##################################################

print("Preparing testing environment ...")
# check if $COPSS_DIR exists
print("--"*1 + "Checking if $COPSS_DIR exists ...")
copss_dir = os.environ.get("COPSS_DIR")
if not copss_dir:
    print("--"*2 + "Error: Unable to locate COPSS directory. Exiting ...")
    print("--"*2 + "Suggestion: add 'export COPSS_DIR=xxxx' to system environment")
    sys.exit(-1)
else:
    print("--"*2 + "COPSS_DIR={0}".format(copss_dir))
    
    
# check if compile tool exists
compile_tool = copss_dir + "/tools/compile.sh"
if not os.path.isfile(compile_tool):
    print("--"*1 + "Error: compile tool: '{0}' does not exist".format(compile_tool))
    sys.exit(-1)
    
    
# check if all test resources are existed
print("--"*1 + "Checking all required files exist for all tests ...")
for key, value in tests.items():
    print("--"*2 + "Checking for test --> {}".format(key))
    for input_file in value["required_inputs"]:
        if not os.path.isfile(copss_dir + "/tests/integration_tests/resources/" \
                              + key + "/" + input_file):
            print("--"*3 + "Error: Required input file '{}' does not exist. Exiting ...".format(input_file))
            sys.exit(-1)
    for output_file in value["validation_outputs"]:
        if not os.path.isfile(copss_dir + "/tests/integration_tests/resources/" \
                              + key + "/output/" + output_file):
            print("--"*3 + "Error: Validation output file '{}' does not exist. Exiting ...".format(output_file))
            sys.exit(-1)

##################################################
# Perform integration tests
##################################################


print("Performing tests ...")
os.chdir(copss_dir+'/tests/integration_tests')
for key, value in tests.items():    
    print("--"*1 + "Performing test --> '{0}'".format(key))
    test_dir = copss_dir + "/tests/integration_tests/resources/" + key
    os.chdir(test_dir)
    
    # recompile package
    compile(value["package"])
    
    # run simulation
    start = time.time()
    returnvalue = subprocess.call("bash run.sh", shell=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    end = time.time()
    if returnvalue != 0:
        print("--"*2 + "Error: COPSS exited with non-zero exit code {}".format(returnvalue))
        subprocess.call('bash zclean.sh', shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sys.exit(returnvalue)
    # compile output files defined in the json file
    is_output_same, diff_file = compare_output(value['validation_outputs'])
    if not is_output_same:
        print("--"*2 + "Error: At least output file '{0}' is not the same with benchmark. Exiting ...".format(diff_file))
        print("--"*2 + "Suggestion: Manually run this test and figure out the differences.")
        subprocess.call('bash zclean.sh', shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sys.exit(-1)
    print("--"*1 + "Timing for test --> '{0}': {1} s".format(key, end-start))
    subprocess.call('bash zclean.sh', shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)

print("All integration tests passed. Exiting ...")
sys.exit(0)
