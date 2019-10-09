from src.validate import Validator
from src.compare import compare_outputs
from src.copss import COPSS
import os


def main():
    # initialize validator
    validator = Validator()
    validator.validate()
    # initialize copss
    copss = COPSS(validator.copss_dir, validator.compile_tool)
    # change directory
    os.chdir(validator.copss_dir + "/tests/integration_tests")
    # loop over all tests
    for test_name, test_dict in validator.tests.items():
        if not test_dict['on']:
            print("####" + "Skipping test '{0}'####".format(test_name))
            continue
        print("####"+ "Testing '{0}'####".format(test_name))
        test_dir = copss.copss_dir + "/tests/integration_tests/resources/" + \
                   test_name

        # recompile package if needed
        print("----recompiling package '{0}' if changes were made in src "
                     "folder  ... ".format(test_dict['package']), end='')
        copss.compile(package=test_dict["package"])

        # run simulation
        print("----running test simulation and generating output file ... ",
              end='')
        copss.run(test_dir=test_dir)

        # loop over all output file defined in required_outputs and validate
        print("----comparing all output files ... ", end='')
        compare_outputs(test_dir, test_dict["required_outputs"])
        print("----test {0} succeed!".format(test_name))

        # clean files generated in test
        COPSS.clean(test_dir)

if __name__ == "__main__":
    main()
