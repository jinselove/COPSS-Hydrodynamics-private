import jsonschema
import json
import os
import sys

class Validator:
    def __init__(self):
        self.json_schema_dir = "json/test_schema.json"
        self.test_schema = None
        self.json_dir = "json/test.json"
        self.tests = None
        self.copss_dir = None
        self.compile_tool = None

    def validate(self):
        """validate all test definitions and environment"""
        # load json schema file for tests
        print("Loading test schema from '{0}' ...".format(
            self.json_schema_dir),  end='')
        self.test_schema = self.load_json(self.json_schema_dir)
        print("Succeed!")

        # load json file for tests
        print("Loading tests from '{0}' ...".format(self.json_dir), end='')
        self.tests = self.load_json(self.json_dir)
        print("Succeed!")

        # validate tests against test_schema
        print("Validating tests against test schema ...")
        self.validate_tests()
        print("Validating tests against test schema Succeed!")

        # validate COPSS
        print("Validating COPSS ...")
        self.validate_copss()
        print("Validating COPSS Succeed!")

        # validate resources
        print("Validating resources ...")
        self.validate_resource()
        print("Validating resources Succeed!")

    def load_json(self, file_dir):
        """ load json file

        :param file_dir: absolute directory to a json file
        :return: dict: a dictionary
        """
        if not os.path.isfile(file_dir):
            print("\nError: {0} does not exist")
            sys.exit(-1)
        else:
            with open(file_dir, 'r') as read_file:
                dict = json.load(read_file)
                return dict

    def validate_tests(self):
        """ validate all tests

        :param tests: a dictionary of tests
        :param test_schema: a dictionary of test schema
        :return: None
        """
        for test_key, test in self.tests.items():
            print("--Validating test {0} ...".format(test_key), end="")
            try:
                jsonschema.validate(test, self.test_schema)
                print("Succeed!")
            except SystemError as err:
                print("Fail!")
                print("Please double check this test in test.json. Error is: {0}".format(err))
                sys.exit(-1)

    def validate_copss(self):
        """ Validate COPSS environement

        :return: None
        """
        print("--"*1 + "Checking if $COPSS_DIR exists ...", end='')
        self.copss_dir = os.environ.get("COPSS_DIR")
        if not self.copss_dir:
            print("\nError: Unable to locate COPSS directory. Exiting ...")
            print("Suggestion: add 'export COPSS_DIR=xxxx' to system environment")
            sys.exit(-1)
        else:
            print("True, COPSS_DIR={0}".format(self.copss_dir))

        # check if compile tool exists
        self.compile_tool = self.copss_dir + "/tools/compile.sh"
        print("--"*1 + "Checking if COPSS compile tool exists ...", end='')
        if not os.path.isfile(self.compile_tool):
            print("--"*1 + "Error: compile tool: '{0}' does not "
                           "exist".format(self.compile_tool))
            sys.exit(-1)
        print("True")

    def validate_resource(self):
        """ Validate if all required resources files exist for tests

        :param tests: tests dictionary loaded from test.json
        :return: None
        """
        # check if all test resources are existed
        for test_name, test_dict in self.tests.items():
            print("--"*1 + "Checking resources for test --> {} ...".format(
                test_name), end='')
            # check 'required_inputs'
            for ifile in test_dict["required_inputs"]:
                if not os.path.isfile(self.copss_dir +
                                      "/tests/integration_tests/resources/" \
                                      + test_name + "/" + ifile):
                    print("--"*3 + "Error: Required input file '{}' does not "
                                   "exist. Exiting ...".format(ifile))
                    sys.exit(-1)
            # check 'validation_outputs'
            for ofile_dict in test_dict["required_outputs"]:
                if not os.path.isfile(self.copss_dir +
                                      "/tests/integration_tests/resources/" \
                                      + test_name + "/output/" + ofile_dict[
                                          "filename"]):
                    print("--"*3 + "Error: Validation output file '{}' does "
                                   "not exist.  Exiting ...".format(
                        ofile_dict["filename"]))
                    sys.exit(-1)
            print("Succeed!")

