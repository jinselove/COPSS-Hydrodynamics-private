import sys
import filecmp
import pandas as pd
import numpy as np
from src.copss import COPSS


def compare_numeric_df(df_file, df_file_bm, tol):
    # drop NA values if there are any
    df_file.dropna(axis=1, inplace=True)
    df_file_bm.dropna(axis=1, inplace=True)
    # check the shapes of df_file and df_file_bm are not zero after dropping na
    if df_file.shape[0] == 0 or df_file.shape[1] == 0 \
            or df_file_bm.shape[0] == 0 or df_file_bm.shape[1] == 0:
        err_msg = "shape of dataframe for after dropping NA is 0".format(file)
        return False, err_msg
    # check if the shape of df_file and df_file_bm are the same
    if not df_file.shape==df_file_bm.shape:
        err_msg = "file and benchmark dataframe shape is different"
        return False, err_msg
    # check if all values are numeric
    if not df_file.apply(lambda s: pd.to_numeric(s, errors='coerce').notnull().all()).all():
        err_msg = "there are non-numeric values in file"
        return False, err_msg
    if not df_file_bm.apply(lambda s: pd.to_numeric(s, errors='coerce').notnull().all()).all():
        err_msg = "there are non-numeric values in benchmark file "
        return False, err_msg
    # check if the maximum absolute difference is within tolerance
    max_diff = np.max(np.abs((df_file - df_file_bm).values))
    if max_diff > tol:
        err_msg = "maximum absolute difference in file and benchmark file is " \
                  "larger than tolerance (" + str(max_diff) + " > " + str(
            tol) + ")"
        return False, err_msg

    return True, None

def compare_output(test_dir, output_dict):
    """compare one output file defined by output_dict

    """
    file = test_dir + "/" + output_dict["filename"]
    file_bm = test_dir + "/output/" + output_dict["filename"]
    cmp_method = output_dict["cmp_method"]
    # compare file directory if cmp_method = "cmp_file"
    if cmp_method == "cmp_file":
        if not filecmp.cmp(file, file_bm):
            return False, file
        return True, None, None
    # read file into dataframe and compare the numeric values if cmp_method = "cmp_df"
    elif cmp_method == "cmp_df":
        # read in file and benchmark file into dataframe
        try:
            df_file = pd.read_csv(file, **output_dict["params"])
        except:
            print("Error: unable to read '{0}' to dataframe using params '{1}'" \
                  .format(file, output_dict["params"]))
            sys.exit(-1)
        try:
            df_file_bm = pd.read_csv(file_bm, **output_dict["params"])
        except:
            print("Error: unable to read benchmark file '{0}' to dataframe using params '{1}'" \
                  .format(file_bm, output_dict["params"]))
            sys.exit(-1)

        # compare numeric dataframe
        is_same, err_msg = compare_numeric_df(df_file, df_file_bm,
                                              output_dict['tol'])
        return is_same, file, err_msg

    else:
        print("Error: Unavailable 'cmp_method specified'")
        sys.exit(-1)


def compare_outputs(test_dir, required_outputs):
    for output_dict in required_outputs:
        is_same, filename, err_msg = compare_output(test_dir,
                                                           output_dict)

        # if not same, output warning and exit
        if not is_same:
            print("--"*2 + "Failed to validate {0}. Exiting ...".format(filename))
            print("--"*2 + "Error Message: {0}".format(err_msg))
            print("--"*2 + "Suggestion: Manually run this test and figure out the differences.")
            # clean files generated in test
            COPSS.clean(test_dir)
            sys.exit(-1)
    print("Succeed!")
