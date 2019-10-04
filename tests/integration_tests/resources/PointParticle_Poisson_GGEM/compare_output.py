
import pandas as pd
import os

def compare_output():
    cwd = os.getcwd()
    outputs = ["output_potential_profile_x_3P.csv", 
                "output_potential_profile_y_3P.csv", 
                "output_potential_profile_z_3P.csv"]
    for output in outputs:
        df = pd.read_csv(cwd + "/" + output)
        df_benchmark = pd.read_csv(cwd + "/output/" + output)
        if not df.equals(df_benchmark):
            print ("Test failed for file: '{0}'".format(output))
            return False, output
            
    return True, None
