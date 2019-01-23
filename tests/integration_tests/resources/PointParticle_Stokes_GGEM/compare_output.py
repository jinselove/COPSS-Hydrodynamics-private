
import pandas as pd
import os

def compare_output():
    cwd = os.getcwd()
    outputs = ["output_velocity_profile_x_3P.txt", 
                "output_velocity_profile_y_3P.txt", 
                "output_velocity_profile_z_3P.txt"]
                
    read_csv_params = {"delimiter": "  ",
                        "header": None,
                        "engine": "python"}
        
    read_csv_params_benchmark = read_csv_params                    
    for output in outputs:
        df = pd.read_csv(cwd + "/" + output, **read_csv_params)
        df_benchmark = pd.read_csv(cwd + "/output/" + output, **read_csv_params_benchmark)
        if not df.equals(df_benchmark):
            print ("Test failed for file: '{0}'".format(output))
            return False, output
            
    return True, None
