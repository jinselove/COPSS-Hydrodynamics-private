import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pandas as pd

box_directions = ['x', 'y', 'z']
sol_types = ["total", "global", "local", "exact"]

f, axes = plt.subplots(nrows=3, ncols=1, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]

for index, box_direction in enumerate(box_directions):
    datafile = "output_potential_profile_" + box_direction + "_3P.csv"
    print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)

    # loop over all solution types
    for i, st in enumerate(sol_types):
        # get solution name, i.e, colume name in df
        sol_name = "phi_" + st
        if i<3:
            ls = '-'
        else:
            ls = '--'
        # plot this solution
        axes[index].plot(df[box_direction], df[sol_name], label=sol_name,
                         color=colorplot[i], linewidth=2, ls=ls)
        axes[index].legend(sol_types)
        axes[index].set_xlabel(box_direction, fontsize=25)
        axes[index].set_ylabel("Potential", fontsize=25)

#plt.show()
plt.savefig('ggem_pointcharge_validation.eps', format='eps', dpi=100)

