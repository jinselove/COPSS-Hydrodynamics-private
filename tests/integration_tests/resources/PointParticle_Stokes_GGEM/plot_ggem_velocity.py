import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pandas as pd

box_directions = ['x', 'y', 'z']
velocity_directions = ['u', 'v', 'w']
sol_types = np.array([ "total", "global", "local", "exact" ])


f, axes = plt.subplots(nrows=3, ncols=3, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]

# loop over three data files
for index, box_direction in enumerate(box_directions):
    datafile = "output_velocity_profile_" + box_direction + "_3P.csv"
    print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)

    # first loop over three directions
    for i, vd in enumerate(velocity_directions):
        # second loop over all solution types
        for j, st in enumerate(sol_types):
            # get solution name, i.e, column name in df
            sol_name = vd + "_" + st
            if j<3:
                ls = '-'
            else:
                ls = '--'
            # plot this solution
            axes[index, i].plot(df[box_direction], df[sol_name], label=sol_name,
                                color=colorplot[j], linewidth=2, ls=ls)
            axes[index, i].legend(sol_types)
            axes[index, i].set_xlabel(box_direction, fontsize=25)
            axes[index, i].set_ylabel(vd, fontsize=25)

#plt.show()
plt.savefig('ggem_pointforce_validation.eps', format='eps', dpi=100)

