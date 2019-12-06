import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pandas as pd

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


box_directions = ['x', 'y', 'z']

f, axes = plt.subplots(nrows=3, ncols=1, dpi=300, figsize=(20,20))
for index, box_direction in enumerate(box_directions):
    datafile = "output_concentration_profile_" + box_direction + ".csv"
    print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)
    cmap = get_cmap(n=df["real_time"].nunique())

    lgd = []
    for step_id, real_time in enumerate(df["real_time"].unique()):
        df_step = df[df["real_time"]==real_time]
        axes[index].plot(df_step[box_direction], df_step["c"], c=cmap(step_id), marker='o', markersize=2, ls="")
        lgd.append("real_time: " + str(real_time) + "(copss)")
        axes[index].plot(df_step[box_direction], df_step["c_exact"], c=cmap(step_id), linewidth=2)
        lgd.append("real_time: " + str(real_time) + "(exact)")

    #     axes[index].legend(lgd)
    axes[index].set_xlabel(box_direction, fontsize=25)
    axes[index].set_ylabel("Concentration", fontsize=25)

axes[0].set_title("Concentration: c (line: exact solution; dot: copss solution)")

# plt.show();
plt.savefig('concentration_comparison.eps', format='eps', dpi=300)