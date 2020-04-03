import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pandas as pd

box_directions = ['x', 'y', 'z']
sol_types = ["total", "global", "local", "exact"]
# sol_types = ["global", "exact"]

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
        axes[index].set_xlabel(box_direction, fontsize=25)
        axes[index].set_ylabel("Potential", fontsize=25)
        axes[index].legend()
        

# plt.show()
plt.savefig('ggem_pointcharge_validation.eps', format='eps', dpi=100)

f, axes = plt.subplots(nrows=3, ncols=1, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]

for index, box_direction in enumerate(box_directions):
    datafile = "output_potential_profile_" + box_direction + "_3P.csv"
    print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)

    # loop over all solution types
    for i, st in enumerate(sol_types):
        # get solution name, i.e, colume name in df
        sol_name = "phi_grad_" + st + "(x)"
        if i<3:
            ls = '-'
        else:
            ls = '--'
        # plot this solution
        axes[index].plot(df[box_direction], df[sol_name], label=sol_name,
                         color=colorplot[i], linewidth=2, ls=ls)
        axes[index].set_xlabel(box_direction, fontsize=25)
        axes[index].set_ylabel("Potential", fontsize=25)
        axes[index].legend()
        

# plt.show()
plt.savefig('ggem_pointcharge_validation_gradient_along_x.eps', format='eps', dpi=100)

f, axes = plt.subplots(nrows=3, ncols=1, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]

for index, box_direction in enumerate(box_directions):
    datafile = "output_potential_profile_" + box_direction + "_3P.csv"
    print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)

    # loop over all solution types
    for i, st in enumerate(sol_types):
        # get solution name, i.e, colume name in df
        sol_name = "phi_grad_" + st + "(y)"
        if i<3:
            ls = '-'
        else:
            ls = '--'
        # plot this solution
        axes[index].plot(df[box_direction], df[sol_name], label=sol_name,
                         color=colorplot[i], linewidth=2, ls=ls)
        axes[index].set_xlabel(box_direction, fontsize=25)
        axes[index].set_ylabel("Potential", fontsize=25)
        axes[index].legend()
        

# plt.show()
plt.savefig('ggem_pointcharge_validation_gradient_along_y.eps', format='eps', dpi=100)

f, axes = plt.subplots(nrows=3, ncols=1, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]

for index, box_direction in enumerate(box_directions):
    datafile = "output_potential_profile_" + box_direction + "_3P.csv"
    print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)

    # loop over all solution types
    for i, st in enumerate(sol_types):
        # get solution name, i.e, colume name in df
        sol_name = "phi_grad_" + st + "(z)"
        if i<3:
            ls = '-'
        else:
            ls = '--'
        # plot this solution
        axes[index].plot(df[box_direction], df[sol_name], label=sol_name,
                         color=colorplot[i], linewidth=2, ls=ls)
        axes[index].set_xlabel(box_direction, fontsize=25)
        axes[index].set_ylabel("Potential", fontsize=25)
        axes[index].legend()
        

# plt.show()
plt.savefig('ggem_pointcharge_validation_gradient_along_z.eps', format='eps', dpi=100)


f, axes = plt.subplots(nrows=3, ncols=1, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]
sol_types = ["total", "exact"]

for index, box_direction in enumerate(box_directions):
    datafile = "output_potential_profile_" + box_direction + "_3P.csv"
    print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)

    # loop over all solution types
    for i, st in enumerate(sol_types):
        # get solution name, i.e, colume name in df
        sol_name = "phi_laplacian_" + st
        if i<1:
            ls = '-'
        else:
            ls = '--'
        # plot this solution
        axes[index].plot(df[box_direction], df[sol_name], label=sol_name,
                         color=colorplot[i], linewidth=2, ls=ls)
        axes[index].set_xlabel(box_direction, fontsize=25)
        axes[index].set_ylabel("Potential", fontsize=25)
        axes[index].legend()


# plt.show()
plt.savefig('ggem_pointcharge_validation_laplacian.eps', format='eps', dpi=100)
