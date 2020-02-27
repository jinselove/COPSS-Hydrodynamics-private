import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pandas as pd

setrun=2;
box_directions = ['x','z']
box_counter=[0,2]
sol_types = ["COPSS", "Comsol"]
# sol_types = ["global", "exact"]
j=0

comsol_dim_fac=10**(-2)
copss_charge_scaling= 1.8*10**(-3)
f, axes = plt.subplots(nrows=len(box_directions), ncols=1, figsize=(20,20))
colorplot = ["red", "green"]

for index, box_direction in enumerate(box_directions):

    datafile = "Neumann_bc_charged_" + box_direction + "_potential.csv"
    #print ("plot data in {0}".format(datafile))
    df = pd.read_csv(datafile)

    datafile_comsol = "data_along_" + box_direction + "_with_discrete_charge_comsol.txt"
    dat_comsol=np.loadtxt(datafile_comsol, skiprows=9,usecols=(box_counter[index],6))
    i=0
    ls = 'x'
        # plot this solution
    axes[index].plot(dat_comsol[:,0], dat_comsol[:,1],color=colorplot[i], marker="x",linestyle = 'None')
    axes[index].set_xlabel(box_direction+ "($\mu$ m)", fontsize=25)
    axes[index].set_ylabel("Potential", fontsize=25)
    axes[index].legend(['COMSOL','COPSS'], fontsize=25)

    sol_name = "phi"

    
    i=1
    ls = '-'       
        # plot this solution
    axes[index].plot(df[box_direction]*comsol_dim_fac, df[sol_name]*copss_charge_scaling, color=colorplot[i], linewidth=2, ls=ls)
    axes[index].set_xlabel(box_direction + "($\mu$ m)", fontsize=25)
    axes[index].set_ylabel("Potential", fontsize=25)
    axes[index].legend(['COMSOL','COPSS'], fontsize=25)
        

# plt.show()
plt.savefig('Poisson_Neumann_BC_validation.eps', format='eps', dpi=100)

