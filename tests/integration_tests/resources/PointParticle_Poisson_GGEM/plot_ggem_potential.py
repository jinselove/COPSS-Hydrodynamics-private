import numpy as np
import matplotlib.pyplot as plt
import math


box_directions = ['x', 'y', 'z']

f, axes = plt.subplots(nrows=3, ncols=1, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]
vtype = np.array([ "Total", "Global", "Local", "Exact" ])
labelplot = [ str( vt ) for vt in vtype ];

for index, box_direction in enumerate(box_directions):
    datafile = "output_potential_profile_" + box_direction + "_3P.txt"
    print ("plot data in {0}".format(datafile))
    tmpA = np.loadtxt(datafile)
    mA = tmpA.shape[0]
    nA = tmpA.shape[1]
    # coordinates of cross section
    x  = tmpA[:,0];      # coordinates
    vx = np.zeros( (mA,4), np.float ); # total | global | local | exact
                                       #           |        |
                                       #          FEM     Green
    for i in range(0,4):
       vx[:,i] = tmpA[:,i+1]
    for i in range(0,4):
       if i<3:
           axes[index].plot(x, vx[:,i], color=colorplot[i], linewidth=2);
           #axes[index].plot(x, vx[:,1], color=colorplot[i], linewidth=2);
       else :
           axes[index].plot(x, vx[:,i], color=colorplot[i], linewidth=2, ls='--');
           #axes[index].plot(x, vx[:,1], color=colorplot[i], linewidth=2, ls='--');
    axes[index].legend(vtype)
    axes[index].set_xlabel(box_direction, fontsize=25)
    axes[index].set_ylabel("Potential", fontsize=25)

axes[0].set_title("Potential: phi")

plt.show()
plt.savefig('ggem_pointcharge_validation.eps', format='eps', dpi=100)

