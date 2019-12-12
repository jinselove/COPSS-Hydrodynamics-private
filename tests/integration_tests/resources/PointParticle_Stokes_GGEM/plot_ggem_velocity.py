import numpy as np
import matplotlib.pyplot as plt
import math


box_directions = ['x', 'y', 'z']

f, axes = plt.subplots(nrows=3, ncols=3, figsize=(20,20))
colorplot = ["red", "blue", "green", "cyan", "black" ]
vtype = np.array([ "Total", "Global", "Local", "Exact" ])
labelplot = [ str( vt ) for vt in vtype ];

for index, box_direction in enumerate(box_directions):
    datafile = "output_velocity_profile_" + box_direction + "_3P.txt"
    print ("plot data in {0}".format(datafile))
    tmpA = np.loadtxt(datafile)
    mA = tmpA.shape[0]
    nA = tmpA.shape[1]    
    # coordinates of cross section
    x  = tmpA[:,0];      # coordinates
    vx = np.zeros( (mA,4), np.float ); # total | global | local | exact
    vy = np.zeros( (mA,4), np.float ); #           |        |
    vz = np.zeros( (mA,4), np.float ); #          FEM      Green
    
    for i in range(0,4):
       vx[:,i] = tmpA[:,i*3+1]
       vy[:,i] = tmpA[:,i*3+2]
       vz[:,i] = tmpA[:,i*3+3]

    for i in range(0, 4):
       if i<3:
           axes[index, 0].plot(x, vx[:,i], color=colorplot[i], linewidth=2);
       else :
           axes[index, 0].plot(x, vx[:,i], color=colorplot[i], linewidth=2, ls='--');        
    axes[index, 0].legend(vtype)
    axes[index, 0].set_xlabel(box_direction, fontsize=25)
    axes[index, 0].set_ylabel("Velocity", fontsize=25)

    for i in range(0, 4):
       if i<3:
           axes[index, 1].plot(x, vy[:,i], label=labelplot[i], color=colorplot[i], linewidth=2);
       else :
           axes[index, 1].plot(x, vy[:,i], label=labelplot[i], color=colorplot[i], linewidth=2, ls='--');        
           
    axes[index, 1].set_xlabel(box_direction, fontsize=25)
    axes[index, 1].set_ylim([-0.05, 0.05])

    for i in range(0, 4):
       if i<3:
           axes[index, 2].plot(x, vz[:,i], label=labelplot[i], color=colorplot[i], linewidth=2);
       else :
           axes[index, 2].plot(x, vz[:,i], label=labelplot[i], color=colorplot[i], linewidth=2, ls='--');        
           
    axes[index, 2].set_xlabel(box_direction, fontsize=25)
    axes[index, 2].set_ylim([-0.05, 0.05])

axes[0, 0].set_title("velocity in x-direction: u")
axes[0, 1].set_title("velocity in y-direction: v")
axes[0, 2].set_title("velocity in z-direction: w")
#plt.show()
plt.savefig('ggem_pointforce_validation.eps', format='eps', dpi=100)

