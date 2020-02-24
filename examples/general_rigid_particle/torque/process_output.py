#!/bin/env python

import glob
import numpy as np
import pandas as pd
from pylab import *
import matplotlib.pyplot as plt
from numpy import linalg as LA
import os
from scipy.spatial import distance
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
from scipy import special, optimize
from scipy.interpolate import InterpolatedUnivariateSpline
import sys

# add cwd to path when script run by slurm
sys.path.append(os.getcwd())

### initial parameters
n_spheres = 9
print ("### number of spheres = {0}".format(n_spheres))
a = 3
r = (2./3.)**(1./3.) * a
top_layer = [12,13,14,15,16,17]
center_layer = [2,4,6,8,10,11]
bottom_layer = [0,1,3,5,7,9]
top_center = 19
bottom_center = 18
color_list = []
for color in colors.cnames:
    color_list.append(color)

# ### 0.2 define helper functions
# 2)preprocess df_time, df_pnode, df_pcenter
def preprocess_dataframe(df_time,df_pcenter,df_pnode):
    # compute total number of nodes among all particles
    num_of_nodes = 0
    for particle_id in df_pnode.particle_id.unique():
        num_of_nodes += df_pnode[df_pnode.particle_id==particle_id].node_id.unique().shape[0]
    num_of_particles = df_pnode.particle_id.unique().shape[0]
    # add "o_step" column to df_pcenter
    df_pcenter["o_step"] = pd.Series(df_pcenter.index.values//num_of_particles,index=df_pcenter.index)
    # add "particle_type_id" column to df_pcenter
    df_pcenter["particle_type_id"] = pd.Series(df_pcenter["particle_id"].values >= n_spheres, index=df_pcenter.index) * 1
    # add "o_step" column to df_node
    df_pnode["o_step"] = pd.Series(df_pnode.index.values // num_of_nodes,index=df_pnode.index)
    # add "particle_type_id" column to df_pnode
    df_pnode["particle_type_id"] = pd.Series(df_pnode["particle_id"].values >= n_spheres, index=df_pnode.index) * 1
    # the minimum final o_step
    final_o_step = np.min([df_pcenter.o_step.values[-1],df_pnode.o_step.values[-1],df_time.o_step.values[-1]])
    print(final_o_step)
    df_time = df_time[df_time.o_step<=final_o_step]
    df_pcenter = df_pcenter[df_pcenter.o_step<=final_o_step]
    df_pnode = df_pnode[df_pnode.o_step<=final_o_step]
    num_of_steps=final_o_step
    print(df_pnode.tail(1))
    print(df_pcenter.tail(1))
    print(df_time.tail(1))
    # add "real_time" column to df_pnode and df_pcenter
    real_time_matrix = df_time.real_time.values
    df_pcenter["real_time"]=df_pcenter["o_step"].apply(lambda x: real_time_matrix[x])
    df_pnode["real_time"] = df_pnode["o_step"].apply(lambda x: real_time_matrix[x])
    return num_of_particles, num_of_nodes, num_of_steps



# ### 1.  Analyze data
# **1. process dataframe**
datapath = os.getcwd()
timefile = "out.time"
pcenterfile = "output_particle_trajectory.csv"
pnodefile = "output_surface_node.csv"
print ("analyze timefile: %s"% timefile)
print ("analyze particle center file: %s"%pcenterfile)
print ("analyze surface node file: %s"%pnodefile)
df_time = pd.read_csv(timefile,delimiter=" ",engine='python')
df_time = df_time.dropna(axis=1)
df_time.columns=["step_id","o_step","real_time"]
df_pcenter = pd.read_csv(pcenterfile,delimiter=" ",engine='python',skiprows=1, header=None,comment="#")
df_pcenter = df_pcenter.dropna(axis=1)
df_pcenter.columns=["particle_id","x_coord","y_coord","z_coord","x_vel","y_vel","z_vel","x_force","y_force","z_force"]
df_pnode = pd.read_csv(pnodefile,delimiter=" ",engine='python',skiprows=1, header=None,comment="#")
df_pnode = df_pnode.dropna(axis=1)
df_pnode.columns=["particle_id","node_id","x_coord","y_coord","z_coord","node_center_distance"]
# process dataframe
num_of_particles, num_of_nodes, num_of_steps = preprocess_dataframe(df_time,df_pcenter,df_pnode)
print ("number_of_particles = %d, num_of_nodes = %d, num_of_steps = %d" %(num_of_particles,num_of_nodes,num_of_steps))
# check node_centerline_distance, error, moment_of_inertia for multiple steps before fit the cylinders
# compute center line unit vector or alltimesteps and append it to df_pcenter
#compute_centerline_unit()
# analyse msd starting from step ...
#o_step = 10
df_pcenter.to_csv(datapath+"/df_pcenter.csv", index=False)
df_pnode.to_csv(datapath+"/df_pnode.csv", index=False)
#df_pcenter[df_pcenter["o_step"]==o_step].to_csv(datapath+"/df_pcenter_step_"+str(o_step)+".csv", index=False)
#df_pnode[df_pnode["o_step"]==o_step].to_csv(datapath+"/df_pnode_step_"+str(o_step)+".csv", index=False)

for o_step in df_pcenter["o_step"].unique():
    df_pcenter[df_pcenter["o_step"]==o_step].to_csv(datapath+"/df_pcenter_step_{}.csv".format(str(o_step).zfill(4)), index=False)
    df_pnode[df_pnode["o_step"]==o_step].to_csv(datapath+"/df_pnode_step_{}.csv".format(str(o_step).zfill(4)), index=False)

print("done!")
