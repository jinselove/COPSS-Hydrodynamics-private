
# coding: utf-8

# In[1]:

import numpy as np


# In[15]:

output_stat = np.loadtxt("output_statistics.dat",skiprows=1)
output_stat.shape
output_stat[1]


# In[11]:

t0 = 3.10458*10e-6
time = output_stat[:,1]*t0


# In[12]:

time


# In[22]:

msd = np.sum(output_stat[:,2:5],1)


# In[24]:

import matplotlib.pyplot as plt


# In[28]:

plt.plot(time,msd)
plt.xlabel("time")
plt.ylabel("msd")
plt.show()


# In[29]:

plt.savefig('msd.png')


# In[ ]:



