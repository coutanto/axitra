#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import numpy as np
from axitra import *
import matplotlib.pyplot as pt


# In[2]:


# 2 sources
# index, lat, lon, depth
sources=np.array([[1, 45.100, 2.000, 5000.000],
                  [2, 45.200, 2.000, 5000.000]])


# In[3]:


# 5 receivers with geographical coordinates
# index, lat, lon, depth
stations=np.array(
        [[1, 45.000, 2.000, 0.000],
         [2, 46.000, 1.000, 0.000],
         [3, 46.000, 3.000, 0.000],
         [4, 44.000, 1.000, 0.000],
         [5, 44.000, 3.000, 0.000]])


# In[4]:


# 2 layers
# thickness (or top), Vp, Vs, rho, Qp, Qs
model = np.array([[1000., 5000., 2886., 2700., 1000., 500.],
                  [0., 6000., 3886., 2700., 1000., 500.]])


# In[5]:


# Fill in the instance of Axitra Class
ap = Axitra(model, stations, sources, fmax=20., duration=50., xl=0., latlon=True)


# In[6]:


# Compute green's function
# fmax = 20Hz
# duration = 50 sec
ap = moment.green(ap)


# In[7]:


# history of source
hist = np.array([[1,7.5e20,148.0,84.0,-47.0,0.,0.,20.0],
                 [2,7.5e20,148.0,84.0,-47.0,0.,0.,20.0]])


# In[8]:


# first convolution example
# source= ricker
# source time width = 3 sec
# output unit = displacement
t, sx, sy, sz = moment.conv(ap,hist,source_type=1,t0=2,unit=1)


# In[9]:


pt.figure(figsize=(18, 9))
ier=pt.plot(t,sx[1,:],t,sx[2,:],t,sx[3,:],t,sx[4,:],t,sx[0,:],)


# In[10]:


# convolution with a user provided source function
# A dirac of unit 1 at 10th sample
sfunc = np.zeros((ap.npt,),dtype='float64')
sfunc[10]=1.
t, sx, sy, sz = moment.conv(ap,hist,source_type=3,t0=3,unit=1,sfunc=sfunc)


# In[11]:


pt.figure(figsize=(18, 9))
ier=pt.plot(t,sx[1,:],t,sx[2,:],t,sx[3,:],t,sx[4,:],t,sx[0,:],)


# In[12]:


# Run a new instance by reading existing axitra input files
ap2=Axitra.read(str(ap.id))
# clean files associtated with previous instance
ap.clean()


# In[13]:


# Run the second instance that should give identical results
moment.green(ap2);
t, sx, sy, sz = moment.conv(ap2,hist,source_type=1,t0=2,unit=1)
pt.figure(figsize=(18, 9))
ier=pt.plot(t,sx[1,:],t,sx[2,:],t,sx[3,:],t,sx[4,:],t,sx[0,:],)


# In[14]:


help(moment.green)


# In[15]:


help(moment.conv)


# In[ ]:




