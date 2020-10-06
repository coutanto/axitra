import numpy as np
import axitraf
import matplotlib.pyplot as plt

# 2 sources
sources=np.array([[1, 45.100, 2.000, 5000.000],
                  [2, 45.200, 2.000, 5000.000]])

# 5 receivers
stations=np.array(
        [[1, 45.000, 2.000, 0.000],
         [2, 46.000, 1.000, 0.000],
         [3, 46.000, 3.000, 0.000],
         [4, 44.000, 1.000, 0.000],
         [5, 44.000, 3.000, 0.000]])

# 2 layers
model = np.array([[1000., 5000., 2886., 2700., 1000., 500.],
                  [0., 6000., 3886., 2700., 1000., 500.]])

ap=axitraf.force_green(model,stations,sources,40.,20.,latlon=True)

# history of source
hist = np.array([[1,1.,0.,1.,10.,10.0],
                 [2,0.,1.,0.,10.,10.0]])

t, sx, sy, sz = axitraf.force_conv(ap,hist,source_type=1,t0=0.5,unit=1)

plt.plot(t,sx[1,:],t,sx[2,:],t,sx[3,:],t,sx[4,:],t,sx[0,:])
plt.pause(0.001)

input("Press [enter] to continue.")

sfunc=np.zeros((ap.npt,),dtype='float64')
sfunc[10]=1
t, s2x, s2y, s2z = axitraf.force_conv(ap,hist,source_type=3,t0=3,unit=1,sfunc=sfunc)

plt.plot(t,s2x[1,:],t,s2x[2,:],t,s2x[3,:],t,s2x[4,:],t,s2x[0,:])
plt.pause(0.001)

input("Press [enter] to continue.")

ap.clean()
