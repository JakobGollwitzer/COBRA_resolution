### import modules
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import cmath
import glob
import re
import Dans_Diffraction as dif

from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #
from lmfit.models import LorentzianModel
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.optimize as opt
from contextlib import contextmanager
import matplotlib.gridspec as gridspec
from scipy.optimize import leastsq, curve_fit
from scipy.interpolate import interp1d
from matplotlib.patches import ConnectionPatch
from operator import add
from mpl_toolkits.mplot3d import Axes3D


##make a reciprocal space plot using the dans diffraction. 

from matplotlib import cm, colors

### open CIF file##

path = "/Users/jakobgo/OneDrive - Cornell University/LAO_lcls/Structural_data_LaAlO3_and_substrates/NSAT/ICSD_NSAT.cif"


xtl = dif.Crystal(path)
params = str(xtl.Cell.latt).strip().split()

#lattice constants#
a = 3.83
b = 3.83
c = 3.83
print(a, b, c)
Energy = 30.0 # in keV SET THE ENERGY

## reciprocal lattice vectors
astar = 2*np.pi/float(a)
bstar = 2*np.pi/float(b)
cstar = 2*np.pi/float(c)


#### ewald sphere surface ##########
wavestar = 2*np.pi/0.4133

u = np.linspace(-(1./8.)*np.pi, (1./8.)*np.pi, 200)
v = np.linspace(np.pi/2.7, np.pi/2., 200)
x = wavestar*np.outer(np.cos(u), np.sin(v))-wavestar
y = wavestar*np.outer(np.sin(u), np.sin(v))
z = wavestar*np.outer(np.ones(np.size(u)), np.cos(v))

#### rotate the ewald sphere 5 degrees around the y axis, corresponds to incidence angle 5 degrees
angleRad = np.deg2rad(-5.0)
Ry = [[np.cos(angleRad), 0, np.sin(angleRad)],[0,1,0],[-np.sin(angleRad), 0, np.cos(angleRad)]] #rotate around y axis 
t = np.transpose(np.array([x,y,z]), (1,2,0))
x,y,z = np.transpose(np.dot(t,Ry), (2,0,1))









# define the Bragg peak here 022
xx=0.0*astar
yy=2.0*bstar
zz=2.0*cstar + 1*(cstar/16.0)## to see resolution add 1/25th of a RLU corresponding to 1 fringe for a 1 nm thick LAO film
### print position in reciprocal space
print(xx, yy, zz)

### initialize array to find minimum. best point on ewald surface closest to bragg reflection point
minarray = np.zeros(5)
minarray[0] =1e8

### array of rotated values. lots of rotation haha
phis = np.linspace(0, 12, 1201)


### for loop of distance minimization routine
for t in range(0, len(phis), 1):
    angleRad = np.deg2rad(phis[t])
    Rz = [[np.cos(angleRad),  -np.sin(angleRad), 0], [np.sin(angleRad),  np.cos(angleRad), 0], [0,0,1]]

    tt = np.transpose(np.array([x,y,z]), (1,2,0))
    xR,yR,zR = np.transpose(np.dot(tt,Rz), (2,0,1))

    if(t % 100 == 0):
        print('.')
    
    
    for f in range(0, len(xR), 1):
        for g in range(0, len(xR[0]), 1):
            minimum = np.sqrt(np.abs(xR[f][g]-xx)**2 + np.abs(yR[f][g]-yy)**2 + np.abs(zR[f][g]-zz)**2)
            if(minarray[0]>minimum):
                minarray[0]=minimum
                minarray[1]=np.rad2deg(angleRad)
                minarray[2]=xR[f][g]
                minarray[3]=yR[f][g]
                minarray[4]=zR[f][g]

#### minarray outputs the angles you need to get to the point in line 75. Depending on how many points you want per fringe, minarray outputs the 
#### phi resolution needed
print(minarray)


plt.show()







