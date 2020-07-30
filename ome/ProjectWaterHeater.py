# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 13:31:58 2019

@author: kdst
"""

from OMPython import ModelicaSystem
from numpy import transpose, ones, linspace

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
degree_sign= u'\N{DEGREE SIGN}'

#line1 = ones(370)

mod = ModelicaSystem("ProjectWaterHeater.mo", "ProjectWaterHeater.SimSlicedWaterHeater1")
#mod.setInputs(["Ti=27", "Vd=13", "up=0.5", "uv=0.5", 'Inversed="False"'])

mod.setSimulationOptions(["startTime=0" ,"stopTime=3600", "stepSize=10"])

mod.simulate()

t = transpose(mod.getSolutions(["time"]))
N = mod.getSolutions("N")
#tarr = range(int(N[0,0]))
#for i in range(int(N[0,0])):
#    Tarr(i) = "T[" + str(i) + "]"

T = mod.getSolutions(["T"])
'''
T = transpose(T)
line1 = ones(len(T))

line = linspace(0,1.5,len(T[1]))
fig = plt.figure(figsize=(12 , 6))
ax = Axes3D(fig)
ax.axis([3600,0,0,1.5])
ax.set_title("Inversed intial temperature and heater is turned on at t=1500s",fontsize=15)
ax.set_ylabel("Heigh position in tank (m)",labelpad=15,fontsize=15)
ax.set_xlabel("Time (s)",labelpad=15,fontsize=15)
ax.set_zlabel("Temperature (" + degree_sign + "C)",labelpad=5,fontsize=15)
ax.plot_surface(t,line,T,cmap='plasma',edgecolor='none')
ax.view_init(30, -45)
plt.show()
'''


