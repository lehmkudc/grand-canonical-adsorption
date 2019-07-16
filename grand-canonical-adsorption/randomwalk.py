# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 17:45:14 2019

@author: Dustin
"""

from random import random, seed
from math import floor, pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os


N_steps = 1000
histx = [0]*N_steps
histy = [0]*N_steps

x = 0
y = 0
for j in range( N_steps ):
    R = random()
    
    if R < 0.25:
        x = x + 1
    elif R < 0.5:
        x = x - 1
    elif R < 0.75:
        y = y + 1
    else:
        y = y - 1
    histx[j] = x
    histy[j] = y

    
plt.plot( histx, histy, 'k-' )