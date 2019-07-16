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
N_runs = 10000
histx = [0]*N_runs
histy = [0]*N_runs

for i in range( N_runs ):
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
            
    histx[i] = x
    histy[i] = y
    
plt.scatter( histx, histy, alpha=0.2,s=2 )