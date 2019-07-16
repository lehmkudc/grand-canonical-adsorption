# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 18:01:42 2019

@author: Dustin
"""
from random import random, seed
from math import floor, pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os



N_points = 1000
xv = [0]*N_points
yv = [0]*N_points
cir = [0]*N_points

for i in range(N_points):
    
    rx = 2*random()-1
    ry = 2*random()-1
    xv[i] = rx
    yv[i] = ry
    
    if (rx*rx + ry*ry < 1):
        cir[i] = 1
    
plt.scatter(xv,yv,c=cir )