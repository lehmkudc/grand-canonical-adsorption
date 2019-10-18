# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 19:56:32 2019

@author: Dustin Lehmkuhl
"""

import numpy as np



N_part = 3
Nm_max = 1000

M = np.zeros(( N_part, 3, Nm_max ))
Nm = [0]*N_part


def move(M, Part, o, x, y, z):
    # Carry out the "move" trial move
    M[Part, 0, o] = x
    M[Part, 1, o] = y
    M[Part, 2, o] = z
    
    return( M )
    
def add(M, Nm, Part, x, y, z):
    # Carry out the "add" trial move
    M[Part, 0, Nm[Part]] = x
    M[Part, 1, Nm[Part]] = y
    M[Part, 2, Nm[Part]] = z
    Nm[Part] = Nm[Part] + 1
    
    return( M, Nm )

def remove(M, Nm, Part, o):
    # Carry out the "remove" trial move
    Nm[Part] = Nm[Part] - 1
    M[Part, 0, o] = M[Part, 0, Nm[Part]]
    M[Part, 1, o] = M[Part, 1, Nm[Part]]
    M[Part, 2, o] = M[Part, 2, Nm[Part]]
    
    return( M, Nm )
    

def load_M( Nm ):
    

print( Nm[0] )
print( M.shape )

M, Nm = remove( M, Nm,1, 1)
M
Nm
