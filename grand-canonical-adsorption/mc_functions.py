from random import random, seed
from math import floor
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os

def load_h(Nh, s_box, Nh_max):
    X = [0]*Nh_max; Y=[0]*Nh_max; Z=[0]*Nh_max;
    for i in range(Nh ):
        X[i] = random()*s_box
        Y[i] = random()*s_box
        Z[i] = random()*s_box
    return X,Y,Z

def load_c(s_box,c_bond):
    Nci = floor( s_box/c_bond )
    Nc = 3*Nci
    Xc = [0]*Nc
    Yc = [0]*Nc
    Zc = [0]*Nc
    for i in range(Nci):
        Xc[i] = c_bond*i + 0.02
        Yc[i] = s_box/2
        Zc[i] = s_box/2
    for j in range(Nci):
        Xc[Nci + j] = s_box/2
        Yc[Nci + j] = c_bond*j + 0.02
        Zc[Nci + j] = s_box/2
    for k in range(Nci):
        Xc[2*Nci + k] = s_box/2
        Yc[2*Nci + k] = s_box/2
        Zc[2*Nci + k] = c_bond*k + 0.02
    return Xc, Yc, Zc, Nc

def move(o, x, y, z):
    X[o] = x
    Y[o] = y
    Z[o] = z
def add(x, y, z):
    X[Nh] = x
    Y[Nh] = y
    Z[Nh] = z
def remove(o):
    X[o] = X[Nh-1]
    Y[o] = Y[Nh-1]
    Z[o] = Z[Nh-1]
def get_h(N):
    return X[:N], Y[:N], Z[:N]


def dist_hi(x,y,z,j):
    dx = x - X[j]
    dy = y - Y[j]
    dz = z - Z[j]
    if (dx > 0.5*s_box):
        dx = dx-s_box
    elif (dx < -0.5*s_box):
        dx = dx + s_box
    if (dy > 0.5*s_box):
        dy = dy-s_box
    elif (dy < -0.5*s_box):
        dy = dy + s_box
    if (dz > 0.5*s_box):
        dz = dz-s_box
    elif (dz < -0.5*s_box):
        dz = dz + s_box
    return dx*dx + dy*dy + dz*dz

def dist_ci(x,y,z,j):
    dx = x - Xc[j]
    dy = y - Yc[j]
    dz = z - Zc[j]
    if (dx > 0.5*s_box):
        dx = dx-s_box
    elif (dx < -0.5*s_box):
        dx = dx + s_box
    if (dy > 0.5*s_box):
        dy = dy-s_box
    elif (dy < -0.5*s_box):
        dy = dy + s_box
    if (dz > 0.5*s_box):
        dz = dz-s_box
    elif (dz < -0.5*s_box):
        dz = dz + s_box
    return dx*dx + dy*dy + dz*dz


def dist_h(x,y,z):
    r2h = [0]*Nh
    for i in range(Nh):
        r2h[i] = dist_hi( x,y,z,i )
    return r2h

def dist_c(x,y,z):
    r2c = [0]*Nc
    for i in range(Nc):
        r2c[i] = dist_ci(x,y,z,i)
    return r2c


def Ui(r2, eps, sig):
    if (r2 <= rc*rc):
        r2i = sig*sig/r2
        r6i = r2i*r2i*r2i
        En = 4*eps*(r6i*r6i-r6i)
        Vir = 48*eps*(r6i*r6i-0.5*r6i)
    else:
        En = 0
        Vir = 0
    return En, Vir


def Up( x, y, z, ia):
    En_move = 0
    Vir_move = 0
    for i in range(Nh):
        if ( i != ia):
            r2 = dist_hi( x,y,z, i)
            ui, viri = Ui( r2, e_hh, s_hh)
            En_move = En_move + ui
            Vir_move = Vir_move + viri
    for ic in range(Nc):
        r2 = dist_ci( x,y,z, ic)
        ui, viri = Ui( r2, e_hc, s_hc)
        En_move = En_move + ui
        Vir_move = Vir_move + viri
    return En_move, Vir_move


def mc_add():
    # Select Random location
    x = random()*s_box
    y = random()*s_box
    z = random()*s_box
    
    #Calculate Energy of Trial Move
    U_move, Vir_move = Up(x,y,z,Nh+1)
    if tailcor:
        U_move = U_move + (Nh+1)*Ucor(rc,(Nh+1)/Vol) - Nh*Ucor(rc, Nh/Vol)
        
    # Probability of accepting trial move
    pa_add = ZZ*Vol*exp( -beta*U_move )/(Nh + 1)
    
    # Accept or Decline the Trial move
    if (random() < pa_add):
        add(x,y,z)
        UT = UT + U_move
        VirT = VirT + Vir_move
        Nh = Nh + 1
        Nacc = Nacc + 1
        if (Nh > Nh_max):
            Panic()
            
            
def mc_remove():
    global Nh
    if (Nh == 0):
        return
    
    # Select Random particle
    o = floor( random()*20 )
    x = X[o]
    y = Y[o]
    z = Z[o]
    
    # Calculate Energy of Trial Move
    U_move, Vir_move = Up( x,y,z, o)
    U_move = -U_move
    Vir_move = -Vir_move
    if tailcor:
        U_move = U_move + ((Nh-1)*Ucor(rc,(Nh-1)/Vol)-Nh*Ucor(rc,Nh/Vol))
    
    # Probability of Accepting Trial Move
    pa_remove = Nh*np.exp( -beta*U_move )/(ZZ*Vol)
    
    # Accept or Decline Trial Move
    if (random() < pa_remove):
        remove(o)
        UT = UT + U_move
        VirT = VirT + Vir_move
        Nh = Nh - 1
        Nacc = Nacc + 1
        if (Nh > Nh_max):
            Panic()
            
            
            
def mc_move():
    if (Nh == 0):
        return
    
    
    
def look():
    plt.scatter( Xc, Yc, c="k")
    plt.scatter( X, Y, c="b")
    
    