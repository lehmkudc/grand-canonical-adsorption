from random import random, seed, randrange
from math import floor, pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
from time import time
import pandas as pd



def load_co( Nco, s_box, N_max ):
    
    Xco = [0]*N_max; Yco =[0]*N_max; Zco =[0]*N_max;
    for i in range(Nco ):
        Xco[i] = random()*s_box #[A]
        Yco[i] = random()*s_box #[A]
        Zco[i] = random()*s_box #[A]
    return Xco, Yco, Zco   

def load_me( Nme, s_box, N_max ):
    
    Xme = [0]*N_max; Yme =[0]*N_max; Zme =[0]*N_max;
    for i in range( Nme ):
        Xme[i] = random()*s_box #[A]
        Yme[i] = random()*s_box #[A]
        Zme[i] = random()*s_box #[A]
    return Xme, Yme, Zme   



def move(spec, o, x, y, z):
    # Carry out the "move" trial move
    if( spec == "co2"): 
        Xco[o] = x
        Yco[o] = y
        Zco[o] = z
        
    elif (spec == "me"):
        Xme[o] = x
        Yme[o] = y
        Zme[o] = z
    
    
def add(spec, x, y, z):
    # Carry out the "add" trial move
    global Nco, Nme
    if( spec == "co2"): 
        Xco[Nco] = x
        Yco[Nco] = y
        Zco[Nco] = z
        Nco = Nco + 1
        
    elif (spec == "me"):
        Xme[Nme] = x
        Yme[Nme] = y
        Zme[Nme] = z
        Nme = Nme + 1
    
    
def remove(o):
    # Carry out the "remove" trial move
    global Nco, Nme
    if( spec == "co2"): 
        Nco = Nco - 1
        Xco[o] = Xco[Nco]
        Yco[o] = Yco[Nco]
        Zco[o] = Zco[Nco]
        
    elif (spec == "me"):
        Nme = Nme - 1
        Xme[o] = Xme[Nme]
        Yme[o] = Yme[Nme]
        Zme[o] = Zme[Nme]
        
    


def dist_hi(spec,x,y,z,j):
    # Distance btw a proposed particle and the ith H2 particle
    if(spec == "co2"):
        dx = x - Xco[j] #[A]
        dy = y - Yco[j] #[A]
        dz = z - Zco[j] #[A]
    elif (spec == "me"):
        dx = x - Xme[j] #[A]
        dy = y - Yme[j] #[A]
        dz = z - Zme[j] #[A]

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
    # Distance btw proposed particle and the ith C particle
    dx = x - Xc[j] #[A]
    dy = y - Yc[j] #[A]
    dz = z - Zc[j] #[A]
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



def Ui(r2, eps, sig):
    # LJ potential btw two particles at r distance away
    if ( r2 < rc*rc ):
        r2i = sig*sig/r2 #[]
        r6i = r2i*r2i*r2i
        En = 4*eps*(r6i*r6i-r6i) #[K]
        Vir = 48*eps*(r6i*r6i-0.5*r6i) #[K]
    else:
        En = 0
        Vir = 0
    return En, Vir

def Up(spec, x, y, z, ia, jb=0):
    # Total LJ potential of proposed particle with all other particles
    # omit the ia'th H2 particle. When not needed, Nh is used.
    # jb used in the UTo operation to avoid overcounting interactions
    En_move = 0
    Vir_move = 0
    
    if(spec == "co2"):
        #Similar Particles
        for i in range(jb, Nco): # H2 particles
            if ( i != ia):
                r2 = dist_hi("co2", x,y,z, i)
                ui, viri = Ui( r2, e_co, s_co)
                En_move = En_move + ui
                Vir_move = Vir_move + viri
        #Other Mobile
        
        for i2 in range(0,Nme):
            r2 = dist_hi("me", x,y,z, i)
            ui, viri = Ui( r2, e_meco, s_meco)
            En_move = En_move + ui
            Vir_move = Vir_move + viri
        
        for ic in range(Nc): # C particles
            r2 = dist_ci( x,y,z, ic)
            ui, viri = Ui( r2, e_cco, s_cco)
            En_move = En_move + ui
            Vir_move = Vir_move + viri
        
    elif (spec == "me"):
        #Similar Particles
        for i in range(jb, Nme): # H2 particles
            if ( i != ia):
                r2 = dist_hi("me", x,y,z, i)
                ui, viri = Ui( r2, e_me, s_me)
                En_move = En_move + ui
                Vir_move = Vir_move + viri
        #Other Mobile
        
        for i2 in range(0,Nco):
            r2 = dist_hi("co2", x,y,z, i)
            ui, viri = Ui( r2, e_meco, s_meco)
            En_move = En_move + ui
            Vir_move = Vir_move + viri
        
        for ic in range(Nc): # C particles
            r2 = dist_ci( x,y,z, ic)
            ui, viri = Ui( r2, e_cme, s_cme)
            En_move = En_move + ui
            Vir_move = Vir_move + viri
                 
    return En_move, Vir_move



e_me = 147.9 # eps over kb[K]
s_me = 3.73 # sigma [A]

e_co = 309
s_co = 3.36

e_c = 28.0 #Activated is 89.44
s_c = 3.4

e_meco = np.sqrt(e_me*e_co)
s_meco = (s_me+s_co)/2

e_cme = np.sqrt(e_me*e_c)
s_cme = (s_me+s_c)/2

e_cco = np.sqrt(e_co*e_c)
s_cco = (s_co + s_c)/2

y_co = 0.5
P_res = 1 #Pa
T = 77 #K





# User Defined Variables
    
N_max = 1000
pi_move = 0.5

ab = seed()

filename = "sbox_updated2.csv"
input_data = pd.read_csv(filename)
N_exp = input_data.shape[0]

rho_red_result = [0]*N_exp
P_red_result = [0]*N_exp
P_var_red_result = [0]*N_exp
rho_result = [0]*N_exp
P_result = [0]*N_exp
P_var_result = [0]*N_exp
time_result = [0]*N_exp
random_seed = [0]*N_exp
print( N_exp, "Runs")

for i in range(N_exp):
    t0 = time()
    
    Pid_red = input_data.Pid_red[i]
    T_red = input_data.T_red[i]
    tailcor = input_data.tailcor[i]
    N_moves = input_data.nmoves[i]
    N_equil = input_data.nequil[i]
    N_prod = input_data.nprod[i]
    s_box = input_data.sbox[i]

    # Define Simulation Properties
    delta = 1
    
    # Define Useful Constants
    s_hh = 3.73 # sigma [A]
    e_hh = 147.9 # eps over kb[K]
    s_me = 3.73 # [A]
    e_me = 147.5 #[K]
    s_hc = 2.74 # sigma [A]
    e_hc = 16.2463 # eps over kb[K]
    kb = 1.3806*10**(7) #[Pa*A^3/K]
    c_bond = 2.24 #[A]
        
    # Calculated Properties
    Pid = e_hh*kb*Pid_red/s_hh**3 # [Pa]
    T = e_hh*T_red # [K]
    Vol = s_box**3 #[A^3]
    beta = 1/T #[K^-1]
    ZZ = beta*Pid
    Nh = floor( ZZ*Vol/kb )
    rc = s_box
    if ( tailcor):
        rc = min([2.5*s_hh,0.5*s_box]) #[A]
        
    # Run Simulation
    random_seed[i] = randrange(10000,99999)
    seed( random_seed[i] )
    rhov, Env, Pv, Nv = mc_run() #[A-3, , MPa, ]
    rho_result[i] = rhov.mean() #[A-3]
    P_result[i] = Pv.mean() #[MPa]
    P_var_result[i] = Pv.std()
    rho_red_result[i] = s_hh**3*rhov.mean()
    P_red_result[i] = s_hh**3/e_hh/kb*Pv.mean()*10**(6)
    P_var_red_result[i] = s_hh**3/e_hh/kb*Pv.std()*10**(6)
    time_result[i] = time() - t0
    print("i =", i,
          "\trho =", round( rho_red_result[i], 5), 
          "\tP =", round(P_red_result[i],5), 
          "\tTime =", round(time_result[i], 5) )
    
    output_data = input_data    
    output_data['rho'] = rho_result
    output_data['P'] = P_result
    output_data['P_var'] = P_var_result
    output_data['rho_red'] = rho_red_result
    output_data['P_red'] = P_red_result
    output_data['P_var_red'] = P_var_red_result
    output_data['time'] = time_result
    output_data['seed'] = random_seed
    output_data.to_csv(filename, index = False)
