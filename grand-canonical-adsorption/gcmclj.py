

from random import random, seed
from math import floor, pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os



def load_h(Nh, s_box, Nh_max):
    # Initializes the Unit Cell with H2 particles
    X = [0]*Nh_max; Y=[0]*Nh_max; Z=[0]*Nh_max;
    for i in range(Nh ):
        X[i] = random()*s_box #[A]
        Y[i] = random()*s_box #[A]
        Z[i] = random()*s_box #[A]
    return X,Y,Z   

def load_c(s_box,c_bond):
    # Initializes the Unit Cell with C particles
    Nci = floor( s_box/c_bond )
    Nc = 3*Nci
    Xc = [0]*Nc
    Yc = [0]*Nc
    Zc = [0]*Nc
    for i in range(Nci):
        Xc[i] = c_bond*i + 0.02 #[A]
        Yc[i] = s_box/2 #[A]
        Zc[i] = s_box/2 #[A]
    for j in range(Nci):
        Xc[Nci + j] = s_box/2 #[A]
        Yc[Nci + j] = c_bond*j + 0.02 #[A]
        Zc[Nci + j] = s_box/2 #[A]
    for k in range(Nci):
        Xc[2*Nci + k] = s_box/2 #[A]
        Yc[2*Nci + k] = s_box/2 #[A]
        Zc[2*Nci + k] = c_bond*k + 0.02 #[A]
    Nc = 0
    Xc=[]; Yc=[]; Zc=[]
    return Xc, Yc, Zc, Nc

def move(o, x, y, z):
    # Carry out the "move" trial move
    X[o] = x
    Y[o] = y
    Z[o] = z
def add(x, y, z):
    # Carry out the "add" trial move
    global Nh
    X[Nh] = x
    Y[Nh] = y
    Z[Nh] = z
    Nh = Nh + 1
def remove(o):
    # Carry out the "remove" trial move
    global Nh
    Nh = Nh - 1
    X[o] = X[Nh]
    Y[o] = Y[Nh]
    Z[o] = Z[Nh]
    
def get_h():
    return X[:Nh], Y[:Nh], Z[:Nh]

def dist_hi(x,y,z,j):
    # Distance btw a proposed particle and the ith H2 particle
    dx = x - X[j] #[A]
    dy = y - Y[j] #[A]
    dz = z - Z[j] #[A]
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

def dist_h(x,y,z):
    # Distance btw proposed particle and all H2 particles
    # Currently not used
    r2h = [0]*Nh 
    for i in range(Nh):
        r2h[i] = dist_hi( x,y,z,i ) #[A^2]
    return r2h

def dist_c(x,y,z):
    # Distance btw proposed particle and all C particles
    # Currently not used
    r2c = [0]*Nc
    for i in range(Nc):
        r2c[i] = dist_ci(x,y,z,i) #[A^2]
    return r2c


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

def Up( x, y, z, ia, jb=0):
    # Total LJ potential of proposed particle with all other particles
    # omit the ia'th H2 particle. When not needed, Nh is used.
    # jb used in the UTo operation to avoid overcounting interactions
    En_move = 0
    Vir_move = 0
    for i in range(jb, Nh): # H2 particles
        if ( i != ia):
            r2 = dist_hi( x,y,z, i)
            ui, viri = Ui( r2, e_hh, s_hh)
            En_move = En_move + ui
            Vir_move = Vir_move + viri
    for ic in range(Nc): # C particles
        r2 = dist_ci( x,y,z, ic)
        ui, viri = Ui( r2, e_hc, s_hc)
        En_move = En_move + ui
        Vir_move = Vir_move + viri
    return En_move, Vir_move

def p_rem( U_move ):
    # Acceptance probability of the "remove" trial move 
    #     given change in LJ potential
    if (U_move*beta > 100):
        return 0
    elif (U_move*beta < -100):
        return 1
    else:
        return Nh*kb*np.exp( -beta*U_move )/(ZZ*Vol)
    
def p_add( U_move ):
    # acceptance probability of the "add" trial move 
    #     given change in LJ potential
    if (U_move*beta > 100):
        return 0
    elif (U_move*beta > 100):
        return 1
    else:
        return ZZ*Vol*np.exp( -beta*U_move )/(Nh + 1)/kb
    
def p_move( U_move ):
    # acceptance probability of the "move" trial move 
    #     given change in LJ potential
    if (U_move*beta > 100):
        return 0
    elif (U_move*beta < -100):
        return 1
    else:
        return np.exp( -beta*U_move)

def mc_add():
    # Attempt an "add" trial move
    global UT, VirT, Nh, Aacc, Aatt
    
    Aatt = Aatt + 1
    # Select Random location
    x = random()*s_box
    y = random()*s_box
    z = random()*s_box
    
    #Calculate Energy of Trial Move
    U_move, Vir_move = Up(x,y,z,Nh+1)
    if tailcor:
        U_move = U_move + (Nh+1)*Ucor(rc,(Nh+1)/Vol) - Nh*Ucor(rc, Nh/Vol)
        
    # Probability of accepting trial move
    pa_add = p_add( U_move )
    
    # Accept or Decline the Trial move
    if (random() < pa_add):
        add(x,y,z)
        UT = UT + U_move
        VirT = VirT + Vir_move
        Aacc = Aacc + 1
        if (Nh > Nh_max):
            Panic()
            
def mc_remove():
    # Attempt a "remove" trial move
    global UT, VirT, Nh, Ratt, Racc
    Ratt = Ratt + 1
    if (Nh == 0):
        return
    
    # Select Random particle
    o = floor( random()*Nh )
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
    pa_remove = p_rem( U_move )
    
    # Accept or Decline Trial Move
    if (random() < pa_remove):
        remove(o)
        UT = UT + U_move
        VirT = VirT + Vir_move
        Racc = Racc + 1
        if (Nh > Nh_max):
            Panic()
            
            
            
def mc_move():
    # Attempt a "move" trial move
    global UT, VirT, Nh, Nacc, Natt
    Natt = Natt + 1
    if (Nh == 0):
        return
    
    # Select Random particle
    o = floor( random()*Nh )
    x = X[o]
    y = Y[o]
    z = Z[o]
    
    # Calculate Energy of current configuration
    U1, V1 = Up( x,y,z, o)
    
    # Calculate new Location
    xn = x + delta*(random()-0.5)
    yn = y + delta*(random()-0.5)
    zn = z + delta*(random()-0.5)
    xn, yn, zn = box_fix( xn, yn, zn )
    
    # Calculate Energy of New Configuration
    U2, V2 = Up( xn,yn,zn, o)
    U_move = U2 - U1
    V_move = V2 - V1
    
    pa_move = p_move( U_move )
        
    # Accept the trial move
    if ( random() < pa_move ):
        Nacc = Nacc + 1
        move(o, xn, yn, zn)
        UT = UT + U_move
        VirT = VirT + V_move
    
def look():
    # Plot a bird's eye view of Unit Cell
    plt.scatter( Xc, Yc, c="k")
    Xh, Yh, Zh = get_h()
    plt.scatter( Xh, Yh, c="b")
        
def box_fix( x, y, z):
    # Correct for unit cell periodicity
    # Currently not used
    if x < 0:
        x = x + s_box
    if x > s_box:
        x = x - s_box
    if y < 0:
        y = y + s_box
    if y > s_box:
        y = y - s_box
    if z < 0:
        z = z + s_box
    if z > s_box:
        z = z - s_box
    return x,y,z

def sample():
    # Extract data after every cycle
    if (Nh == 0):
        rho = 0
        Enp = 0
        Pressure = 0
    else:
        rho = Nh/Vol #[A^-3]
        Enp = UT/Nh #[K]
        Pressure = ( Nh/Vol/beta + VirT/(3*Vol) )*kb*10**(-6) #[MPa] 
        if tailcor:
            Pressure = Pressure + Pcor(rc, Nh/Vol)
    return rho, Enp, Pressure, Nh

def UTo():
    # Total potential energy of the system (C-C interactions omitted)
    Ut = 0
    Virt = 0
    for i in range(Nh):
        xn = X[i]
        yn = Y[i]
        zn = Z[i]
        
        U_p, Vir_p = Up( xn,yn,zn,Nh, i+1)
        Ut = Ut + U_p
        Virt = Virt + Vir_p
    
    if tailcor:
        Ut = Ut + Nh*Ucor( rc, Nh/Vol)
        
    return Ut, Virt

def mc_step():
    # Perform one step of simulation
    if ( random() < pi_move ):
        mc_move()
    else:
        if ( random()< 0.5):
            mc_add()
        else:
            mc_remove()
            
def mc_cycle():
    # perform one cycle of simulation
    global Natt, Nacc, Aatt, Aacc, Ratt, Racc, delta
    rhow=np.zeros(N_moves); Enw=np.zeros(N_moves); Pw=np.zeros(N_moves);
    Nw=np.zeros(N_moves)
    for i in range(N_moves):
        mc_step()
        rhow[i], Enw[i], Pw[i], Nw[i] = sample()
    adjust()
    return rhow.mean(), Enw.mean(), Pw.mean(), Nw.mean()

def mc_run():
    # Perform simulation
    global X,Y,Z,Xc,Yc,Zc,Nc, UT, VirT, Natt, Nacc, Aatt, Aacc, Ratt, Racc
    
    # Initialize Unit Cell
    X,Y,Z = load_h( Nh, s_box, Nh_max)
    Xc, Yc, Zc, Nc = load_c(s_box,c_bond)
    UT, VirT = UTo()
    
    
    # Equilibration Step
    print( "Equlibration")
    Natt = 0; Nacc=0; Aatt=0; Aacc=0; Ratt=0;Racc=0;Pacc = 0; Patt = 0
    for j in range(N_equil):
        mc_cycle()
        if( (j)%floor(N_equil/10) == 0 ):
            print( str(100*j/N_equil) + "% Completed"  )
            print( "\tDelta: " + str( delta ))
            print( "\tMove acceptance: "  + str(Nacc) + " out of " + str(Natt) + " attempts." )
            print( "\tAdd acceptance:  " + str(Aacc) + " out of " + str(Aatt) + " attempts." )
            print( "\tRem acceptance:  " + str(Racc) + " out of " + str(Ratt) + " attempts." )
            Natt = 0; Nacc=0; Aatt=0; Aacc=0; Ratt=0;Racc=0 
            

    #Production Step
    rhov=np.zeros(N_prod); Env=np.zeros(N_prod); Pv=np.zeros(N_prod);
    Nv=np.zeros(N_prod); 
    print( "Production")
    for j in range(N_prod):
        rhov[j], Env[j], Pv[j], Nv[j] = mc_cycle()
        if( (j)%floor(N_prod/10) == 0 ):
            print( str(100*j/N_prod) + "% Completed"  )
    #plt.plot( Pv )
    return (rhov,Env,Pv,Nv)

def adjust():
    # Change average "move" distance by acceptance history
    # Currently not being used
    global delta, Natt, Nacc
    if (Natt == 0 or Nacc > Natt):
        return
    else:
        frac = Nacc/Natt
        dro = delta
        delta = delta*abs(frac/0.5)
        if (delta/dro > 1.5 ):
            delta = dro*1.5
        if (delta/dro < 0.5):
            delta = dro*0.5
        if (delta > s_box*0.25):
            delta = s_box*0.25
    return

def Ucor(r, rho):
    # Tail correction for LJ Potential
    sr = s_hh/r
    return (8/3)*pi*e_hh*rho*s_hh**3*(1/3*sr**9 - sr**3)

def Pcor(r, rho):
    # Tail correction for Pressure
    sr = s_hh/r
    return (16/3)*pi*e_hh*rho**2*s_hh**3*( 2/3*sr**9 - sr**3)*kb*10**(-6) #[MPa] 

def move_test():
    o = floor( random()*Nh )
    x = X[o]
    y = Y[o]
    z = Z[o]
    print("x: ", x, " y: ", y," z: ", z," o: ", o )

    # Calculate Energy of current configuration
    U1, V1 = Up( x,y,z, o)
    print("U1: ", U1, "V1: ", V1 )

    # Calculate new Location
    xn = x + delta*(random()-0.5)
    yn = y + delta*(random()-0.5)
    zn = z + delta*(random()-0.5)
    print("xn: ", xn," yn: ", yn,"zn: ", zn)
    xn, yn, zn = box_fix( xn, yn, zn )
    print("xn: ", xn," yn: ", yn,"zn: ", zn)

    # Calculate Energy of New Configuration
    U2, V2 = Up( xn,yn,zn, o)
    print("U2: ", U2, " V2: ", V2 )
    U_move = U2 - U1
    V_move = V2 - V1
    print("Umove: ", U_move, " Vmove: ", V_move )

    if (U_move*beta > 100):
        pa_move = 0
    elif (U_move*beta < -100):
        pa_move = 1
    else:
        pa_move = p_move( U_move )
    print("pa_move: ", pa_move )
    
def add_test(Pid_red, T_red, Nh,x=-1.0, y=-1.0, z=-1.0):
    if ( x<0.0 and y<0.0 and z<0.0):
        x = random()*s_box
        y = random()*s_box
        z = random()*s_box

    print( "x: ", x, " y: ", y, " z: ", z)
    #Calculate Energy of Trial Move
    U_move, Vir_move = Up(x,y,z,Nh+1)
    print("Umove: ", U_move, " Vmove: ", Vir_move )
    print("Pcorr: ", ZZ*Vol/(Nh + 1)/kb )

    # Probability of accepting trial move
    if (U_move*beta > 100):
        pa_add = 0
    elif (U_move*beta < -100):
        pa_add = 1
    else:
        pa_add = p_add( U_move )
    print("pa_add: ", pa_add )
    return( pa_add, U_move )
    
def remove_test(Pid_red, T_red, Nh, o=-1):
    if (o==-1):
        o = floor( random()*Nh )
        
    x = X[o]
    y = Y[o]
    z = Z[o]
    print("x: ", x, " y: ", y," z: ", z," o: ", o )

    # Calculate Energy of Trial Move
    U_move, Vir_move = Up( x,y,z, o)
    U_move = -U_move
    Vir_move = -Vir_move
    
    Pid = e_hh*kb*Pid_red/s_hh**3 # [Pa]
    T = e_hh*T_red # [K]
    Vol = s_box**3 #[A^3]
    beta = 1/T #[K^-1]
    ZZ = beta*Pid

    print("Umove: ", U_move, " Vmove: ", Vir_move )
    print("Pcorr: ", Nh*kb/(ZZ*Vol))

    # Probability of Accepting Trial Move
    if (U_move*beta > 100):
        pa_remove = 0
    elif (U_move*beta < -100):
        pa_remove = 1
    else:
        pa_remove = p_rem(U_move)

    print("pa_remove: ", pa_remove)
    return( pa_remove, U_move )
    
    
def mc_grand( Pid_red, T_red ):
    
    # Computed Properties
    Pid = e_hh*kb*Pid_red/s_hh**3 # [Pa]
    T = e_hh*T_red # [K]
    Vol = s_box**3 #[A^3]
    beta = 1/T #[K^-1]
    ZZ = beta*Pid
    Nh = floor( ZZ*Vol/kb )
    print( Nh )
    rc = s_box
    if ( tailcor):
        rc = min([2.5*s_hh,0.5*s_box]) #[A]
    rhov, Env, Pv, Nv = mc_run(Vol, beta, ZZ, Nh, rc)
    
    rho_red = s_hh**3*rhov.mean()
    P_red = s_hh**3/e_hh/kb*Pv.mean()*10**(6)
    print("rho*: ", rho_red, " P*: ", P_red)
    return rho_red, P_red


# User Defined Variables
Pid_red = 2
T_red = 2

# Define Simulation Properties
Nh_max = 1000
random_seed = 7
tailcor = False
N_moves = 100
N_equil = 500
N_prod = 1000
pi_move = 0.5
s_box = 5.64*3.73 #[A]
delta = 1

# Define Useful Constants
s_hh = 3.73 # sigma [A]
e_hh = 147.5 # eps over kb[K]
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
print( Nh )
rc = s_box
if ( tailcor):
    rc = min([2.5*s_hh,0.5*s_box]) #[A]
    
# Run Simulation
seed(random_seed)
rhov, Env, Pv, Nv = mc_run() #[A-3, , MPa, ]

rho_red = s_hh**3*rhov.mean()
P_red = s_hh**3/e_hh/kb*Pv.mean()*10**(6)
print("rho*: ", rho_red, " P*: ", P_red)


A_vdw = 6.832*(10**10)
B_vdw = 43.206
P_vdw = kb*T/(1/rhov.mean() - B_vdw) - A_vdw*rhov.mean()*rhov.mean()
P_vdw_red = P_vdw*s_hh**3/e_hh/kb


s_hh**3/e_hh/kb*6*10**6
P_id = rhov.mean()*kb*T #Pa
P_id_red = s_hh**3/e_hh/kb*P_id

Nv.mean()/Vol*s_hh**3


def VDW( rho_red, T_red ):
    rho = (rho_red/s_hh**3)*10**30 #[Part/m3]
    R = 8.3144598 # [m3 Pa/K mol]
    T = T_red*e_hh #[K]
    a = 0.2283 #[Pa m6/mol2]
    b = 4.278*10**(-5) #[m3/mol]
    Vm = s_hh**3/rho_red*6.02*10**(-7)
    P =  R*T/( Vm - b  ) - a*(Vm)**(-2) #[Pa]
    P_red = s_hh**3/e_hh/kb*P
    return( P_red )
    
VDW(0.5311, 2)

P = Pv.mean()*10**6 #[Pa]
rho = rhov.mean()
Z = P/rho/kb/T

P_id
Pid
Pv.mean()/rhov.mean()/kb/T

8.314*T/P_id
