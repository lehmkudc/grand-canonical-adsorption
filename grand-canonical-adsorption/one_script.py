from random import random, seed
from math import floor, pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
%matplotlib inline
import numpy as np
import os


def load_h(Nh, s_box, Nh_max):
    X = [0]*Nh_max; Y=[0]*Nh_max; Z=[0]*Nh_max;
    for i in range(Nh ):
        X[i] = random()*s_box #[A]
        Y[i] = random()*s_box #[A]
        Z[i] = random()*s_box #[A]
    return X,Y,Z   

def load_c(s_box,c_bond):
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
    return Xc, Yc, Zc, Nc

def move(o, x, y, z):
    X[o] = x
    Y[o] = y
    Z[o] = z
def add(x, y, z):
    global Nh
    X[Nh] = x
    Y[Nh] = y
    Z[Nh] = z
    Nh = Nh + 1
def remove(o):
    global Nh
    X[o] = X[Nh-1]
    Y[o] = Y[Nh-1]
    Z[o] = Z[Nh-1]
    Nh = Nh - 1
def get_h():
    return X[:Nh], Y[:Nh], Z[:Nh]


def dist_hi(x,y,z,j):
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
    r2h = [0]*Nh 
    for i in range(Nh):
        r2h[i] = dist_hi( x,y,z,i ) #[A^2]
    return r2h

def dist_c(x,y,z):
    r2c = [0]*Nc
    for i in range(Nc):
        r2c[i] = dist_ci(x,y,z,i) #[A^2]
    return r2c


def Ui(r2, eps, sig):
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
    En_move = 0
    Vir_move = 0
    for i in range(jb, Nh):
        if ( i != ia):
            r2 = dist_hi( x,y,z, i)
            ui, viri = Ui( r2, e_hh, s_hh)
            En_move = En_move + ui
            Vir_move = Vir_move + viri
    #print( En_move, Vir_move)
    for ic in range(Nc):
        r2 = dist_ci( x,y,z, ic)
        ui, viri = Ui( r2, e_hc, s_hc)
        En_move = En_move + ui
        Vir_move = Vir_move + viri
    return En_move, Vir_move

def p_rem( U_move ):
    return Nh*kb*np.exp( -beta*U_move )/(ZZ*Vol)
def p_add( U_move ):
    return ZZ*Vol*np.exp( -beta*U_move )/(Nh + 1)/kb
def p_move( U_move ):
    return np.exp( -beta*U_move)

def mc_add():
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
    if (U_move*beta > 100):
        pa_add = 0
    elif (U_move*beta < -100):
        pa_add = 1
    else:
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
    if (U_move*beta > 100):
        pa_remove = 0
    elif (U_move*beta < -100):
        pa_remove = 1
    else:
        pa_remove = p_rem(U_move)
    
    # Accept or Decline Trial Move
    if (random() < pa_remove):
        remove(o)
        UT = UT + U_move
        VirT = VirT + Vir_move
        Racc = Racc + 1
        if (Nh > Nh_max):
            Panic()
            
            
            
def mc_move():
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
    
    if (U_move*beta > 100):
        pa_move = 0
    elif (U_move*beta < -100):
        pa_move = 1
    else:
        pa_move = p_move( U_move )
        
        
    if ( random() < pa_move ):
        Nacc = Nacc + 1
        move(o, xn, yn, zn)
        UT = UT + U_move
        VirT = VirT + V_move
    
def look():
    plt.scatter( Xc, Yc, c="k")
    Xh, Yh, Zh = get_h()
    plt.scatter( Xh, Yh, c="b")
        
def box_fix( x, y, z):
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
    return rho, Enp, Pressure

def UTo():
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
    if ( random() < pi_move ):
        mc_move()
    else:
        if ( random()< 0.5):
            mc_add()
        else:
            mc_remove()
            
def mc_cycle():
    global Natt, Nacc, Aatt, Aacc, Ratt, Racc
    Natt = 0; Nacc=0; Aatt=0; Aacc=0; Ratt=0;Racc=0 
    for i in range(N_moves):
        mc_step()
    adjust()

def mc_run():
    global X,Y,Z,Xc,Yc,Zc,Nc, UT, VirT
    X,Y,Z = load_h( Nh, s_box, Nh_max)
    Xc, Yc, Zc, Nc = load_c(s_box,c_bond)
    UT, VirT = UTo()
    prod = False
    print( "Equlibration")
    for j in range(N_equil):
        mc_cycle()
        if( (j)%floor(N_equil/10) == 0 ):
            print( str(100*j/N_equil) + "% Completed"  )
            print( "\tMove acceptance: "  + str(Nacc/Natt) )
            print( "\tAdd acceptance:  " + str(Aacc/Aatt) )
            print( "\tRem acceptance:  " + str(Racc/Ratt) )
    prod = True
    rhov=np.zeros(N_prod); Env=np.zeros(N_prod); Pv=np.zeros(N_prod);
    Virv=np.zeros(N_prod); 
    print( "Production")
    for j in range(N_prod):
        mc_cycle()
        rhov[j], Env[j], Pv[j] = sample()
        if( (j)%floor(N_prod/10) == 0 ):
            print( str(100*j/N_prod) + "% Completed"  )
    plt.plot( Pv )
    return (rhov,Env,Pv)

def adjust():
    return
    global delta
    if (Natt == 0 or Nacc > Natt):
        naccp = Nacc
        attempp = Natt
    else:
        frac = (Nacc-naccp)/(Attemp-attempp)
        dro = delta
        delta = delta*abs(frac/0.5)
        if (delta/dro > 1.5 ):
            delta = dro*1.5
        if (delta/dro < 0.5):
            delta = dro*0.5
        if (delta > s_box*0.25):
            delta = s_box*0.25
        naccp = Nacc
        attempp = Attemp
    return

def Ucor(r, rho):
    sr = s_hh/r
    return (8/3)*pi*e_hh*rho*s_hh**3*(1/3*sr**9 - sr**3)


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
    
def add_test():
    x = random()*s_box
    y = random()*s_box
    z = random()*s_box

    print( "x: ", x, " y: ", y, " z: ", z)
    #Calculate Energy of Trial Move
    U_move, Vir_move = Up(x,y,z,Nh+1)
    print("Umove: ", U_move, " Vmove: ", Vir_move )

    # Probability of accepting trial move
    if (U_move*beta > 100):
        pa_add = 0
    elif (U_move*beta < -100):
        pa_add = 1
    else:
        pa_add = p_add( U_move )
    print("pa_add: ", pa_add )
    

def remove_test():
    o = floor( random()*Nh )
    x = X[o]
    y = Y[o]
    z = Z[o]
    print("x: ", x, " y: ", y," z: ", z," o: ", o )

    # Calculate Energy of Trial Move
    U_move, Vir_move = Up( x,y,z, o)
    U_move = -U_move
    Vir_move = -Vir_move

    print("Umove: ", U_move, " Vmove: ", V_move )

    # Probability of Accepting Trial Move
    if (U_move*beta > 100):
        pa_remove = 0
    elif (U_move*beta < -100):
        pa_remove = 1
    else:
        pa_remove = p_rem(U_move)

    print("pa_remove: ", pa_remove)
    
    

# Define Simulation Properties
Nh_max = 1000
seed = 1731
tailcor = False
N_moves = 1000 # Moves per cycle
N_equil = 500 #Cycles of equilibration step
N_prod = 500 #Cycles of poduction step
pi_move = 0.5 #Probability of initiating move step

# Define Useful Constants
s_hh = 2.958 # sigma [A]
e_hh = 36.7 # eps over kb[K]
s_me = 3.73 # [A]
e_me = 147.5 #[K]
s_hc = 3.216 # sigma [A]
e_hc = 41.924 # eps over kb[K]
kb = 1.3806*10**(7) #[Pa*A^3/K]
c_bond = 1.54 #[A]


# Define System Properties
Pid = 19.56*10**(6) #[Pa]
T = 36.7 #[K]
s_box = 10 #[A]
delta = 2


# Computed Properties
Vol = s_box**3 #[A^3]
beta = 1/T #[K^-1]
Nh = floor( Pid*Vol*beta/kb )
rc = s_box
#rc = min([2.5*s_hh,0.5*s_box]) #[A]
ZZ = beta*Pid
UT = 0
VirT = 0

# Run the whole thing
#rhov, Env, Pv = mc_run()
#print( Pv.mean() )
#print( rhov.mean() )