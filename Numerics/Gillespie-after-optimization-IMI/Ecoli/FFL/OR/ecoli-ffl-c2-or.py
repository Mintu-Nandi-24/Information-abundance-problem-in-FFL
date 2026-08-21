# =============================================================================
# FILE: Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-c2-or.py
# PAPER: An information-theoretic perspective on feed-forward loop abundances in transcriptional networks
# AUTHORS: Mintu Nandi, Sudip Chattopadhyay, and Suman K Banik
# CONTACT: mintunandi@ubi.s.u-tokyo.ac.jp; sudip@chem.iiests.ac.in;
#          skbanik@jcbose.ac.in
#
# PURPOSE
#   Runs a direct Gillespie simulation of the C2 feed-forward-loop (FFL) in E. coli with an OR output gate and estimates total Gaussian MI I(X;Z).
#   The MI is calculated from simulated means, variances, and covariances
#   using the Gaussian second-moment formula; it is not a histogram-based MI.
#
# MODEL CONDITION
#   Motif: C2; edge signs (X->Y, X->Z, Y->Z) = (-,-,+).
#   X->Y: repression; X->Z: repression; Y->Z: activation; gate: OR.
#   Simulation time: 1000000.
#   Hard-coded parameter set: <x>,<y>,<z>=(10,17.25,99.95); beta=(0.013426,0.093002,0.27098); K=(64.61,41.42,26.9).
#   Production-rate scales are calculated from the operating-point conditions.
#
# INPUTS AND EXECUTION
#   No command-line arguments or external input files are used.
#   Run from the organism root so the output is collected there:
#   cd Numerics/Gillespie-after-optimization-IMI/Ecoli
#   python FFL/OR/ecoli-ffl-c2-or.py
#
# OUTPUT
#   total-MI-abund-c2ffl-or-num.dat in the current working directory.
#   File format: <motif label><tab><MI value>, with six decimal places.
#
# CODE-TO-MANUSCRIPT NOTATION
#   a_x,a_y,a_z -> alpha_X,alpha_Y,alpha_Z; b_x,b_y,b_z -> beta_X,beta_Y,beta_Z.
#   Kxy,Kxz,Kyz -> K_XY,K_XZ,K_YZ.
#   xav_param,yav_param,zav_param -> <x>,<y>,<z>.
#   xcv2,zcv2,xzcv2 -> eta_X^2, eta_Z^2, zeta_XZ.
#   Ixz -> total mutual information I(X;Z).
# =============================================================================

import numpy as np

# Random seed
RANDOM_SEED = 1202
np.random.seed(RANDOM_SEED)
# import matplotlib.pyplot as plt

# Production function and derivatives
def cal_func_f(Kxy,Kxz,Kyz,x,y,z):
    fx = 1
    fy = Kxy/(Kxy+x)
    fzx = Kxz/(Kxz+x)
    fzy = y/(Kyz+y)
    
    return fx, fy, fzx, fzy

# Function to calculate propensities
def calculate_propensities(a_x,b_x,a_y,b_y,a_z,b_z,Kxy,Kxz,Kyz,x,y,z):
    fx, fy, fzx, fzy = cal_func_f(Kxy, Kxz, Kyz, x, y, z)
    
    psx = a_x * fx
    pdx = b_x * x 
    psy = a_y * fy 
    pdy = b_y * y 
    pszx = a_z * fzx
    pszy = a_z * fzy
    pdz = b_z * z           
    return [psx,pdx,psy,pdy,pszx,pszy,pdz]

#Common parameters
simulation_time = 1000000

Ixz_array = []

motif = "C2"

xav_param = 10.00
yav_param = 17.25
zav_param = 99.95

b_x = 0.013426
b_y = 0.093002
b_z = 0.27098

Kxy = 64.61
Kxz = 41.42
Kyz = 26.90

fx, fy, fzx, fzy = cal_func_f(Kxy,Kxz,Kyz,xav_param,yav_param,zav_param)

a_x = b_x * xav_param / fx  
a_y = b_y * yav_param / fy
a_z = b_z * zav_param / (fzx + fzy)

x = xav_param
y = yav_param
z = zav_param

time = 0
nu = 7   

n = 0
tot_x = 0.0
tot_z = 0.0

tot_x2 = 0.0   # sum of squared deviations for x
tot_z2 = 0.0   # sum of squared deviations for z
tot_xz = 0.0   # sum for covariance between x and z

while time < simulation_time:
    
    propensities = calculate_propensities(a_x,b_x,a_y,b_y,a_z,b_z,Kxy,Kxz,Kyz,x,y,z)
    total_propensity = sum(propensities)
    
    if total_propensity <= 0.0:
        break

    # Calculate time until the next event
    delta_t = -np.log(np.random.random()) / total_propensity
    time += delta_t

    # Choose which event occurs
    r2 = np.random.random()
    r2a0 = r2 * total_propensity
    sum2 = 0.0
    event = -1

    for j in range(nu):
        sum2 += propensities[j]
        if sum2 >= r2a0:
            event = j
            break

    if event == 0:  
        x += 1
    elif event == 1:
        x -= 1
    elif event == 2:
        y += 1
    elif event == 3:
        y -= 1
    elif event == 4:  
        z += 1
    elif event == 5:  
        z += 1
    elif event == 6:
        z -= 1
    else:
        continue
    
    # Upadte the steps
    n += 1
    
    tot_x += x
    tot_z += z
    
    tot_x2 += x*x
    tot_z2 += z*z
    
    tot_xz += x*z
    
xav = tot_x/n
zav = tot_z/n

xvar = (tot_x2/n) - (xav*xav)
zvar = (tot_z2/n) - (zav*zav)

xzcov = (tot_xz/n) - (xav*zav)

xcv2 = xvar/(xav**2)
zcv2 = zvar/(zav**2)
xzcv2 = xzcov / (xav*zav)

Ixz = 0.5 * np.log2(xcv2*zcv2/(xcv2*zcv2 - xzcv2**2))
Ixz_array.append(Ixz)

with open('total-MI-abund-c2ffl-or-num.dat', 'w') as f:
    f.write(f"{motif}\t{Ixz_array[0]:.6f}\n")
