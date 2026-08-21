# =============================================================================
# FILE: Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-i2-and.py
# PAPER: An information-theoretic perspective on feed-forward loop abundances in transcriptional networks
# AUTHORS: Mintu Nandi, Sudip Chattopadhyay, and Suman K Banik
# CONTACT: mintunandi@ubi.s.u-tokyo.ac.jp; sudip@chem.iiests.ac.in;
#          skbanik@jcbose.ac.in
#
# PURPOSE
#   Runs a Gillespie simulation of the I2 matched open-loop (OL) reference in S. cerevisiae (yeast) with an AND output gate and estimates I(X, X_tilde; Z).
#   The MI is calculated from simulated means, variances, and covariances
#   using the Gaussian second-moment formula; it is not a histogram-based MI.
#
# MODEL CONDITION
#   Motif: I2; edge signs (X->Y, X->Z, Y->Z) = (-,-,-).
#   X->Y: repression; X->Z: repression; Y->Z: repression; gate: AND.
#   Simulation time: 5000000.
#   Hard-coded parameter set: <x1>,<x2>,<y>,<z>=(50.01,50.01,446.66,306.79); beta=(0.0012494,0.0012494,0.062408,0.02536); K=(217.3,157.6,25.69).
#   Production-rate scales are calculated from the operating-point conditions.
#
# INPUTS AND EXECUTION
#   No command-line arguments or external input files are used.
#   Run from the organism root so the output is collected there:
#   cd Numerics/Gillespie-after-optimization-IMI/Yeast
#   python OL/AND/yeast-ol-i2-and.py
#
# OUTPUT
#   total-MI-abund-i2ol-and-num.dat in the current working directory.
#   File format: <motif label><tab><MI value>, with six decimal places.
#
# CODE-TO-MANUSCRIPT NOTATION
#   x1 -> X, the direct-path input; x2 -> X_tilde, the independent indirect-path input.
#   a_x1,a_x2,a_y,a_z -> alpha_X,alpha_Xtilde,alpha_Y,alpha_Z.
#   b_x1,b_x2,b_y,b_z -> beta_X,beta_Xtilde,beta_Y,beta_Z.
#   Kx1z -> K_XZ; Kx2y -> the matched K_XY; Kyz -> K_YZ.
#   x1cv2,x2cv2,zcv2 -> eta_X^2, eta_Xtilde^2, eta_Z^2.
#   x1zcv2,x2zcv2 -> normalized covariances with Z.
#   zcx1x2 -> eta_Z|X,Xtilde^2; Ix1x2z -> I(X,X_tilde;Z)=I_path(X;Z).
# =============================================================================

import numpy as np

# Random seed
RANDOM_SEED = 2312
np.random.seed(RANDOM_SEED)
# import matplotlib.pyplot as plt

# Production function and derivatives
def cal_func_f(Kx2y,Kx1z,Kyz,x1,x2,y,z):
    fx1 = 1
    fx2 = 1
    fy = Kx2y/(Kx2y+x2)
    fz = (Kyz/(Kyz+y)) * (Kx1z/(Kx1z+x1)) 
    
    return fx1, fx2, fy, fz

# Function to calculate propensities
def calculate_propensities(a_x1,a_x2,b_x1,b_x2,a_y,b_y,a_z,b_z,Kx2y,Kx1z,Kyz,x1,x2,y,z):
    fx1, fx2, fy, fz= cal_func_f(Kx2y, Kx1z, Kyz, x1, x2, y, z)
    
    psx1 = a_x1 * fx1
    pdx1 = b_x1 * x1 
    psx2 = a_x2 * fx2
    pdx2 = b_x2 * x2
    psy = a_y * fy 
    pdy = b_y * y 
    psz = a_z * fz 
    pdz = b_z * z           
    return [psx1,pdx1,psx2,pdx2,psy,pdy,psz,pdz]

#Common parameters
simulation_time = 1000000*5

Ix1x2z_array = []

motif = "I2"

x1av_param = 50.01
x2av_param = 50.01
yav_param = 446.66
zav_param = 306.79

b_x1 = 0.0012494
b_x2 = 0.0012494
b_y = 0.062408
b_z = 0.02536

Kx2y = 217.30
Kx1z = 157.60
Kyz = 25.69

fx1, fx2, fy, fz = cal_func_f(Kx2y,Kx1z,Kyz,x1av_param,x2av_param,yav_param,zav_param)

a_x1 = b_x1 * x1av_param / fx1
a_x2 = b_x2 * x2av_param / fx2
a_y = b_y * yav_param / fy
a_z = b_z * zav_param / fz

x1 = x1av_param
x2 = x2av_param
y = yav_param
z = zav_param

time = 0
nu = 8   

n = 0
tot_x1 = 0.0
tot_x2 = 0.0
tot_z = 0.0

tot_x1_2 = 0.0   # sum of squared deviations for x
tot_x2_2 = 0.0   # sum of squared deviations for x
tot_z2 = 0.0   # sum of squared deviations for z
tot_x1z = 0.0   # sum for covariance between x and z
tot_x2z = 0.0   # sum for covariance between x and z

while time < simulation_time:
    
    propensities = calculate_propensities(a_x1,a_x2,b_x1,b_x2,a_y,b_y,a_z,b_z,Kx2y,Kx1z,Kyz,x1,x2,y,z)
    total_propensity = sum(propensities)

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
        x1 += 1
    elif event == 1:
        x1 -= 1
    elif event == 2:
        x2 += 1
    elif event == 3:
        x2 -= 1
    elif event == 4:
        y += 1
    elif event == 5:
        y -= 1
    elif event == 6:  
        z += 1
    elif event == 7:
        z -= 1
    else:
        continue
    
    # Upadte the steps
    n += 1
    
    tot_x1 += x1
    tot_x2 += x2
    tot_z += z
    
    tot_x1_2 += x1*x1
    tot_x2_2 += x2*x2
    tot_z2 += z*z
    
    tot_x1z += x1*z
    tot_x2z += x2*z
    
x1av = tot_x1/n
x2av = tot_x2/n
zav = tot_z/n

x1var = (tot_x1_2/n) - (x1av*x1av)
x2var = (tot_x2_2/n) - (x2av*x2av)
zvar = (tot_z2/n) - (zav*zav)

x1zcov = (tot_x1z/n) - (x1av*zav)
x2zcov = (tot_x2z/n) - (x2av*zav)

x1cv2 = x1var/(x1av**2)
x2cv2 = x2var/(x2av**2)
zcv2 = zvar/(zav**2)
x1zcv2 = x1zcov / (x1av*zav)
x2zcv2 = x2zcov / (x2av*zav)

zcx1x2 = zcv2 - ((x1zcv2**2/x1cv2) + (x2zcv2**2/x2cv2))
Ix1x2z = 0.5 * np.log2(zcv2/zcx1x2)
Ix1x2z_array.append(Ix1x2z)

with open('total-MI-abund-i2ol-and-num.dat', 'w') as f:
    f.write(f"{motif}\t{Ix1x2z_array[0]:.6f}\n")
