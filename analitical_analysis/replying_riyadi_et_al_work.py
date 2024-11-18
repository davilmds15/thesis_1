"""
Filename: replying_riyadi_et_al_work.py
Author: Davi LIMA MENDES DOS SANTOS
Institution: Universidade Federal de Santa Maria
Date Created: 09/11/2024
Last Modified: 09/11/2024
Description:
            This code analyzes thermal stress at the interface between ceramic and metal substrates, 
            considering the impact of various parameters, including:
            1. Young's modulus of the ceramic material.
            2. Coefficient of thermal expansion of the materials.
            3. Thickness of the ceramic layer.
            4. Thickness of the metal substrate.
            5. Deposition temperature of the ceramic material.
            The code calculates and compares thermal stress variations based on these parameters.

License:
    TO DO
"""

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

# sputtered ceramic (TiN)
e_c = 600e9 # yungs modulus [GPa]
v_c = 0.25 # poisson coef      
t_c = 2.5e-6 # thickness [m]  

# metalic substrate (Fe)
e_s = 200e9 # yungs modulus [GPa]
v_s = 0.3 # poisson coef      
t_s = 3e-3 # thickness [m]     

# temperatures [ºC]
T_d = 300 # deposition      TO DO
T_r = 25 # room            

# data from figure 3: 
T = np.array([100, 250, 500, 650, 1000]) # [ºC] deposition temperature
alpha_c = np.array([0.67, 0.77, 0.86, 1]) * 10**-5 # [ºC^-1]
alpha_s = np.array([1.2, 1.5, 1.6, 1.85, 2.4]) * 10**-5 # [ºC^-1]

# interpolation of the temperature
T_alpha_c = T[:len(alpha_c)]  
interp_alpha_c = interp1d(T_alpha_c, alpha_c, fill_value="extrapolate")
interp_alpha_s = interp1d(T, alpha_s, fill_value="extrapolate")

def integrand(T):
    return interp_alpha_s(T) - interp_alpha_c(T)

# Thermal stresses of thin coatings, σc for simple planar geometry
def thermal_stress():
    e_efc = e_c / (1 - v_c)
    e_efs = e_s / (1 - v_s)
    result, error = quad(integrand, T_r, T_d) # if T_r is outside the interpolation 
    #                                           range of alpha_c/alpha_s, imprecision 
    #                                           can be generated
    integral = result #              TO DO
    numerator = e_efc * integral
    denominator = 1 + 4*(e_efc / e_efs)*(t_c / t_s)
    stress_c = numerator / denominator
    return stress_c

def thermal_stress_temp_room(T_d):
    T = T_d
    e_efc = e_c / (1 - v_c)
    e_efs = e_s / (1 - v_s)
    result, error = quad(integrand, T_r, T) # if T_r is outside the interpolation 
    #                                           range of alpha_c/alpha_s, imprecision 
    #                                           can be generated
    #print("Integration error =", error)
    integral = result #              TO DO
    numerator = e_efc * integral
    denominator = 1 + 4*(e_efc / e_efs)*(t_c / t_s)
    stress_c = numerator / denominator
    return stress_c

stress_c = thermal_stress()
print("Thermal stress =", stress_c*10**-9, "[GPa]")

T_d = np.array([100, 200, 300, 400, 500])
N = np.size(T_d)
stress = np.zeros(N)

for i in range(N):
    stress[i] = thermal_stress_temp_room(T_d[i])

print("deposition temp.", T_d, "ºC")
print("stresses", stress*10**-6, "MPa")
print("stress paper ~[200 500 800 1050 1350] MPa")