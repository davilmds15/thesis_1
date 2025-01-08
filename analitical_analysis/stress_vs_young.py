"""
Filename: stress_vs_young.py
Author: Davi LIMA MENDES DOS SANTOS
Institution: Universidade Federal de Santa Maria
Date Created: 09/11/2024
Last Modified: 19/12/2024
Description:
            This code finds the thermal stress at the interface between a ceramic coating and metalic substrate, 
            considering the impact the parameter indicated bellow:

            1. Deposition temperature of the ceramic material.
            2. Young's modulus of the ceramic material.            <---
            3. Coefficient of thermal expansion of the materials.   
            4. Thickness of the ceramic layer.                      
            5. Thickness of the metal substrate. 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

plt.close('all')

# Given data
T = np.array([25, 100, 200, 300, 400, 500])  # Temperature [ºC]
alpha_c = np.array([7.7, 7.7, 8.4, 9, 9.2, 9.4]) * 10**-6  # Coefficient of thermal expansion [ºC^-1]
#alpha_s = np.array([1.2, 1.22, 1.42, 1.54, 1.56, 1.6]) * 10**-5 # [ºC^-1]
alpha_s = np.array([1.22, 1.22, 1.22, 1.22, 1.22, 1.22]) * 10**-5 # [ºC^-1]

# sputtered ceramic (TiN)
N = 1000 #points
e_c = np.linspace(250e9, 600e9, N) # yungs modulus [Pa]
v_c = 0.25 # poisson coef      
t_c = 2.5e-6 # thickness [m]  

# metalic substrate (Fe)
e_s = 200e9 # yungs modulus [Pa] 193-200e9 Pa
v_s = 0.3 # poisson coef      
t_s = 3e-3 # thickness [m]  

# Interpolation of thermal expansion coefficients alpha_c and alpha_s
interp_alpha_c = interp1d(T, alpha_c, kind='linear', fill_value="extrapolate")
interp_alpha_s = interp1d(T, alpha_s, kind='linear', fill_value="extrapolate")

# Define a finer range of deposition temperatures for plotting
n = np.size(T)
T_fine = np.linspace(T[0], T[n-1], 475)  # Fine grid between 25°C and 500°C with 475 points
alpha_c_fine = interp_alpha_c(T_fine)  # Interpolated values
alpha_s_fine = interp_alpha_s(T_fine)  # Interpolated values
n = np.size(e_c)

def thermal_stress(T_r, T_d, *data):

    e_c, v_c, t_c, e_s, v_s, t_s = data
    e_efc = e_c / (1 - v_c)
    e_efs = e_s / (1 - v_s)

    # Integration using scipy.integrate.quad
    integral_c, error_c = quad(interp_alpha_c, T_r, T_d)
    integral_s, error_s = quad(interp_alpha_s, T_r, T_d)
    integral = integral_s - integral_c
    
    numerator = e_efc * integral
    denominator = 1 + 4*(e_efc/e_efs)*(t_c/t_s)

    stress = numerator / denominator
    return stress, error_c, error_s

# Defining sputtering temperature
T_d = 500 # ºC
T_r = 25 # ºC
stress = np.zeros(N)

for i in range(N):
    data = [e_c[i], v_c, t_c, e_s, v_s, t_s]
    stress[i], error_c, error_s = thermal_stress(T_r, T_d, *data)

stress_paper = np.linspace(567, 1355, N) # [MPa]

print("Inegration error of alpha_c:", error_c, "%")
print("Inegration error of alpha_s:", error_s, "%")

plt.figure(figsize=(8, 6))
plt.plot(T_fine, alpha_c_fine, label='Interpolated $\u03B1_c$', color='blue')
plt.plot(T_fine, alpha_s_fine, label='Interpolated $\u03B1_s$', color='red')
plt.scatter(T, alpha_s, color='green', label='Original Data $\u03B1_s$', zorder=5)
plt.scatter(T, alpha_c, color='purple', label='Original Data $\u03B1_c$', zorder=5)
plt.title('Coefficient of Thermal Expansion ($\u03B1$) vs. Temperature', fontsize=14)
plt.xlabel('Temperature [ºC]', fontsize=12)
plt.ylabel('Coefficient of Thermal Expansion [$ºC^{-1}$]', fontsize=12)
plt.legend()
plt.grid(True)

plt.figure(figsize=(8, 6))
plt.plot(e_c/10**9, stress/10**6, label='Calculated stress', color='blue')
plt.plot(e_c/10**9, stress_paper, label='Original paper stress', color='red')
plt.title('Thermal Stress as a function of Coating Young\'s Modulus', fontsize=14)
plt.xlabel('Coating Young\'s Modulus $E_c$ [GPa]', fontsize=12)
plt.ylabel('Thermal Stress $\\sigma_c$ [$MPa$]', fontsize=12)
plt.legend()
plt.grid(True)

plt.show()
