# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 14:43:48 2022

@authors: Cesar Diaz, Ben Furst, Jorge Muñoz, Takuro Daimaru, Scott Roberts

cesar.dc1509@gmail.com, ben.i.furst@jpl.nasa.gov, jamunoz@utep.edu, 
takuro.daimaru@jpl.nasa.gov, scott.n.roberts@jpl.nasa.gov

The University of Texas at El Paso
NASA Jet Propulsion Laboratory, California Institute of Technology

The following program was developed to calculate the operational limits of
Oscilating Heat Pipes (OHP), based on the analytical limts derived by Drolen 
and Smoot (Performance Limits of Oscillating Heat Pipes: Theory and 
Validation), being the collaboration a joint research effort between the JPL 
Thermal Additive Group and the Research Group of Prof. Jorge Muñoz at the 
University of Texas at El Paso.

The code accounts for calculating the vapor inertia limit, critical heat flux
limit, sonic limit, viscous limit, evaporation fraction, swept length limit, 
and modified bond number limit, as outlined in Drolen and Smoot's paper. 
Property tables are downloaded from the NIST Refprop Data Base and imported 
into Python as csv files. All properties are taken to be at the adiabatic 
temperature on turn as saturation temperature.

All properties must be in "fundamental SI Units" when inserted into the 
defined calculation functions. Proper conversions must be performed in the 
calculations or plotting sections before calling the methods. Inputs to the 
declared functions can be numpy arrays of any amount of elements, as required 
per used needs.

References:
- Drolen, Bruce L., and Christopher D. Smoot. "Performance limits of oscillating 
heat pipes: Theory and validation." Journal of Thermophysics and Heat Transfer 
31.4 (2017): 920-936.
- Diaz-Caraveo, C., Wolk, K., Miesner, S., Montemayor, M., Rodriguez, A., Kumar, V., 
Muñoz, J.A., Daimaru, T., Furst, B.I., & Roberts, S. N. Performance-Dryout Limits 
of Oscillating Heat Pipes: A Comprehensive Theoretical and Experimental Determination. 
Journal of Thermophysics and Heat Transfer, 1-11 (2024).

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve

#%% Parameters declaration section
# Variables are just declared here for purposes of being used across different 
# files.

# Important Note: "Values do not neet to be modified by the user here, all of 
# them are declared empty"

T_cold_plate = None

# Geomertical Features

# Number of channels present on the evaporator
N = None
# Number of condensors on the OHP
Cn = None
# Channel Diameter
diameter = None # m
# Length per channel
L_channel = None # m
# Condensor section length per channel
L_cond = None # m
# Heater section length per channel
L_evap = None # m


# Fill fraction, volume, and mass

# Total OHP channel volume
V_total_OHP = None
# Total working fluid mass
mass_total = None
# Fill Fraction (decimal)
fill_fraction = None
#Average density 
avg_density = None


# Contact angles

# Static contact angle
theta_s = None
# Receding contact angle 
theta_r = None


# Other manual imported parameters

# Evaporation Fraction
evap_fraction_value_manual = None
# Average Channel Velocity
average_channel_velocity = None
# Heat Transfer Coefficient
heat_transfer_coeff_manual = None
# Temperature increment Evaporator Condensor
temp_increment_evap_cond = None
# OHP Name of Identification Data
name = None

#%% Nomenclature

# q - Heat Flux in W/m^2
# Q - Heating Power in W
# f - Fill Fraction
# st - Surface Tension
# pv - Vapor Density
# pl - Liquid Density
# hfg - Vaporization Enthalpy
# evap_fraction - Evaporation Fraction
# theta_r - Receding Contact Angle
# theta_a - Advancing Contact Angle

# All equation numbers refer to Drolen and Smooth (2017)


#%% 1.- Modified Bond Number Limit

# Equation for calculating the critical diameter of the OHP (eqn 6)
def dh_critical(surface_tension, liquid_density, vapor_density, gravity):
    
    return 2.74 * np.sqrt(surface_tension / (gravity*(liquid_density - vapor_density)))

# Equation for calulcating the square root of the Bond Number (eotvos number), using eqn 6
def eotvos_number_calc(surface_tension, liquid_density, vapor_density, gravity, 
                    diameter):
    
    return diameter * np.sqrt(gravity*(liquid_density - vapor_density) / 
                              surface_tension)


#%% 2.- Vapor Inertia Limit Equations

# Defining Q max vapor per channel laminar (eqn 17)
def vapor_inertia_Q_max_laminar_per_channel(st, pv, hfg, r0, theta_r, theta_a,
                                        evap_fraction):
    
    return (1.225/evap_fraction)*math.pi*r0*hfg*np.sqrt(st*pv*r0*
                                    (np.cos(theta_r) - np.cos(theta_a)))


# Defining q max laminar (eqn 18)
def vapor_inertia_q_max_laminar(st, pv, hfg, r0, theta_r, theta_a, evap_fraction):
    
    return (1.225/evap_fraction) * (np.sqrt(st*pv*(hfg**2)/r0) * 
                                   np.sqrt(np.cos(theta_r) - np.cos(theta_a)))


# Defining Q max vapor per channel turbulent (eqn 21)
def vapor_inertia_Q_max_turbulent_per_channel(st, vapor_density, hfg, r0,
                                      theta_r, theta_a, pv, evap_fraction):
    
    return (1.145/evap_fraction)*math.pi*r0*hfg*np.sqrt(st*pv*r0*(np.cos(theta_r) - 
                                                            np.cos(theta_a)))


# Defining q max laminar (eqn 22)
def vapor_inertia_q_max_turbulent(st, pv, hfg, r0, theta_r, theta_a, evap_fraction):
    
    return (1.145/evap_fraction)*(np.sqrt(st*pv*(hfg**2)/r0) * 
                          np.sqrt(np.cos(theta_r) - np.cos(theta_a)))


# Defining Q max laminar, net heat transport capability (eqn 20)
def vapor_inertia_Q_max_laminar(st, pv, hfg, r0, N, Cn,
                                theta_r, theta_a, evap_fraction):
    
    return ((1.924*N*Cn*r0*hfg)/evap_fraction * 
            np.sqrt(st*pv*r0*(np.cos(theta_r) - np.cos(theta_a))))


# Defining Q max turbulent, net heat transport capability (eqn 23)
def vapor_inertia_Q_max_turbulent(st, pv, hfg, r0, N, Cn, theta_r, theta_a, evap_fraction):
    
    return ((1.8*N*Cn*r0*hfg)/evap_fraction * 
            np.sqrt(st*pv*r0*(np.cos(theta_r) - np.cos(theta_a))))

def vapor_inertia_Q_max_Yin_et_al(visc_l, hfg, L_channel, fill_fraction):
    
    return 3*math.pi*visc_l*hfg*L_channel*fill_fraction

#%% 3.- Critical Heat Flux Limit

# Defining maximum q in the channel for Critical Heat Flux Limit (eqn 33)
def chf_limit_q_channel(Lh, r0, pv, pl, st, hfg, evap_fraction):
    
    return (40147/evap_fraction)*((Lh/r0)**1.02)*((pv/pl)**0.76)*(st*pv*(hfg**2)/r0)**0.5


# Defining maximum Q in the channel for CHF (eqn 34)
def chf_limit_Q_total(Cn, N, Lh, r0, pv, pl, st, hfg, evap_fraction):
    
    return ((126128/evap_fraction)*(Cn*N/2)*((Lh/r0)**1.02)*((pv/pl)**0.76)*
            (st*pv*(hfg**2)*(r0**3))**0.5)


#%% 4.- Sonic Limit

# Defining maximum Q for the whole OHP for Sonic Limit (eqn 47)
def sonic_limit_Q_total(Cn, N, Dh, hfg, f, pv, pl, T_sat, cv, evap_fraction):
    
    return (((Cn*N*math.pi*(Dh**2)) / (8*evap_fraction*(f*(pl/pv)+(1-f))*((1/pv)-(1/pl)))) * 
            ((hfg**2) / np.sqrt(cv*T_sat)))


# Defining maximum q per channel (eqn 48)
def sonic_limit_q_channel(hfg, evap_fraction, f, pl, pv, cv, T_sat):
    
    return (1/(evap_fraction*(f*(pl/pv)+(1-f))*((1/pv)-(1/pl)))) * ((hfg**2) / np.sqrt(cv*T_sat))


#%% 5.- Viscous Limit

# Defining maximum Q for the whole OHP for Viscous Limit (eqn 55)
def viscous_limit_Q_total(r0, N, Cn, L_channel, evap_fraction, f, pv, pl, hfg, T_sat, 
                              visc_l, T_evap, T_cond):
    
    return ((math.pi*(r0**4)*N*(Cn**2) / (32*evap_fraction*L_channel*f)) *
            ((pv**2)*pl*(hfg**2) / (visc_l*T_sat*(pl-pv))) * (T_evap - T_cond))


# Defining maximum q per channel (eqn 56)
def viscous_limit_q_channel_max(r0, Cn, L_channel, evap_fraction, f, pv, pl, hfg, visc_l, 
                              T_sat, T_evap, T_cond):
    
    return (((r0**2)*Cn / (16*evap_fraction*L_channel*f)) * 
            (((pv**2)*pl*(hfg**2)) / (visc_l*T_sat*(pl-pv))) * (T_evap-T_cond))


#%% 6.- Swept Length Limit (eqn 79)

def swept_length_limit_q_reduced(q_max, l_s, l_htr):
    
    return q_max * (l_s / l_htr) ** 1.35


#%% 7.- Evaporation fraction calculation (eqn 64)

def evaporation_fraction_calc(pl, cpl, pv, hfg, T_adiab, T_cond, r0, h, L_cond, m_dot_l):
    
    return (1 / (((pl * cpl) / (pv * hfg)) * (T_adiab - T_cond) * 
                 (1 - np.exp((-2*math.pi*r0*h*L_cond)/(m_dot_l*cpl))) + 1))


#%% 8.- Heat Transfer Coefficient calculation (eqn 62)

# Calculation for the value of the heat transfer coefficient for sensible 
# heating of the liquid slug in turbulent flow. Authors assume this to be as 
# the coefficient of heat transfer to be used in the evaporation fraction 
# equation. 

def h_sens_calculator(kl, u, pl, visc_l, cp_l, diameter):
    
    pr_l = prandtl_number_calc(cp_l, visc_l, kl)
    re_l = reynolds_number_calc(pl, visc_l, diameter, u)
    
    return 0.023 * (kl / diameter) * re_l**0.8 * pr_l**0.4


#%% 9.- Fill fraction calculation

# Fill fraction is typically a known/measured property of an OHP, but it is 
# referenced to room temeprature. The change of density of the working fluid 
# with temperature leads to changes in fill fraction inside the OHP, with can 
# be easily incorporated into the OHP model with the following equation.

def fill_fraction_calc(pv, pl):
    
    # Average density
    p = avg_density
    
    return (p - pv)/(pl - pv)


#%% 10.- Average velocity calculators

# Calculations for the average velocity on channel pipe as defined in the 
# Vapor Inertia Limit Section

# Defining average velocity for laminar flow (eqn 13)
def average_velocity_laminar(st, pv, r0, theta_r, theta_a):
    
    return 1.225 * np.sqrt((st/(r0*pv)) * (np.cos(theta_r) - np.cos(theta_a)))


# Defining average velocity for turbulent flow (eqn 15)
def average_velocity_turbulent(st, pv, r0, theta_r, theta_a):
    
    return 1.145 * np.sqrt((st/(r0*pv)) * (np.cos(theta_r) - np.cos(theta_a)))


#%% 11.- Dynamic contact angle calculation

# Calculation of the dynamic contact angle, which according to Drolen and Smooth,
# it is taken as the advancing contact angle (theta_a)

# Note: as per math library requirement, angle inputs for the function must be 
# in radians, and the result will also be returned in radians. Proper 
# conversions must be performed before inputing values into the equation. 

def dynamic_contact_angle_calc(theta_s, ca):
    
    return np.arccos(np.cos(theta_s) - ((np.cos(theta_s) + 1) * 
                                          np.tanh(4.96*(ca**0.702))))


#%% 12.- Dynamic contact angle and average velocity solver

# Taking the the equations for the dynamic contact angle and plugging inside 
# the average velocity equations, we get an expression to solve for the dynamic
# contact angle, and once we get the dynamic contact angle we calculate the 
# average velocity

# Equations are solved using the SciPy fsolve function 

def solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l):
    
    
    if(isinstance(st, (np.floating, float, int)) and isinstance(r0, (np.floating, float, int))
        and isinstance(pv, (np.floating, float, int)) and isinstance(theta_r, (np.floating, float, int))
        and isinstance(theta_s, (np.floating, float, int)) and isinstance(visc_l, (np.floating, float, int))):
        
        
        function = lambda theta_a :  (np.arccos(np.cos(theta_s) - 
                                     np.tanh(4.96*((1.225*visc_l/st) * 
                                     np.sqrt((st/(r0*pv)) * (np.cos(theta_r) - 
                                         np.cos(theta_a))))**0.702) * (np.cos(theta_s) + 1)) -
                                     theta_a)
        
        theta_a_first_guess = math.pi/4
        theta_a_solution = fsolve(function, theta_a_first_guess)
        
        avg_velocity = average_velocity_laminar(st, pv, r0, theta_r, 
                                                theta_a_solution)
        
        return avg_velocity
    
    
    if (len(st) != len(r0) or len(st) != len(pv)):
       error_message = ("If not a single value, property arrays passed into" + 
                         "the functions must have the same length")
        
       raise Exception(error_message)
    
    
    else:
        
        function = lambda theta_a :  (np.arccos(np.cos(theta_s) - 
                                     np.tanh(4.96*((1.225*visc_l/st) * 
                                     np.sqrt((st/(r0*pv)) * (np.cos(theta_r) - 
                                         np.cos(theta_a))))**0.702) * (np.cos(theta_s) + 1)) -
                                     theta_a)
        
        theta_a_first_guess = np.full(len(st), math.pi/4)
        theta_a_solution = fsolve(function, theta_a_first_guess)
        
        avg_velocity = average_velocity_laminar(st, pv, r0, theta_r, 
                                                theta_a_solution)
        
        return avg_velocity


def solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l):
    
    if(isinstance(st, (np.floating, float, int)) 
       and isinstance(r0, (np.floating, float, int))
        and isinstance(pv, (np.floating, float, int)) and isinstance(theta_r, (np.floating, float, int))
        and isinstance(theta_s, (np.floating, float, int)) and isinstance(visc_l, (np.floating, float, int))):
        
        
        function = lambda theta_a :  np.arccos(np.cos(theta_s) - 
                                     np.tanh(4.96*((1.145*visc_l/st) * 
                                     np.sqrt((st/(r0*pv)) * (np.cos(theta_r) - 
                                         np.cos(theta_a))))**0.702) * (np.cos(theta_s) + 1) -
                                     theta_a)
        
        theta_a_first_guess = math.pi/4
        theta_a_solution = fsolve(function, theta_a_first_guess)
        
        avg_velocity = average_velocity_laminar(st, pv, r0, theta_r, 
                                                theta_a_solution)
        
        return avg_velocity
    
    
    if (len(st) != len(r0) or len(st) != len(pv)):
        error_message = ("If not a single value, property arrays passed into" + 
                         "the functions must have the same length")
        
        raise Exception(error_message)
    
    
    else:
        
        function = lambda theta_a :  np.arccos(np.cos(theta_s) - 
                                     np.tanh(4.96*((1.145*visc_l/st) * 
                                     np.sqrt((st/(r0*pv)) * (np.cos(theta_r) - 
                                         np.cos(theta_a))))**0.702) * (np.cos(theta_s) + 1) -
                                     theta_a)
        
        theta_a_first_guess = np.full(len(st), math.pi/4)
        theta_a_solution = fsolve(function, theta_a_first_guess)
        
        avg_velocity = average_velocity_laminar(st, pv, r0, theta_r, 
                                                theta_a_solution)
        
        return avg_velocity


#%% 13.- OHP Transport Factors

def inertial_factor_calc(hfg, st, pv):
    
    return hfg * np.sqrt(st*pv)


def chf_factor_calc(hfg, st, pv, pl):
    
    return (pv/pl)**0.76 * hfg * (st*pv)**0.5


def viscous_factor_calc(hfg, pv, pl, visc_l, T_sat):
    
    return (pv**2 * pl * hfg**2) / (visc_l * T_sat * (pl - pv))


def swept_length_factor_calc(pv, pl, st, g):
    
    return ((pl/pv)**0.5) * (st/(g*(pl - pv)))**0.25
    
    # Note for swept length - Gravity cannot be zero because of the equations
    # Take into account when calculating performance in orbit. 


#%% 14.- OHP Transport Factors Plotting Functions

def plotting_OHP_inertial_factor(file_name, refrigerant, color, linestyle):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    inertial_factors_array = np.zeros([len(saturation_table[:, 1])])
    
    i = 0
    
    for temperature in saturation_table[:, 1]:
        
        hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
        st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
        pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
        
        #print('Temperature =', temperature)
        #print('Heat of Vapor =', hfg)
        #print('Surf. Tension =', st)
        #print('Liquid Density (kg/m3) =', pv)
        
        inertial_factor = inertial_factor_calc(hfg, st, pv)
        inertial_factors_array[i] = inertial_factor
        
        i += 1
    
    #print('\n')
    
    plt.plot(saturation_table[:, 1], inertial_factors_array, color = color, 
             linestyle = linestyle, label = refrigerant)


def plotting_OHP_chf_factor(file_name, refrigerant, color, linestyle):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    chf_factors_array = np.zeros([len(saturation_table[:, 1])])
    
    i = 0
    
    for temperature in saturation_table[:, 1]:
        
        hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
        st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
        pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
        pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
        
        chf_factor = chf_factor_calc(hfg, st, pv, pl)
        chf_factors_array[i] = chf_factor
        
        i += 1
    
    plt.plot(saturation_table[:, 1], chf_factors_array, color = color, 
             linestyle = linestyle, label = refrigerant)


def plotting_OHP_viscous_factor(file_name, refrigerant, color, linestyle):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    viscous_factors_array = np.zeros([len(saturation_table[:, 1])])
    
    i = 0
    
    for temperature in saturation_table[:, 1]:
        
        hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
        pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
        pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
        visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
        
        temperature += 273.15 # Converstion from Celsius to Kelvin 
        
        viscous_factor = viscous_factor_calc(hfg, pv, pl, visc_l, temperature)
        viscous_factors_array[i] = viscous_factor
        
        i += 1
        
    plt.plot(saturation_table[:, 1], viscous_factors_array, color = color, 
             linestyle = linestyle, label = refrigerant)
    
    
def plotting_OHP_viscous_factor_logaritmic(file_name, refrigerant, color, 
                                           linestyle):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    viscous_factors_array = np.zeros([len(saturation_table[:, 1])])
    
    i = 0
    
    for temperature in saturation_table[:, 1]:
        
        hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
        pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
        pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
        visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
        
        temperature += 273.15 # Converstion from Celsius to Kelvin 
        
        viscous_factor = viscous_factor_calc(hfg, pv, pl, visc_l, temperature)
        viscous_factors_array[i] = viscous_factor
        
        i += 1
        
    plt.plot(saturation_table[:, 1], viscous_factors_array, color = color, 
             linestyle = linestyle, label = refrigerant)


def plotting_OHP_swept_length_factor(file_name, refrigerant, color, linestyle):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    swept_length_factors_array = np.zeros([len(saturation_table[:, 1])])
    
    i = 0
    
    for temperature in saturation_table[:, 1]:
        
        st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
        pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
        pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
        
        swept_length_factor = swept_length_factor_calc(pv, pl, st, 9.81)
        swept_length_factors_array[i] = swept_length_factor
        
        i += 1
        
    plt.plot(saturation_table[:, 1], swept_length_factors_array, color = color, 
             linestyle = linestyle, label = refrigerant)
    
    
def plotting_OHP_swept_length_factor_logaritmic(file_name, refrigerant, color, 
                                                linestyle):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    swept_length_factors_array = np.zeros([len(saturation_table[:, 1])])
    
    i = 0
    
    for temperature in saturation_table[:, 1]:
        
        st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
        pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
        pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
        
        swept_length_factor = swept_length_factor_calc(pv, pl, st, 9.81)
        swept_length_factors_array[i] = swept_length_factor
        
        i += 1
        
    plt.plot(saturation_table[:, 1], swept_length_factors_array, color = color, 
             linestyle = linestyle, label = refrigerant)
    plt.yscale('log')


#%% 15.- Max Channel Heat Flux Plotting Functions

# Main set of functions to be called to calculate the performance limits of an 
# Oscillating Heat Pipes. Callings to the methods from the paper above to
# perform the calculations are implemented here.

# Notes:
#   - On Matplotlib Documentation, alpha refers to the transparency value of the 
#     line being plotted, and it is a parameter that can be passed to the 
#     plt.plot function. A value of alpha = 1 means a solid line and alpha = 0 
#     means a fully transparent line. 


def plotting_bond_number_limit(file_name, limit_label, color, linestyle, alpha, 
                               ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    i = 0
    
    global diameter
    g = 9.81
    
    for temperature in saturation_table[:, 1]:
        
        st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
        pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
        pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
        
        eotvos_number = eotvos_number_calc(st, pl, pv, g, diameter)
        
        if(eotvos_number >= 2.74):
            break
    
    ax.axvline(temperature, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)


def plotting_max_heat_flux_inertial(file_name, limit_label, color, linestyle, 
                                    theta_r, theta_s, alpha, T_cold_plate, fluid, 
                                    ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    inertial_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    theta_a_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature.
        if(temperature > (T_cold_plate)): 
            
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Conversion into Kelvin
            temperature += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            
            ca = capilarity_number_calc(u, visc_l, st)
            theta_a = dynamic_contact_angle_calc(theta_s, ca)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                ca = capilarity_number_calc(u, visc_l, st)
                theta_a = dynamic_contact_angle_calc(theta_s, ca)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, 
                                                          temperature, T_cond, r0, 
                                                          h, L_cond, mass_flow_rate_l)
                
                inertial_limit = vapor_inertia_q_max_turbulent(st, pv, hfg, r0, 
                                                               theta_r, theta_a, 
                                                               evap_fraction)
            else:
                inertial_limit = vapor_inertia_q_max_laminar(st, pv, hfg, r0, 
                                                             theta_r, theta_a, 
                                                             evap_fraction)
            
            
            inertial_limits_array[i] = inertial_limit
            
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            theta_a_array[i] = theta_a
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    inertial_limit_data = np.hstack([temperatures, inertial_limits_array, 
                                     velocities_array, mass_flow_rates_array, 
                                     theta_a_array, h_array, 
                                     evaporation_fractions_array])
    
    
    df_inertial_limits_array = pd.DataFrame(inertial_limit_data, 
                                            columns = ('Temperature', 
                                                       'Inertial Limits (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       'Advancing Angle (rad)', 
                                                       'Heat Transfer Coefficient (w/m2K)', 
                                                       'Evaporation Fraction'))
    df_inertial_limits_array.to_excel('Inertial Limit Data {} - CDO.xlsx'.format(name))
    
    ax.plot(saturation_table[:, 1], inertial_limits_array, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)


def plotting_max_heat_flux_CHF(file_name, limit_label, color, linestyle, alpha, 
                               T_cold_plate, fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    chf_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature. Single conversion back into celsius.
        if(temperature > (T_cold_plate)): 
            
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Conversion into Kelvin
            temperature += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)

            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                          T_cond, r0, h, L_cond, mass_flow_rate_l)
                
                chf_limit = chf_limit_q_channel(L_evap, r0, pv, pl, st, hfg, evap_fraction)
                
            else:
                chf_limit = chf_limit_q_channel(L_evap, r0, pv, pl, st, hfg, evap_fraction)
            
            chf_limits_array[i] = chf_limit
            
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    chf_limit_data = np.hstack([temperatures, chf_limits_array, 
                                velocities_array, mass_flow_rates_array, 
                                h_array, evaporation_fractions_array])
    
    df_chf_limit_data = pd.DataFrame(chf_limit_data, columns = ('Temperature', 
                                                       'CHF Limit (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       'Heat Transfer Coefficient (w/m2)', 
                                                       'Evaporation Fraction'))
    df_chf_limit_data.to_excel('CHF Limit Data {} - CDO.xlsx'.format(name))
    
    ax.plot(saturation_table[:, 1], chf_limits_array, color = color,
             linestyle = linestyle, label = limit_label, alpha = alpha)


def plotting_max_heat_flux_sonic(file_name, limit_label, color, linestyle, alpha, 
                                 T_cold_plate, fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    sonic_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature. Single conversion back into celsius.
        if(temperature > (T_cold_plate)): 
        
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            cv_vapor = df_saturation_table.loc[i, "Vapor Cv (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Conversion into Kelvin
            temperature += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            
            #h = 10000
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            sonic_limit = sonic_limit_q_channel(hfg, evap_fraction, fill_fraction, 
                                                pl, pv, cv_vapor, temperature)
            
            if(re > 2200):
                
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                          T_cond, r0, h, L_cond, mass_flow_rate_l)
                
                sonic_limit = sonic_limit_q_channel(hfg, evap_fraction, fill_fraction, 
                                                    pl, pv, cv_vapor, temperature)
                
            else:
                sonic_limit = sonic_limit_q_channel(hfg, evap_fraction, fill_fraction, 
                                                    pl, pv, cv_vapor, temperature)
            
            
            sonic_limits_array[i] = sonic_limit
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    sonic_limit_data = np.hstack([temperatures, sonic_limits_array, 
                                     velocities_array, mass_flow_rates_array, 
                                     h_array, evaporation_fractions_array])
    
    df_sonic_limit_data = pd.DataFrame(sonic_limit_data, columns = ('Temperature', 
                                                       r'Sonic Limit (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       r'Heat Transfer Coefficient (w/m2)', 
                                                       'Evaporation Fraction'))
    df_sonic_limit_data.to_excel('Sonic Limit Data {} - CDO.xlsx'.format(name))
    
    ax.plot(saturation_table[:, 1], sonic_limits_array, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)
    

def plotting_max_heat_flux_viscous(file_name, limit_label, color, linestyle, 
                                   alpha, T_cold_plate, fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    viscous_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0    
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature. Single conversion back into celsius.
        if(temperature > (T_cold_plate)):
            
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Assumption of evaporator and condensor temperatures being 1 degrees
            # above and below the adiabatic temperature, according to Drolen and 
            # Smoot
            T_evap_viscous = temperature + 1
            T_cond_viscous = temperature - 1
            
            # Conversion into Kelvin
            temperature += 273.15
            T_evap_viscous += 273.15
            T_cond_viscous += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                          T_cond, r0, h, L_cond, mass_flow_rate_l)
                
                viscous_limit = viscous_limit_q_channel_max(r0, Cn, L_channel, evap_fraction, 
                                                            fill_fraction, pv, pl, hfg, visc_l, 
                                                            temperature, T_evap_viscous, T_cond_viscous)
                
            else:
                viscous_limit = viscous_limit_q_channel_max(r0, Cn, L_channel, evap_fraction, 
                                                            fill_fraction, pv, pl, hfg, visc_l, 
                                                            temperature, T_evap_viscous, T_cond_viscous)
            
            viscous_limits_array[i] = viscous_limit
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
        
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    viscous_limit_data = np.hstack([temperatures, viscous_limits_array, 
                                     velocities_array, mass_flow_rates_array, 
                                     h_array, evaporation_fractions_array])
    
    df_viscous_limit_data = pd.DataFrame(viscous_limit_data, columns = ('Temperature', 
                                                       'Viscous Limit (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       'Heat Transfer Coefficient (w/m2)', 
                                                       'Evaporation Fraction'))
    df_viscous_limit_data.to_excel('Viscous Limit Data {} - CDO.xlsx'.format(name))
    
    ax.plot(saturation_table[:, 1], viscous_limits_array, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)


#%% 16.- Max Channel Heat Load Plotting Functions

# Main set of functions to be called to calculate the performance limits of an 
# Oscillating Heat Pipes. Callings to the methods from the paper above to
# perform the calculations are implemented here.

# Notes:
#   - On Matplotlib Documentation, alpha refers to the transparency value of the 
#     line being plotted, and it is a parameter that can be passed to the 
#     plt.plot function. A value of alpha = 1 means a solid line and alpha = 0 
#     means a fully transparent line.


def plotting_max_heat_load_inertial(file_name, limit_label, color, linestyle, 
                                    theta_r, theta_s, alpha, T_cold_plate, 
                                    fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    inertial_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    theta_a_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature.
        if(temperature > (T_cold_plate)): 
            
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Conversion into Kelvin
            temperature += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            
            ca = capilarity_number_calc(u, visc_l, st)
            theta_a = dynamic_contact_angle_calc(theta_s, ca)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                ca = capilarity_number_calc(u, visc_l, st)
                theta_a = dynamic_contact_angle_calc(theta_s, ca)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, 
                                                          temperature, T_cond, r0, 
                                                          h, L_cond, mass_flow_rate_l)
                
                inertial_limit = vapor_inertia_Q_max_turbulent(st, pv, hfg, r0, 
                                                               N, Cn, theta_r, 
                                                               theta_a, evap_fraction)
            else:
                inertial_limit = vapor_inertia_Q_max_laminar(st, pv, hfg, r0, 
                                                               N, Cn, theta_r, 
                                                               theta_a, evap_fraction)
            
            
            inertial_limits_array[i] = inertial_limit
            
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            theta_a_array[i] = theta_a
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    inertial_limit_data = np.hstack([temperatures, inertial_limits_array, 
                                     velocities_array, mass_flow_rates_array, 
                                     theta_a_array, h_array, 
                                     evaporation_fractions_array])
    
    
    df_inertial_limits_array = pd.DataFrame(inertial_limit_data, 
                                            columns = ('Temperature', 
                                                       'Inertial Limits (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       'Advancing Angle (rad)', 
                                                       'Heat Transfer Coefficient (w/m2K)', 
                                                       'Evaporation Fraction'))
    df_inertial_limits_array.to_excel('Inertial Limit Data {} {} C Cold Plate- CDO.xlsx'.format(name, T_cold_plate))
    
    ax.plot(saturation_table[:, 1], inertial_limits_array, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)
    
    return inertial_limits_array


def plotting_max_heat_load_inertial_Yin_et_al(file_name, limit_label, color, linestyle, 
                                              theta_r, theta_s, alpha, T_cold_plate, 
                                              fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    inertial_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature.
        if(temperature > (T_cold_plate + 5)): 
            
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Conversion into Kelvin
            temperature += 273.15
            T_cond += 273.15
            
            fill_fraction_var = fill_fraction_calc(pv, pl)
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            
            ca = capilarity_number_calc(u, visc_l, st)
            theta_a = dynamic_contact_angle_calc(theta_s, ca)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                ca = capilarity_number_calc(u, visc_l, st)
                theta_a = dynamic_contact_angle_calc(theta_s, ca)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, 
                                                          temperature, T_cond, r0, 
                                                          h, L_cond, mass_flow_rate_l)
                
                inertial_limit = vapor_inertia_Q_max_Yin_et_al(visc_l, hfg, L_channel, fill_fraction_var)
            
            else:
                inertial_limit = vapor_inertia_Q_max_Yin_et_al(visc_l, hfg, L_channel, fill_fraction_var)
                

            inertial_limits_array[i] = inertial_limit
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    inertial_limits_data = np.hstack([temperatures, inertial_limits_array])
    
    df_inertial_limits_Yin_et_al_array = pd.DataFrame(inertial_limits_data, 
                                                      columns = ('Temperature', 
                                                       'Inertial Limits (w/m2)'))
    df_inertial_limits_Yin_et_al_array.to_excel('Inertial Limit Data Yin et al. {} - CDO.xlsx'.format(name))
    
    ax.plot(saturation_table[:, 1], inertial_limits_array, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)


def plotting_max_heat_load_CHF(file_name, limit_label, color, linestyle, alpha, 
                               T_cold_plate, fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    chf_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature. Single conversion back into celsius.
        if(temperature > (T_cold_plate + 5)): 
            
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Conversion into Kelvin
            temperature += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                          T_cond, r0, h, L_cond, mass_flow_rate_l)
                
                chf_limit = chf_limit_Q_total(Cn, N, L_evap, r0, pv, pl, st,
                                              hfg, evap_fraction)
                
            else:
                chf_limit = chf_limit_Q_total(Cn, N, L_evap, r0, pv, pl, st, 
                                              hfg, evap_fraction)
            
            chf_limits_array[i] = chf_limit
            
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    chf_limit_data = np.hstack([temperatures, chf_limits_array, 
                                velocities_array, mass_flow_rates_array, 
                                h_array, evaporation_fractions_array])
    
    df_chf_limit_data = pd.DataFrame(chf_limit_data, columns = ('Temperature', 
                                                       'CHF Limit (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       'Heat Transfer Coefficient (w/m2)', 
                                                       'Evaporation Fraction'))
    df_chf_limit_data.to_excel('CHF Limit Data {} - CDO.xlsx'.format(name))
    
    ax.plot(saturation_table[:, 1], chf_limits_array, color = color,
             linestyle = linestyle, label = limit_label, alpha = alpha)


def plotting_max_heat_load_sonic(file_name, limit_label, color, linestyle, alpha, 
                                 T_cold_plate, fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    sonic_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0    
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature. Single conversion back into celsius.
        if(temperature > (T_cold_plate + 5)): 
        
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            cv_vapor = df_saturation_table.loc[i, "Vapor Cv (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Conversion into Kelvin
            temperature += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            fill_fraction_var = fill_fraction_calc(pv, pl)
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                          T_cond, r0, h, L_cond, mass_flow_rate_l)
                
                sonic_limit = sonic_limit_Q_total(Cn, N, diameter, hfg, 
                                                  fill_fraction_var, pv, pl, 
                                                  temperature, cv_vapor, evap_fraction)
                
            else:
                sonic_limit = sonic_limit_Q_total(Cn, N, diameter, hfg, 
                                                  fill_fraction_var, pv, pl, 
                                                  temperature, cv_vapor, evap_fraction)
            
            
            sonic_limits_array[i] = sonic_limit
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    sonic_limit_data = np.hstack([temperatures, sonic_limits_array, 
                                     velocities_array, mass_flow_rates_array, 
                                     h_array, evaporation_fractions_array])
    
    df_sonic_limit_data = pd.DataFrame(sonic_limit_data, columns = ('Temperature', 
                                                       r'Sonic Limit (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       r'Heat Transfer Coefficient (w/m2)', 
                                                       'Evaporation Fraction'))
    df_sonic_limit_data.to_excel('Sonic Limit Data {} - CDO.xlsx'.format(name))
    
    ax.plot(saturation_table[:, 1], sonic_limits_array, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)


def plotting_max_heat_load_viscous(file_name, limit_label, color, linestyle, 
                                   alpha, T_cold_plate, fluid, ax):
    
    df_saturation_table = pd.read_csv(file_name)
    df_saturation_table = pd.DataFrame(df_saturation_table)
    saturation_table = np.array(df_saturation_table)
    
    viscous_limits_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    #Arrays Section for outputting data
    velocities_array = np.zeros((len(saturation_table[:, 1]), 1))
    mass_flow_rates_array = np.zeros((len(saturation_table[:, 1]), 1))
    h_array = np.zeros((len(saturation_table[:, 1]), 1))
    evaporation_fractions_array = np.zeros((len(saturation_table[:, 1]), 1))
    
    i = 0
    r0 = diameter / 2
    
    for temperature in saturation_table[:, 1]:
        
        # Comparing for ensuring temperatures in the system are greater than 
        # the condenser temperature. Single conversion back into celsius.
        if(temperature > (T_cold_plate)):
            
            hfg = df_saturation_table.loc[i, "Heat of Vapor (kJ/kg)"] * 1000 # Converting kJ/kg to J/kg
            pv = df_saturation_table.loc[i, "Vapor Density (kg/m3)"]
            pl = df_saturation_table.loc[i, "Liquid Density (kg/m3)"]
            st = df_saturation_table.loc[i, "Surf. Tension (N/m)"]
            kl = df_saturation_table.loc[i, "Liquid Therm. Cond. (W/m-K)"]
            cp_liquid = df_saturation_table.loc[i, "Liquid Cp (kJ/kg-K)"] * 1000 # Converting kJ/kg-K to J/kg-K
            visc_l = df_saturation_table.loc[i, "Liquid Viscosity (Pa-s)"]
            
            # Calculating cold plate temperature from experimental fit
            T_cond = condenser_temperature_calc(temperature, T_cold_plate, fluid)
            
            # Assumption of evaporator and condensor temperatures being 1 degrees
            # above and below the adiabatic temperature, according to Drolen and 
            # Smoot
            T_evap_viscous = temperature + 1
            T_cond_viscous = temperature - 1
            
            # Conversion into Kelvin
            temperature += 273.15
            T_evap_viscous += 273.15
            T_cond_viscous += 273.15
            T_cond += 273.15
            
            u = solver_laminar_flow(st, r0, pv, theta_r, theta_s, visc_l)
            
            mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
            h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
            
            re = reynolds_number_calc(pl, visc_l, diameter, u)
            fill_fraction_var = fill_fraction_calc(pv, pl)
            
            evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                      T_cond, r0, h, L_cond, mass_flow_rate_l)
            
            if(re > 2200):
                
                u = solver_turbulent_flow(st, r0, pv, theta_r, theta_s, visc_l)
                mass_flow_rate_l = mass_flow_rate_calc(pl, u, diameter)
                
                h = h_sens_calculator(kl, u, pl, visc_l, cp_liquid, diameter)
                
                evap_fraction = evaporation_fraction_calc(pl, cp_liquid, pv, hfg, temperature, 
                                                          T_cond, r0, h, L_cond, mass_flow_rate_l)
                
                viscous_limit = viscous_limit_Q_total(r0, N, Cn, L_channel,
                                                      evap_fraction, fill_fraction_var, 
                                                      pv, pl, hfg, temperature, 
                                                      visc_l, T_evap_viscous, 
                                                      T_cond_viscous)
            
            else:
                viscous_limit = viscous_limit_Q_total(r0, N, Cn, L_channel,
                                                      evap_fraction, fill_fraction_var, 
                                                      pv, pl, hfg, temperature, 
                                                      visc_l, T_evap_viscous, 
                                                      T_cond_viscous)
            
            viscous_limits_array[i] = viscous_limit
            velocities_array[i] = u
            mass_flow_rates_array[i] = mass_flow_rate_l
            h_array[i] = h
            evaporation_fractions_array[i] = evap_fraction
            
        i += 1
    
    
    temperatures = saturation_table[:, 1].reshape((len(saturation_table[:, 1]), 1))
    
    viscous_limit_data = np.hstack([temperatures, viscous_limits_array, 
                                     velocities_array, mass_flow_rates_array, 
                                     h_array, evaporation_fractions_array])
    
    df_viscous_limit_data = pd.DataFrame(viscous_limit_data, columns = ('Temperature', 
                                                       'Viscous Limit (w/m2)', 
                                                       'Velocitiy (m/s)', 
                                                       'Mass flow rates (kg/s)', 
                                                       'Heat Transfer Coefficient (w/m2)', 
                                                       'Evaporation Fraction'))
    df_viscous_limit_data.to_excel('Viscous Limit Data {} {} C Cold Plate- CDO.xlsx'.format(name, T_cold_plate))
    
    ax.plot(saturation_table[:, 1], viscous_limits_array, color = color, 
             linestyle = linestyle, label = limit_label, alpha = alpha)


#%% 17.- Defining features for annotations in graphs (matplotlib stuff)

props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.5)


#%% 18.- Additional/Supportive Methods

'''Condenser temperatures calculation based on the experimental results obtained
   with R-134a and R-123. Calculation of a condenser temperature value is based on a 
   polynomial fit of order 3 from the average condenser temperature vs the 
   adiabatic temperature. At temperatures above 65 C for R-134a and 110 for R-123, 
   it is changed for a linear fit, as the few data available makes the polynomial 
   fit fail. Predictions for R-123 at 10 C chiller below 60 C are also replaced
   by the linear fit due to the unreliable data problem described on Diaz-Caraveo et al. 
   (data was discarted)'''

def condenser_temperature_calc(temperature, T_cold_plate, fluid):
    
    if(fluid == 'R-134a' or fluid == 'r-134a'):
        if (T_cold_plate == 10):
            if(temperature <= 60):
                return ((-8.68625172e-05*(temperature**3)) + 
                        (1.18424157e-02*(temperature**2)) + 1.68320952e-01*temperature + 6.91341458)
            else:
                return (0.64709184*temperature + 1.46582525)
        elif (T_cold_plate == 30):
            if(temperature <= 60):
                return ((-2.25136039e-04)*(temperature ** 3) + (3.76745875e-02 * (temperature ** 2)) -
                        (1.44659874*temperature) + 45.0742129)
            else:
                return (0.59167*temperature + 9.64488)
        elif (T_cold_plate == 50):
            if(temperature <= 60):
                return ((-1.64031998e-03)*(temperature ** 3) + (0.275121734 * (temperature ** 2)) -
                        (14.9246462 * temperature) + 311.144008)
            else:
                return (0.41795*temperature + 26.67707447)
                
        else:
            error_message = ("Specified cold plate temperature does not have an " +
                             "available experimental fit for determining condenser " +
                             "temperature (Temperatures must be in Celsius)")
            
            raise Exception(error_message)
    
    elif(fluid == 'R-123' or fluid == 'r-123'):
        if (T_cold_plate == 10):
            #if(temperature <= 80 or temperature >= 100):
            return (0.71833215*temperature + 1.94532823)
        elif (T_cold_plate == 30):
            if(temperature >= 100):
                return (0.66698906*temperature + 10.38928906)
            else:
                return ((-5.67423408e-05)*(temperature ** 3) + (1.22968875e-02 * (temperature ** 2)) +
                        (-1.57393761e-01 * temperature) + 2.71358419e+01)
        elif (T_cold_plate == 50):
            if(temperature >= 100):
                return (0.64014731*temperature + 19.33641085)
            else:
                return ((-6.09366749e-05)*(temperature ** 3) + (1.40791178e-02 * (temperature ** 2)) +
                        (-3.70166280e-01 * temperature) + 4.15785601e+01)
        else:
            error_message = ("Specified cold plate temperature does not have an " +
                             "available experimental fit for determining condenser " +
                             "temperature (Temperatures must be in Celsius)")
            
            raise Exception(error_message)
            
    else:
        if (T_cold_plate == 10): # Average fit between R-123 and R-134a
            return (0.682711995*temperature + 1.70557674)
        else:
            error_message = ("Specified fluid has not been tested at JPL yet nor " + 
                             "experimental fits for the adiabatic and condenser temperature " +
                             "relationship exists")
            
            raise Exception(error_message)


def reynolds_number_calc(p, visc, d, u):
    
    return (p * u * d) / visc

def prandtl_number_calc(cp, visc, k):
    
    return cp * visc / k

def capilarity_number_calc(u, visc, st):
    
    return (u * visc) / st

def mass_flow_rate_calc(p, v, d):
    
    a = math.pi * d**2 / 4
    
    return p * v * a
