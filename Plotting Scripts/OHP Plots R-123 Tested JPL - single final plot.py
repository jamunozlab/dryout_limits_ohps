# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:36:47 2022

@author: Cesar Diaz-Caraveo
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import math

import OHP_Operational_Limits_Calculator_Methods_Constant_Cold_Plate_Temp as Methods

#%% Parameters input section 
# The following are constants that must be introduced by the users depending 
# on the specifications of their OHP device

# Geometrical Features

# Number of channels present on the evaporator
Methods.N = 42
# Number of condensors on the OHP
Methods.Cn = 1
# Channel Diameter
Methods.diameter = 0.001 # m
# Length per channel
Methods.L_channel = 0.184 # m, 184 mm
# Condensor section length per channel
Methods.L_cond = 0.08675 # m, 86.75 mm
# Heater section length per channel
Methods.L_heater = 0.0762 # m, 76.2 mm


# Fill fraction, volume, and mass

# Total OHP channel volume
Methods.V_total_OHP = -1
# Total working fluid mass
Methods.mass_total = -1
# Fill Fraction (decimal)
Methods.fill_fraction = 0.501
# Average density
Methods.avg_density = 757.69


# Contact angles

# Static contact angle
Methods.theta_s = 63.19 * (math.pi/180) # 45 deg radians
# Receding contact angle 
Methods.theta_r = 0.0


# Other manual imported parameters

# Name of OHP or Test
Methods.name = "R-123 Tested at JPL Summer 2022"
# Calculating total cross sectional area
total_area = (math.pi*(Methods.diameter ** 2) / 4) * Methods.N #/ 2


#%% Props definition section. Matplotlib stuff

props = dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=1)  
# Props for bbox parameter in the text box of the graphs

#%% Data definition section

adiabatic_temperature_10_C = 113 # C
dry_out_power_10_C = (560 - 19.34) # W
power_density_10_C = dry_out_power_10_C / (total_area)

adiabatic_temperature_10_C_cond_dryout = 79.9
dry_out_power_10_C_cond_dryout = (360 - 10.37) # W
power_density_10_C_cond_dryout = dry_out_power_10_C_cond_dryout / (total_area)


adiabatic_temperature_30_C = 121 # C
dry_out_power_30_C = (500 - 19.22) # W
power_density_30_C = dry_out_power_30_C / total_area

adiabatic_temperature_30_C_cond_dryout = 88 # C
dry_out_power_30_C_cond_dryout = (310 - 10.45) # W
power_density_30_C_cond_dryout = dry_out_power_30_C_cond_dryout / total_area


adiabatic_temperature_50_C = 125
dry_out_power_50_C = (410 - 18.86) # W
power_density_50_C = dry_out_power_50_C / total_area

adiabatic_temperature_50_C_cond_dryout = 100.9 # C
dry_out_power_50_C_cond_dryout = (280 - 11.65) # W
power_density_50_C_cond_dryout = dry_out_power_50_C_cond_dryout / total_area


dryout_limit_10_C = 266.81
dryout_limit_30_C = 218.32
dryout_limit_50_C = 176.88

#%% Plots of all temperatures together - 10 C, 30 C, and 50 C - heat flux

Methods.T_cold_plate = 10

'''region_plot, ax2 = plt.subplots(figsize = (11, 7))

Methods.plotting_bond_number_limit(file_name = 'R-123 Saturation Table.csv', 
                                   limit_label = 'Bond Number Limit', color = 'black',
                                   linestyle = '--', alpha = 1, ax = ax2)
Methods.plotting_max_heat_flux_inertial(file_name = 'R-123 Saturation Table.csv', 
                                limit_label = 'Vapor Inertia Limit', color = 'red', 
                                linestyle = '-', theta_r = Methods.theta_r, 
                                theta_s = Methods.theta_s, alpha = 1, 
                                T_cold_plate = Methods.T_cold_plate, fluid = 'R-123', ax = ax2)
Methods.plotting_max_heat_flux_viscous(file_name = 'R-123 Saturation Table.csv', 
                               limit_label = 'Viscous Limit', linestyle = '--', 
                               color = 'gray', alpha = 1, T_cold_plate = Methods.T_cold_plate, 
                               fluid = 'R-123', ax = ax2)

ax2.scatter(adiabatic_temperature_10_C, power_density_10_C, label = '10 C Chiller R-123 Dryout Point', 
            color = 'tab:blue', marker = '^')
ax2.scatter(adiabatic_temperature_10_C_cond_dryout, power_density_10_C_cond_dryout,
            label = '10 C Chiller R-123 Insipient Dryout Point', color = 'tab:blue', marker = 'o')
ax2.scatter(adiabatic_temperature_30_C, power_density_30_C, 
            label = '30 C Chiller R-123 Dryout Point', color = 'tab:purple', marker = '^')
ax2.scatter(adiabatic_temperature_30_C_cond_dryout, power_density_30_C_cond_dryout,
            label = '30 C Chiller R-123 Insipient Dryout Point', color = 'tab:purple', marker = 'o')
ax2.scatter(adiabatic_temperature_50_C, power_density_50_C, label = '50 C Chiller R-123 Dryout Point', 
            color = 'tab:red', marker = '^')
ax2.scatter(adiabatic_temperature_50_C_cond_dryout, power_density_50_C_cond_dryout,
            label = '50 C Chiller R-123 Insipient Dryout Point', color = 'tab:red', marker = 'o')

ax2.set_title("Heat Flux Limits for R-123 Al10SiMg OHP Tested at JPL")
ax2.set_xlabel("Adiabatic Temperature (°C)")
ax2.set_ylabel(r"Channel Heat Flux (W/$m^{2}$)")
ax2.set_xlim([0, 170])
#ax2.set_ylim([1, 2e7])
ax2.set_yscale('log')
ax2.grid()

ax2.legend()
plt.show()'''


#%% Plots of all temperatures together - 10 C, 30 C, and 50 C - heat load

Methods.T_cold_plate = 10
region_plot, ax2 = plt.subplots(figsize = (11, 7))

Methods.plotting_bond_number_limit(file_name = 'R-123 Saturation Table.csv', 
                                   limit_label = 'Bond Number Limit', color = 'black',
                                   linestyle = '--', alpha = 1, ax = ax2)
Methods.plotting_max_heat_load_inertial(file_name = 'R-123 Saturation Table.csv', 
                                limit_label = 'Vapor Inertia Limit', color = 'red', 
                                linestyle = '-', theta_r = Methods.theta_r, 
                                theta_s = Methods.theta_s, alpha = 1, 
                                T_cold_plate = Methods.T_cold_plate, fluid = 'R-123', ax = ax2)
Methods.plotting_max_heat_load_viscous(file_name = 'R-123 Saturation Table.csv', 
                               limit_label = 'Viscous Limit', linestyle = '--', 
                               color = 'gray', alpha = 1, T_cold_plate = Methods.T_cold_plate, 
                               fluid = 'R-123', ax = ax2)


ax2.scatter(adiabatic_temperature_10_C_cond_dryout, dry_out_power_10_C_cond_dryout,
            label = '10 C Cold Plate Incipient Dryout Point', color = 'tab:blue', marker = 'o')
ax2.scatter(adiabatic_temperature_30_C_cond_dryout, dry_out_power_30_C_cond_dryout,
            label = '30 C Cold Plate Incipient Dryout Point', color = 'tab:purple', marker = 'o')
ax2.scatter(adiabatic_temperature_50_C_cond_dryout, dry_out_power_50_C_cond_dryout,
            label = '50 C Cold Plate Incipient Dryout Point', color = 'tab:red', marker = 'o')
ax2.scatter(adiabatic_temperature_10_C, dry_out_power_10_C, label = '10 C Cold Plate Dryout Point', 
            color = 'tab:blue', marker = '^')
ax2.scatter(adiabatic_temperature_30_C, dry_out_power_30_C, 
            label = '30 C Cold Plate Dryout Point', color = 'tab:purple', marker = '^')
ax2.scatter(adiabatic_temperature_50_C, dry_out_power_50_C, label = '50 C Cold Plate Dryout Point', 
            color = 'tab:red', marker = '^')

# Matplotlib fluid anotation in the corner
ax2.text(0.92, 0.96, 'R-123', fontsize = 13, transform=ax2.transAxes, verticalalignment='top', bbox=props)


#ax2.set_title("Heat Flux Limits for R-123 Al10SiMg OHP Tested at JPL")
ax2.set_xlabel("Adiabatic Temperature (°C)", fontsize = 13)
ax2.set_ylabel(r"System Heat Load (W))", fontsize = 13)
#ax2.set_title("Heat Load Limits for R-123 Al10SiMg OHP Tested at JPL")
ax2.set_xlim([0, 190])
ax2.set_ylim([0, 700])
ax2.set_yscale('linear')
ax2.grid()
ax2.tick_params(labelsize=13)

ax2.legend(loc='upper left', fontsize = 12)
plt.savefig("R_123_predictions_experimental_data_JPL.jpg", dpi = 700, bbox_inches = 'tight')

plt.show()

#%% Plots different temperatures together

region_plot, ax2 = plt.subplots(figsize = (11, 8))
Methods.T_cold_plate = 10

Methods.plotting_bond_number_limit(file_name = 'R-123 Saturation Table.csv', 
                                   limit_label = 'Bond Number Limit', color = 'black',
                                   linestyle = '--', alpha = 1, ax = ax2)
Methods.plotting_max_heat_load_inertial(file_name = 'R-123 Saturation Table.csv', 
                                limit_label = 'Vapor Inertia Limit 10 C', color = 'red', 
                                linestyle = '-', theta_r = Methods.theta_r, 
                                theta_s = Methods.theta_s, alpha = 1, 
                                T_cold_plate = Methods.T_cold_plate, fluid = 'R-123', ax = ax2)
Methods.plotting_max_heat_load_viscous(file_name = 'R-123 Saturation Table.csv', 
                               limit_label = 'Viscous Limit 10 C', linestyle = '-', 
                               color = 'gray', alpha = 1, T_cold_plate = Methods.T_cold_plate, 
                               fluid = 'R-123', ax = ax2)

Methods.T_cold_plate = 30

'''Methods.plotting_bond_number_limit(file_name = 'R-123 Saturation Table.csv', 
                                   limit_label = 'Bond Number Limit', color = 'black',
                                   linestyle = '--', alpha = 1, ax = ax2)'''
Methods.plotting_max_heat_load_inertial(file_name = 'R-123 Saturation Table.csv', 
                                limit_label = 'Vapor Inertia Limit 30 C', color = 'red', 
                                linestyle = '--', theta_r = Methods.theta_r, 
                                theta_s = Methods.theta_s, alpha = 1, 
                                T_cold_plate = Methods.T_cold_plate, fluid = 'R-123', ax = ax2)
Methods.plotting_max_heat_load_viscous(file_name = 'R-123 Saturation Table.csv', 
                               limit_label = 'Viscous Limit 30 C', linestyle = '--', 
                               color = 'gray', alpha = 1, T_cold_plate = Methods.T_cold_plate, 
                               fluid = 'R-123', ax = ax2)

Methods.T_cold_plate = 50

'''Methods.plotting_bond_number_limit(file_name = 'R-123 Saturation Table.csv', 
                                   limit_label = 'Bond Number Limit', color = 'black',
                                   linestyle = '--', alpha = 1, ax = ax2)'''
Methods.plotting_max_heat_load_inertial(file_name = 'R-123 Saturation Table.csv', 
                                limit_label = 'Vapor Inertia Limit 50 C', color = 'red', 
                                linestyle = ':', theta_r = Methods.theta_r, 
                                theta_s = Methods.theta_s, alpha = 1, 
                                T_cold_plate = Methods.T_cold_plate, fluid = 'R-123', ax = ax2)
Methods.plotting_max_heat_load_viscous(file_name = 'R-123 Saturation Table.csv', 
                               limit_label = 'Viscous Limit 50 C', linestyle = ':', 
                               color = 'gray', alpha = 1, T_cold_plate = Methods.T_cold_plate, 
                               fluid = 'R-123', ax = ax2)


ax2.scatter(adiabatic_temperature_10_C_cond_dryout, dry_out_power_10_C_cond_dryout,
            label = 'Incipient Dryout Point 10 C Cold Plate', color = 'tab:blue', marker = 'o')
ax2.scatter(adiabatic_temperature_30_C_cond_dryout, dry_out_power_30_C_cond_dryout,
            label = 'Incipient Dryout Point 30 C Cold Plate', color = 'tab:purple', marker = 'o')
ax2.scatter(adiabatic_temperature_50_C_cond_dryout, dry_out_power_50_C_cond_dryout,
            label = 'Incipient Dryout Point 50 C Cold Plate', color = 'tab:red', marker = 'o')
ax2.scatter(adiabatic_temperature_10_C, dry_out_power_10_C, label = 'Final Dryout Point 10 C Cold Plate', 
            color = 'tab:blue', marker = '^')
ax2.scatter(adiabatic_temperature_30_C, dry_out_power_30_C, 
            label = 'Final Dryout Point 30 C Cold Plate', color = 'tab:purple', marker = '^')
ax2.scatter(adiabatic_temperature_50_C, dry_out_power_50_C, label = 'Final Dryout Point 50 C Cold Plate', 
            color = 'tab:red', marker = '^')

ax2.scatter(adiabatic_temperature_10_C_cond_dryout, dryout_limit_10_C,
            label = 'Prediction Incipient Dryout Point 10 C', marker = 'o', 
            facecolors='none', edgecolors='tab:blue')
ax2.scatter(adiabatic_temperature_30_C_cond_dryout, dryout_limit_30_C,
            label = 'Prediction Incipient Dryout Point 30 C', marker = 'o', 
            facecolors='none', edgecolors='tab:purple')
ax2.scatter(adiabatic_temperature_50_C_cond_dryout, dryout_limit_50_C,
            label = 'Prediction Incipient Dryout Point 50 C', marker = 'o', 
            facecolors='none', edgecolors='tab:red')

# Matplotlib fluid anotation in the corner
ax2.text(0.92, 0.96, 'R-123', fontsize = 13, transform=ax2.transAxes, verticalalignment='top', bbox=props)


#ax2.set_title("Heat Flux Limits for R-123 Al10SiMg OHP Tested at JPL")
ax2.set_xlabel("Adiabatic Temperature (°C)", fontsize = 13)
ax2.set_ylabel(r"System Heat Load (W)", fontsize = 13)
#ax2.set_title("Heat Load Limits for R-123 Al10SiMg OHP Tested at JPL")
ax2.set_xlim([0, 190])
ax2.set_ylim([0, 900])
ax2.set_yscale('linear')
ax2.grid()
ax2.tick_params(labelsize=13)

ax2.legend(loc='upper left', fontsize = 11, ncol = 2)
plt.savefig("R_123_predictions_experimental_data_JPL.jpg", dpi = 700, bbox_inches = 'tight')

plt.show()
