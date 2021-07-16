#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 14:40:33 2021

@author: jliao
"""

import numpy as np
import pandas as pd
import math 
import statistics
import itertools as it

MMS_ENERGY_BIN = np.array([1.8099999, 3.5100000,6.7700000, 13.120000, 25.410000, 49.200001, 95.239998, 184.38000,356.97000, 691.10999,1338.0400,2590.4900, 5015.2900, 9709.7900, 18798.590, 32741.160])
MMS_ENERGY_HIGH = np.array([2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79, 41598.37])
MMS_ENERGY_LOW = np.array([1.27, 2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79])
O_MASS = 15.89 #u
AVOGADRO_CONSTANT = 6.0221410e+23 # mol-1
ELECTRON_CHARGE = 1.6021766e-19 # eV
EARTH_RADIUS = 6371. #km

def find_nearest_value(value, np_array):
    index = (np.abs(np_array - value)).argmin()
    output_bin = np_array[index]
    return(output_bin)

def find_value_bin(value, np_array, high_boundary, low_boundary): 
    index = (value > low_boundary) & (value <= high_boundary)
    output_bin = None
    if len(np_array[index]) > 0:
        output_bin = np_array[index]
    return(output_bin)
    
def find_nearest_energy_bin(value):
    return(find_nearest_value(value, MMS_ENERGY_BIN))
   
def find_energy_bin(value):
    return(find_value_bin(value, MMS_ENERGY_BIN, MMS_ENERGY_HIGH, MMS_ENERGY_LOW))
    
def find_energy_bin_df(energy_df):
    energy_binned = (np.apply_along_axis(find_energy_bin, 0, [energy_df]))
    if energy_binned.shape[0] == 1:
        energy_binned_df = energy_binned[0]
    else:
        energy_binned_df = energy_binned
    
    return(energy_binned_df)
    
def convert_Re_to_km(r_re):
    r_km = r_re * EARTH_RADIUS
    return(r_km)
         
def calcualte_energy_from_velocity(velocity, ion = 'O+'):   
    if ion == 'O+':
        ion_mass = O_MASS
        
    energy = np.array(pow(velocity*1000.,2)/2./ELECTRON_CHARGE*(ion_mass/AVOGADRO_CONSTANT/1000))
    return(energy)
    
def calculate_flux(itime, energy, result, distance, avg_time):
    st = itime*avg_time
    et = (itime+1)*avg_time
    index = (result['time'] >= st) & (result['time'] < et) & (result['distance_Re'] == distance) & (result['energy_binned'] == energy)
    
    return(len(result.loc[index,'energy_binned']))
    