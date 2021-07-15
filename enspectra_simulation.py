#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 14:40:33 2021

@author: jliao
"""


MMS_ENERGY_BIN = np.array([1.8099999, 3.5100000,6.7700000, 13.120000, 25.410000, 49.200001, 95.239998, 184.38000,356.97000, 691.10999,1338.0400,2590.4900, 5015.2900, 9709.7900, 18798.590, 32741.160])

AVOGADRO_CONSTANT = 6.0221410e+23 # 

O_MASS = 15.89 #
ELECTRON_CHARGE = 1.6021766e-19 # eV

EARTH_RADIUS = 6371. #KM

def find_nearest_value(value, np_array):
    value_binned = np_array[(np.abs(np_array - value)).argmin()]
    return(value_binned)
    
    
def convert_Re_to_km(r_re):
    EARTH_RADIUS = 6371.
    r_km = r_re * EARTH_RADIUS
    return(r_km)
    
    
def bin_energy_in_mms:
    
    
    
    