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
MMS_DENERGY = np.array([0.5,1,1.91,3.71,7.18,13.92,26.93,52.15,100.97,195.42,378.36,732.53,1418.24,2745.8,5315.99,6713.96])

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
   
def calcualte_velocity_from_energy(energy, ion = 'O+'):   
    if ion == 'O+':
        ion_mass = O_MASS    
    velocity = math.sqrt((energy*ELECTRON_CHARGE/(ion_mass/AVOGADRO_CONSTANT/1000)*2.))/1000.
    return(velocity)
    
def calculate_flux(itime, energy, result, distance, avg_time):
    st = itime*avg_time
    et = (itime+1)*avg_time
    index = (result['time'] >= st) & (result['time'] < et) & (result['distance_Re'] == distance) & (result['energy_binned'] == energy)
    return(len(result.loc[index,'energy_binned']))

def find_denergy(energy):
    index = MMS_ENERGY_BIN == energy
    return(MMS_DENERGY[index])

def calculate_energy_spectra_mean(result, target_distance, avg_time):
    index = (result['distance_Re'] == target_distance)
    ntime = round(max(result.loc[index,'time'])/avg_time)
    time_loop = np.arange(0,ntime,1)
    energy_spectra = pd.DataFrame(data = np.empty((ntime,3)), columns =['time','energy','flux'])

    for itime in time_loop:
        st = itime * avg_time
        et = (itime+1) * avg_time
        index = (result['time'] >= st) & (result['time'] < et) & (result['distance_Re'] == target_distance)
        energy_spectra.loc[itime,'flux'] = len(result.loc[index,'energy_binned'])
        if energy_spectra.loc[itime,'flux'] > 0:
            energy_spectra.loc[itime,'energy'] = statistics.mean(result.loc[index,'energy_binned'])
        energy_spectra.loc[itime,'time'] = (st + et)/2
    
    energy_spectra['energy_binned'] = find_energy_bin_df(energy_spectra['energy'])
    return(energy_spectra)


def calculate_energy_spectra(result, target_distance, avg_time):
    index = (result['distance_Re'] == target_distance)
    ntime = round(max(result.loc[index,'time'])/avg_time)
    time_loop = np.arange(0,ntime,1)
    energy_spectra = pd.DataFrame(np.array(np.meshgrid(time_loop,MMS_ENERGY_BIN)).T.reshape(-1,2), columns=["itime","energy_binned"])

    energy_spectra['time'] = (energy_spectra['itime']+0.5)*avg_time
    energy_spectra['flux'] = energy_spectra.apply(
        lambda x: calculate_flux(x['itime'], x['energy_binned'], result, target_distance, avg_time), axis=1)
    
    index_valid = (energy_spectra['flux'] > 0)
    return(energy_spectra.loc[index_valid,:])

def identify_beam(energy_spectra):
    ntime = len(np.unique(energy_spectra['time']))
    beam = pd.DataFrame(np.empty((ntime,2)), columns=['time','energy_binned'])
    i = 0
    for itime in np.unique(energy_spectra['time']):
        index = energy_spectra['time'] == itime
        beam.loc[i,'time'] = itime
        if energy_spectra.loc[index,'flux'].max() == 0:
            beam.loc[i, 'energy_binned'] = None
        else:
            beam.loc[i, 'energy_binned'] = energy_spectra.loc[(energy_spectra.loc[index,'flux']).idxmax(),'energy_binned']
        i = i + 1   
    
    return(beam)
       
def identify_dispersion(beam):
    dispersion = beam.copy()
    dispersion['decrease'] = (-dispersion['energy_binned'].diff()) > 0
    dispersion['equal'] = (-dispersion['energy_binned'].diff()) ==  0
    dispersion['energy_binned'] = None
#    dispersion['ndispersion'] = None
    
    index = dispersion.loc[:,'decrease']
    ndisperison = 0
    start_index = -1
    end_index = -1
    length = 0
    nequal = 0
    for iindex in dispersion['time'].index:
        if dispersion.loc[iindex,'decrease'] == True:
            if length == 0:
                length = 2
            else:
                length = length + 1
            if start_index == -1:
                start_index = iindex - 1
            nequal = 0
        else:
            if (dispersion.loc[iindex,'equal'] == True) & (nequal < 2):
                if length == 0:
                    length = 2
                else:
                    length = length + 1
                if start_index == -1:
                    start_index = iindex - 1                 
                if nequal == 0:
                    nequal = 2
                else:
                    nequal = nequal + 1
            else:      
                if length > 4:
                    ndisperison = ndisperison + 1
                    end_index = iindex - 1
                    dispersion.loc[start_index:end_index, 'dispersion_length'] = length
                    dispersion.loc[start_index:end_index, 'energy_binned'] = beam.loc[start_index:end_index, 'energy_binned']
                    dispersion.loc[start_index:end_index, 'denergy'] = dispersion.loc[start_index:end_index, 'energy_binned'].apply(find_denergy)
                    dispersion.loc[start_index:end_index, 'ndisperison'] = ndisperison
                    
                    dispersion.loc[start_index:end_index, 'inverse_v'] = 1/(dispersion.loc[start_index:end_index, 'energy_binned'].apply(calcualte_velocity_from_energy))
                    
                    dispersion.loc[start_index:end_index, 'd_inverse_v'] = dispersion.loc[start_index:end_index, 'inverse_v'] -  1/(dispersion.loc[start_index:end_index, 'energy_binned'] + dispersion.loc[start_index:end_index, 'denergy']).apply(calcualte_velocity_from_energy)
                    
                    start_index = -1
                    end_index = -1
                    length = 0
                    nequal = 0
                else:
                    length = 0
                    start_index = -1
                    end_index = -1
                    nequal = 0
    
    if start_index != -1 & length > 4:
        ndisperison = ndisperison + 1
        end_index = iindex
        dispersion.loc[start_index:end_index, 'dispersion_length'] = length

        dispersion.loc[start_index:end_index, 'energy_binned'] = beam.loc[start_index:end_index, 'energy_binned']
        dispersion.loc[start_index:end_index, 'denergy'] = dispersion.loc[start_index:end_index, 'energy_binned'].apply(find_denergy)
        dispersion.loc[start_index:end_index, 'ndisperison'] = ndisperison

        dispersion.loc[start_index:end_index, 'inverse_v'] = 1/(dispersion.loc[start_index:end_index, 'energy_binned'].apply(calcualte_velocity_from_energy))

        dispersion.loc[start_index:end_index, 'd_inverse_v'] = dispersion.loc[start_index:end_index, 'inverse_v'] -  1/(dispersion.loc[start_index:end_index, 'energy_binned'] + dispersion.loc[start_index:end_index, 'denergy']).apply(calcualte_velocity_from_energy)    
    return(dispersion)

                        