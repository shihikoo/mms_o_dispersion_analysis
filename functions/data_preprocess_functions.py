#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 20:08:54 2021

@author: jliao
"""

#import pandas as pd
from scipy import stats
import math
import datetime as datetime
import pandas as pd
import numpy as np
import statistics


MMS_ENERGY_BIN = np.array([1.8099999, 3.5100000,6.7700000, 13.120000, 25.410000, 49.200001, 95.239998, 184.38000,356.97000, 691.10999,1338.0400,2590.4900, 5015.2900, 9709.7900, 18798.590, 32741.160])
MMS_ENERGY_BIN_INT = np.round(MMS_ENERGY_BIN)
np.array([2, 4, 7, 13, 25, 49, 95, 184,357, 691,1338,2590, 5015, 9710, 18799, 32741])
MMS_ENERGY_HIGH = np.array([2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79, 41598.37])
MMS_ENERGY_LOW = np.array([1.27, 2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79])
MMS_DENERGY = np.array([0.5,1,1.91,3.71,7.18,13.92,26.93,52.15,100.97,195.42,378.36,732.53,1418.24,2745.8,5315.99,6713.96])

PA_BIN_SIZE = 11.25
PI = 3.1415926

def remove_outside_data(df):
    new_df = df.loc[df['r'] <= 20,:].reindex()
    return new_df

def read_beam_csv(beam_filenames):
    nbatch = 800
    
    if len(beam_filenames) < nbatch:
        df_beam = read_multiple_csv(beam_filenames[0:nbatch])
        ind = (df_beam['Flux_para'] > 0) | (df_beam['Flux_anti'] > 0)
        df_beam = df_beam.loc[ind,:]
        
    else:
        df_beam1 = read_multiple_csv(beam_filenames[0:nbatch])
        ind = (df_beam1['Flux_para'] > 0) | (df_beam1['Flux_anti'] > 0)
        df_beam1 = df_beam1.loc[ind,:]

        df_beam2 = read_multiple_csv(beam_filenames[800:len(beam_filenames)])
        ind = (df_beam2['Flux_para'] > 0) | (df_beam2['Flux_anti'] > 0)
        df_beam2 = df_beam2.loc[ind,:]

        df_beam = df_beam1.append(df_beam2, ignore_index = True) 
    
    df_beam = df_beam.rename(columns=str.upper)
    
    df_beam['EN_PARA'] = df_beam['EN']
    df_beam['EN_ANTI'] = df_beam['EN']
            
    return df_beam

def read_from_raw_csv(beam_filenames, ext_filenames):
    
    df_beam = read_beam_csv(beam_filenames)
    df_ext = read_external_csv(ext_filenames)
    
    df0 = pd.merge(df_beam, df_ext, on = 'TIME')

    df = preprocess_data(df0)
    return df

def read_external_csv(external_filenames):
    df_ext = read_multiple_csv(external_filenames)
    
    df_ext = df_ext.rename(columns=str.upper)
    df_ext['r'] = (df_ext['GSM_Y']**2 + df_ext['GSM_Z']**2).apply(math.sqrt)

    df_ext = remove_outside_data(df_ext)

    return df_ext

def read_multiple_csv(filenames):    
    li = []
    for filename in filenames:
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)

    df_beam = pd.concat(li, axis=0, ignore_index=True)
    return df_beam

def extract_date(input_datetime_obj):
    date = input_datetime_obj.strftime("%Y-%m-%d")
    return(date)

def extract_time(input_datetime_obj):
    time = input_datetime_obj.strftime("%H-%M-%S")
    return(time)
    
def extract_dispersion_list(mydata, direction_name = 'PARA'):
    estimated_distance_name = 'ESTIMATED_DISTANCE_' + direction_name
    energy_name = 'EN_' + direction_name
    chisq_name = 'DIS_FITTING_CHISQ_' + direction_name
    dof_name = 'DIS_FITTING_DOF_' + direction_name
    rsquare_name = 'DIS_FITTING_RSQUARE_' + direction_name
    n_dispersion_name = 'N_DISPERSION_' + direction_name
    model_field_length_name = 'FLLEN_' + direction_name
    index = mydata.loc[:,estimated_distance_name].notna()
    mydata = mydata.loc[index,:]
    
    dispersion = mydata.groupby([estimated_distance_name,'date',  n_dispersion_name]).agg({'xgse':'count'
                               , chisq_name:'mean', rsquare_name:'mean', dof_name:'mean',  energy_name:'mean', model_field_length_name:'mean'
                               , 'time':'mean', 'xgsm':'mean', 'ygsm':'mean', 'zgsm':'mean', 'MLT':'median', 'L':'mean',  'STORM_PHASE':'max', 'bx':'mean'
                               , 'dist':'mean', 'beta':'mean', 'datetime_str':'min', 'kp':'mean', 'swp':'mean', 'dst':'mean', 'IMF_BY':'mean', 'IMF_BZ':'mean'
                               }).reset_index()
    
    dispersion = dispersion.rename(columns={estimated_distance_name:'estimated_distance', n_dispersion_name:'n_dispersion', 'GSE_X':'dispersion_length',  chisq_name:'chisq', rsquare_name:'rsquare', dof_name:'dof', energy_name:'energy',model_field_length_name:'model_field_line_length_idl'})

    dispersion["direction"] = direction_name

    return(dispersion)

def calculate_pvalue(dispersion):
    p = 1-stats.chi2.cdf(dispersion['chisq'], dispersion['dof'])
    return(p)

def identify_location(onedata):
    if onedata.xgsm < -15:
        region = 'Middle Tail'
    else:
        region = 'Near Earth'
    return(region)

def identify_region(onedata):
    if (onedata['mlt'] >= 8.) & (onedata['mlt'] < 16.):
        region = 'Dayside'
    else:
        if onedata.beta < 0.05:
            region = 'Lobe'
        elif onedata.beta < 1.:
            region = 'BL'
        else:
            region = 'PS'
    return(region)

def estimate_density(onedata):
    Avogadro_constant = 6.02214086e23 # moe-1
    electron_charge = 1.60217662e-19  # coulombs
    mass_o = 15.89
    
    velocity = math.sqrt(2.*onedata.energy*electron_charge/(mass_o/Avogadro_constant/1e3))/1e3
    
    density = onedata.intergrated_flux / velocity / 1e5
    
    return(density)

def closest(lst, K):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

def calculate_B(onedata):
    return(math.sqrt(onedata['bx']**2 + onedata['BY_GSM']**2 + onedata['BZ_GSM']**2))

def calculate_velocity(energy, ion_mass = 15.89):
    Avogadro_constant = 6.02214086e23 # moe-1
    electron_charge = 1.60217662e-19  #coulombs
    velocity = math.sqrt(2.*energy*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3 # 4.577e7*sqrt(data_energy/1e3/AMU) 
    return(velocity)

def find_denergy(energy_int):
    energy_int = closest(MMS_ENERGY_BIN_INT, energy_int)
    
    if energy_int is not None:
        de = MMS_DENERGY[MMS_ENERGY_BIN_INT == energy_int]
        if len(de) == 0:
            print(energy_int)
        return(de[0])
    else:
        return(None)
    
def preprocess_data(data):
    cooked_data = data
    
    cooked_data.rename(columns={'TIME':'time','BX_GSM':'bx','DIST':'dist', 'GSE_X':'xgse', 'GSE_Y':'ygse','GSE_Z':'zgse','GSM_X':'xgsm', 'GSM_Y':'ygsm','GSM_Z':'zgsm', 'BETA':'beta', 'O_V':'velocity_o_all', 'O_VPAR':'v_par_all', 'O_VPERP':'v_perp_all', 'O_N':'density_o_all', 'O_P':'pressure_o_all', 'H_N':'density_h_all', 'H_V':'velocity_h_all', 'KP':'kp', 'DST':'dst','MLT':'mlt'}, inplace = True)

    cooked_data['datetime_str'] = data.loc[:,'time'].apply(datetime.datetime.utcfromtimestamp)
    cooked_data['date'] = data.loc[:,'datetime_str'].apply(extract_date)   
                                
    cooked_data['year'] = cooked_data['datetime_str'].dt.to_period('Y')
    cooked_data['kp_gt_2'] = cooked_data['kp'] > 2 
    cooked_data['storm'] = cooked_data['STORM_PHASE'] > 0

    cooked_data['storm_phase'] = pd.Categorical(cooked_data['STORM_PHASE']).rename_categories({0:'nonstorm',1:'prestorm',2:'main phase',3:'fast recovery', 4:'long recovery'})
    
    cooked_data['region'] = cooked_data.apply(identify_region, axis=1)

    cooked_data['compression_mode'] = (cooked_data['datetime_str'] < pd.Timestamp('2019-4-16')) | (cooked_data['datetime_str'] > pd.Timestamp('2019-8-17'))

    cooked_data['start_time'] = (((cooked_data['datetime_str'].dt.hour/4).apply(int)))*4
    cooked_data['end_time'] = cooked_data['start_time'] + 4
    cooked_data['start_time_dt'] = cooked_data['datetime_str'].apply(datetime.datetime.combine,time=datetime.time.min) + cooked_data['start_time'].apply(pd.Timedelta,unit="h")
    cooked_data['end_time_dt'] = cooked_data['datetime_str'].apply(datetime.datetime.combine,time=datetime.time.min) + cooked_data['end_time'].apply(pd.Timedelta,unit="h")

    cooked_data['o_beam_filepath'] = 'obeam_day/'+cooked_data['start_time_dt'].apply(pd.Timestamp.strftime,format='%Y') +'/o_beam' + cooked_data['start_time_dt'].apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') +'_to_' + cooked_data['end_time_dt'].apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') + '_plasma_condition_short.png'

    index = (cooked_data['dist'] >= 7) & (cooked_data['dist'] < 9)
    cooked_data.loc[index,'dist_region'] = 'near'
    index = cooked_data['dist'] >= 9
    cooked_data.loc[index,'dist_region'] = 'tail'

    index = ((cooked_data['xgsm'] > -1) & (cooked_data['zgsm'] < 0)) | ((cooked_data['xgsm'] < -1) & (cooked_data['bx'] < 0))
    cooked_data.loc[index,'hemi'] = 'south'
    index = ((cooked_data['xgsm'] > -1) & (cooked_data['zgsm'] > 0)) | ((cooked_data['xgsm'] < -1) & (cooked_data['bx'] > 0))
    cooked_data.loc[index,'hemi'] = 'north'

    cooked_data.loc[:, 'flag'] = 0
    index = ((cooked_data['hemi'] == 'south') & (data['FLAG_PARA'] == 1))
    cooked_data.loc[index, 'flag'] = 1
    cooked_data.loc[index, 'flux'] = data['FLUX_PARA']
    cooked_data.loc[index, 'int_flux'] = data['INT_FLUX_PARA']
    cooked_data.loc[index, 'energy'] = data['EN_PARA']
    cooked_data.loc[index, 'eflux'] = data['EFLUX_PARA']
    cooked_data.loc[index, 'imfBy'] = data['IMF_BY_PARA']
    cooked_data.loc[index, 'imfBz'] = data['IMF_BZ_PARA']
    cooked_data.loc[index, 'pa'] = data['PA_PARA']
    cooked_data.loc[index, 'pa_range'] = data['PA_RANGE_PARA']

    cooked_data.loc[index, 'imfBy'] = data['IMF_BY_PARA']
    cooked_data.loc[index, 'imfBz'] = data['IMF_BZ_PARA']
    cooked_data.loc[index, 'swp'] = data['SW_P_PARA']
    cooked_data.loc[index, 'swv'] = data['SW_V_PARA']    

    index = ((cooked_data['hemi'] == 'north') & (data['FLAG_ANTI'] == 1))
    cooked_data.loc[index, 'flag'] = -1
    cooked_data.loc[index, 'flux'] = data['FLUX_ANTI']
    cooked_data.loc[index, 'int_flux'] = data['INT_FLUX_ANTI']
    cooked_data.loc[index, 'energy'] = data['EN_ANTI']
    cooked_data.loc[index, 'eflux'] = data['EFLUX_ANTI']
    cooked_data.loc[index, 'imfBy'] = data['IMF_BY_ANTI']
    cooked_data.loc[index, 'imfBz'] = data['IMF_BZ_ANTI']
    cooked_data.loc[index, 'pa'] = data['PA_ANTI']
    cooked_data.loc[index, 'pa_range'] = data['PA_RANGE_ANTI']

    cooked_data.loc[index, 'imfBy'] = data['IMF_BY_ANTI']
    cooked_data.loc[index, 'imfBz'] = data['IMF_BZ_ANTI']
    cooked_data.loc[index, 'swp'] = data['SW_P_ANTI']
    cooked_data.loc[index, 'swv'] = data['SW_V_ANTI']
    
    cooked_data['energy_int'] = round(cooked_data['energy'])    
    cooked_data['denergy'] = cooked_data['energy_int'].apply(find_denergy)
    cooked_data['intergrated_flux'] =  cooked_data['int_flux'] * cooked_data['denergy']
    
    cooked_data['density_est'] =  cooked_data.apply(estimate_density, axis=1)    
    
    cooked_data['r'] =  (cooked_data['ygsm']**2 + cooked_data['zgsm']**2).apply(math.sqrt)

    cooked_data = remove_outside_data(cooked_data)
    
    cooked_data = cooked_data.sort_values(by=['datetime_str'])
    
    identify_location
    
    return(cooked_data)

def aggregate_angle(df):      
    df = df.loc[(df.loc[:,'pa']).apply(np.isfinite),:]
    
    agg_data = df.groupby(['time','energy']).agg({'xgse':'count' ,'date':'min', 'datetime_str':'min'
                                                  , 'pa':'mean','pa_range':'mean','int_flux':'mean'
                                                  , 'flux':'sum', 'eflux':'sum', 'intergrated_flux':'sum', 'density_est':'sum'
                                                  , 'denergy':'mean','r':'mean'
                                                  , 'xgsm':'min', 'ygsm':'min', 'zgsm':'min'
                                                  , 'ygse':'min', 'zgse':'min', 'mlt':'min', 'L':'min'
                                                  ,  'bx':'min' , 'BY_GSM':'min','BZ_GSM':'min'
                                                  , 'dist':'min', 'beta':'min'
                                                  , 'kp':'min', 'swp':'min', 'swv':'min', 'dst':'min', 'imfBy':'min'
                                                  , 'imfBz':'min', 'storm_phase':'min'
                                                  ,'density_o_all':'mean', 'velocity_o_all':'mean','pressure_o_all':'mean'
                                                  , 'density_h_all':'mean', 'velocity_h_all':'mean'}).reset_index()
    agg_data.rename(columns={'xgse':'nbeam'}, inplace = True)
    
    return(agg_data)

def aggregate_energy(df):      
    
    df = df.loc[(df.loc[:,'energy']).apply(np.isfinite),:]
    
    agg_data = df.groupby(['time']).agg({'nbeam':'sum' ,'date':'min', 'flux':'sum', 'eflux':'sum','int_flux':'mean'
                                         , 'energy':'mean',  'denergy':'mean', 'intergrated_flux':'sum', 'density_est':'sum'
                                         , 'xgsm':'min', 'ygsm':'min', 'zgsm':'min', 'ygse':'min', 'zgse':'min'
                                         , 'pa':'mean', 'pa_range':'mean','r':'mean'
                                         , 'mlt':'min', 'L':'min',  'bx':'min' , 'BY_GSM':'min' 
                                         , 'BZ_GSM':'min', 'dist':'min', 'beta':'min', 'datetime_str':'min'
                                         , 'kp':'min', 'swp':'min','swv':'min', 'dst':'min'
                                         , 'imfBy':'min', 'imfBz':'min', 'storm_phase':'min'
                                         ,'density_o_all':'mean', 'velocity_o_all':'mean','pressure_o_all':'mean'
                                         , 'density_h_all':'mean', 'velocity_h_all':'mean'}).reset_index()
    agg_data['location'] = agg_data.apply(identify_location, axis=1)
    agg_data['region'] = agg_data.apply(identify_region, axis=1)
    agg_data['B'] = agg_data.apply(calculate_B, axis=1)
    
    return(agg_data)

def calculate_density_heatmap(varx, vary, varz, xedges, yedges):
    xstep = abs(xedges[1] - xedges[0])
    ystep = abs(yedges[1] - yedges[0])
    output = np.array([[np.nan for _ in range(len(xedges))] for _ in range(len(yedges))])
    for ix in range(len(xedges)):
        for iy in range(len(yedges)):
            ixedge = xedges[ix]
            iyedge = yedges[iy]
            ind = (varx > ixedge ) & (varx <= (ixedge+xstep)) & (vary > iyedge) & (vary <= (iyedge+ystep))
            if sum(ind) == 0:
                continue
            output[ix,iy] = statistics.mean(varz[ind])
    output = output.T

    return output     

def calculate_occurrence_rate(varx, vary, xedges, yedges):
    xstep = abs(xedges[1] - xedges[0])
    ystep = abs(yedges[1] - yedges[0])
    output = np.array([[np.nan for _ in range(len(xedges))] for _ in range(len(yedges))])
    for ix in range(len(xedges)):
        for iy in range(len(yedges)):
            ixedge = xedges[ix]
            iyedge = yedges[iy]
            ind = (varx > ixedge ) & (varx <= (ixedge+xstep)) & (vary > iyedge) & (vary <= (iyedge+ystep))
            output[ix,iy] = sum(ind)
            
    output = output.T
    return output    

def preprocess_dispersion_list(dispersion_list):
    dispersion_list['p_value'] = dispersion_list.apply(calculate_pvalue,axis = 1)
    dispersion_list['region'] = dispersion_list.apply(identify_region, axis=1)

    dispersion_list['index'] = dispersion_list.index

    dispersion_list.loc[ (dispersion_list['GSM_Z'] < 0), 'location'] = 'south'
    dispersion_list.loc[ (dispersion_list['GSM_Z'] > 0), 'location'] = 'north'

    dispersion_list.loc[((dispersion_list['GSM_X'] < -1) & (((dispersion_list['direction'] == 'PARA') & (dispersion_list['BX_GSM'] > 0)) | ((dispersion_list['direction'] == 'ANTI') & (dispersion_list['BX_GSM'] < 0)))) | ((dispersion_list['GSM_X'] > -1) & (((dispersion_list['direction'] == 'ANTI') & (dispersion_list['GSM_Z'] < 0)) | ((dispersion_list['direction'] == 'PARA') & (dispersion_list['GSM_Z'] > 0)))), 'direction_et'] = 'earthward'

    dispersion_list.loc[((dispersion_list['GSM_X'] < -1) & (((dispersion_list['direction'] == 'ANTI') & (dispersion_list['BX_GSM'] > 0)) | ((dispersion_list['direction'] == 'PARA') & (dispersion_list['BX_GSM'] < 0)))) | ((dispersion_list['GSM_X'] > -1) & (((dispersion_list['direction'] == 'PARA') & (dispersion_list['GSM_Z'] < 0)) | ((dispersion_list['direction'] == 'ANTI') & (dispersion_list['GSM_Z'] > 0)))) , 'direction_et'] = 'outward'

    dispersion_list['dispersion_time'] = 2. * (dispersion_list['dof']+2)

    return(dispersion_list)

def extract_dispersions(data, save_to_filename = 'output/dispersion_list.csv'):
    dispersion_para = extract_dispersion_list(data, direction_name = 'PARA')
    dispersion_anti = extract_dispersion_list(data, direction_name = 'ANTI')
    dispersion_list = pd.concat([dispersion_para,dispersion_anti],ignore_index = True)  
    
    dispersion_list = preprocess_dispersion_list(dispersion_list)
    
    dispersion_list.to_csv(save_to_filename)
    
    return(dispersion_list)

