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

MMS_ENERGY_BIN = np.array([1.8099999, 3.5100000,6.7700000, 13.120000, 25.410000, 49.200001, 95.239998, 184.38000,356.97000, 691.10999,1338.0400,2590.4900, 5015.2900, 9709.7900, 18798.590, 32741.160])
MMS_ENERGY_BIN_INT = np.round(MMS_ENERGY_BIN)
np.array([2, 4, 7, 13, 25, 49, 95, 184,357, 691,1338,2590, 5015, 9710, 18799, 32741])
MMS_ENERGY_HIGH = np.array([2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79, 41598.37])
MMS_ENERGY_LOW = np.array([1.27, 2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79])
MMS_DENERGY = np.array([0.5,1,1.91,3.71,7.18,13.92,26.93,52.15,100.97,195.42,378.36,732.53,1418.24,2745.8,5315.99,6713.96])

PA_BIN_SIZE = 22.5
PI = 3.1415926


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
       
    dispersion = mydata.groupby([estimated_distance_name,'date',  n_dispersion_name]).agg({'GSE_X':'count'
                               , chisq_name:'mean', rsquare_name:'mean', dof_name:'mean',  energy_name:'mean', model_field_length_name:'mean'
                               , 'TIME':'mean', 'GSM_X':'mean', 'GSM_Y':'mean', 'GSM_Z':'mean', 'MLT':'median', 'L':'mean',  'STORM_PHASE':'max', 'BX_GSM':'mean'
                               , 'DIST':'mean', 'BETA':'mean', 'datetime_str':'min', 'KP':'mean', 'SW_P':'mean', 'DST':'mean', 'IMF_BY':'mean', 'IMF_BZ':'mean'
                               }).reset_index()
    
    dispersion = dispersion.rename(columns={estimated_distance_name:'estimated_distance', n_dispersion_name:'n_dispersion', 'GSE_X':'dispersion_length',  chisq_name:'chisq', rsquare_name:'rsquare', dof_name:'dof', energy_name:'energy',model_field_length_name:'model_field_line_length_idl'})

    dispersion["direction"] = direction_name

    return(dispersion)

def calculate_cooked_data(dispersion):
    p = 1-stats.chi2.cdf(dispersion['chisq'], dispersion['dof'])
    return(p)
    
def identify_region(onedata):
    if (onedata['MLT'] >= 8.) & (onedata['MLT'] < 16.):
        region = 'Dayside'
    else:
        if onedata.BETA < 0.05:
            region = 'Lobe'
        elif onedata.BETA < 1.:
            region = 'BL'
        else:
            region = 'PS'
            
    return(region)

def closest(lst, K):
     
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

def calculate_B(onedata):
    return(math.sqrt(onedata['BX_GSM']**2 + onedata['BY_GSM']**2 + onedata['BZ_GSM']**2))

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
    cooked_data['datetime_str'] = data.loc[:,'TIME'].apply(datetime.datetime.utcfromtimestamp)
    cooked_data['date'] = data.loc[:,'datetime_str'].apply(extract_date)   
    cooked_data['year'] = cooked_data['datetime_str'].dt.to_period('Y')   
    cooked_data['pitch angle'] = pd.concat([cooked_data['PA_PARA'] ,cooked_data['PA_ANTI']],ignore_index = True)

    cooked_data['bx'] = pd.concat([cooked_data['BX_GSM'] ,cooked_data['BX_GSM']],ignore_index = True)
    
    cooked_data['dist'] = pd.concat([data['DIST'] ,data['DIST']],ignore_index = True)
    cooked_data['xgsm'] = pd.concat([data['GSM_X'] ,data['GSM_X']],ignore_index = True)
    cooked_data['ygsm'] = pd.concat([data['GSM_Y'] ,data['GSM_Y']],ignore_index = True)
    cooked_data['zgsm'] = pd.concat([data['GSM_Z'] ,data['GSM_Z']],ignore_index = True)
    cooked_data['beta'] = pd.concat([data['BETA'] ,data['BETA']],ignore_index = True)
    cooked_data['imfBy'] = pd.concat([data['IMF_BY_PARA'] ,data['IMF_BY_ANTI']],ignore_index = True)
    cooked_data['imfBz'] = pd.concat([data['IMF_BZ_PARA'] ,data['IMF_BZ_ANTI']],ignore_index = True)
    cooked_data['storm_phase'] = pd.concat([data['STORM_PHASE'] ,data['STORM_PHASE']],ignore_index = True)
    cooked_data['velocity_o_all'] = pd.concat([data['O_V'] ,data['O_V']],ignore_index = True)
    cooked_data['v_par_all'] = pd.concat([data['O_VPAR'] ,data['O_VPAR']],ignore_index = True)
    cooked_data['v_perp_all'] = pd.concat([data['O_VPERP'] ,data['O_VPERP']],ignore_index = True)
    cooked_data['density_o_all'] = pd.concat([data['O_N'] ,data['O_N']],ignore_index = True)
    cooked_data['pressure_o_all'] = pd.concat([data['O_P'] ,data['O_P']],ignore_index = True)
    cooked_data['density_h_all'] = pd.concat([data['H_N'] ,data['H_N']],ignore_index = True)
    cooked_data['swp'] = pd.concat([data['SW_P_PARA'] ,data['SW_P_ANTI']],ignore_index = True)
    cooked_data['swv'] = pd.concat([data['SW_V_PARA'] ,data['SW_V_ANTI']],ignore_index = True)
    cooked_data['kp'] = pd.concat([data['KP'] ,data['KP']],ignore_index = True)
    cooked_data['dst'] = pd.concat([data['DST'] ,data['DST']],ignore_index = True)
    cooked_data['year'] = cooked_data['datetime_str'].dt.to_period('Y')
    
    cooked_data['storm'] = cooked_data['STORM_PHASE'] > 0
    cooked_data['kp_gt_2'] = cooked_data['KP'] > 2 
    cooked_data['storm_phase'] = pd.Categorical(cooked_data['STORM_PHASE']).rename_categories({0:'nonstorm',1:'prestorm',2:'main phase',3:'fast recovery', 4:'long recovery'})
    
    cooked_data['region'] = cooked_data.apply(identify_region, axis=1)

#    cooked_data.loc[cooked_data['BETA'] < 0.05,'region'] = 'lobe'
#    cooked_data.loc[(cooked_data['BETA'] < 1) & (cooked_data['BETA'] >= 0.05),'region'] = 'bl'
#    cooked_data.loc[cooked_data['BETA'] >= 1,'region'] = 'ps'

    cooked_data['compression_mode'] = (cooked_data['datetime_str'] < pd.Timestamp('2019-4-16')) | (cooked_data['datetime_str'] > pd.Timestamp('2019-8-17'))

    cooked_data['start_time'] = (((cooked_data['datetime_str'].dt.hour/4).apply(int)))*4
    cooked_data['end_time'] = cooked_data['start_time'] + 4
    cooked_data['start_time_dt'] = cooked_data['datetime_str'].apply(datetime.datetime.combine,time=datetime.time.min) + cooked_data['start_time'].apply(pd.Timedelta,unit="h")
    cooked_data['end_time_dt'] = cooked_data['datetime_str'].apply(datetime.datetime.combine,time=datetime.time.min) + cooked_data['end_time'].apply(pd.Timedelta,unit="h")

    cooked_data['o_beam_filepath'] = 'obeam_day/'+cooked_data['start_time_dt'].apply(pd.Timestamp.strftime,format='%Y') +'/o_beam' + cooked_data['start_time_dt'].apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') +'_to_' + cooked_data['end_time_dt'].apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') + '_plasma_condition_short.png'

    index = (cooked_data['DIST'] >= 7) & (cooked_data['DIST'] < 9)
    cooked_data.loc[index,'dist_region'] = 'near'
    index = cooked_data['DIST'] >= 9
    cooked_data.loc[index,'dist_region'] = 'tail'

    index = ((cooked_data['GSM_X'] > -1) & (cooked_data['GSM_Z'] < 0)) | ((cooked_data['GSM_X'] < -1) & (cooked_data['BX_GSM'] < 0))
    cooked_data.loc[index,'hemi'] = 'south'
    index = ((cooked_data['GSM_X'] > -1) & (cooked_data['GSM_Z'] > 0)) | ((cooked_data['GSM_X'] < -1) & (cooked_data['BX_GSM'] > 0))
    cooked_data.loc[index,'hemi'] = 'north'

    cooked_data.loc[:, 'flag'] = 0
    index = ((cooked_data['hemi'] == 'south') & (data['FLAG_PARA'] == 1))
    cooked_data.loc[index, 'flag'] = 1
    cooked_data.loc[index, 'flux'] = data['FLUX_PARA']
    cooked_data.loc[index, 'energy'] = data['EN_PARA']
    cooked_data.loc[index, 'eflux'] = data['EFLUX_PARA']
    cooked_data.loc[index, 'imfBy'] = data['IMF_BY_PARA']
    cooked_data.loc[index, 'imfBz'] = data['IMF_BZ_PARA']
    cooked_data.loc[index, 'pa'] = data['PA_PARA']


    index = ((cooked_data['hemi'] == 'north') & (data['FLAG_ANTI'] == 1))
    cooked_data.loc[index, 'flag'] = -1
    cooked_data.loc[index, 'flux'] = data['FLUX_ANTI']
    cooked_data.loc[index, 'energy'] = data['EN_ANTI']
    cooked_data.loc[index, 'eflux'] = data['EFLUX_ANTI']
    cooked_data.loc[index, 'imfBy'] = data['IMF_BY_ANTI']
    cooked_data.loc[index, 'imfBz'] = data['IMF_BZ_ANTI']
    cooked_data.loc[index, 'pa'] = data['PA_ANTI']

    cooked_data['energy_int'] = round(cooked_data['energy'])
    
    cooked_data['denergy'] = cooked_data['energy_int'].apply(find_denergy)

    cooked_data['intergrated_flux'] =  cooked_data['flux'] * cooked_data['denergy']* 2*PI * (PI/8)
    
    cooked_data = cooked_data.sort_values(by=['datetime_str'])
    
    return(cooked_data)

def aggregate_energy(df):      
    
    df = df.loc[(df.loc[:,'energy']).apply(np.isfinite),:]
    
    agg_data = df.groupby(['TIME']).agg({'GSE_X':'count' ,'date':'min', 'flux':'sum', 'eflux':'sum'
                                         , 'energy':'mean',  'denergy':'mean', 'intergrated_flux':'sum'
                                         , 'xgsm':'min', 'ygsm':'min', 'zgsm':'min', 'pa':'mean'
                                         , 'MLT':'min', 'L':'min',  'BX_GSM':'min' , 'BY_GSM':'min' 
                                         ,'BZ_GSM':'min', 'DIST':'min', 'BETA':'min', 'datetime_str':'min'
                                         , 'KP':'min', 'SW_P':'min', 'DST':'min', 'imfBy':'min', 'imfBz':'min'
                                         , 'BX_GSM':'min', 'BY_GSM':'min', 'BZ_GSM':'min'}).reset_index()

    agg_data['region'] = agg_data.apply(identify_region, axis=1)
    agg_data['B'] = agg_data.apply(calculate_B, axis=1)

    agg_data.rename(columns={'GSE_X':'nbeam'}, inplace = True)
    
    return(agg_data)

def preprocess_dispersion_list(dispersion_list):
    dispersion_list['p_value'] = dispersion_list.apply(calculate_cooked_data,axis = 1)
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

