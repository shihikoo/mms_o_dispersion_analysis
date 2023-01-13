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

def calculate_ccooked_data(dispersion):
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
    
    
def calculate_velocity(energy, ion_mass = 15.89):
    Avogadro_constant = 6.02214086e23 # moe-1
    electron_charge = 1.60217662e-19  #coulombs
    velocity = math.sqrt(2.*energy*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3 # 4.577e7*sqrt(data_energy/1e3/AMU) 
    return(velocity)
        
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

    cooked_data.loc[cooked_data['BETA'] < 0.05,'region'] = 'lobe'
    cooked_data.loc[(cooked_data['BETA'] < 1) & (cooked_data['BETA'] >= 0.05),'region'] = 'bl'
    cooked_data.loc[cooked_data['BETA'] >= 1,'region'] = 'ps'

    cooked_data['compression_mode'] = (cooked_data['datetime_str'] < pd.Timestamp('2019-4-16')) | (cooked_data['datetime_str'] > pd.Timestamp('2019-8-17'))

    cooked_data['start_time'] = (((cooked_data['datetime_str'].dt.hour/4).apply(int)))*4
    cooked_data['end_time'] = cooked_data['start_time'] + 4
    cooked_data['start_time_dt'] = cooked_data['datetime_str'].apply(datetime.datetime.combine,time=datetime.time.min) + cooked_data['start_time'].apply(pd.Timedelta,unit="h")
    cooked_data['end_time_dt'] = cooked_data['datetime_str'].apply(datetime.datetime.combine,time=datetime.time.min) + cooked_data['end_time'].apply(pd.Timedelta,unit="h")

    cooked_data['o_beam_filepath'] = 'idl_plots/obeam_day/'+cooked_data['start_time_dt'].apply(pd.Timestamp.strftime,format='%Y') +'/o_beam' + cooked_data['start_time_dt'].apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') +'_to_' + cooked_data['end_time_dt'].apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') + '_plasma_condition_short.png'

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

    index = ((cooked_data['hemi'] == 'north') & (data['FLAG_ANTI'] == 1))
    cooked_data.loc[index, 'flag'] = -1
    cooked_data.loc[index, 'flux'] = data['FLUX_ANTI']
    cooked_data.loc[index, 'energy'] = data['EN_ANTI']
    cooked_data.loc[index, 'eflux'] = data['EFLUX_ANTI']
    cooked_data.loc[index, 'imfBy'] = data['IMF_BY_ANTI']
    cooked_data.loc[index, 'imfBz'] = data['IMF_BZ_ANTI']

    cooked_data['energy_int'] = round(cooked_data['energy'])

    cooked_data = cooked_data.sort_values(by=['datetime_str'])

    return(cooked_data)

def preprocess_dispersion_list(dispersion_list):
    dispersion_list['p_value'] = dispersion_list.apply(calculate_ccooked_data,axis = 1)
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

