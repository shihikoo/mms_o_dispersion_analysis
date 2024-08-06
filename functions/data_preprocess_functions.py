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
from functions import geopack_wrapper
 

MMS_ENERGY_BIN = np.array([1.8099999, 3.5100000,6.7700000, 13.120000, 25.410000, 49.200001, 95.239998, 184.38000,356.97000, 691.10999,1338.0400,2590.4900, 5015.2900, 9709.7900, 18798.590, 32741.160])
MMS_ENERGY_BIN_INT = np.round(MMS_ENERGY_BIN)
np.array([2, 4, 7, 13, 25, 49, 95, 184,357, 691,1338,2590, 5015, 9710, 18799, 32741])
MMS_ENERGY_HIGH = np.array([2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79, 41598.37])
MMS_ENERGY_LOW = np.array([1.27, 2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79])
MMS_DENERGY = np.array([0.5,1,1.91,3.71,7.18,13.92,26.93,52.15,100.97,195.42,378.36,732.53,1418.24,2745.8,5315.99,6713.96])

PA_BIN_SIZE = 11.25
PI = 3.1415926



def read_multiple_csv(filenames, lowcase_colname = True):    
    li = []
    for filename in filenames:
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)

    df_output = pd.concat(li, axis=0, ignore_index=True)
    if lowcase_colname:
        df_output = df_output.rename(columns=str.lower)
    
    return df_output

def read_dispersion_csv(dispersion_filenames):
    df_dis = read_multiple_csv(dispersion_filenames, lowcase_colname = True)
    
    return df_dis

def read_beam_csv(beam_filenames):
    batch_size = 800
    li = []
    n_batch = math.ceil(len(beam_filenames) / batch_size)

    for ibatch in range(n_batch):
        df_beam = read_multiple_csv(beam_filenames[(ibatch*batch_size):min((ibatch+1)*batch_size, len(beam_filenames))], lowcase_colname = True)
        ind = (df_beam['flux_para'] > 0) | (df_beam['flux_anti'] > 0)
        li.append(df_beam.loc[ind,:])
    df_output = pd.concat(li, axis=0, ignore_index=True)
    
#     if len(beam_filenames) < nbatch:
#         df_beam = read_multiple_csv(beam_filenames[0:nbatch], lowcase_colname = True)
#         ind = (df_beam['flux_para'] > 0) | (df_beam['flux_anti'] > 0)
#         df_beam = df_beam.loc[ind,:]
        
#     else:
#         df_beam1 = read_multiple_csv(beam_filenames[0:nbatch], lowcase_colname = True)
#         ind = (df_beam1['flux_para'] > 0) | (df_beam1['flux_anti'] > 0)
#         df_beam1 = df_beam1.loc[ind,:]

#         df_beam2 = read_multiple_csv(beam_filenames[800:len(beam_filenames)], lowcase_colname = True)
#         ind = (df_beam2['flux_para'] > 0) | (df_beam2['flux_anti'] > 0)
#         df_beam2 = df_beam2.loc[ind,:]

#         df_output = pd.concat([df_beam1, df_beam2], axis=0, ignore_index=True)
    return df_output

def read_external_csv(external_filenames):
    df_ext = read_multiple_csv(external_filenames, lowcase_colname = True)
    
    return df_ext

def read_beam_external_from_csv(beam_filenames, ext_filenames):
    
    df_beam = read_beam_csv(beam_filenames)
    df_ext = read_external_csv(ext_filenames)
    
    df0 = pd.merge(df_beam, df_ext, on = 'time', how='outer')

    df = preprocess_data(df0)
    return df

def extract_date(input_datetime_obj):
    date = input_datetime_obj.strftime("%Y-%m-%d")
    return(date)

def extract_time(input_datetime_obj):
    time = input_datetime_obj.strftime("%H-%M-%S")
    return(time)

def remove_large_y(df):
    new_df = df.loc[df['ygsm'] <= 20,:].reindex()
    return new_df

def remove_outside_magnetosphere(df):
    new_df = df.loc[ (df['region'] == 'Lobe') | (df['region'] == 'BL') | (df['region'] == 'PS') | (df['region'] == 'Dayside') ,:].reindex()
    return new_df

def calculate_pvalue(dispersion):
    p = 1-stats.chi2.cdf(dispersion['chisq'], dispersion['dof'])
    return(p)

def identify_region(onedata):
    if(onedata['region'] >= 10):
        region = 'outside'
    elif(onedata['mlt'] >= 8.) & (onedata['mlt'] < 16.):
        region = 'Dayside'
    elif(onedata['region'] == 1):
        region = 'Lobe'
    elif(onedata['region'] == 2):
        region = 'BL'
    elif(onedata['region'] == 3):
        region = 'PS'
    else:
        region = 'outside'   
    
#     if (onedata['mlt'] >= 8.) & (onedata['mlt'] < 16.):
#         region = 'Dayside'
#     else:
#         if onedata.beta < 0.05:
#             region = 'Lobe'
#         elif onedata.r >= 15:
#             if onedata.beta < 1.:
#                 region = 'BL'
#             else:
#                 region = 'PS'
#         else:
#             if onedata.beta > math.exp(0.14*onedata.r)-2.1:
#                 region = 'BL'
#             else:
#                 region = 'PS'             
            
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
    return(math.sqrt(onedata['bx']**2 + onedata['by_gsm']**2 + onedata['bz_gsm']**2))

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
    
def extract_beam_info(mydf, direction):
    if direction == 'para':
        flag = 1
        hemi = 'south'
    else:
        flag = -1
        hemi = 'north'
    
    index = ((mydf['hemi'] == hemi) & (mydf['flag_'+direction] == 1))
    mydf.loc[index, 'flag'] = flag
    mydf.loc[index, 'flux'] = mydf.loc[index,'flux_'+direction]
    mydf.loc[index, 'int_flux'] = mydf.loc[index,'int_flux_'+direction]
    mydf.loc[index, 'energy'] = mydf.loc[index,'en_'+direction]
    mydf.loc[index, 'eflux'] = mydf.loc[index,'eflux_'+direction]
    mydf.loc[index, 'pa'] = mydf.loc[index,'pa_'+direction]
    mydf.loc[index, 'pa_range'] = mydf.loc[index,'pa_range_'+direction]
    mydf.loc[index, 'n'] = mydf.loc[index,'n_'+direction]
    mydf.loc[index, 't'] = mydf.loc[index,'t_'+direction]
    mydf.loc[index, 'p'] = mydf.loc[index,'p_'+direction]
    
    mydf.loc[index, 'n_weighted'] = mydf.loc[index,'n_'+direction+'_weighted']
    
    mydf.loc[index, 'imfBy'] = mydf.loc[index,'imf_by_'+direction+'_1h']
    mydf.loc[index, 'imfBz'] = mydf.loc[index,'imf_bz_'+direction+'_1h']
    mydf.loc[index, 'swp'] = mydf.loc[index,'sw_p_'+direction+'_1h']
    mydf.loc[index, 'swv'] = mydf.loc[index,'sw_v_'+direction+'_1h']
    
    mydf.loc[index, 'fllen'] = mydf.loc[index,'fllen_'+direction]

    return mydf

def find_filepath(datetime_str, dir_name='', file_append_name = 'identification', avg_hour = 6):
    
    start_time = (((datetime_str.dt.hour/avg_hour).apply(int)))*avg_hour
    end_time = start_time + avg_hour
    start_time_dt = datetime_str.apply(datetime.datetime.combine,time=datetime.time.min) + start_time.apply(pd.Timedelta,unit="h")
    end_time_dt = datetime_str.apply(datetime.datetime.combine,time=datetime.time.min) + end_time.apply(pd.Timedelta,unit="h")
    
    o_beam_filepath = dir_name+'o_beam' + start_time_dt.apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') +'_to_' + end_time_dt.apply(pd.Timestamp.strftime,format='%Y%m%d_%H%M%S') + '_'+file_append_name+'.png'
    return o_beam_filepath
    
def extract_hemisphere(cooked_data):
    index = ((cooked_data['xgsm'] > -1) & (cooked_data['zgsm'] < 0)) | ((cooked_data['xgsm'] < -1) & (cooked_data['bx'] < 0))
    cooked_data.loc[index,'hemi'] = 'south'
    index = ((cooked_data['xgsm'] > -1) & (cooked_data['zgsm'] > 0)) | ((cooked_data['xgsm'] < -1) & (cooked_data['bx'] > 0))
    cooked_data.loc[index,'hemi'] = 'north'
    return cooked_data

# this function is used to 
def preprocess_data(data, remove_large_y = False, avg_hour = 6):
    cooked_data = data
    
    cooked_data.rename(columns={'bx_gsm':'bx','gse_x':'xgse', 'gse_y':'ygse','gse_z':'zgse','gsm_x':'xgsm', 'gsm_y':'ygsm','gsm_z':'zgsm', 'o_vpar':'v_par_o_all', 
'o_vperp':'v_perp_o_all', 'h_vperp':'v_perp_h_all','h_vpar':'v_par_h_all',
'o_n':'density_o_all', 'o_v':'velocity_o_all','o_p':'pressure_o_all', 'h_n':'density_h_all', 'h_v':'velocity_h_all','h_p':'pressure_h_all',}, inplace = True)
    
    if 'b_model' not in cooked_data.columns:
        cooked_data['b_model'] = 0
    if 'fllen_anti' not in cooked_data.columns:
        cooked_data['fllen_anti'] = 0
    if 'fllen_para' not in cooked_data.columns:
        cooked_data['fllen_para'] = 0    
    
    if 'en' in cooked_data.columns:
        cooked_data['en_para'] = cooked_data['en']
        cooked_data['en_anti'] = cooked_data['en']
    
    # datetime extraction
    cooked_data['datetime_str'] = cooked_data.loc[:,'time'].apply(datetime.datetime.utcfromtimestamp)
    cooked_data['date'] = cooked_data.loc[:,'datetime_str'].apply(extract_date)   
    cooked_data['year'] = cooked_data['datetime_str'].dt.to_period('Y')
                            
#     cooked_data['kp_gt_2'] = cooked_data['kp'] > 2 
    cooked_data['storm'] = cooked_data['storm_phase'] > 0
    cooked_data['storm_phase'] = pd.Categorical(cooked_data['storm_phase']).rename_categories({0:'nonstorm',1:'prestorm',2:'main phase',3:'fast recovery', 4:'long recovery'})
    
    # extra infomation
    cooked_data['r'] = (cooked_data['ygsm']**2 + cooked_data['zgsm']**2).apply(math.sqrt)

    cooked_data['region'] = cooked_data.apply(identify_region, axis=1)

    cooked_data['compression_mode'] = (cooked_data['datetime_str'] < pd.Timestamp('2019-4-16')) | (cooked_data['datetime_str'] > pd.Timestamp('2019-8-17'))

    cooked_data['b'] = cooked_data.apply(calculate_B, axis=1)
    
    # parameters for extracting the tplot images filepath
    cooked_data['o_beam_filepath'] = find_filepath(cooked_data['datetime_str'], dir_name='plots/obeam_day/identification/', file_append_name = 'identification', avg_hour = 6)
    # extract south/north hemisphere according to definitions in the identification
    cooked_data = extract_hemisphere(cooked_data)
    
    cooked_data['flag'] = 0
    cooked_data = extract_beam_info(cooked_data,'para')
    cooked_data = extract_beam_info(cooked_data,'anti')
   
    cooked_data['energy_int'] = round(cooked_data['energy'])    
    cooked_data['denergy'] = cooked_data['energy_int'].apply(find_denergy)
    cooked_data['intergrated_flux'] =  cooked_data['int_flux'] * cooked_data['denergy']
    
    cooked_data['density_est'] =  cooked_data.apply(estimate_density, axis=1)    
       
    cooked_data = cooked_data.sort_values(by=['datetime_str'])
       
    if remove_large_y:
        cooked_data = remove_large_y(cooked_data)

    remove_outside_magnetosphere(cooked_data)
        
    #     index = (cooked_data['dist'] >= 7) & (cooked_data['dist'] < 9)
#     cooked_data.loc[index,'dist_region'] = 'near'
#     index = cooked_data['dist'] >= 9
#     cooked_data.loc[index,'dist_region'] = 'tail'

    return cooked_data

def aggregate_angle(df):      
    df = df.loc[(df.loc[:,'pa']).apply(np.isfinite),:]
    
    agg_data = df.groupby(['time','energy']).agg({ 'xgse':'count' , 'flag':'mean'
                                                  , 'date':'first', 'datetime_str':'first', 'year':'first'
                                                  , 'xgsm':'first', 'ygsm':'first', 'zgsm':'first'
                                                  , 'ygse':'first', 'zgse':'first', 'mlt':'first', 'l':'first'
                                                  ,  'bx':'first' , 'by_gsm':'first','bz_gsm':'first', 'b':'first'
                                                  , 'dist':'first', 'beta':'first'
                                                  , 'kp':'first', 'swp':'first', 'swv':'first', 'dst':'first','f107':'first'
                                                  , 'imfBy':'first', 'imfBz':'first' 
                                                  , 'storm_phase':'first', 'compression_mode':'first'
                                                  , 'density_o_all':'first', 'velocity_o_all':'first','pressure_o_all':'first'
                                                  , 'density_h_all':'first', 'velocity_h_all':'first','pressure_h_all':'first'
                                                  , 'v_par_o_all':'mean','v_perp_o_all':'mean'
                                                  , 'v_par_h_all':'mean', 'v_perp_h_all':'mean'  
                                                  , 'r':'first', 'region':'first'
                                                  , 'o_beam_filepath':'first'
                                                  , 'pa':'mean','pa_range':'mean','int_flux':'mean'
                                                  , 'flux':'mean', 'eflux':'sum', 'intergrated_flux':'sum', 'density_est':'sum'
                                                  , 't':'mean', 'n':'sum','p':'sum' , 'denergy':'mean','n_weighted':'sum'
                                                  , 'b_model':'mean', 'fllen':'mean'}).reset_index()
    agg_data.rename(columns={'xgse':'nbeam'}, inplace = True)
    
    return(agg_data)

def aggregate_energy(df):      
    
    df = df.loc[(df.loc[:,'energy']).apply(np.isfinite),:]
    
    agg_data = df.groupby(['time']).agg({'nbeam':'sum' , 'energy':'mean', 'flag':'mean'
                                         , 'date':'first', 'datetime_str':'first', 'year':'first'
                                         , 'xgsm':'first', 'ygsm':'first', 'zgsm':'first'
                                          , 'ygse':'first', 'zgse':'first', 'mlt':'first', 'l':'first'
                                          , 'bx':'first' , 'by_gsm':'first','bz_gsm':'first', 'b':'first'
                                          , 'dist':'first', 'beta':'first'
                                          , 'kp':'first', 'swp':'first', 'swv':'first', 'dst':'first','f107':'first'
                                          , 'imfBy':'first', 'imfBz':'first' 
                                          , 'storm_phase':'first', 'compression_mode':'first'
                                          , 'density_o_all':'first', 'velocity_o_all':'first','pressure_o_all':'first'
                                          , 'density_h_all':'first', 'velocity_h_all':'first','pressure_h_all':'first'
                                          , 'v_par_o_all':'mean', 'v_perp_o_all':'mean'
                                          , 'v_par_h_all':'mean', 'v_perp_h_all':'mean'  
                                          , 'denergy':'mean','r':'first', 'region':'first'
                                          , 'o_beam_filepath':'first'
                                          , 'pa':'mean','pa_range':'mean','int_flux':'mean'
                                          , 'flux':'mean', 'eflux':'sum', 'intergrated_flux':'sum', 'density_est':'sum'
                                          , 't':'mean', 'n':'sum','p':'sum','n_weighted':'sum'
                                        , 'b_model':'mean', 'fllen':'mean'}).reset_index()   
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

# this function is for old data when dispersion is saved in beam data
def preprocess_dispersion_list(dispersion_list, model = 't89'):   
    dispersion_list.rename(columns={'bx_gsm':'bx','gse_x':'xgse', 'gse_y':'ygse','gse_z':'zgse','gsm_x':'xgsm', 'gsm_y':'ygsm','gsm_z':'zgsm'}, inplace = True)
    
    if 'b_model' not in dispersion_list.columns:
        dispersion_list['b_model'] = 0
    if 'fllen_anti' not in dispersion_list.columns:
        dispersion_list['fllen_anti'] = 0
    if 'fllen_para' not in dispersion_list.columns:
        dispersion_list['fllen_para'] = 0

    dispersion_list['model_field_line_length_python'] = dispersion_list.apply(geopack_wrapper.get_magnetic_model, model = model, axis=1)

    dispersion_list['dispersion_time'] = 2. * (dispersion_list['dis_fitting_dof']+2)
    dispersion_list['bias'] = abs((dispersion_list['estimated_distance']-dispersion_list['model_field_line_length_python']) / dispersion_list['model_field_line_length_python'])
    
    dispersion_list['temporal_dispersion'] = dispersion_list['bias'] <= 0.3
    
    dispersion_list = dispersion_list.sort_values(by=['datetime_str'])
    
    return(dispersion_list)


def extract_dispersion_list_old(mydata, direction_name = 'para'):
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

def extract_dispersions_old(data, save_to_filename = 'output/dispersion_list.csv'):
    dispersion_para = extract_dispersion_list_old(data, direction_name = 'para')
    dispersion_anti = extract_dispersion_list_old(data, direction_name = 'anti')
    dispersion_list = pd.concat([dispersion_para,dispersion_anti],ignore_index = True)  
    
    dispersion_list = preprocess_dispersion_list(dispersion_list)
    
    dispersion_list.to_csv(save_to_filename)
    
    return(dispersion_list)
