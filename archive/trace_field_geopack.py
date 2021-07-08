#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trace model field lines with package geopack

Created on Wed Jun  9 17:57:28 2021

@author: jliao
"""

# In[1]: 
from functions import geopack_wrapper
from functions import data_preprocess_functions
import pandas as pd

model = "t89"
filename = 'data/fulldata_20160101_to_20171231.csv'

# In[2]: 
data = data_preprocess_functions.preprocess_data(pd.read_csv(filename))
dispersion_list = data_preprocess_functions.extract_dispersions(data)

# In[3]: 
dispersion_list['model_field_line_length_python'] = None

# In[4]:
print(model)
dispersion_list['FlLen'] = dispersion_list.apply(geopack_wrapper.get_magnetic_model, model = model, axis=1)
# In[4]:

#dispersion_index = dispersion_list.index[dispersion_list.loc[:, 'model_field_line_length_python'].apply(lambda x: x is None)]
#for ii in dispersion_index:
#    if (dispersion_list.loc[ii,'GSM_X'] < -15):
#        continue
#    print(ii)
#    dispersion_list.loc[ii,'model_field_line_length_python'] = geopack_wrapper.get_magnetic_model(dispersion_list.loc[ii,:], model = model)

# In[8]:

dispersion_list.to_csv('output/'+model+'_dispersion_list.csv')