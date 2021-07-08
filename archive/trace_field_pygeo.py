#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis for dispersion observations of O+

Created on Wed Jun  9 17:57:28 2021

@author: jliao
"""

# In[1]: 
import pandas as pd
from functions import pygeo_wrapper
from functions import data_preprocess_functions

model = "t89"
filename = 'data/fulldata_20160101_to_20171231.csv'
print(model)

# In[2]: 
data = data_preprocess_functions.preprocess_data(pd.read_csv(filename))
dispersion_list = data_preprocess_functions.extract_dispersions(data)

# In[4]:
dispersion_list['model'] = dispersion_list.apply(pygeo_wrapper.get_magnetic_model,axis = 1)
dispersion_list['model_field_line_length_python_geopack'] = dispersion_list['model'].apply(pygeo_wrapper.extract_field_line_length)
dispersion_list['model_x'] = dispersion_list['model'].apply(pygeo_wrapper.extract_x)
dispersion_list['model_y'] = dispersion_list['model'].apply(pygeo_wrapper.extract_y)
dispersion_list['model_z'] = dispersion_list['model'].apply(pygeo_wrapper.extract_z)

# In[8]:
dispersion_list.drop(['model'], axis=1)

dispersion_list.to_csv('output/'+model+'_dispersion_list.csv')