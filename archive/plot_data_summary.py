#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 16:28:51 2021

@author: jliao
"""

# In[1]: 
import pandas as pd
import plotly.express as px
import plotly.io as pio
pio.renderers.default='browser'

# In[2]: 
filename = 'data/fulldata_20160101_to_20171231.csv'

data = pd.read_csv(filename)

foo = pd.concat([data['En_para'] ,data['En_anti']],ignore_index = True)

# In[3]: 
fig = px.histogram(foo,  nbins=100)
fig.show()