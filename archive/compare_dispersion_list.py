#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 22:27:50 2021

@author: jliao
"""

# In[1]: 
from collections import Counter
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import functions

pio.renderers.default='browser'

goodness_of_fit_threshold = 0.6
model = "t89"
radius_earth = 6378.14 #km
# In[2]: 
dispersion_filename = 'output/'+model+'_dispersion_list.csv'
dispersion_list = pd.read_csv(dispersion_filename)

dispersion_m_filename = 'data/dispersion list - mms.csv'
dispersion_list_m = pd.read_csv(dispersion_m_filename)
# In[3]:
dispersion_list_m['max_v'] = dispersion_list_m['Max E'].apply(functions.calculate_velocity)
dispersion_list_m['min_v'] = dispersion_list_m['Min E'].apply(functions.calculate_velocity)
dispersion_list_m['mean_E'] = (dispersion_list_m['Max E'] + dispersion_list_m['Min E'])/2
dispersion_list_m['estimated_distance_m'] = dispersion_list_m['Duration'] * 60. / (1./dispersion_list_m['min_v'] - 1./dispersion_list_m['max_v']) / radius_earth

# In[3]:
dispersion_list_m_merged = pd.merge(left=dispersion_list_m, right=dispersion_list, how='left', left_on='Full list Index', right_on='index')
dispersion_list_m_merged['good_fit'] = dispersion_list_m_merged['p_value'] > goodness_of_fit_threshold

# In[9]:
index = dispersion_list_m_merged['p_value'] > 0
fig = px.scatter(dispersion_list_m_merged.loc[index,:], x="model_field_line_length_idl", y="estimated_distance_m"
                 , size="dispersion_time"                 , hover_name="Event Time"
                 , hover_data = ["index", "BETA", "direction", "GSM_Z", "p_value", 'GSM_X', 'BX_GSM']
#                , trendline="ols"
   #             # , symbol = "direction"
          #       , facet_col="direction_et"
         #        , facet_row="region"
                , color = np.log10(dispersion_list_m_merged.loc[index,"mean_E"]), range_color=[3,5]
                 )
reference_line = go.Scatter(x=[0, 35], y=[0, 35], mode="lines", line=go.scatter.Line(color="gray"), showlegend=False)

fig.add_trace(reference_line)
#fig.add_trace(reference_line, row=1, col=1)
#fig.add_trace(reference_line, row=1, col=2)

#fig.update_coloraxes(cmin=-2)
#fig.update_coloraxes(cmax=1)
fig.update_layout(#title = ' ' + ', goodness of fit > ' + str(goodness_of_fit_threshold),
 font=dict(
        family="Courier New, monospace",
        size=18    ), legend_x = 0, legend_y=1)

fig.update_layout(coloraxis_colorbar=dict(
    title="energy",
    tickvals=[1,2,3,4,5],
    ticktext=["10", "100", "1k","100k","100k"],
))

fig.show()

# In[10]:
#fig = px.histogram(dispersion_list, x="goodness_of_fit", nbins=13)
#fig.show()