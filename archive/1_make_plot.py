#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 20:21:47 2021
@author: jliao
"""

# In[1]: 
from collections import Counter
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'

model = "t89"
dispersion_filename = 'output/'+model+'_dispersion_list.csv'
dispersion_list = pd.read_csv(dispersion_filename)
# In[3]:
goodness_of_fit_threshold = 0.6
index = (dispersion_list.loc[:,'p_value'] > goodness_of_fit_threshold)

# In[9]:
fig = px.scatter(dispersion_list.loc[index,:], x="model_field_line_length_idl", y="estimated_distance"
                 , size="dispersion_time", hover_name="datetime_str"
                 , hover_data = ["index", "BETA", "direction", "GSM_Z", "p_value", 'GSM_X', 'BX_GSM']
                , trendline="ols"
                # , symbol = "direction"
          #       , facet_col="direction_et"
         #        , facet_row="region"
               , color = "energy" #np.log10(dispersion_list.loc[index,"Beta"]) 
                 )
reference_line = go.Scatter(x=[0, 70], y=[0, 70], mode="lines", line=go.scatter.Line(color="gray"), showlegend=False)

fig.add_trace(reference_line)
#fig.add_trace(reference_line, row=1, col=1)
#fig.add_trace(reference_line, row=1, col=2)

#fig.update_coloraxes(cmin=-2)
#fig.update_coloraxes(cmax=1)
fig.update_layout(title = 'Dispersion analysis, ' + model+ ', goodness of fit > ' + str(goodness_of_fit_threshold), font=dict(
        family="Courier New, monospace", size=18), legend_x = 0, legend_y=1)
#fig.update_layout(coloraxis_colorbar=dict(
#    title="Beta",
#    tickvals=[-2,-1,1,2],
#    ticktext=["0.01", "0.1", "1", "10"],
#))

fig.show()
# In[10]:
#fig = px.histogram(dispersion_list, x="goodness_of_fit", nbins=13)
#fig.show()