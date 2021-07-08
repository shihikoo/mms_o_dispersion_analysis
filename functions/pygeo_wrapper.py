#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 23:05:43 2021
@author: jliao
"""

import PyGeopack as gp
#gp.UpdateParameters(SkipWParameters=True)

def get_field(onedata, model = 'T89'):
    x = onedata['GSM_X']
    y = onedata['GSM_Y']
    z = onedata['GSM_Z']
        
    date = float(onedata['datetime_str'].strftime("%Y%m%d"))
    ut = float(onedata['datetime_str'].strftime("%H")) + float(onedata['datetime_str'].strftime("%M"))/60. + float(onedata['datetime_str'].strftime("%S"))/3600.

    Bx,By,Bz = gp.ModelField(x,y,z,date,ut,Model=model, CoordIn='GSM', CoordOut='GSM')
    return(Bx,By,Bz)

def get_magnetic_model(onedata):
    x = onedata['GSM_X']
    y = onedata['GSM_Y']
    z = onedata['GSM_Z']
    date = int(onedata['datetime_str'].strftime("%Y%m%d"))
    ut = float(onedata['datetime_str'].strftime("%H")) + float(onedata['datetime_str'].strftime("%M"))/60. + float(onedata['datetime_str'].strftime("%S"))/3600.
    
    T = gp.TraceField(x,y,z,date,ut
                      , Model='T96',CoordIn='GSM',CoordOut='GSM', alt=100.0, MaxLen=1000, DSMax=1.0,FlattenSingleTraces=True,Verbose=True
                      , Kp = onedata['KP'], Pdyn = onedata['SW_P'], Symh = onedata['DST'], By = onedata['IMF_BY'], Bz = onedata['IMF_BZ'])
    return(T)
    
def extract_field_line_length(T):
    field_line_length = T.FlLen
    return(field_line_length)
    
def extract_x(T):
    field_line_length = T.x
    return(field_line_length)
    
def extract_y(T):
    field_line_length = T.y
    return(field_line_length)
    
def extract_z(T):
    field_line_length = T.z
    return(field_line_length)
    