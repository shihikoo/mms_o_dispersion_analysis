#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 22:44:32 2021

@author: jliao
"""

from matplotlib import pyplot as plt
from matplotlib.patches import Wedge, Circle
from geopack import geopack, t89, t96, t01, t04
import numpy as np
import math
import plotly.express as px


def dual_half_circle(center=(0,0), radius=1, angle=90, ax=None, colors=('w','k','k'),
                     **kwargs):
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[1], **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[0], **kwargs)
   
    cr = Circle(center, radius, fc=colors[2], fill=False, **kwargs)
    for wedge in [w1, w2, cr]:
        ax.add_artist(wedge)
    return [w1, w2, cr]

def setup_fig(xlim=(10,-30),ylim=(-20,20),xlabel='X GSM [Re]',ylabel='Z GSM [Re]'):

    fig = plt.figure(figsize=(15,10))
    ax  = fig.add_subplot(111)
    ax.axvline(0,ls=':',color='k')
    ax.axhline(0,ls=':',color='k')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.set_aspect('equal')
    w1,w2,cr = dual_half_circle(ax=ax)
    
    return ax

def test_t89(data,i):
    idata = data.loc[i,:]
    ut = idata["time"]
    ps = geopack.recalc(ut)
    xgsm = idata['xgsm']
    ygsm = idata['ygsm']
    zgsm = idata['zgsm']
    kp = idata['kp']
    
    b0xgsm,b0ygsm,b0zgsm = geopack.igrf_gsm(xgsm,ygsm,zgsm)
    dbxgsm,dbygsm,dbzgsm = t89.t89(int(kp), ps, xgsm,ygsm,zgsm) 
    
    bxgsm,bygsm,bzgsm = [b0xgsm+dbxgsm,b0ygsm+dbygsm,b0zgsm+dbzgsm]
    return(bxgsm,bygsm,bzgsm )
    
def test_t96(data,i):
    idata = data.loc[i,:]
    ut = idata["Time"]
    ps = geopack.recalc(ut)
    xgsm = idata['GSM_X']
    ygsm = idata['GSM_Y']
    zgsm = idata['GSM_Z']
    par = [idata['SW_p'], idata['Dst'], idata['IMF_Bz'],idata['IMF_By'],0.,0., ps, xgsm,ygsm,zgsm]
    
    b0xgsm,b0ygsm,b0zgsm = geopack.igrf_gsm(xgsm,ygsm,zgsm)
    dbxgsm,dbygsm,dbzgsm = t96.t96(par, ps, xgsm,ygsm,zgsm)
                                   
    bxgsm,bygsm,bzgsm = [b0xgsm+dbxgsm,b0ygsm+dbygsm,b0zgsm+dbzgsm]
    return(bxgsm,bygsm,bzgsm )
    
def test_t01(data,i):
    idata = data.loc[i,:]
    ut = idata["Time"]
    ps = geopack.recalc(ut)
    xgsm = idata['GSM_X']
    ygsm = idata['GSM_Y']
    zgsm = idata['GSM_Z']
    par = [idata['SW_p'], idata['Dst'], idata['IMF_Bz'],idata['IMF_By'],0.,0., ps, xgsm,ygsm,zgsm]

    b0xgsm,b0ygsm,b0zgsm = geopack.igrf_gsm(xgsm,ygsm,zgsm)
    dbxgsm,dbygsm,dbzgsm = t01.t01(par, ps,xgsm,ygsm,zgsm)  

    bxgsm,bygsm,bzgsm = [b0xgsm+dbxgsm,b0ygsm+dbygsm,b0zgsm+dbzgsm]
    return(bxgsm,bygsm,bzgsm )
    
def test_t04(data,i):
    idata = data.loc[i,:]
    ut = idata["Time"]
    ps = geopack.recalc(ut)
    xgsm = idata['GSM_X']
    ygsm = idata['GSM_Y']
    zgsm = idata['GSM_Z']
    par = [idata['SW_p'], idata['Dst'], idata['IMF_Bz'],idata['IMF_By'],0.,0., ps, xgsm,ygsm,zgsm] 

    b0xgsm,b0ygsm,b0zgsm = geopack.igrf_gsm(xgsm,ygsm,zgsm)
    dbxgsm,dbygsm,dbzgsm = t04.t04(par, ps,xgsm,ygsm,zgsm)  

    bxgsm,bygsm,bzgsm = [b0xgsm+dbxgsm,b0ygsm+dbygsm,b0zgsm+dbzgsm]
    return(bxgsm,bygsm,bzgsm )

def test_trace(data, i, model):
    idata = data.loc[i,:]
    ut = idata["time"]
    ps = geopack.recalc(ut)
    xgsm = idata['xgsm']
    ygsm = idata['ygsm']
    zgsm = idata['zgsm']
    par = [idata['swp'], idata['dst'], idata['imfbz'],idata['imfby'],0.,0., ps, xgsm,ygsm,zgsm]
    dir = 1
    
    x,y,z,xx,yy,zz = geopack.trace(xgsm,ygsm,zgsm,dir=dir,rlim=100,r0=0.99999,parmod=par,exname = model,inname='igrf',maxloop=10000)
    
    ax=setup_fig()
    ax.plot(xx,zz)
    plt.show()
    return x,y,z,xx,yy,zz

def calculate_fieldlineLen(time, xgsm,ygsm,zgsm, dir, kp, pdyn, dst, imfBz, imfBy, model):
    ps = geopack.recalc(time)
    if model == 't89': 
        parmod = int(kp)
    else:
        parmod=[pdyn,dst,imfBz,imfBy,0.,0.,ps, xgsm,ygsm,zgsm]
        if math.isnan(parmod[0]):
            return None, None, None, None

    x,y,z,xx,yy,zz = geopack.trace(xgsm,ygsm,zgsm,dir= dir,rlim=100,r0=0.99999, parmod=parmod, exname = model, inname='igrf', maxloop=10000)    
    n = len(xx)    
    flLen = np.sum(np.sqrt(np.square(xx[1:n] - xx[0:n-1]) + np.square(yy[1:n] - yy[0:n-1]) + np.square(zz[1:n] - zz[0:n-1])))
    return flLen,xx,yy,zz
        
    
def get_magnetic_model(onedata, plot = False, model = "t89"):
    xgsm = onedata['xgsm']
    if (math.isnan(xgsm)):
        return None 
        
    ygsm = onedata['ygsm']
    zgsm = onedata['zgsm']
    time = onedata['time']
    
    flLen,xx,yy,zz = calculate_fieldlineLen(time,xgsm,ygsm,zgsm, dir =  onedata["flag"], kp = onedata['kp'], pdyn = onedata['swp'], dst = onedata['dst'], imfBz = onedata['imfBz'], imfBy =onedata['imfBy'], model = model)

    if plot:
        ax=setup_fig()
        ax.plot(xx,zz)
        #plt.show()
        plt.savefig('plots/'+str(int(time)) + '.png')
    
    return(flLen)


def plot_mag_model(xx,yy,zz,filename = None): 
    fig = px.line(x = xx, y=zz)
    fig.update_xaxes(range=[10,-70])
    fig.update_yaxes(range=[-20,20])
    fig.show()
       
    return 0
    