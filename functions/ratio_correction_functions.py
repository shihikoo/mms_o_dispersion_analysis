#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 15:11:54 2023

@author: jliao
"""

import numpy as np
import statistics
import pandas as pd
import matplotlib.pyplot as plt

def compare_correction_figure(com_bl,ratio_bl, correction_x_bl, avg_ratio_bl, std_ratio_bl, correction_x_bl_2, avg_ratio_bl_2, std_ratio_bl_2,region, height, width):

    fig, ax = plt.subplots()
    fig.set_figheight(height)
    fig.set_figwidth(width)

    ax.scatter(com_bl, ratio_bl,label='HPCA observation',s=10)
    ax.errorbar(correction_x_bl, avg_ratio_bl, std_ratio_bl, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='HPCA observation')
    ax.errorbar(correction_x_bl_2, avg_ratio_bl_2, std_ratio_bl_2, marker='s', color='red', mfc='black',mec='black', ms=3, mew=1, label='Median Value')

    ax.legend()
    ax.set_title(region)
    ax.set_xlabel('Occurrence frequency with compression')
    ax.set_ylabel('Occurrence frequency correction factor')
    
    return 1

def draw_correction_figure(com_all, ratio_all,correction_x_all, avg_ratio_all, std_ratio_all,region):

    fig, ax = plt.subplots()
    fig.set_figheight(4)
    fig.set_figwidth(12)

    ax.scatter(com_all, ratio_all,label='HPCA observation',s=10)
    ax.errorbar(correction_x_all, avg_ratio_all, std_ratio_all, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='HPCA observation')

    ax.legend()

    ax.set_title(region)

    ax.set_xlabel('Occurrence frequency with compression')

    ax.set_ylabel('Occurrence frequency correction factor')
    
    return 1


def draw_correction_figures(com_lobe, ratio_lobe, correction_x_lobe,avg_ratio_lobe,  std_ratio_lobe, com_bl,ratio_bl, correction_x_bl, avg_ratio_bl, std_ratio_bl, com_ps,ratio_ps,correction_x_ps, avg_ratio_ps, std_ratio_ps ,  com_all, ratio_all,correction_x_all, avg_ratio_all, std_ratio_all , filename, yrange=[0,0]):
    
    # Draw the figure 
    fig, axs = plt.subplots(1, 4,sharex=True, sharey=True)
    fig.set_figheight(4)
    fig.set_figwidth(12)
    axs[0].set_xlim([0, 1])
    
    if yrange[0] != yrange[1]:
        axs[0].set_ylim(yrange)
    
    axs[0].scatter(com_lobe, ratio_lobe,label='HPCA observation',s=10)
    axs[0].errorbar(correction_x_lobe, avg_ratio_lobe, std_ratio_lobe, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='Median Value')

    axs[1].scatter(com_bl, ratio_bl,label='HPCA observation',s=10)
    axs[1].errorbar(correction_x_bl, avg_ratio_bl, std_ratio_bl, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='Median Value')
#    axs[1].errorbar(correction_x_bl_2, avg_ratio_bl_2, std_ratio_bl_2, marker='s', color='red', mfc='black',mec='black', ms=3, mew=1, label='Median Value')

    axs[2].scatter(com_ps, ratio_ps,label='HPCA observation', s=10)
    axs[2].errorbar(correction_x_ps, avg_ratio_ps, std_ratio_ps, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='Median Value')

    axs[3].scatter(com_all, ratio_all,label='HPCA observation', s=10)
    axs[3].errorbar(correction_x_all, avg_ratio_all, std_ratio_all, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='Median Value')
    
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[3].legend()
    
    axs[0].set_title('Lobe')
    axs[1].set_title('PSBL')
    axs[2].set_title('Plasma Sheet')
    axs[3].set_title('All Region')

    axs[0].set_xlabel('Occurrence frequency')
    axs[1].set_xlabel('Occurrence frequency')
    axs[2].set_xlabel('Occurrence frequency')
    axs[3].set_xlabel('Occurrence frequency')

    axs[0].set_ylabel('Occurrence frequency correction factor')
    axs[1].set_ylabel('Occurrence frequency correction factor')
    axs[2].set_ylabel('Occurrence frequency correction factor')
    axs[3].set_ylabel('Occurrence frequency correction factor')

    plt.savefig(filename + '.svg', format='svg', dpi=300)
    plt.savefig(filename + '.png', format='png', dpi=300)
    
    return 1

def draw_correction_figures_3(com_lobe, ratio_lobe, correction_x_lobe,avg_ratio_lobe,  std_ratio_lobe, com_bl,ratio_bl, correction_x_bl, avg_ratio_bl, std_ratio_bl, com_ps,ratio_ps,correction_x_ps, avg_ratio_ps, std_ratio_ps ,  filename, yrange=[0,0]):
    
    # Draw the figure 
    fig, axs = plt.subplots(1, 3,sharex=True, sharey=True)
    fig.set_figheight(4)
    fig.set_figwidth(12)
    axs[0].set_xlim([0, 1])
    
    if yrange[0] != yrange[1]:
        axs[0].set_ylim(yrange)
    
    axs[0].scatter(com_lobe, ratio_lobe,label='HPCA observation',s=10)
    axs[0].errorbar(correction_x_lobe, avg_ratio_lobe, std_ratio_lobe, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='Median Value')

    axs[1].scatter(com_bl, ratio_bl,label='HPCA observation',s=10)
    axs[1].errorbar(correction_x_bl, avg_ratio_bl, std_ratio_bl, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='Median Value')
#    axs[1].errorbar(correction_x_bl_2, avg_ratio_bl_2, std_ratio_bl_2, marker='s', color='red', mfc='black',mec='black', ms=3, mew=1, label='Median Value')

    axs[2].scatter(com_ps, ratio_ps,label='HPCA observation', s=10)
    axs[2].errorbar(correction_x_ps, avg_ratio_ps, std_ratio_ps, marker='s', color='black', mfc='black',mec='black', ms=3, mew=1, label='Median Value')
    
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    
    axs[0].set_title('Lobe')
    axs[1].set_title('PSBL')
    axs[2].set_title('Plasma Sheet')

    axs[0].set_xlabel('Occurrence frequency after compression')
    axs[1].set_xlabel('Occurrence frequency after compression')
    axs[2].set_xlabel('Occurrence frequency after compression')

    axs[0].set_ylabel('Occurrence frequency correction factor')
    axs[1].set_ylabel('Occurrence frequency correction factor')
    axs[2].set_ylabel('Occurrence frequency correction factor')

    plt.savefig(filename + '.svg', format='svg', dpi=300)
    plt.savefig(filename + '.png', format='png', dpi=300)
    
    return 1

def read_data(path1, path2, grid_str, region, direction_str, start_time, end_time, energy = ''):
    file_1 = path1 + "/" + grid_str + "/" + region + "_" + direction_str + "/events/2d/" + region + "_" + direction_str + energy + "__ratio_" + start_time + "_to_" + end_time + "_X_GSM_vs_Y_GSM.csv"

    data_1 = np.array(pd.read_csv(file_1))
    data_1 = data_1[:,1:21]
    data_1 =  data_1.reshape((400))  

    file_2 = path2 + "/" + grid_str + "/" + region + "_" + direction_str + "/events/2d/" + region + "_" + direction_str + energy + "__ratio_" + start_time + "_to_" + end_time + "_X_GSM_vs_Y_GSM.csv"

    data_2 = np.array(pd.read_csv(file_2))
    data_2 = data_2[:,1:21]
    data_2 = data_2.reshape((400))  

    index = (~np.isnan(data_1)) & (~np.isnan(data_2)) & (data_2 != 0)
    nocom = data_1[index]
    com = data_2[index]

    index = np.argsort(com)
    com = com[index]
    nocom = nocom[index]

    ratio = nocom/com
    
    return com, nocom, ratio

def calculate_fit(x,y):
    slope = np.nansum((x - np.nanmean(x)) * (y - np.nanmean(y))) / np.nansum(((x - np.nanmean(x)) * (x - np.nanmean(x))))
    intercept = np.nanmean(y) - slope * np.nanmean(x)
    return(slope,intercept)

def calculate_avg(x,y, x_low, x_high, ):
    lower_edge = min(x_low)
    nbins = len(x_low)
    
    x_center = (x_low + x_high)/2

    y_avg = np.zeros(nbins)
    y_std = np.zeros([2, nbins])

    for i in range(0, x_low.shape[0]):
        index = (x > x_low[i]) & (x <= x_high[i])
        if y[index].shape[0] > 0:
            y_avg[i] = statistics.mean(y[index])
        if y[index].shape[0] > 2:
            y_std[:,i] = abs(statistics.quantiles(y[index],n=3)-y_avg[i])
#            y_std[:,i] = [statistics.stdev(y[index],y_avg[i]),statistics.stdev(y[index],y_avg[i])]

    return y_avg, y_std, x_center

def add_ratio_1(r1,r2,r3):
    
    return np.append(r1,1),np.array([np.append(r2[0,:],0), np.append(r2[1,:],0)]),np.append(r3,1)
    
def apply_correction_old(x, correction):
    lower_edge = 0
    nbins = correction.shape[0]
    step = 1/nbins
    x_low = step*np.array(list(range(lower_edge,nbins)))
    x_high = x_low + step
    
    
    x_corrected = np.zeros(x.shape[0])
    for i in range(lower_edge,nbins):
        index = (x > x_low[i]) & (x <= x_high[i])
        if x[index].shape[0] > 0:
            x_corrected[index] = x[index] * correction[i]
            
    return x_corrected

def apply_correction(x, ratio, x_low, x_high):
    lowest_edge = 0
    highest_edge = 0
    
    lower_edge = min(x_low)
    nbins = len(x_low)
    x_center = (x_low + x_high)/2
    x_center = np.append(x_center,1)
    ratio = np.append(ratio,1)
    
    x_corrected = np.zeros(x.shape[0])
    correction = np.zeros(x.shape[0])
    
    k = np.zeros(x_center.shape[0])
    m = np.zeros(x_center.shape[0])
    
    for i in range(int(lower_edge), nbins):
        index = (x >= x_center[i]) & (x <= x_center[i+1])
        if x[index].shape[0] > 0:
            k[i] = (ratio[i+1]-ratio[i])/(x_center[i+1]-x_center[i])
            m[i] = (x_center[i]*ratio[i+1]-x_center[i+1]*ratio[i])/(x_center[i]-x_center[i+1])
            correction[index] = x[index]*k[i]+m[i]
            x_corrected[index] = x[index]*correction[index]
    index = (x < x_center[0])
    if x[index].shape[0] > 0: 
        correction[index] = x[index]*k[0]+m[0]
        x_corrected[index] = x[index]*correction[index]
        
    return x_corrected

def inverse_func(x, a, b):
    return a/x+b

def log_func(x, a, b,c):
    return a*np.log(x)**2+b*np.log(x)+c