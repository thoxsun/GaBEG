#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:34:11 2023

@author: tom
"""




#%% Init

import csv
import numpy as np
# import os
# import glob
import matplotlib.pyplot as plt
# from zero_point import zpt

# from scipy.stats import kde

import GAIA_plot as gp

#%% ### ------------ Set default plot style------------------------------------
#https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure

# Note: fig_width is fixed in plotStyle.py for MasterThesis specifications
#       fig_height can and should be adjusted to special needs of individual plots
#       fig_height of 5inch is approximatelly half of height on A4 MasterThesis      
fig_height = 5
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)



#%% ### ------------ Set default plot colors --------------------------------------

t00 = 'lightgrey'
t01 = 'grey'
t02 = 'black' 

t_fit = 'black'
t_lost = 'firebrick'
# pm 0.5
t03 = 'silver'
t03 = 'darkgreen'
# pm 1.0 & zpc_pm1.0 & pm 1.0 Ag
t04 = 'whitesmoke'
t04 = 'purple'
# Simbad
t05 = 'forestgreen'
# parallax
t06 = 'lightseagreen'

# dataDR3 - viaux13 convolution
t07 = 'orangered'

# Av extinction models
t08 = 'lightgrey'

# RGB
t09 = 'teal'
t_RGB = t07
# AGB
t10 = 'purple'
t_AGB = t09

# AGB RGB
t_ARGB = t10

#%% ### ------------ Load /Store data settings ---------------------------------------


## ------------ load data into variables --------------------------------------
import os
import glob

# Epoch photometry can be downloaded by selecting SingleObject at:
# https://gea.esac.esa.int/archive/
folder_name = '/data/singleObjectAnalysis/'
# search_string = '*_EPOCH_PHOTOMETRY-Gaia DR3*'
search_string = '*.csv'
# Get current work directory
work_dir = os.getcwd()
file_path = work_dir + folder_name + search_string 

file_paths = glob.glob(file_path)

star_epoch_photometry = []
for i, el in enumerate(file_paths):
      file_data = np.genfromtxt(el,skip_header = 0, encoding = None, delimiter = ',' ,dtype=None, names = True)
      star_epoch_photometry.append(file_data)


#---------Store GAIA data --------------------------------

file_sub_path = 'data/trgbEvalData/'
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### ----------- Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/trgb_eval/'

#%% calculate weighted mean
def weighted_mean_1(timestamps,measurements):
  
  # Calculate time intervals (differences between consecutive timestamps)
  time_intervals = np.diff(timestamps)
  
  # Calculate weights (inverse of time intervals)
  weights = 1.0 / time_intervals
  
  # Same calculation
  # weighted_mean_np = np.average(measurements[:-1],weights=weights)
  # print(weighted_mean_np) 
  
  # Calculate the weighted sum of measurements
  weighted_sum = np.sum(measurements[:-1] * weights)
  
  # Calculate the total weight
  total_weight = np.sum(weights)
  
  # Calculate the weighted mean
  weighted_mean = weighted_sum / total_weight
  
  print("Weighted Mean:", weighted_mean)
  
  return weighted_mean


#%% Single Object Analysis plot
# mbol_fig, mbol_ax = plt.subplots(1,len(star_epoch_photometry),
mbol_fig, mbol_ax = plt.subplots(2,2,
                                 figsize = (fig_width,fig_height*0.6),
                                 sharex=(True),
                                 sharey=True)
if ready_to_publish is False:
  mbol_fig.suptitle('Single Object Analysis')
        

# mbol_ax[0].set_ylabel('$ M [mag] $')  
mbol_ax[0,0].set_ylabel(r'$ \Delta $[mag]')  
mbol_ax[1,0].set_ylabel(r'$ \Delta $[mag]')  
# mbol_ax[1,0].set_xlabel('$ Time $')
# mbol_ax[1,1].set_xlabel('$ Time $')
# mbol_ax.set_xlabel('$ Time $')
mbol_ax[0,0].invert_yaxis()

# mbol_ax_stars = [mbol_ax[0,0],mbol_ax[0,1],mbol_ax[1,0],mbol_ax[1,1]]

# Gaia order
# 3 4 1 2
mbol_ax_stars = [mbol_ax[1,0],mbol_ax[1,1],mbol_ax[0,0],mbol_ax[0,1]]

# Viaux order to Gaia 1,   2   ,  3 ,   4
# viaux_markers = [' (2)',' (3) ','  ',' (1)']
# 2,3,x,1

# convert to 
# 3 4 1 2
# 0 1 2 3

viaux_markers = [r' ',r' (V13-1) ',r' (V13-2) ',r' (V13-3)']
label_viaux ='V13'

# gaia_markers = [r' (G3-3)',r' (G3-4)',r' (G3-1)',r' (G3-2)']
gaia_markers = [r' (3)',r' (4)',r' (1)',r' (2)']
label_gaia ='G3'


markers_stars = gaia_markers

a = np.arange(-.5,.6,0.1,dtype = float)


from astropy.time import Time
import datetime
# times = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']
# t = Time(times, format='isot', scale='utc')
# t = Time(times, format='jd', scale='tcb')


for i, star in enumerate(star_epoch_photometry): 
 
 
  star_id = star['source_id'][0]
  print()
  print(star_id)
  print(file_paths[i])
  print()
  mask_quality = (star['rejected_by_photometry'] == False) & \
                 (star['rejected_by_variability'] == False)
  
  mask_G =   (star['band'] == 'G')  & mask_quality
  mask_RP =  (star['band'] == 'RP') & mask_quality
  mask_BP =  (star['band'] == 'BP') & mask_quality
  
  
  g = star['mag'][mask_G]
  # g = g - np.mean(g)
  print(g.size)
  g_time = star['time'][mask_G]
  # g_time = Time(g_time, format='jd', scale='tcb')
  # g_time.to_datetime()
  
  g_w_mean = weighted_mean_1(g_time,g)
  
  
  #  source for OBMT time format in Gaia:
  # https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_photometry/ssec_dm_epoch_photometry.html
  offset_days = 2455197.5
  # unit = days, format = julian date (jd), scale = barycentric coordinate time (franc tcb)
  g_time = Time(offset_days + g_time, format='jd', scale='tcb')
  g_time = g_time.to_datetime()
  
  
  same_day_measurement_filter = False
  if same_day_measurement_filter is True:
    g_time_dif =g_time[1:] -g_time[:-1]
    
    mask_same_time = (g_time_dif < datetime.timedelta(days=0.5))
    mask_same_time = np.append(mask_same_time,False)
    mask_unic_time = np.logical_not(mask_same_time)
    
    g_time = g_time[mask_unic_time] 
    g = g[mask_unic_time] 
  
  
  
  
  
  
  rp = star['mag'][mask_RP]
  rp_time = star['time'][mask_RP]
  rp_w_mean = weighted_mean_1(rp_time,rp)
  # rp_time = Time(offset_days + rp_time, format='jd', scale='tcb')
  # rp_time = rp_time.to_datetime()
  
  bp = star['mag'][mask_BP]
  bp_time = star['time'][mask_BP]
  bp_w_mean = weighted_mean_1(bp_time,bp)
  # bp_time = Time(offset_days + bp_time, format='jd', scale='tcb')
  # bp_time = bp_time.to_datetime()
  
  # Structure rp and bp measuerments by transit ------------------------
  # Note: store only rp bp measuerement pairs to calc bp-rp
  rp_star = star[mask_RP]
  bp_star = star[mask_BP]
  
  # Base array is rp_star
  mask_RP_BP_transite = np.isin(rp_star['transit_id'], bp_star['transit_id'])
  # rp_star[mask_RP_BP_transite]
  
  # Base array is bp_star
  mask_BP_RP_transite = np.isin(bp_star['transit_id'], rp_star['transit_id'])
  # bp_star[mask_BP_RP_transite]
  
  # Test structure --------------------------------------------
  if False:
    test_array_entry_same_transit_id = \
      np.all(bp_star['transit_id'][mask_BP_RP_transite] == \
             rp_star['transit_id'][mask_RP_BP_transite])
   
    print('Test: Transit IDs per entry in both data sets equal')
    print('Result: '+str(test_array_entry_same_transit_id))
   
  # Derive bp-rp and adjoined time ----------------------------------------
  # Get transite times of pairs
  # bp_rp_time = bp_star['time'][mask_BP_RP_transite]
  bp_rp_time = bp_time[mask_BP_RP_transite]
  
  # Calc bp-rp for transite pairs
  bp_rp      = bp_star['mag'][mask_BP_RP_transite] - \
               rp_star['mag'][mask_RP_BP_transite]
  
  bp_rp_w_mean = weighted_mean_1(bp_rp_time,bp_rp)
  bp_rp_time = Time(offset_days + bp_rp_time, format='jd', scale='tcb')
  bp_rp_time = bp_rp_time.to_datetime()
  # bp_rp = bp_rp - np.mean(bp_rp)
  
  # plot data -------------------------------------------------------------
  
  # g_mean_dev = g- np.mean(g)
  g_mean_dev = g- g_w_mean
  bp_rp_mean_dev = bp_rp - np.mean(bp_rp)
  bp_rp_w_dev = (bp_rp.max() - bp_rp.min())/2
  bp_rp_w_dev = np.abs(bp_rp_mean_dev).max()
  
  g_w_dev = (g.max()-g.min())/2
  g_w_dev = np.abs(g_mean_dev).max()
  # a = np.arange(g_mean_dev.min()-0.1,g_mean_dev.max()+0.1,0.1,dtype = float)
  
  # label_g = r'$ m_{g} $'+'   = '+r'({0:5.2f}  $ \pm $ {1:5.2f}) mag'.format(np.mean(g),np.std(g))
  # label_bp_rp = r' bp-rp ='+r'({0:1.2f}  $ \pm $ {1:1.2f}) mag'.format(np.mean(bp_rp),np.std(bp_rp))
  label_g = r'$ m_{g} $'+'   = '+r'({0:5.2f}  $ \pm $ {1:5.2f}) mag'.format(g_w_mean,g_w_dev)
  label_bp_rp = r' bp-rp ='+r'({0:1.2f}  $ \pm $ {1:1.2f}) mag'.format(bp_rp_w_mean,bp_rp_w_dev)
  
  
  mbol_ax_stars[i].axhline(0, color=t_fit,linewidth=0.4,linestyle='-',
                           # label= star_id
                           )
  
  for el in a: 
    mbol_ax_stars[i].axhline(el, color=t01,linewidth=0.2,linestyle='-') 
    
  
  # bp_rp_time.plot_date
  
  mbol_ax_stars[i].plot(bp_rp_time,bp_rp_mean_dev,
                  # label = label_bp_rp,
                  color=t01,
                  linestyle ='-',
                  linewidth = 0.5,
                  marker='o',
                  markersize=2
                  ) 
  
  mbol_ax_stars[i].plot(g_time,g_mean_dev,
                  # label = label_g,
                  color=t_lost,
                  linestyle ='-',
                  linewidth = 0.5,
                  marker='+',
                  markersize=4)
  
  
  
  # mbol_ax_stars[i].legend(title=star_id)
  # mbol_ax_stars[i].set_title( 'DR3 '+ str(star_id)+markers_stars[i], loc='left',fontsize=9)  
  mbol_ax_stars[i].set_title( label_g +markers_stars[i]+'\n'+label_bp_rp, loc='left',fontsize=9)    
  
  
  # Rotate the x-axis labels for better readability (optional)
  # mbol_ax_stars[i].tick_params(axis='x', labelrotation=35)
  min_x = g_time.min()
  # min_x.day = 1
  max_x = g_time.max()
  x_tick_entries = np.arange(min_x, max_x, (max_x-min_x)/5)
  mbol_ax_stars[i].xaxis.set_ticks(x_tick_entries)
  mbol_ax_stars[i].set_xticklabels(x_tick_entries, rotation=45, horizontalalignment = 'right')
  # mbol_ax_stars[i].locator_params(nbins=4) 
  # Format the x-axis as date
  plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d'))

if savefig is True:
  mbol_fig.savefig(fig_sub_path+'SingleObjectAnalysis')


