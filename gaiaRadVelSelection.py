#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 13:21:47 2023

@author: Thomas Suntinger

Input:
  - load: Gaia query data set as csv file (e.g. from GAIA_query_circle.py)

Dependencies:
  - numpy
  - csv
  - GAIA_load.py
  - GAIA_save.py
  - curve_of_growth.py
  
Dependencies for visualization:
  - matplotlib
  - plotStyle.py
  - GAIA_plot.py
  
Functionality:
  - Using a Gaussian Kernel Density Estimator (KDE) included in 
    the propability_density_estimation() function included in 
    curve_of_growth.py, a probability density estimation distribution in 
    1D phase space (radial velocity in km/s) and an adjoined curve of growth
    is created. The coresponding isopleth to the highest probability density 
    encloses the radial velocity cluster member selection.
  - Since radial velocity measurements are rare, csv tables of all sets can be
    stored (see 'Output' for more details)
  - A variety of plots can be made and saved, if adjoined parameters 
    are set to True. (see 'Output plots' for more details)

Structure ( & parameters):
  - Init 
  - Plot settings 
      (figure height, default plot colors)  
  - Load /Store data settings 
      (path-, file- & data set names)
  - Toggle plots on / off
  - Store figures settings
  
  - Radial velocity cluster member selection
    - First iteration and curve of growth function call
    - Second iteration unfinished
  - Store data of single iteration
  
  - Histogram
  - CMD  
  - Histogram + CMD comparison  
  
Output:
  - store: csv Gaia_M5_r-9arcmin_radVel-center.csv  (radial velocity center in km/s)
  - store: csv Gaia_M5_r-9arcmin_excluded_radVel (stars without rad vel measurement)
  - store: csv Gaia_M5_r-9arcmin_valid_radVel    (stars with    rad vel measurement)
  - store: csv Gaia_M5_r-9arcmin_radVel_iso_1    (rad vel  cluster member selection)
  
  - print: number of stars in data sample
  - print: number of stars without rad vel
  - print: number of stars with    rad vel
  - print: number of stars with    rad vel enclosed by isopleth 1.0 prop.max 
  - print: number of stars with    rad vel outside     isopleth 1.0 prop.max 
  - print: rad vel KDE peak in km/s

Output plots:
  - radial velocity 1D KDE curve of growth plot
  - radial velocity histogram and 1D KDE function plot
  - radial velocity CMD  
  - radial velocity histogram + CMD comparison
 
Improvements:
  - Clean unnecessary comments and lines of code
  - Write documentation and improve in code documentation.
  - Bonus: Use 'Structure:Second Iteration' to create iterative selection 
      (see as example proper_motion_evaluation_loop.py)

"""

#%% Init
import csv
import numpy as np

import matplotlib.pyplot as plt
import foundation.GAIA_plot as gp

#%% ### ------------ Set default plot style------------------------------------
#https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure

# Note: fig_width is fixed in plotStyle.py for MasterThesis specifications
#       fig_height can and should be adjusted to special needs of individual plots
#       fig_height of 5inch is approximatelly half of height on A4 MasterThesis      
fig_height = 5
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)



### ------------ Set default plot colors --------------------------------------

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

# rad velocity
t_radVel = 'blue'
# dataDR3 - Starhorse convolution
t07 = 'orangered'

# Av extinction models
t08 = 'lightgrey'

# RGB
t09 = 'teal'
# AGB
t10 = 'purple'

#%% ### ------------ Load /Store data settings --------------------------------

folder_name = 'data/queryData/'

## -------------- GAIA data ---------------------------------------------------
# data_set_selected = 'Gaia_M3_r-18arcmin.csv'
data_set_selected = 'Gaia_M5_r-9arcmin.csv'
# data_set_selected = 'Gaia_M67_r-25arcmin.csv'

file_name = data_set_selected[:-4]

## ------------ load data into variables --------------------------------------

from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

# gaia_dr3_array = file_data
# label_gaia_dr3_array = 'gaia dr3 data set'
# label_gaia_dr3_array_n = 'gaia dr3 data set ('+str(gaia_dr3_array.size)+' stars)'

gaia_dr3_array = file_data
label_gaia_dr3_array = 'data sample'
label_gaia_dr3_array_n = 'data sample ('+str(gaia_dr3_array.size)+' stars)'
print('number of stars in data sample:')
print(gaia_dr3_array.size)

#---------Store GAIA proper motion data --------------------------------

file_sub_path = 'data/radVelData/'
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### -------------Toggle plots on / off  -----------------------------------

radVel_curve_of_growth_plot       = True
radVel_sel_plot_hist_cmd_on       = True

radVel_single_plots               = False



if radVel_single_plots is True:
  radVel_hist_KDE_plot            = True
  radVel_sel_cmd_plot             = True
else:
  radVel_hist_KDE_plot            = False
  radVel_sel_cmd_plot             = False



#%% ### ----------- Store figures settings ------------------------------------

savefig = True
fig_sub_path = 'gfx/rad_vel_eval/'


# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)





#%% ### rad vel cluster member selection ----------------- --------------------

#-------------------------- curve of growth .------------------------  
# import curve of growth  
from curve_of_growth import probability_KDE_1D as p_KDE_1D

# prepare sample without rad vel measurement
mask_rad_vel_is_NaN = np.isnan(gaia_dr3_array['radial_velocity'] )
rad_vel_is_NaN_array = gaia_dr3_array[mask_rad_vel_is_NaN]
label_rad_vel_is_NaN_array = 'no RV'
label_rad_vel_is_NaN_array_n = label_rad_vel_is_NaN_array+' ('+str(rad_vel_is_NaN_array.size)+' stars)'
# print(label_rad_vel_is_NaN_array_n)
print('number of stars without rad vel: ')
print(rad_vel_is_NaN_array.size)



# prepare sample with rad vel measurement
mask_rad_vel_no_NaN = np.logical_not(mask_rad_vel_is_NaN)
rad_vel_no_NaN_array = gaia_dr3_array[mask_rad_vel_no_NaN]
label_rad_vel_no_NaN_array = 'RV sample'
label_rad_vel_no_NaN_array_N = label_rad_vel_no_NaN_array+' ('+str(rad_vel_no_NaN_array.size)+' stars)'
#print(label_rad_vel_no_NaN_array_N)
print('number of stars with rad vel: ') 
print(rad_vel_no_NaN_array.size)


# Set input and parameters for curve of growth function
  
radVel_x = rad_vel_no_NaN_array['radial_velocity']
radVel_x_array = rad_vel_no_NaN_array

radVel_N_bins = 10000
radVel_N_levels = 2000
radVel_levels_scaling = 'lin'
fig_name = file_name + '_round1' # Added in function:  + '_curve of growth'
radVel_label_x_data = label_rad_vel_no_NaN_array
radVel_label_x_quantity = 'stars'
radVel_result_unit = 'km/s'


# call curve of growth function
[radVel_levels,radVel_density_sums,radVel_normed_density_sums, \
        radVel_mask_sample_in_levels_of_interesst, \
          radVel_levels_of_interesst,\
          radVel_contour_labels, \
            radVel_x_center,
            radVel_x_fit, radVel_prob_fit            
            ] = \
  p_KDE_1D(radVel_x, 
          # levels_of_interesst = levels_of_interesst, 
          N_bins = radVel_N_bins,
          levels = radVel_N_levels, 
          levels_scaling = radVel_levels_scaling,
          savefig = savefig,
          fig_name = fig_name,
          fig_sub_path = fig_sub_path,
          label_x_data = radVel_label_x_data,
          label_x_quantity = radVel_label_x_quantity,
          curve_of_growth_plot = radVel_curve_of_growth_plot
          )

# First and only iteration ----------------------------------------------------

# Rad vel selection with iso 1.0 prop. max
radVel_kde_1D = radVel_x[radVel_mask_sample_in_levels_of_interesst[-1]]
radVel_kde_1D_arry = rad_vel_no_NaN_array[radVel_mask_sample_in_levels_of_interesst[-1]]
label_radVel_kde_1D = radVel_contour_labels[-1]
#
gaia_dr3_array_radVel_iso_1       = radVel_kde_1D_arry
label_gaia_dr3_array_radVel_iso_1 = label_radVel_kde_1D
#
label_gaia_dr3_array_radVel_iso_1_N = label_gaia_dr3_array_radVel_iso_1+' ('+str(gaia_dr3_array_radVel_iso_1.size)+' stars)'
#print(label_gaia_dr3_array_radVel_iso_1_N)
print('number of stars with rad vel enclosed by isopleth 1.0 prop.max: ')
print(gaia_dr3_array_radVel_iso_1.size)

# residue between (valid-) and (enclosed by isopleth 1.0 prop max) rad vel 
equal_mask = np.isin(rad_vel_no_NaN_array['source_id'],gaia_dr3_array_radVel_iso_1['source_id'])
gaia_dr3_array_not_in_radVel_selection       = rad_vel_no_NaN_array[np.logical_not(equal_mask)]
#print('Residue Iso100 to RV sample ('+str(gaia_dr3_array_not_in_radVel_selection.size)+')')
print('number of stars with rad vel outside isopleth 1.0 prop.max: ')
print(gaia_dr3_array_not_in_radVel_selection.size)

# KDE peak of radial velocity distribution
label_radVel_kde_1D_val = '{0:s}:{1:10.3g}{2:s}'.format(r'KDE peak',radVel_x_center[0],radVel_result_unit)
print('rad vel '+label_radVel_kde_1D_val)

#%% Second iteration (unfinished) ---------------------------------------------
# Note: check 'first and only iteration' and 'Store data of single iteration'
#       for style, format, and naming (above lines outdated!) 

## adjust set of data (comparable to proper motion loop)--------
# fig_name = file_name + '_round2' #  '_curve of growth_' added in same named function!
# radVel_x = radVel_x[radVel_mask_sample_in_levels_of_interesst[0]]
# radVel_x_array = rad_vel_no_NaN_array[radVel_mask_sample_in_levels_of_interesst[0]]
# label_radVel_x_array = radVel_contour_labels[0] + 'derived from round 1'
# radVel_N_levels = 2000
# radVel_levels_scaling = 'lin'
# -------------------------------------------------------------

# ## Adjust levels from round 1 and use them (set of data unchanged) -----
# fig_name = file_name + '_round2_levelAdjust'
# # Note radVel_N_levels is input parameter levels and can also be array
# # -------------------------------------------------------------
# radVel_N_levels = radVel_levels[2:]
# radVel_levels_scaling = 'custom'
# # -------------------------------------------------------------
# # Note: did not yield nice results
# # add_levels = np.linspace(radVel_levels[1],radVel_levels[2],100)
# # radVel_N_levels = np.concatenate((radVel_levels[0:1],add_levels,radVel_levels[2:]))
# # -------------------------------------------------------------

# radVel_N_bins = 10000
# radVel_label_x_quantity = 'stars'
# radVel_result_unit = 'mas'

# [radVel_levels,radVel_density_sums,radVel_normed_density_sums, \
#         radVel_mask_sample_in_levels_of_interesst, \
#           radVel_levels_of_interesst,\
#           radVel_contour_labels, \
#             radVel_x_center,
#             radVel_x_fit, radVel_prob_fit
            
#             ] = \
#   p_KDE_1D(radVel_x, 
#           # levels_of_interesst = levels_of_interesst, 
#           N_bins = radVel_N_bins,
#           levels = radVel_N_levels, 
#           levels_scaling = radVel_levels_scaling,
#           savefig = savefig,
#           fig_name = fig_name,
#           fig_sub_path = fig_sub_path,
#           label_x_data = radVel_label_x_data,
#           label_x_quantity = radVel_label_x_quantity,
#           curve_of_growth_plot = radVel_curve_of_growth_plot
#           )

# label_radVel_kde_1D_val = '{0:s}:{1:10.3g}{2:s}'.format(r'1D KDE peak',radVel_x_center[0],radVel_result_unit)


# radVel_kde_1D = radVel_x[radVel_mask_sample_in_levels_of_interesst[-1]]
# radVel_kde_1D_arry = radVel_x_array[radVel_mask_sample_in_levels_of_interesst[-1]]
# # radVel_kde_1D_arry = radVel_nar_window_array[mask_sample_in_levels_of_interesst[-1]]
# label_radVel_kde_1D = radVel_contour_labels[-1]
#
# # Store data of 2 iterations ------------------------------------------------
# iteration_label = '_round 2_'
# radVel_eval_name = file_name +iteration_label+ 'radVel-center.csv'



#%% Store data of single iteration---------------------------------------------

iteration_label = '_'
radVel_eval_name = file_name +iteration_label+ 'radVel-center.csv'

# create new csv file and store radial cluster velocity in km/s
with open(file_sub_path+radVel_eval_name, 'w', newline='') as csvfile:
  # create new csv object
  csv_write = csv.writer(csvfile, delimiter=',')
  # store data in new csv file
  csv_write.writerows(label_radVel_kde_1D_val)
  
# store csv data sample with valid radial velocity
rad_vel_no_NaN_array_name = '_valid_radVel'
save_gaia(rad_vel_no_NaN_array, file_name + rad_vel_no_NaN_array_name,file_sub_path=file_sub_path)

# store csv data sample with rad vel enclosed by isopleth with 1.0 prob.max
gaia_dr3_array_radVel_iso_1_name = '_radVel_iso_1'
save_gaia(gaia_dr3_array_radVel_iso_1,file_name +gaia_dr3_array_radVel_iso_1_name,file_sub_path=file_sub_path)

# store residue between (valid-) and (enclosed by isopleth 1.0 prop max) rad vel 
gaia_dr3_array_not_in_radVel_selection_name = '_excluded_radVel'
save_gaia(gaia_dr3_array_not_in_radVel_selection, file_name + gaia_dr3_array_not_in_radVel_selection_name,file_sub_path=file_sub_path)

# If multiple iterations: store csv radial velocity data selection 
# radVel_selection_array      = gaia_dr3_array_radVel_iso_1
# radVel_selection_array_name = file_name +iteration_label+ 'radVel_iso_1'
# save_gaia(radVel_selection_array,radVel_selection_array_name,file_sub_path=file_sub_path)



#%% Histogram
#-----------------------------------------------------------------------------

if radVel_hist_KDE_plot is True:

  radVel_N_bins = 100 
  
  hist_radVel_sel_fig,hist_radVel_sel_ax = plt.subplots(1,1, sharey=True)
  
  radVel_n,radVel_bins,radVel_patches = hist_radVel_sel_ax.hist(radVel_x,color = t00, bins = radVel_N_bins, label=radVel_label_x_data)
  
  hist_radVel_sel_ax.hist(radVel_kde_1D,color = t_radVel, bins = radVel_bins,label=label_radVel_kde_1D)
  
  hist_radVel_sel_ax.axvline(radVel_x_center,
                        label = label_radVel_kde_1D_val,\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '-' \
                            )
  
  hist_radVel_sel_ax2 = hist_radVel_sel_ax.twinx()
  hist_radVel_sel_ax2.set_ylabel('probability')
  hist_radVel_sel_ax2.plot(radVel_x_fit, radVel_prob_fit,label='1D KDE',color = 'k')
  
  
  hist_radVel_sel_ax.set_xlabel('radial velocity [km/s]')
  hist_radVel_sel_ax.set_ylabel('count of stars [N]')
  hist_radVel_sel_ax.legend(loc='upper left')
  
  if savefig is True:
    hist_radVel_sel_fig.savefig(fig_sub_path+file_name+'_Hist_radVel_1D_KDE')

#%% CMD 
  
if radVel_sel_cmd_plot is True:
    list_comparison = []
    comparison_labels = [radVel_label_x_data,label_radVel_kde_1D]
    comparison_colors = [t02,t_radVel]
    comparison_markers = ['']
    
    list_comparison.append(rad_vel_no_NaN_array)
    list_comparison.append(radVel_kde_1D_arry)
   
    
    
    
    # Note: cmd plot savefig in function!!!
    [cmd_fig,cmd_ax,all_in_one] = gp.cmd_plot_full(list_comparison, \
                labels=comparison_labels, 
                 colors=comparison_colors, \
                    # data_background = file_data,\
                    # color_data_background = t00,\
                    # label_data_background = label_file_data,\
                    super_title=file_name+'_CMD_radVel_KDE_1D',
                     markers = comparison_markers,
                      all_in_one=True,
                        ready_to_publish=True,
                          savefig = savefig,
                            fig_sub_path = fig_sub_path
                          )



#%% Histogram + CMD comparison

if radVel_sel_plot_hist_cmd_on is True:
  rgb_agb_lower_bound = 12.3
  mask = (gaia_dr3_array_not_in_radVel_selection['phot_g_mean_mag'] < rgb_agb_lower_bound )
  gaia_dr3_array_not_in_radVel_selection_trgb = gaia_dr3_array_not_in_radVel_selection[mask]
  
  mask = (rad_vel_is_NaN_array['phot_g_mean_mag'] < rgb_agb_lower_bound )
  rad_vel_is_NaN_array_trgb = rad_vel_is_NaN_array[mask]
  
  
  marker_size = 0.4
  marker_alpha = 0.4
  
  marker_size_lost = 30
  marker_alpha_lost = 1
  marker = '+'
  
  radVel_N_bins = 100
  
  
  
  pos_fig, pos_ax = plt.subplots(1,2,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[9,8]},sharex=False,sharey=False)
  # gs = gridspec.GridSpec(1, 2,width_ratios=[2,1])
  
  
  hist_radVel_sel_ax = pos_ax[0]
  
  radVel_n,radVel_bins,radVel_patches = hist_radVel_sel_ax.hist(radVel_x,color = t00, bins = radVel_N_bins, label=radVel_label_x_data)
  
  hist_radVel_sel_ax.hist(radVel_kde_1D,color = t_radVel, alpha = 0.4, bins = radVel_bins,label=label_radVel_kde_1D)
  
  hist_radVel_sel_ax.axvline(radVel_x_center,
                        label = label_radVel_kde_1D_val,\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '-' \
                            )
  delta = 60
  hist_radVel_sel_ax.set_xlim(radVel_x_center-delta,radVel_x_center+delta)
  hist_radVel_sel_ax.set_xlabel('radial velocity [km/s]')
  hist_radVel_sel_ax.set_ylabel('count of stars [N]')
  # hist_radVel_sel_ax.legend(loc='upper left')
  hist_radVel_sel_ax.legend()
  
  
  # hist_radVel_sel_ax = hist_radVel_sel_ax2.twinx()
  hist_radVel_sel_ax2 = hist_radVel_sel_ax.twinx()
  
  
  
  hist_radVel_sel_ax2.set_ylabel('probability')
  hist_radVel_sel_ax2.plot(radVel_x_fit, radVel_prob_fit,label='1D KDE',color = 'k')
  
  hist_radVel_sel_ax2.set_ylim(0,radVel_prob_fit.max()+0.0003)
  
  
  # hist_radVel_sel_ax2.legend()
  # hist_radVel_sel_ax2.set_aspect('equal', adjustable='box')
  
  
  
  
  
  
  pos_ax[1].scatter(rad_vel_is_NaN_array_trgb['bp_rp'],rad_vel_is_NaN_array_trgb['phot_g_mean_mag'],label=label_rad_vel_is_NaN_array,marker = '<',s=10,alpha=marker_alpha_lost,color=t02)  
  pos_ax[1].scatter(gaia_dr3_array_not_in_radVel_selection_trgb['bp_rp'],gaia_dr3_array_not_in_radVel_selection_trgb['phot_g_mean_mag'],label='lost '+label_rad_vel_no_NaN_array,marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t02)  
  
  pos_ax[1].axhline(rgb_agb_lower_bound, label= 'bound: '+str(rgb_agb_lower_bound) +' mag', linestyle = 'dashed',  color = 'k',linewidth=0.5, zorder = 1.5)  
    
  pos_ax[1].scatter(gaia_dr3_array_radVel_iso_1['bp_rp'],gaia_dr3_array_radVel_iso_1['phot_g_mean_mag'],label= label_gaia_dr3_array_radVel_iso_1, s=marker_size,alpha=marker_alpha,color=t_radVel)  
  
  pos_ax[1].set_xlabel('BP-RP [mag]')
  pos_ax[1].set_ylabel('phot_g_mean_mag [mag]')
  pos_ax[1].invert_yaxis()
  pos_ax[1].legend(loc = 'lower center')
  
    
  # pos_ax[0,i].set_aspect('equal', adjustable='box')
  # pos_ax[1].set_aspect('equal')
  
  fig_name = file_name
  
  if savefig is True:
    pos_fig.savefig(fig_sub_path+fig_name+'_radVel_cmd_comp.jpg')


