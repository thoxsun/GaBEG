# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 13:21:47 2023

@author: Thomas Suntinger

Input:
  - load: Gaia query data set as csv file (e.g. from GAIA_query_circle.py)
  - type: search cone radius of loaded Gaia query

Dependencies:
  - numpy
  - csv
  - zero_point (Parallax zero point recipe for (e)DR3. From Lindegren et al. 2020)
  - GAIA_load.py
  - GAIA_save.py
  - curve_of_growth.py
  - calcSearchRadius.py
  
Dependencies for visualization:
  - matplotlib
  - plotStyle.py
  - GAIA_plot.py
  
Functionality:
  - Gaia parallax measurements need to be zero point corrected (zpc). The algorithm
    provided by Lindegren et al. 2020 is used for zpc of the Gaia (e)DR3 data.
  - Using a Gaussian Kernel Density Estimator (KDE) included in 
    the propability_density_estimation() function included in 
    curve_of_growth.py, a probability density estimation distribution in 
    1D phase space (parallaxes in mas) and an adjoined curve of growth
    is created. The coresponding isopleth to the highest probability density 
    encloses the parallax cluster member selection.
  - A variety of plots can be made and saved, if adjoined parameters 
    are set to True. (see 'Output plots' for more details)

Structure ( & parameters):
  - Init 
  - Plot settings 
      (figure height, default plot colors)  
  - Load /Store data settings 
      (path-, file- & data set names, search cone radius)
  - Toggle plots on / off
  - Store figures settings
  - Data preparation (data set label, parallax quality filters)
  - zero point correction for parallax (par_zpc)
  
  - parallax cluster member selection (par_zpc_sel)
    - curve of growth 1
    - curve of growth 2 level adjusted
    
  - Search radius converted to parallax 
  - Store data (of second iteration)
  
  - Histogram
  - CMD  
  - Histogram + CMD comparison  
  
Output:
  - print: zpc mean 
  - print: zpc std  
  - print: Number of lost stars due to zero point correction
  - print: 1D KDE parallax peak in mas
  - print: lower and upper parallax search cone bound in mas 
  
  - store: csv Gaia_M5_r-9arcmin_par-center.csv  (parallax center in mas)
  
  - store: csv [file_name]_zpc_invalid (data sample with invalid par zpc)
  - store: csv [file_name]_zpc         (data sample with par zpc)
  - store: csv [file_name]_round2_par_iso_1   
           (data sample with par enclosed by isopleth with 1.0 prob.max)
  - print: number of lost stars in parallax selection
  
Output plots:
  - radial velocity 1D KDE curve of growth plot
  - radial velocity histogram and 1D KDE function plot
  - radial velocity CMD  
  - radial velocity histogram + CMD comparison
 
Improvements:
  - Clean unnecessary comments and lines of code
  - Write documentation and improve in code documentation.
    based on Gaußian KDE (see also gaiaDistEval.py)
  - Bonus: Change iterative selection in 'Structure:curve of growth' 
      (see as example proper_motion_evaluation_loop.py)

"""

#%% Init

import csv
import numpy as np
from zero_point import zpt
import calcSearchRadius

import matplotlib.pyplot as plt
import GAIA_plot as gp

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

# dataDR3 - Starhorse convolution
t07 = 'orangered'

# Av extinction models
t08 = 'lightgrey'

# RGB
t09 = 'teal'
# AGB
t10 = 'purple'


#%% ### ------------ Load /Store data settings ---------------------------------------

folder_name = 'data/queryData/'

## -------------- GAIA data ---------------------------------------------------
# data_set_selected = 'Gaia_M3_r-18arcmin.csv'
data_set_selected = 'Gaia_M5_r-9arcmin.csv'
# data_set_selected = 'Gaia_M67_r-25arcmin.csv'

#Search radius angle in arcmin
SR_angle = 9
SR_angle = calcSearchRadius.convertArcMinToRad(SR_angle)

file_name = data_set_selected[:-4]

## ------------ load data into variables --------------------------------------

from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

gaia_dr3_array = file_data
label_gaia_dr3_array = 'gaia dr3 data set'
label_gaia_dr3_array_n = 'gaia dr3 data set ('+str(gaia_dr3_array.size)+' stars)'



#---------Store GAIA proper motion data --------------------------------

file_sub_path = 'data/parData/'
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### -------------Toggle plots on / off  -----------------------------------

par_curve_of_growth_plot          = True
par_sel_plot_hist_cmd_on          = True

par_single_plots                  = True


if par_single_plots is True:
  par_hist_KDE_plot               = True
  par_sel_cmd_plot                = True
else:
  par_hist_KDE_plot               = False
  par_sel_cmd_plot                = False
  
  

#%% ### ----------- Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/par_sel/'


# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)



#%% ### ------------- Data preparation ----------------------------------------


gaia_dr3_array = file_data


### Testing quality of parallax data

# Filter stars with positiv parallax 
# mask_positive_parallax = (gaia_dr3_array['parallax'] > 0)

# Filter stars with parallax error < 20% of respectiv parallax
# mask_parallax_error_selection = \
#   np.less( np.abs( \
#                   gaia_dr3_array['parallax_error'] /    \
#                   gaia_dr3_array['parallax']) , 0.2 )   
# Filter stars with parallax error < 0.2mas
# mask_parallax_error_selection = np.less(gaia_dr3_array['parallax_error'], 0.1 )

# Choose above masks to select parallax data set
# mask_parallax_selection = mask_positive_parallax  & mask_parallax_range & mask_parallax_error_selection
# mask_parallax_selection = mask_positive_parallax  & mask_parallax_error_selection
# mask_parallax_selection = mask_positive_parallax  & mask_parallax_range 
# mask_parallax_selection = mask_positive_parallax  
# mask_parallax_selection = mask_parallax_error_selection
# mask_parallax_selection = mask_parallax_range 

# gaia_dr3_array = gaia_dr3_array[mask_parallax_selection] 



label_gaia_dr3_array = 'data sample'
label_gaia_dr3_array_n = 'data sample ('+str(gaia_dr3_array.size)+' stars)'

#%% ### zero point correction for parallax (par_zpc) --------------------------
#Source: https://gitlab.com/icc-ub/public/gaiadr3_zeropoint
# NOTE: if any of the inputs values is NaN, the output will be NaN. Also, if the warnings are turned off and the
# source probided falls outside the valid range specified in describtion of get_zpt(),
# its zero-point will be NaN! 
# Note: Parallax error is relative value and therefore without zpc
# Note: For future use, NaN values have to be filtered after the procedure 
# Note: Check, if filtered stars are of interesst for you with cmd plot. 

zpt.load_tables()
# array of corrected values
par_zpc = zpt.get_zpt(gaia_dr3_array['phot_g_mean_mag'], \
                      gaia_dr3_array['nu_eff_used_in_astrometry'], \
                        gaia_dr3_array['pseudocolour'], \
                          gaia_dr3_array['ecl_lat'], \
                            gaia_dr3_array['astrometric_params_solved'],\
                              _warnings=False)


  
# Calculate corrected parallax
parallax_data_zpc = gaia_dr3_array['parallax'] - par_zpc
# Create copy of original data and replace parallax with zpc parallax
gaia_dr3_array_par_zpc = gaia_dr3_array.copy()
gaia_dr3_array_par_zpc['parallax'] =  parallax_data_zpc
# label_gaia_dr3_array_par_zpc  = str(gaia_dr3_array_par_zpc.size)+' zpc stars'

# Filter NaN values from corrected parallax (NaN = par_zpc not valid)
# Note: Following arrays contain less stars compared to original array!
#       Lost stars are printed beneath
mask_non_nan_after_correction = np.isfinite(parallax_data_zpc)
# Derive mean and st.dev from zpc
par_zpc_no_nan = par_zpc[mask_non_nan_after_correction]
print('zpc mean = '+ str(par_zpc_no_nan.mean()))
print('zpc std  = '+ str(par_zpc_no_nan.std()))

parallax_data_zpc_no_nan = parallax_data_zpc[mask_non_nan_after_correction]
gaia_dr3_array_par_zpc_no_nan = gaia_dr3_array_par_zpc[mask_non_nan_after_correction] 
label_gaia_dr3_array_par_zpc_no_nan    = 'M5 r9 zpc sample'
label_gaia_dr3_array_par_zpc_no_nan_N  = str(gaia_dr3_array_par_zpc.size)+' M5 r9 zpc sample'

mask_nan_after_correction = np.invert(mask_non_nan_after_correction)
filtered_star_array = gaia_dr3_array_par_zpc[mask_nan_after_correction]
label_filtered_star_array = 'zpc invalid'

print('')
print('Lost stars in zero point correction: '+str(filtered_star_array.size))
print('')


#%% ### parallax cluster member selection (par_zpc_sel) -----------------
  
#-------------------------- curve of growth .------------------------

  
from curve_of_growth import probability_KDE_1D as p_KDE_1D


par_x = gaia_dr3_array_par_zpc_no_nan['parallax']

# Note: Include iteration_label in following code?
# iteration_label = '_round 1_'

#-------------------------- curve of growth 1 ------------------------
par_N_bins = 10000
par_N_levels = 2000
par_levels_scaling = 'lin'
fig_name = file_name + '_round1' # Added in function:  + '_curve of growth'
par_label_x_data = label_gaia_dr3_array_par_zpc_no_nan
par_label_x_quantity = 'stars'
par_result_unit = 'mas'

[par_levels,par_density_sums,par_normed_density_sums, \
        par_mask_sample_in_levels_of_interesst, \
          par_levels_of_interesst,\
          par_contour_labels, \
            par_x_center,
            par_x_fit, par_prob_fit
            
            ] = \
  p_KDE_1D(par_x, 
          # levels_of_interesst = levels_of_interesst, 
          N_bins = par_N_bins,
          levels = par_N_levels, 
          levels_scaling = par_levels_scaling,
          savefig = savefig,
          fig_name = fig_name,
          fig_sub_path = fig_sub_path,
          label_x_data = par_label_x_data,
          label_x_quantity = par_label_x_quantity,
          curve_of_growth_plot = par_curve_of_growth_plot
          )

label_par_kde_1D_val = '{0:s} at:{1:10.3g}{2:s}'.format(r'Gaußian KDE parallax peak',par_x_center[0],par_result_unit)

# For single iteration uncooment ----------------------
# par_kde_1D = par_x[par_mask_sample_in_levels_of_interesst[-1]]
# par_kde_1D_arry = gaia_dr3_array_par_zpc_no_nan[par_mask_sample_in_levels_of_interesst[-1]]
# # par_kde_1D_arry = par_nar_window_array[mask_sample_in_levels_of_interesst[-1]]
# label_par_kde_1D = par_contour_labels[-1]
# -----------------------------------------------------------

# For single iteration comment ----------------------
par_kde_1D = par_x[par_mask_sample_in_levels_of_interesst[-1]]
par_kde_1D_arry = gaia_dr3_array_par_zpc_no_nan[par_mask_sample_in_levels_of_interesst[-1]]
# par_kde_1D_arry = par_nar_window_array[mask_sample_in_levels_of_interesst[-1]]
label_par_kde_1D = par_contour_labels[-1]


#-------------------------- curve of growth 2 level adjusted ------------------
### round 2----------------------------------------------------------------

## adjust set of data (comparable to proper motion loop)--------
# fig_name = file_name + '_round2' #  '_curve of growth_' added in same named function!
# par_x = par_x[par_mask_sample_in_levels_of_interesst[0]]
# par_x_array = gaia_dr3_array_par_zpc_no_nan[par_mask_sample_in_levels_of_interesst[0]]
# label_par_x_array = par_contour_labels[0] + 'derived from round 1'
# par_N_levels = 2000
# par_levels_scaling = 'lin'
# -------------------------------------------------------------

## Adjust levels from round 1 and use them (set of data unchanged) -----
par_x =       gaia_dr3_array_par_zpc_no_nan['parallax']
par_x_array = gaia_dr3_array_par_zpc_no_nan
par_label_x_data = label_gaia_dr3_array_par_zpc_no_nan
fig_name = file_name + '_round2_levelAdjust'
# Note par_N_levels is input parameter levels and can also be array
# -------------------------------------------------------------
par_N_levels = par_levels[2:]
par_levels_scaling = 'custom'
# -------------------------------------------------------------
# Note: did not yield nice results
# add_levels = np.linspace(par_levels[1],par_levels[2],100)
# par_N_levels = np.concatenate((par_levels[0:1],add_levels,par_levels[2:]))
# -------------------------------------------------------------

par_N_bins = 10000
par_label_x_quantity = 'stars'
par_result_unit = 'mas'

[par_levels,par_density_sums,par_normed_density_sums, \
        par_mask_sample_in_levels_of_interesst, \
          par_levels_of_interesst,\
          par_contour_labels, \
            par_x_center,
            par_x_fit, par_prob_fit
            
            ] = \
  p_KDE_1D(par_x, 
          # levels_of_interesst = levels_of_interesst, 
          N_bins = par_N_bins,
          levels = par_N_levels, 
          levels_scaling = par_levels_scaling,
          savefig = savefig,
          fig_name = fig_name,
          fig_sub_path = fig_sub_path,
          label_x_data = par_label_x_data,
          label_x_quantity = par_label_x_quantity,
          curve_of_growth_plot = par_curve_of_growth_plot
          )

label_par_kde_1D_val = '{0:s}:{1:10.3g}{2:s}'.format(r'1D KDE peak',par_x_center[0],par_result_unit)


par_kde_1D = par_x[par_mask_sample_in_levels_of_interesst[-1]]
par_kde_1D_arry = par_x_array[par_mask_sample_in_levels_of_interesst[-1]]
# par_kde_1D_arry = par_nar_window_array[mask_sample_in_levels_of_interesst[-1]]
label_par_kde_1D = par_contour_labels[-1]



print_par = True
if print_par is True:
  print(label_par_kde_1D_val)
  # print('Above value was derived with rad.velocity data larger 19.0km/h')


#%% Search radius converted to distance

# parallax center in [mas] from gaiaParallaxSelection for M5 with 9arcmins
par_center = par_x_center[0]


# print('(1D KDE peak = par center) at: '+str(par_center))

# distance in kPc
d = 1/par_center

# search radius in kPc
SR = calcSearchRadius.calcSearchRadius(d, SR_angle)
print('search radius in kPc: '+str(SR))

d_min = d - SR
par_max = 1/d_min
 
delta_par_max = par_max - par_center


d_max = d + SR
par_min = 1/d_max

delta_par_min = par_center - par_min

par_bound_min_name = 'converted search cone radius to lower parallax bound yields: '+str(delta_par_min)+' mas'
par_bound_max_name = 'converted search cone radius to upper parallax bound yields: '+str(delta_par_max)+' mas'

par_bound_mask = ((par_x > par_min)&(par_x <  par_max))
par_bound_array = gaia_dr3_array_par_zpc_no_nan[par_bound_mask]

print(par_bound_min_name)
print(par_bound_max_name)






#%% Store data -----------------------------------------------
par_eval_name = file_name +'_round 2'+ '_par-center.csv'

par_store_list = [label_par_kde_1D_val,par_bound_min_name,par_bound_max_name]


# create new csv file
with open(file_sub_path+par_eval_name, 'w', newline='') as csvfile:
  # create new csv object
  csv_write = csv.writer(csvfile, delimiter=',')
  # store data in new csv file
  csv_write.writerows(par_store_list)
  
# data sample with invalid par zpc 
save_gaia(filtered_star_array, file_name +'_zpc_invalid' ,file_sub_path=file_sub_path)

# data sample with par zpc 
save_gaia(gaia_dr3_array_par_zpc_no_nan, file_name +'_zpc' ,file_sub_path=file_sub_path)

# data sample with par enclosed by isopleth with 1.0 prob.max
gaia_dr3_array_par_iso_1 =  par_kde_1D_arry
label_gaia_dr3_array_par_iso_1 = label_par_kde_1D

par_selection_array = gaia_dr3_array_par_iso_1
par_selection_array_name = file_name +'_round2_'+ 'par_iso_1'


save_gaia(gaia_dr3_array_par_iso_1,par_selection_array_name,file_sub_path=file_sub_path)


equal_mask = np.isin(gaia_dr3_array_par_zpc_no_nan['source_id'],gaia_dr3_array_par_iso_1['source_id'])
gaia_dr3_array_not_in_par_selection = gaia_dr3_array_par_zpc_no_nan[np.logical_not(equal_mask)]
print('gaia_dr3_array_not_in_par_selection: '+str(gaia_dr3_array_not_in_par_selection.size))


#%% Histogram
#-----------------------------------------------------------------------------

if par_hist_KDE_plot is True:

  par_N_bins = 100 
  
  hist_par_sel_fig,hist_par_sel_ax = plt.subplots(1,1, sharey=True)
  
  par_n,par_bins,par_patches = hist_par_sel_ax.hist(par_x,color = t00, bins = par_N_bins, label=par_label_x_data)
  
  hist_par_sel_ax.hist(par_kde_1D,color = t06, bins = par_bins,label=label_par_kde_1D)
  
  hist_par_sel_ax.axvline(par_x_center,
                        label = label_par_kde_1D_val,\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '-' \
                            )
  
  hist_par_sel_ax2 = hist_par_sel_ax.twinx()
  hist_par_sel_ax2.set_ylabel('probability')
  hist_par_sel_ax2.plot(par_x_fit, par_prob_fit,label='1D KDE',color = 'k')
  
  
  hist_par_sel_ax.set_xlabel('parallax [mas]')
  hist_par_sel_ax.set_ylabel('count of stars [N]')
  hist_par_sel_ax.legend(loc='upper left')
  
  if savefig is True:
    hist_par_sel_fig.savefig(fig_sub_path+file_name+'_Hist_par_1D_KDE')
  
#%% CMD plot
  
if par_sel_cmd_plot is True:
    list_comparison = []
    comparison_labels = [par_label_x_data,label_par_kde_1D]
    comparison_colors = [t02,t06]
    comparison_markers = ['']
    
    list_comparison.append(gaia_dr3_array_par_zpc_no_nan)
    list_comparison.append(par_kde_1D_arry)
   
    
    import GAIA_plot as gp
    
    # Note: cmd plot savefig in function!!!
    [cmd_fig,cmd_ax,all_in_one] = gp.cmd_plot_full(list_comparison, \
                labels=comparison_labels, 
                 colors=comparison_colors, \
                    # data_background = file_data,\
                    # color_data_background = t00,\
                    # label_data_background = label_file_data,\
                    super_title=file_name+'_CMD_par_KDE_1D',
                     markers = comparison_markers,
                      all_in_one=True,
                        ready_to_publish=True,
                          savefig = savefig,
                            fig_sub_path = fig_sub_path
                          )



#%% Histogram + CMD comparison

if par_sel_plot_hist_cmd_on is True:
  rgb_agb_lower_bound = 12.3
  mask = (gaia_dr3_array_not_in_par_selection['phot_g_mean_mag'] < rgb_agb_lower_bound )
  gaia_dr3_array_not_in_par_selection_trgb = gaia_dr3_array_not_in_par_selection[mask]
  
  marker_size = 0.4
  marker_alpha = 0.4
  
  marker_size_lost = 30
  marker_alpha_lost = 1
  marker = '+'
  
  par_N_bins = 1000 
  
  pos_fig, pos_ax = plt.subplots(1,2,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[9,8]},sharex=False,sharey=False)
  # gs = gridspec.GridSpec(1, 2,width_ratios=[2,1])
  
  hist_par_sel_ax2 = pos_ax[0]
  # hist_par_sel_ax = pos_ax[0]
  
  hist_par_sel_ax2.set_ylabel('probability')
  hist_par_sel_ax2.set_xlabel('parallax [mas/year]')
  hist_par_sel_ax2.plot(par_x_fit, par_prob_fit,label='1D KDE',color = 'k')
  
  hist_par_sel_ax2.set_ylim(0,par_prob_fit.max()+0.001)
  
  
  hist_par_sel_ax = hist_par_sel_ax2.twinx()
  # hist_par_sel_ax2 = hist_par_sel_ax.twinx()
  
  par_n,par_bins,par_patches = hist_par_sel_ax.hist(par_x,color = t00, bins = par_N_bins, label=par_label_x_data)
  
  hist_par_sel_ax.hist(par_kde_1D,color = t06, bins = par_bins,label=label_par_kde_1D)
  
  # hist_par_sel_ax.hist(par_bound_array['parallax'],color = t07, bins = par_bins,\
  #                      # label=par_bound_name
  #                      )
  
  
  hist_par_sel_ax.axvline(par_x_center,
                        label = label_par_kde_1D_val,\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '-' \
                            )
  
  hist_par_sel_ax.axvline(par_max,
                        # label = par_bound_name,\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '--' \
                            )
  
  hist_par_sel_ax.axvline(par_min,
                        # label = par_bound_name,\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '--' \
                            )
  
  # hist_par_sel_ax.scatter(par_x_center,par_prob_fit[(par_x_center == par_x_fit)],label='KDE max '+str(par_x_center),marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t02)    
  
  hist_par_sel_ax2.axhline(par_prob_fit[(par_x_center == par_x_fit)],
                          label='KDE max '+str(par_x_center),
                          color = 'k',\
                            linewidth=0.4, \
                              zorder = 4,\
                            linestyle = '-' \
                              )
    
  delta = 1.5
  hist_par_sel_ax.set_xlim(par_x_center-delta,par_x_center+delta)
  hist_par_sel_ax.set_xlabel('parallax [mas]')
  hist_par_sel_ax.set_ylabel('count of stars [N]')
  hist_par_sel_ax.legend(loc='upper left')
  # hist_par_sel_ax.legend()
  # hist_par_sel_ax2.legend()
  # hist_par_sel_ax2.set_aspect('equal', adjustable='box')
  
  
  
  
  
  
  mask = (filtered_star_array['phot_g_mean_mag'] < rgb_agb_lower_bound )
  pos_ax[1].scatter(filtered_star_array['bp_rp'][mask],filtered_star_array['phot_g_mean_mag'][mask],label='lost '+label_filtered_star_array+' stars',marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t_lost)  
  
  pos_ax[1].scatter(gaia_dr3_array_not_in_par_selection_trgb['bp_rp'],gaia_dr3_array_not_in_par_selection_trgb['phot_g_mean_mag'],label='lost '+label_gaia_dr3_array_par_zpc_no_nan,marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t02)  
  
  pos_ax[1].axhline(rgb_agb_lower_bound, linestyle = 'dashed', label= 'bound: '+str(rgb_agb_lower_bound)+' mag', color = 'k',linewidth=0.5, zorder = 1.5)  
  
  pos_ax[1].scatter(gaia_dr3_array_par_iso_1['bp_rp'],gaia_dr3_array_par_iso_1['phot_g_mean_mag'],label= label_gaia_dr3_array_par_iso_1, s=marker_size,alpha=marker_alpha,color=t06)
  
  # pos_ax[1].scatter(par_bound_array['bp_rp'],par_bound_array['phot_g_mean_mag'],label= '?', s=marker_size,alpha=marker_alpha,color=t07)
  
  pos_ax[1].set_xlabel('BP-RP [mag]')
  pos_ax[1].set_ylabel('phot_g_mean_mag [mag]')
  pos_ax[1].invert_yaxis()
  # pos_ax[1].legend(loc = 'lower center')
  
    
  # pos_ax[0,i].set_aspect('equal', adjustable='box')
  # pos_ax[1].set_aspect('equal')
  
  
  fig_sub_path = 'gfx/par_sel/'
  fig_name = data_set_selected
  
  if savefig is True:
    pos_fig.savefig(fig_sub_path+fig_name+'_par_cmd_comp.jpg')
  
  
  
  
  
 














