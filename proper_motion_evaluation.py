#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 21:35:28 2022

@author: Thomas Suntinger

Import recommendation:
  import proper_motion_evaluation as pme
  
Input:
  - load: Gaia query data set as csv file (e.g. from GAIA_query_circle.py)
  - load: Only for testing ([pmra,pmdec].npy selection of M5 candidates)

Dependencies:
  - numpy
  - scipy
  - matplotlib
  - csv
  - GAIA_load.py
  - GAIA_save.py
  
Dependencies for visualization:
  - plotStyle.py
  - GAIA_plot.py
  
Functionality:
  - propability_density_estimation() function estimates probability density
    distribution in 2D phase space (proper motion right ascension (pmra) and
    declination (pmdec)) with a Gaussian Kernel Density Estimator (KDE).
    
  - Can be executed as script, since in Structure:main, 
    proper_motion_evaluation() is called. Results can be stored in test folders.
  - To use this file as a script, one can also uncomment marked areas below.
    Additionally, function definition and return commant must be commented.
    This enables access to variable values in Spyder Variable Explorer 
    meanwhile testing. 
    
  - A variety of plots can be made and saved, if adjoined parameters 
    are set to True.

Structure ( & parameters):
  - Init 
  - Plot settings 
      (figure height)
  - Load /Store data settings 
      (path-, file- & data set names)
  - Test KDE with different data distributions
      (different data sets available to uncomment)
  - Default values for function and if used as script
      (resolution of meshgrid, 
       scaling and resolution of isopleths,
       isopleths of interesst)
  - Default proper motion evaluation plots
  
  - propability_density_estimation function
    - KDE (f_scaling scales smallest probability z != 0.0 for log scale)
    - Create isopleths (levels)
      - levels_spacing_plot
    - Calculate curve of growth
      - Derive isopleths / levels of interest
    - Derive pm center & contour plot of interest
    - proper_motion_surface_plot
      
  - main routine and function call 
      (isopleth of interesst labels and adjoined file names to store data as csv)

Output:
  - store: csv files of pm center values ([ra,dec] in mas)
  - print: pm center values ([ra,dec] in mas)
  - store: csv tables of stars enclosed by isopleths of interest 

Output plots:
  - KDE curve of growth 
  - KDE pm phase space (contour_plot_of_interesst)
  
  - levels_spacing_plot (shows all available isopleth scalings )             
  - contour_plot ( shows isopleths to chosen scaling on pm phase space )
  - proper_motion_surface_plot (3D plot with 2D pm phase space map 
                                and probability density on z-axis)

Improvements:
  - Clean unnecessary comments and lines of code
    - Plots in Script mode plotted twice -> change that! 
    - Get rid of N_Simbad
    - Change 'interesst' in grammatically correct 'interest'
  - Write documentation and improve in code documentation.
  
  - Check if extended area of interesst has influence on result and delete if not
  - Make f_scaling function parameter and / or automize determination
  - Create a generalized version of probability_density_estimation() to use it
    for two arbitrary parameters.
    
  - Include the generalized version in curve_of_growth.py as KDE_2D()
  

"""
# %% This is a cell. You can run its content using 
#       Shift+Enter (run cell and advance to the next one) 
#       Ctrl+Enter (run cell and remain on it).

#%% Init

import numpy as np
from scipy.stats import kde
import matplotlib as mlt
import matplotlib.pyplot as plt
import csv

#%% ### Plot settings -------------------------------------------

fig_height = 4
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)

savefig = True
fig_name = 'proper motion estimation'
fig_sub_path = 'gfx/pm_eval/test/'

# if savefig is True:
#   cmd_fig.savefig(fig_sub_path+file_name)

#%% ### ------------ Load /Store data settings ---------------------------------------

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

gaia_dr3_array = file_data
label_gaia_dr3_array = 'gaia dr3 data set'
label_gaia_dr3_array_n = 'gaia dr3 data set ('+str(gaia_dr3_array.size)+' stars)'



#---------Store GAIA proper motion data --------------------------------

file_sub_path = 'data/pmData/test/'
from GAIA_save import save_gaia as save_gaia
save_data = True


# %% Test KDE with different distributions


# x = np.array( [4.,2.,4.,2.] )
# y =  np.array( [2.,4.,4.,2.] )

# x = np.array( [0,  1  ,1.2  ,4.6  ,8.7  ,3.1])
# y = np.array([0,  2  ,4.3  ,2.8  ,1.7  ,7.6])

# x = np.array([0,2,1])
# y = np.array([0,1,2])


# Generating random numbers and random numbers from special distribution
# Source: https://numpy.org/doc/stable/reference/random/generator.html

# Generate BitGenerator with seed (random numbers reproducable!)
# rng = np.random.default_rng(seed=42)

# Use random() to generate 3x3 matrix with float values [0.0-1.0)
# x = rng.random((25),dtype=np.float64)
# y = rng.random((25),dtype=np.float64)



#Generate a normal distribution in 3D with mean at [0,0], covariant matrix [[1, 0.5], [0.5, 3]] and 2000 points 
# Source: https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.multivariate_normal.html#numpy.random.Generator.multivariate_normal
# data = rng.multivariate_normal([0, 0], [[1, 0.5], [0.5, 3]], 1000)
# x, y = data.T


# Test array derived from M5
# [x,y] = np.load('selected_pm_array.npy')
# save_data = False


# %% Default values for function and if used as script

[x,y] = gaia_dr3_array['pmra'],gaia_dr3_array['pmdec']


# For M5 = 6336
N_simbad = 6336


# Default values for functioin:
N_bins = 250

N_levels = 100
levels = N_levels

extention_value = 5
extended_area = True

levels_of_interesst = 0
# levels_of_interesst = [0.2,0.66,1.]

# default scaling is 'lin'
levels_scaling = 'log'

# Change, only for testing function
# levels_scaling = 'lin'


# %% Default proper motion evaluation plots

# Interessting for data evaluation
levels_spacing_plot             = True
contour_plot                    = True
proper_motion_surface_plot_on   = True

# Results for proper motion evaluation
curve_of_growth_plot            = True
contour_plot_of_interesst       = True





# %% propability_density_estimation

# as script
# if True:
  
# as function
def propability_density_estimation(x,y, 
          N_simbad                    = 0,
          levels_of_interesst         = levels_of_interesst, 
          N_bins                      = N_bins,
          levels                      = N_levels, 
          levels_scaling              = levels_scaling,
          ready_to_publish            = ready_to_publish,
          extended_area               = extended_area, 
          extention_value             = extention_value, 
          
          levels_spacing_plot         = levels_spacing_plot ,
          contour_plot                = contour_plot ,
            
          curve_of_growth_plot        = curve_of_growth_plot,
          contour_plot_of_interesst   = contour_plot_of_interesst,
          
          proper_motion_surface_plot_on = proper_motion_surface_plot_on,
          savefig                     = savefig,
          fig_name                    = fig_name,
          fig_sub_path                = fig_sub_path
          ):


  if extended_area == True:
    xx, yy = np.mgrid[x.min()-extention_value:x.max()+extention_value:N_bins*1j,\
                      y.min()-extention_value:y.max()+extention_value:N_bins*1j]
  else:
    # xx, yy, zz are 2D arrays 
    # 1j is complex number
    xx, yy = np.mgrid[x.min():x.max():N_bins*1j, y.min():y.max():N_bins*1j]
  
  #%% KDE
  
  # Area of one bin
  delta_x = xx[1,0]-xx[0,0]
  delta_y = yy[0,1]-yy[0,0]
  bin_A = delta_x * delta_y
  
  
  k = kde.gaussian_kde([x,y])
  z = k(np.vstack([xx.flatten(), yy.flatten()]))
  # multiply propability density with bin area to get propability 
  # to find data point in bin
  z = z * bin_A 
  # Scale data points to avoid float inaccuracy 
  # z = z*10**3 
  zz = z.reshape(xx.shape)
  
  
  
  z_min = z.min()
  z_max = z.max() 
  # For logarithmic scale, find smallest z excluding 0.0
  z_log_min = z[(z > 0.0)].min()
  z_log_min = np.log10(z_log_min) 
  # scaling factor of z_min, break down for 0.05  < f_scaling < 0.006
  # f_scaling = 0.006
  f_scaling = 0.018
  z_log_min = z_log_min*f_scaling
  z_max_log = np.log10(z_max)
  
    
  #%% Create isopleths
  # create linearly spaced levels
  (levels_lin,step) = np.linspace(z_min,z_max,num = N_levels+1,\
                                  endpoint=True, retstep = True)
  levels_lin = np.append(0.,levels_lin)
  
  # create log spaced levels
  levels_log_scale = np.logspace(z_log_min,z_max_log,num = N_levels+1,\
                                 endpoint=True)
  levels_log_scale = np.append(0.,levels_log_scale)
  
  # create inverted log spaced levels
  levels_invert_log_scale = z_max - levels_log_scale 
  levels_sortet_invert_log_scale= np.flip(levels_invert_log_scale)
  
  
  # levels_spacing_plot
  # Create plot to show scaling
  if levels_spacing_plot is True:
      
    levels_index = np.arange(levels_log_scale.size)
    
    fig_level,ax_level=plt.subplots()
    ax_level.scatter(levels_index,levels_lin,s=3,label='lin')
    ax_level.scatter(levels_index,levels_log_scale,s=3,label='log')
    ax_level.scatter(levels_index,levels_sortet_invert_log_scale,\
                     s=3,label='orderd invers log')  
    ax_level.legend()
    
    if savefig is True:
      fig_level.savefig(fig_sub_path+fig_name+'_scaling')
    
    
  # Get rid of double entries
  levels_lin, lin_mask        = np.unique(levels_lin,       return_index = True)
  levels_log_scale,log_mask   = np.unique(levels_log_scale, return_index = True)
  levels_sortet_invert_log_scale,inv_log_mask = \
    np.unique(levels_sortet_invert_log_scale, return_index = True)
  
  # Choose levels:
    # 'lin','log','log-1'    
  if levels_scaling == 'log':
    levels = levels_log_scale
  elif levels_scaling == 'loginv':
    levels = levels_sortet_invert_log_scale
  elif levels_scaling == 'lin':
    levels = levels_lin
  else:
    print('Invalid input, default is linear scaling')
    levels = levels_lin
  
  #%%  Create scatter plot including contours plot:  
  if contour_plot is True:
    fig_contour, ax_contour = plt.subplots()
    if ready_to_publish == False:
      ax_contour.set_title('Levels used to calculate curve of growth')
    ax_contour.scatter(x,y,s=3,alpha = 0.5)
    density_QuadContourSet = ax_contour.contour(xx, yy, zz, 
                                                 levels= levels , alpha = 1 )
    if savefig is True:
      fig_contour.savefig(fig_sub_path+fig_name+'_scaling_contours')
  
  #%% Calculate curve of growth
  
  z_sample = k([x,y]) * bin_A
  
  density_sums = np.zeros(levels.size)
  normed_density_sums = density_sums.copy()
  bin_sums = np.zeros(levels.size,dtype=np.int64)
  sample_sums = np.zeros(levels.size,dtype=np.int64)
  
  sample_level_masks = []
  z_sample_in_levels = np.zeros((levels.size,2))
  
  for i, l in enumerate(levels):
    mask = (zz >= l) 
    
    mask_sample = (z_sample >= l)
    sample_level_masks.append( mask_sample )
    sample_sums[i] = z_sample[mask_sample].size
    
    density_sums[i] = np.sum(zz[mask])
    bin_sums[i] = zz[mask].size
    
  
  # average propability per pin in area enclosed by level
  # normed_density_sums = density_sums / bin_sums
  # propability per squareroot bin, why does it work, what does it mean? 
  #      -> Correction for enclosed Area!
  bin_mask = bin_sums != 0.
  normed_density_sums = density_sums[bin_mask] / np.sqrt(bin_sums[bin_mask])
  scaled_sample_sums = sample_sums / z_sample.size
  normed_sample_sums = scaled_sample_sums[bin_mask] /np.sqrt(bin_sums[bin_mask])
  # ---------------------------------------------------------------------------
  
  
  # ---------- Derive (isopleths / levels) of interest --------------------------
  mask_max = normed_density_sums.max() == normed_density_sums 
  level_max = levels[bin_mask][mask_max]
  
  # Check documentation for false positives in case of values < 1
  # https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
  mask_half_max = np.isclose(normed_density_sums, \
                             normed_density_sums.max()/2, atol = 0., rtol= 0.02) 
  levels_half_max = levels[bin_mask][mask_half_max]
  
  if levels_of_interesst == 0:    
    # levels_of_interesst = np.zeros(1+levels_half_max.size)
    levels_of_interesst = np.zeros(2)
    
    levels_of_interesst[0] = level_max[0]
    # levels_of_interesst[1:] = levels_half_max 
    levels_of_interesst[1] = levels_half_max[-1] 
    levels_of_interesst= np.sort(levels_of_interesst)
    contour_infield_labels = '0.5','1.0'
    contour_labels = 'Iso50','Iso100'  
    
  else:
    levels_of_interesst =  np.array(levels_of_interesst) * levels.max()
    levels_of_interesst= np.sort(levels_of_interesst)
    
    contour_infield_labels = '0.5','1.0'
    contour_labels = 'Iso50','Iso100' 

  z_sample_in_levels_of_interesst = []
  mask_sample_in_levels_of_interesst = []
  
  z_sample_sums_in_levels_of_interesst = np.zeros_like(levels_of_interesst)
  sample_bin_sums_of_interesst = np.zeros_like(levels_of_interesst)
  
  for i, l in enumerate (levels_of_interesst):
     mask = (z_sample >= l)
     mask_sample_in_levels_of_interesst.append( mask )
     z_sample_in_levels_of_interesst.append([ x[mask], y[mask] ]   )
     z_sample_sums_in_levels_of_interesst[i] = x[mask].size
     sample_bin_sums_of_interesst[i] = zz[ (zz >= l) ].size
  
  scaled_z_sample_sums_in_levels_of_interesst = \
    z_sample_sums_in_levels_of_interesst / z_sample.size
  normed_z_sample_sums_in_levels_of_interesst = \
    scaled_z_sample_sums_in_levels_of_interesst  \
      / np.sqrt(sample_bin_sums_of_interesst)
    
  
  
  
  # %% Derive pm center & contour plot of interest
  if curve_of_growth_plot is True:
  
    fig_density, ax_density = \
      plt.subplots(2,1,figsize=(fig_width, fig_height*0.7),sharex=True)
      # plt.subplots(2,1,sharex=True,figsize=(fig_width, 10))
      
    
    ax_density[0].plot(levels,density_sums, label = 'KDE from data sample', \
                       color = 'grey', alpha = 1)
    # ax_density[0].scatter(levels,density_sums, label = 'estimated propability', s= 1)
    # ax_density[1].plot(levels[bin_mask],normed_density_sums,label= 'normed propability')  
    
    
    ax_density[0].scatter(levels,scaled_sample_sums, label = 'data sample', \
                          color = 'darkred', s=1,zorder=2) 
    # ax_density[0].scatter(levels_of_interesst, scaled_z_sample_sums_in_levels_of_interesst , \
    #                       # label = 'scaled sample data points of interesst', \
    #                         color = 'black', s= 15, marker = 'x')
    
    
      
    # ax_density[1].scatter(levels[bin_mask],normed_density_sums, s = 1, \
    #                       label= 'normed estimated propability')  
    ax_density[1].plot(levels[bin_mask],normed_density_sums,  \
                          label= 'normed estimated propability',\
                            color = 'grey', alpha = 1)
      
    # ax_density[1].scatter(levels[bin_mask][1:],normed_diff, s = 1, label= 'normed diff')  
    ax_density[1].scatter(levels[bin_mask],normed_sample_sums, \
                          label = 'normed sample data points', \
                            color = 'darkred', s=1,zorder=2)   
    
    # ax_density[1].scatter(levels_of_interesst,normed_z_sample_sums_in_levels_of_interesst, \
    #                       label = 'normed sample data points', \
    #                         color = 'black', s=15, marker = 'x')    
      
    
    # contour_infield_labels = '0.5','1.0'
    # contour_labels = '0.5 - 50% normed prop. max','1.0 - normed prop. max'     
    linewidths = (0.6,0.6)
    linestyles = ['--','-.']
    
    for i, level in enumerate(levels_of_interesst):
      ax_density[0].axvline(level, label = contour_labels [i],\
                            color = 'k',linewidth=linewidths[i], zorder = 0.5,\
                              linestyle = linestyles[i])
      ax_density[1].axvline(level, label = contour_labels [i],\
                            color = 'k',linewidth=linewidths[i], zorder = 0.5,\
                              linestyle = linestyles[i])
       
    # for i, l in enumerate(levels_of_interesst):
    #     ax_density[0].annotate( '{0:s}: {1:1.2g}'.format(r'N',z_sample_sums_in_levels_of_interesst[i]), \
    #                   xy=(l, scaled_z_sample_sums_in_levels_of_interesst[i]), \
    #                   textcoords='offset pixels', \
    #                   xytext=(25,-40) \
    #                 )
    #     ax_density[1].annotate( '{0:s}: {1:1.2e}'.format(r'Isopleth',l), \
    #                   xy=(l, normed_z_sample_sums_in_levels_of_interesst[i]), \
    #                   textcoords='offset pixels', \
    #                   xytext=(70,-10) \
    #                 )
    
    ax_density[0].invert_xaxis()
    ax_density[0].set_xscale('log')
    ax_density[0].set_ylabel('propability')
    ax_density_0_1 = ax_density[0].twinx()
    ax_density_0_1.set_ylim(0,x.size )
    ax_density_0_1.set_ylabel('stars')
    
    # ax_density_0_1.axhline(N_simbad, \
    if N_simbad is not 0:
      ax_density[0].axhline(N_simbad/z_sample.size, \
                           label = 'N=6336 verified by SIMBAD',\
                             color = 'black', linestyle = 'dashed') 
      # No meaning, therefore not included in plot!
      # ax_density[1].axhline(N_simbad / (z_sample.size*np.sqrt(N_bins)), \
      #                         label = 'Verified by SIMBAD', \
      #                           color = 'black', linestyle = 'dashed') 
    
    ax_density[1].set_xlabel('isopleth [propability]')
    ax_density[1].set_ylabel('normed propability')
    ax_density_1_1 = ax_density[1].twinx()
    ax_density_1_1.set_ylim(0,x.size / np.sqrt(N_bins))
    ax_density_1_1.set_ylabel('normed stars')
    
    
    
    ax_density[0].legend(framealpha=1.)
    # ax_density[1].legend()
    
    if savefig is True:
      fig_density.savefig(fig_sub_path+fig_name+'_curve of growth')
    
  # ---------------------------------------------------------------------------  
  
  # center
  mask_max = (zz.max() == zz)
  x_center = xx[mask_max]
  y_center = yy[mask_max]

  pm_ra_kde_center  = x_center
  pm_dec_kde_center = y_center
    
  if contour_plot_of_interesst is True:
    prop_mot_dens_fig, hist2d_prop_mot_dens_axs = \
      plt.subplots(1,1, figsize = (fig_width,5), sharey=False,sharex=False)
    
    # plot 2D Density with shading and contour
    if ready_to_publish == False:
      hist2d_prop_mot_dens_axs.set_title('Contours of interesst')
    
    
    first_legend_handles = []   
    second_legend_handles = []
    
    
    
    # density (Colorbar replaces legend)
    im_mesh = \
      hist2d_prop_mot_dens_axs.pcolormesh(xx, yy, zz, \
                                          shading='gouraud', \
                                            cmap=plt.cm.BuGn_r,\
                                              zorder = 0)
    
    

    center_2D_KDE = \
    hist2d_prop_mot_dens_axs.scatter(xx[mask_max], yy[mask_max], \
                                     label = '2D KDE center', \
                                     marker = 'x', color = 'black', s = 55,\
                                       zorder = 1.9)
    # legend entry  
    first_legend_handles.append(center_2D_KDE)
      
    
    x_delta = np.min( [np.abs(x.min() - xx[mask_max]), \
                       np.abs(x.max() - xx[mask_max]) ])
    y_delta = np.min([ np.abs(y.min() - yy[mask_max]), \
                      np.abs(y.max() - yy[mask_max]) ])
    
    delta = np.min([x_delta,y_delta])  
    # delta = np.max([x_delta,y_delta])  
    
    x_lim_cont_plot = (x_center - delta,x_center + delta)
    y_lim_cont_plot = (y_center - delta,y_center + delta)

    hist2d_prop_mot_dens_axs.set(xlim = x_lim_cont_plot,\
                                  ylim= y_lim_cont_plot)
    
    hist2d_prop_mot_dens_axs.set_aspect('equal')
    
    # data sample
    col = 'lightgrey'
    label = 'data sample'
    data_sample_handle = \
      hist2d_prop_mot_dens_axs.scatter(x,y,color = col, \
                                       s=0.2, marker = '.',\
                                         zorder = 0.2,\
                                           label=label)
    # entry = mlt.collections.PathCollection([], [], color= col , marker='o',\
    #                     markersize=5, label=label)  
    first_legend_handles.append(plt.scatter([],[],color= col , marker='.',\
                        s = 10, label=label))
    
    cont_colors = ('silver','whitesmoke')
    #contour
    im = hist2d_prop_mot_dens_axs.contour(xx, yy, zz, \
                                          levels = levels_of_interesst, \
                                            colors = cont_colors,\
                                              zorder = 0.5)
    
    # contour_infield_labels = '0.5','1.0'
    # contour_labels = '0.5 - 50% normed prop. max','1.0 - normed prop. max'   
    
    # Does copy contour lines into legend, which exceeds legend
    # Use proxy for legend instead!!!
    # for i in range(len(contour_labels)):
    #   im.collections[i].set_label(contour_labels[i])
    #   first_legend_handles.append(im.collections[i])
    
    for i, c in enumerate(cont_colors):
      cont_entry = mlt.lines.Line2D([], [], color= c , marker='',\
                          markersize=15, label=contour_labels[i])  
      first_legend_handles.append(cont_entry)

    fmt = {}
    contour_infield_labels = '0.5','1.0'
    for i, inf_label in zip(im.levels,contour_infield_labels):
      fmt[i] = inf_label
    hist2d_prop_mot_dens_axs.clabel(im,fmt=fmt,inline = 1)
    
    
    zz_max = np.max(zz)
    zz_max_dez = np.abs(np.log10(zz_max))
    zz_max_dez = int(zz_max_dez)+1
    colorbar_max = np.round(zz_max,zz_max_dez)
    col_ticks = np.arange(colorbar_max,0.0,-10**-zz_max_dez)      
    col_bar = prop_mot_dens_fig.colorbar(im_mesh, ax = hist2d_prop_mot_dens_axs,\
                                ticks = col_ticks,\
                                  fraction = 0.046,
                                  pad = 0.04
                                  )
    # col_bar = prop_mot_dens_fig.colorbar(im_mesh, ax = hist2d_prop_mot_dens_axs, \
                                # ticks = im.levels[::5])
    # Duplicate in information
    # col_bar = prop_mot_dens_fig.colorbar(im, ax = hist2d_prop_mot_dens_axs)
    
    col_bar.set_label(('{0:s} {1:3.1g} {2:s}').format(  \
                       'propability to find a star in an area of',bin_A,'mas²'))
    
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Set aspect of the main axes.
    # hist2d_prop_mot_dens_axs.set_aspect(1.)
    
    
    
    # create new axes on the right and on the top of the current axes
    divider = make_axes_locatable(hist2d_prop_mot_dens_axs)
    # below height and pad are in inches
    ax_histx = divider.append_axes("top", 0.5, pad=0.1, \
                                   sharex=hist2d_prop_mot_dens_axs)
    ax_histy = divider.append_axes("right", 0.5, pad=0.1, \
                                   sharey=hist2d_prop_mot_dens_axs)
    
    # make some labels invisible
    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histy.yaxis.set_tick_params(labelleft=False)
    
    # Number of samples per data point in histogram 1D
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html
    bins_1D_hist = 1000
    # Returns [n, bin_edges, patches] 
    hist_1D_x = \
      ax_histx.hist(x,bins = bins_1D_hist, density = False, stacked = True,\
                    label = 'data sample',\
                      color = 'mediumseagreen', alpha = .8)
    
    hist_1D_y = \
      ax_histy.hist(y,bins = bins_1D_hist, density = False, stacked = True , \
                    label = 'data sample pm dec', \
                  color = 'mediumseagreen', alpha = .8, orientation='horizontal')
    
    # Only one entry necessary
    second_legend_handles.append(hist_1D_x[-1][0])
    
    
    hist_1D_x_delta = hist_1D_x[1][1]-hist_1D_x[1][0]
    hist_1D_y_delta = hist_1D_y[1][1]-hist_1D_y[1][0]
    
    k_x = kde.gaussian_kde(x)
    # x_1D = xx[:,0]
    x_1D = hist_1D_x[1]
    x_gauss = k_x(x_1D)
    
    k_y = kde.gaussian_kde(y)
    # y_1D = yy[0,:]
    y_1D = hist_1D_y[1]
    y_gauss = k_y(y_1D)
    
    # propability density function * bin_widt yields propability 
    # np.sum(y_gauss)*hist_1D_y_delta = 1 = np.sum(y_gauss*hist_1D_y_delta)
    # multiply propability with number of data points to get absolute value
    gauss_1D_x = \
    ax_histx.plot(x_1D, x_gauss*hist_1D_x_delta*x.size, \
                  label = '1D KDE',color = 'black' )
    gauss_1D_y = \
      ax_histy.plot(y_gauss*hist_1D_y_delta*y.size,y_1D, \
                  label = 'pm ra gauss',color = 'black' )
    # Only one entry necessary
    second_legend_handles.append(gauss_1D_x[0])
    
    
    mask = (x_gauss == np.max(x_gauss) ) 
    index_max = np.nonzero(mask)
    gauss_1D_x_max = (x_1D[index_max][0]+x_1D[index_max[0]+1][0])/2
    pm_ra_gauss_center = gauss_1D_x_max 
    
    ax_histx.axvline(gauss_1D_x_max , label = 'pm ra_max gauss',\
                     color = 'k',linewidth=0.5)
      
    gauss_1D_x_max_line = \
      hist2d_prop_mot_dens_axs.axvline(gauss_1D_x_max ,label = '1D KDE max',\
                                     color = 'k',linewidth=0.5, zorder = 1.5)
   
    mask = (y_gauss == np.max(y_gauss) ) 
    index_max = np.nonzero(mask)
    gauss_1D_y_max = (y_1D[index_max][0]+y_1D[index_max[0]+1][0])/2
    pm_dec_gauss_center = gauss_1D_y_max
    ax_histy.axhline(gauss_1D_y_max , label = '1D KDE max',\
                     color = 'k',linewidth=0.5)
      
    gauss_1D_y_max_line = \
      hist2d_prop_mot_dens_axs.axhline(gauss_1D_y_max, label = 'pm dec max gauss',\
                                     color = 'k',linewidth=0.5, zorder = 1.5)
    
    # Only one entry necessary
    second_legend_handles.append(gauss_1D_x_max_line)
    
    hist2d_prop_mot_dens_axs.set_ylabel('pm dec [mas/year]')
    hist2d_prop_mot_dens_axs.set_xlabel('pm ra [mas/year]')
    
    ax_histx.set_ylabel('N')
    ax_histy.set_xlabel('N')
    
    first_legend = \
      hist2d_prop_mot_dens_axs.legend(handles = first_legend_handles,\
                                      loc = 'upper left',\
                                        framealpha = 1)
    hist2d_prop_mot_dens_axs.add_artist(first_legend)
    hist2d_prop_mot_dens_axs.legend(handles = second_legend_handles, \
                                    loc='lower right',\
                                      framealpha = 1)
    # ax_histx.legend(handles = second_legend_handles, \
    #                                 bbox_to_anchor=(1.05, 1),\
    #                                   loc='upper left', borderaxespad=0.)
    if savefig is True:
      prop_mot_dens_fig.savefig(fig_sub_path+fig_name+'_KDE')
  else:
    print('pm_ra_gauss_center & pm_dec_gauss_center cannot be calculated with contour_plot_of_interesst = False')
    print('values set to 0.0')
    pm_ra_gauss_center = 0.0
    pm_dec_gauss_center = 0.0
  # ---------------------------------------------------------------------------   
  # %%  proper_motion_surface_plot
  if proper_motion_surface_plot_on is True:
    fig_3d = plt.figure()
    ax_3d = fig_3d.add_subplot(projection='3d')
    surf = ax_3d.plot_surface(xx, yy, zz, rstride=1, cstride=1,
                           linewidth=0, antialiased=False)
    ax_3d.contour3D(xx, yy, zz)
    if savefig is True:
      fig_3d.savefig(fig_sub_path+fig_name+'_3D_KDE')
    
  # ---------------------------------------------------------------------------  
  # uncomment to use as script  
  return [levels,density_sums,normed_density_sums, \
          mask_sample_in_levels_of_interesst, \
            levels_of_interesst,\
            contour_labels, \
              [pm_ra_kde_center,pm_dec_kde_center], \
                [pm_ra_gauss_center,pm_dec_gauss_center] ]
    

# %% main routine and function call 

# import proper_motion_evaluation as pme

[levels, density_sums, normed_density_sums, \
  mask_sample_in_levels_of_interesst, \
  levels_of_interesst, contour_labels, \
    [pm_ra_kde_center,pm_dec_kde_center], \
    [pm_ra_gauss_center,pm_dec_gauss_center] ] = \
propability_density_estimation(x,y, \
                      levels_scaling             = levels_scaling, \
                      levels_of_interesst        = levels_of_interesst, \
                      curve_of_growth_plot       = curve_of_growth_plot ,\
                      contour_plot_of_interesst  = contour_plot_of_interesst,\
                    
                      levels_spacing_plot        = levels_spacing_plot ,\
                      contour_plot               = contour_plot,\
                      proper_motion_surface_plot_on = proper_motion_surface_plot_on ,\
                     
                      ready_to_publish           = ready_to_publish,
                      savefig                    = savefig,
                      fig_name                   = fig_name,
                      fig_sub_path               = fig_sub_path
                      )

print ("pm kde center [ra,dec]:    [%0.3f , %0.3f] mas" % (pm_ra_kde_center, pm_dec_kde_center))
print ("pm gauss center [ra,dec]:  [%0.3f , %0.3f] mas" % (pm_ra_gauss_center, pm_dec_gauss_center))
  
if save_data is True:  
  pm_eval_name = file_name +'_pm-center.csv'
  # save_gaia(iso,iso_name,file_sub_path=file_sub_path)
  
  pm_list= []
  pm_list.append(["pm kde center; [ra,dec];    [%0.4f , %0.4f] mas" % (pm_ra_kde_center, pm_dec_kde_center),\
                  "pm gauss center [ra,dec];  [%0.4f , %0.4f] mas" % (pm_ra_gauss_center, pm_dec_gauss_center) ] )
  
  # create new csv file
  with open(file_sub_path+pm_eval_name, 'w', newline='') as csvfile:
    # create new csv object
    csv_write = csv.writer(csvfile, delimiter=',')
    # store data in new csv file
    csv_write.writerows(pm_list)
  
  # data sample with pm enclosed by isopleth with 1.0 prob.max
  gaia_dr3_array_pm_iso_05 =  gaia_dr3_array[mask_sample_in_levels_of_interesst[0]]  
  gaia_dr3_array_pm_iso_1 =  gaia_dr3_array[mask_sample_in_levels_of_interesst[1]]  
  
  isos_labels = ['pm_iso_05','pm_iso_1']
  isos = [gaia_dr3_array_pm_iso_05,gaia_dr3_array_pm_iso_1]
  
  for j, iso in enumerate(isos):
    iso_name = file_name +'_'+isos_labels[j]
    save_gaia(iso,iso_name,file_sub_path=file_sub_path)

