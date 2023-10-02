#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 20:59:28 2023

@author: Thomas Suntinger

Import recommendation:
  from curve_of_growth import probability_KDE_1D as p_KDE_1D
  
  
Functionality:
 
  - probability_KDE_1D() function estimates probability density distribution 
    in 1D phase space with a Gaussian Kernel Density Estimator (KDE).
  - To use this file as a script, one can uncomment marked areas below. 
    Additionally, function definition and return commant must be commented.
  - A variety of plots can be made and saved, if adjoined parameters 
    are set to True.

Improvements:
  - Clean unnecessary comments and lines of code
  - Write documentation and improve in code documentation.
  - For Script mode: Import load_query from Gaia_load.py and load data
  - Create a generalized version of probability_density_estimation() to use it
    for two arbitrary parameters.
  - Include the generalized version of proper_motion_evaluation.py 2D KDE function.
  
  
"""


import numpy as np
from scipy.stats import kde
import matplotlib as mlt
import matplotlib.pyplot as plt


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
# data = rng.multivariate_normal([0, 0], [[1, 0.5], [0.5, 3]], 10000)
# x, y = data.T



# %% Default values for file used as script




# [x,y] = np.load('selected_pm_array.npy')


# Default values for functioin:
N_bins = 250

N_levels = 100
levels = N_levels

extention_value = 5
extended_area = True



# Interessting for data evaluation
levels_spacing_plot   = True
contour_plot          = True
proper_motion_surface_plot_on = True

# Results for proper motion evaluation
curve_of_growth_plot  = True
contour_plot_of_interesst = True


levels_of_interesst = 0
# levels_of_interesst = [0.2,0.66,1.]
# default scaling is 'lin'
levels_scaling = 'lin'

# Change, only for testing function
# levels_scaling = 'lin'


#%%
fig_height = 4
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)


savefig = False
fig_name = '1D_KDE_analysis'
fig_sub_path = 'gfx/'

label_x_quantity = 'stars [N]'
label_x_data = 'GAIA sample'

label_x_quantity = 'data points [N]'
label_x_data = 'data sample'


# %% main routine

# as script
# if True:



def probability_KDE_1D(x, 
          levels_of_interesst = levels_of_interesst, 
          N_bins = N_bins,
          levels = N_levels, 
          levels_scaling = levels_scaling,
          savefig = savefig,
          fig_name = fig_name,
          fig_sub_path = fig_sub_path,
          label_x_data = label_x_data,
          label_x_quantity = label_x_quantity,
          curve_of_growth_plot = curve_of_growth_plot
          ):
  
  # xx, yy = np.mgrid[x.min():x.max():N_bins*1j, y.min():y.max():N_bins*1j]
  xx = np.linspace(x.min(),x.max(),N_bins)
  
  # Area of one bin
  delta_x = xx[1]-xx[0]
  # delta_y = yy[0,1]-yy[0,0]
  bin_A = delta_x
  # bin_A = delta_x * delta_y
  
  
  k = kde.gaussian_kde(x)
  z = k(xx)
  # multiply propability density with bin area to get propability 
  # to find data point in bin
  z = z * bin_A 
  # Scale data points to avoid float inaccuracy 
  # z = z*10**3 
  zz = z.reshape(xx.shape)
  
  z_min = z.min()
  z_max = z.max() 
  
  #Get x for KDE maximum 
  x_center = xx[np.equal(z_max,z)]
  
  
  # For logarithmic scale, find smallest z excluding 0.0
  z_log_min = z[(z > 0.0)].min()
  z_log_min = np.log10(z_log_min) 
  # scaling factor of z_min, break down for 0.05  < f_scaling < 0.006
  # z_log_min = z_log_min*0.006
  z_log_min = z_log_min*0.018
  z_max_log = np.log10(z_max)

  

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
  
  
  #%%
   
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
  elif levels_scaling == 'custom':
    levels = levels
  else:
    print('Invalid input, default is linear scaling')
    levels = levels_lin
  
  

  #%%
  z_sample = k(x) * bin_A
  
  density_sums = np.zeros(levels.size)
  # normed_density_sums = density_sums.copy()
  bin_sums = np.zeros(levels.size,dtype=np.int64)
  sample_sums = np.zeros(levels.size,dtype=np.int64)
  
  sample_level_masks = []
  
  
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
  bin_mask = bin_sums != 0.
  normed_density_sums = density_sums[bin_mask] / np.sqrt(bin_sums[bin_mask])
  scaled_sample_sums = sample_sums / z_sample.size
  normed_sample_sums = scaled_sample_sums[bin_mask] /np.sqrt(bin_sums[bin_mask])
  # ---------------------------------------------------------------------------
  
  
  mask_max = np.equal(normed_density_sums.max() , normed_density_sums )
  # Chose the smallest level, if multiple levels have the same norm.dens.max
  level_max = np.min( levels[bin_mask][mask_max] )
  # Adjust mask_max to have only one entry
  mask_max = np.equal(level_max , levels[bin_mask] )
  
 

   
  # Data interessting between max and norm_dens_min
  norm_dens_max = normed_density_sums.max()
  
  # highest level = levels[-1] corresponds to norm_dens_max[-1]
  # found before norm_dens_max
  # This data is not interessting for curve of growth selection criterion
  # It is too close to center!
  
  # lowest level = levels[0] corresponds to norm_dens_min
  # found after norm_dens_max
  norm_dens_min_low_level = normed_density_sums[0]
  
  if norm_dens_min_low_level > norm_dens_max / 2:
  
    norm_delta = norm_dens_max - norm_dens_min_low_level 
    norm_half_delta = norm_dens_max - norm_delta/2
  else:
    norm_half_delta  = norm_dens_max / 2
  
  mask_search_half_range = levels[bin_mask] <= level_max 
  
  mask_half_range = norm_half_delta  > normed_density_sums 
  mask_half_range &= mask_search_half_range
  
  if mask_half_range.any():
    level_half_range = levels[bin_mask][mask_half_range][-1]
    mask_half_range = np.equal(level_half_range, levels[bin_mask] )
    dens_half_range = normed_density_sums[bin_mask]
    dens_half_range = dens_half_range[mask_half_range]
    dens_prozent_half_range = dens_half_range/norm_dens_max
    only_one_level_of_interesst = False
  else:
    only_one_level_of_interesst = True
    
 
   
  
  
  # First data point next to half max
  # mask_half_max = normed_density_sums.max()/2 < normed_density_sums
  
  # dens_larger_half_max = levels[bin_mask][mask_half_max]
  # levels_half_max = [dens_larger_half_max[-1]  ]
  
  # mask_half_max = levels_half_max == levels[bin_mask]
  
  
  if only_one_level_of_interesst is False:
  
    if levels_of_interesst == 0:    
      # levels_of_interesst = np.zeros(1+levels_half_max.size)
      levels_of_interesst = np.zeros(2)
      
      levels_of_interesst[0] = level_max
      # levels_of_interesst[1:] = levels_half_max 
      # levels_of_interesst[1] = levels_half_max[-1] 
      levels_of_interesst[1] = level_half_range
      levels_of_interesst= np.sort(levels_of_interesst)
      # contour_infield_labels = '0.5','1.0'
      
      # contour_label_half_range = "%0.2f" % (dens_prozent_half_range)
      contour_label_half_range = "%2.0f" % (dens_prozent_half_range*100)
      contour_infield_labels = contour_label_half_range,'100'
      contour_labels = 'Iso'+ contour_label_half_range,'Iso100'  
      
    else:
      levels_of_interesst =  np.array(levels_of_interesst) * levels.max()
      levels_of_interesst= np.sort(levels_of_interesst)
      
      contour_infield_labels = '50','100'
      # contour_infield_labels = '0.5','1.0'
      contour_labels = 'Iso50','Iso100' 

  else:
    levels_of_interesst = [level_max]
    contour_labels = ['Iso100']



  z_sample_in_levels_of_interesst = []
  mask_sample_in_levels_of_interesst = []
  
  z_sample_sums_in_levels_of_interesst = np.zeros_like(levels_of_interesst)
  sample_bin_sums_of_interesst = np.zeros_like(levels_of_interesst)
  
  
  
  for i, l in enumerate (levels_of_interesst):
     mask = (z_sample >= l)
     mask_sample_in_levels_of_interesst.append( mask )
     z_sample_in_levels_of_interesst.append( x[mask]   )
     z_sample_sums_in_levels_of_interesst[i] = x[mask].size
     sample_bin_sums_of_interesst[i] = zz[ (zz >= l) ].size
  
  scaled_z_sample_sums_in_levels_of_interesst = \
    z_sample_sums_in_levels_of_interesst / z_sample.size
  normed_z_sample_sums_in_levels_of_interesst = \
    scaled_z_sample_sums_in_levels_of_interesst  \
      / np.sqrt(sample_bin_sums_of_interesst)
  
    
  #%%
  if curve_of_growth_plot is True:
  
    fig_density, ax_density = \
      plt.subplots(2,1,figsize=(fig_width,fig_height*.7),sharex=True)

      
    
    ax_density[0].plot(levels,density_sums, label = label_x_data+' KDE', \
                       color = 'grey', alpha = 1)

  
    
    ax_density[0].scatter(levels,scaled_sample_sums, label = label_x_data, \
                          color = 'darkred', s=1,zorder=2) 

    

    ax_density[1].plot(levels[bin_mask],normed_density_sums,  \
                          label= 'normed estimated propability',\
                            color = 'grey', alpha = 1)
      
    # ax_density[1].scatter(levels[bin_mask][1:],normed_diff, s = 1, label= 'normed diff')  
    ax_density[1].scatter(levels[bin_mask],normed_sample_sums, \
                          label = 'normed sample data points', \
                            color = 'darkred', s=1,zorder=2)   
    

  
    linewidths = (0.6,0.6)
    linestyles = ['--','-.']
    
    for i, level in enumerate(levels_of_interesst):
      ax_density[0].axvline(level, label = contour_labels [i],\
                            color = 'k',linewidth=linewidths[i], zorder = 0.5,\
                              linestyle = linestyles[i])
      ax_density[1].axvline(level, label = contour_labels [i],\
                            color = 'k',linewidth=linewidths[i], zorder = 0.5,\
                              linestyle = linestyles[i])
       

    ax_density[0].invert_xaxis()
    ax_density[0].set_xscale('log')
    ax_density[0].set_ylabel('propability')
    ax_density_0_1 = ax_density[0].twinx()
    ax_density_0_1.set_ylim(0,x.size )
    ax_density_0_1.set_ylabel(label_x_quantity)
    
   

    ax_density[1].set_xlabel('isopleth [propability]')
    ax_density[1].set_ylabel('normed propability')
    ax_density_1_1 = ax_density[1].twinx()
    ax_density_1_1.set_ylim(0,x.size / np.sqrt(N_bins))
    ax_density_1_1.set_ylabel('normed '+label_x_quantity)
    
    
    
    ax_density[0].legend(framealpha=1., loc = 'lower right')
    
    if savefig is True:
      fig_density.savefig(fig_sub_path+fig_name+'_curve of growth')
     
       
  return [levels,density_sums,normed_density_sums, \
          mask_sample_in_levels_of_interesst, \
            levels_of_interesst,\
            contour_labels, \
              x_center, xx, zz  ]