#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 13:21:47 2023

@author: Thomas Suntinger

Import recommendation:

  
Input:
  - load: data set of a cluster as csv file (e.g. from gaiaClusterMemberSelection.py)
  
Dependencies:
  - GAIA_load.py
  - GAIA_save.py
  - parallax_evaluation.py
  - calcSearchRadius.py

Dependencies for visualization:
  - matplotlib
  - plotStyle.py
  
Functionality:
  - Script loads a cluster member data set and evaluates parallax and distance
    to the clusters center.
  - For the evaluation, a variety of parameters can be finetunned to get 
    the most accurate results (see 'Structure:Parameter configuration for evaluation').
  - The error distribution and a histogram with the parallax/distance including
     the Lorentz fit evaluation are visualized (see Output plots)
  
Structure ( & parameters):
  - Init 
  - Plot settings 
          (figure height)
          (default plot colors)
  - Store figure settings
  - Distance and distance error calculation
  - Store distance and error into csv file
          (file path and name)  
  - Parameter configuration for evaluation
          (plot and print output on/off)
          (Search radius on/off, value)
          (error selection on/off, colors)
          (Toggle parallax/distance for Lorentz evaluation)
          (outlier filter settings for parallax/distance and adjoined error)
      - Outlier filter  
  - Error filter for distance evaluation 
      -Prepare parallax/distance error filters
          (filter bounds for errors)
      - Prepare data sets for error_selection evaluation
  - Work around negative parallaxes (not working, do not waist time with it!)
  - Lorentz evaluation
          (number of bins for histogram)
  - Structure of structured numpy arrays
  - Visualization of chosen data sets distance and adjoint error distribution
          (see Output plots)
    
Output:
  - store: csv file cluster member selection with distance and dist. error parameters
  - print: outlier filter for unrealistic data
  - print: Lorentz function results for error selections (see parallax_evaluation.py)
  
Output plots:
  - parallax error distribution & parallax/distance evaluation histogram 
      (includes parallax/distance result of evaluation)
  
Improvements:
  - Clean unnecessary comments and lines of code
  - Write documentation and improve in code documentation.

  - Bonus: Work around negative parallaxes (not working, do not waist time with it!)
            - Try other clusters and check, if M5 is special 
            - Wait for DR4/DR5 and redo evaluation
            - If resolved, implement distance error selection 

"""
#%% Init

import numpy as np
import parallax_evaluation as p_e
import calcSearchRadius

import matplotlib.pyplot as plt

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

# dataDR3 - Starhorse convolution
t07 = 'orangered'

# Av extinction models
t08 = 'lightgrey'

# RGB
t09 = 'teal'
# AGB
t10 = 'purple'


#%% ### ----------- Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/clus_sel/'


# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)


#%% ### ------------ Load /Store data settings ---------------------------------------

folder_name = 'data/clusterMemberSelection/'

## -------------- GAIA data ---------------------------------------------------
# data_set_selected = 'Gaia_M3_r-18arcmin.csv'
data_set_selected = 'Gaia_M5_r-9arcmin_zpc_pm_clusterMemberSelection'
# data_set_selected = 'Gaia_M67_r-25arcmin.csv'

file_name = data_set_selected

## ------------ load data into variables --------------------------------------

from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

gaia_cluster_members = file_data
label_gaia_cluster_members = 'cluster members'
label_gaia_cluster_members_n = 'cluster members ('+str(gaia_cluster_members.size)+' stars)'

print()
print(label_gaia_cluster_members_n)
print()
print()



#---------Store data --------------------------------

file_sub_path = folder_name
from GAIA_save import save_gaia as save_gaia
save_data = True

#%% ### Distance and distance error calculation -------------------------------

# Store parallax
par_x = gaia_cluster_members['parallax']
par_x_error = gaia_cluster_members['parallax_error']

#Store distance
# Calc distance 
dist_x = 1/par_x
# Calc distance error
error_x = par_x_error / (par_x**2 - par_x_error**2)

#%% ### Store distance and error into csv file
# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'distance', dist_x, dtypes='<f8')
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'distance_error', error_x, dtypes='<f8')

save_gaia(gaia_cluster_members, file_name +'_distEval' ,file_sub_path=file_sub_path)
#%% ### Parameter configuration for evaluation

# Turn on / off output
par_Lor_sel_plot_on = True
lorentz_plot_on = False
lorentz_text_on = False

# Enable search radius comparison 
search_radius_on = True

if search_radius_on is True:
  # Search radius of Gaia query in arcmin:
  SR_angle = 9
  SR_angle = calcSearchRadius.convertArcMinToRad(SR_angle)
  

# Turn error selection on or off
error_selection = True
# Chose colors for error selection
color_selection = [t06,t07,t10]


# Use parallax or distance for Lorentz fit
parallax_eval_on = True



if parallax_eval_on is False: 
  # distance selection
  data_x_label = 'distance'
  data_x_unit  = '[kPc]'
  data_x = dist_x
  data_error_x = error_x
  
  # Unrealistic data filter in kPc
  unreal_data_min = -1000
  unreal_data_max = 1000
  # Unrealistic data error filter kPc
  unreal_data_error_min = -1000
  unreal_data_error_max = 1000

else:
  # parallax selection
  data_x_label = 'parallax'
  data_x_unit  = '[mas]'
  data_x = par_x
  data_error_x = par_x_error
  
  # Unrealistic data filter in mas
  unreal_data_min = -1000
  unreal_data_max = 1000
  # Unrealistic data error filter mas
  unreal_data_error_min = 0
  unreal_data_error_max = 1000

print('Lorentz fit applied to: '+data_x_label)
print()

# Outlier filter ---------------------------------------------------
mask_dist_valid = np.logical_and(data_x < unreal_data_max, data_x > unreal_data_min)
mask_dist_outliers = np.logical_not(mask_dist_valid)

# print lost stars from parallax/distance filter
print('unrealistic outlier filter constrains data eval to:')
print('['+str(unreal_data_min)+' , '+ str(unreal_data_max)+']'+str(data_x_unit))
print('Number of outliers: ' + str(data_error_x[mask_dist_outliers].size))
print()
mask_dist_error_valid = np.logical_and(data_error_x < unreal_data_error_max, data_error_x > unreal_data_error_min)
mask_dist_error_outliers = np.logical_not(mask_dist_error_valid)

#print lost stars from parallax/distance error filter
print('unrealistic error outlier filter constrains data eval to:')
print('['+str(unreal_data_error_min)+' , '+ str(unreal_data_error_max)+']'+str(data_x_unit))
print('Number of outliers: ' + str(data_x[mask_dist_error_outliers ].size))
print()


#%% ### Error filter for distance evaluation

# Prepare parallax/distance error filters ---------------
if parallax_eval_on is True:
  
  # mask for errors smaller 0.02 mas
  error_upper_bound = 0.02
  mask_dist_error_01kpc  = (data_error_x < error_upper_bound)
  error_upper_bound_label = ' $<$ '+str(error_upper_bound)+''+data_x_unit
  
  # mask for errors smaller 20% 
  # Note: Negative values for error possible! 
  #         - Set unrealistic_data_error_min = 0  to exclude negative errors
  error_ratio = 0.2
  mask_dist_error_01frac = (data_error_x / np.abs(data_x)) < error_ratio  
  error_ratio_label = ' $<$ '+str(error_ratio*100)+'\%'

else:
  print('No implementation for distance, until fix of negative parallax values.')
  print('error_selection is turned off!')
  error_selection = False



# Prepare data sets for error_selection evaluation-----------
print('error_selection: '+str(error_selection) )

if error_selection is True:
  # data_sub_sets_label = [' full',error_ratio_label,error_upper_bound_label]
  data_sub_sets_label = ['',error_ratio_label,error_upper_bound_label]
  masks               = [ mask_dist_error_valid & mask_dist_valid,
                          mask_dist_error_valid & mask_dist_valid & mask_dist_error_01kpc,
                          mask_dist_error_valid & mask_dist_valid & mask_dist_error_01frac]
else:
    data_sub_sets_label = [' full']
    masks               = [mask_dist_error_valid & mask_dist_valid ]












#%% ### Work around negative parallaxes (not working, do not waist time with it!)
# Using math to shift parallaxes far from 0, calc pseudo distance and do 
# histogram. Maybe back convertable to distance without a gap at 0 ?!
# Note: parallax = 0 yields star infinitely far away, negative par behind infinity ;) 
# Note: Investigation in parallax determination from Gaia DPAC necessary!!!
# Note: Do not waist more time on that and wait for DR4 or try other cluster!
 
workaround_on = False

if workaround_on is True: 
  par_workaround = par_x[masks[0]]
  a = 1000
  c = 1/(par_workaround + a)
  d_workaround = c/(1-c*a)
  
  dist_fig, dist_axs = plt.subplots(1,2)
  dist_axs[0].hist( \
                              c,\
                            bins = 1000,
                            label = 'test'                          
                                )
  
  dist_axs[1].hist( \
                              d_workaround,\
                            bins = 1000,
                            label = 'test'                          
                                )
  
    
  # Normal distribution inverted------------------------------------
  # Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
  from scipy.stats import norm
  
  mean, var, skew, kurt = norm.stats(moments='mvsk')
  x = np.linspace(norm.ppf(0.01), norm.ppf(0.99), 100)
  
  fig, ax = plt.subplots(1, 1)
  
  ax.plot(x, norm.pdf(x),'r-', lw=5, alpha=0.6, label='norm pdf')
  
  rv = norm()
  ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
  
  y = 1/rv.pdf(x)
  # ax.plot(x, y, 'k-', lw=2, label='inverse')
  
  
  vals = norm.ppf([0.001, 0.5, 0.999])
  np.allclose([0.001, 0.5, 0.999], norm.cdf(vals))
  r = norm.rvs(size=1000)
  
  ax.hist(1/r, density=True, bins='auto', histtype='stepfilled', alpha=0.2)
  # ax.hist(r, density=True, bins='auto', histtype='stepfilled', alpha=0.2)
  ax.set_xlim([x[0], x[-1]])
  ax.legend(loc='best', frameon=False)
  plt.show()

  
  
  
#%% ### Prepare storage for results

#Create empty structured numpy array
# Note: numpy arrays have specific size at initialization!!!
# https://www.datasciencelearner.com/numpy-structured-array-example/
np_dtype = [                         
            ('data_set',  'U1000'      ), 
            ('color',     'U1000'      ), 
            ('lor_amp',   '<f8'         ),
            ('lor_cen',   '<f8'         ),
            ('lor_width', '<f8'         )
            ]

# Structured numpy array including evaluation data
# lor_width fwhm / 2 (error)
dist_eval_array       = np.zeros(3,dtype = np_dtype)

# Dictionary including plot data
# keys = ['data_set','lor_curve','hist_data']
# dist_eval_plot_dict = dict.fromkeys(keys)
dist_eval_plot_dict_list = []




#%% ### Lorentz evaluation

# # for comparison parallax center in [mas] from gaiaParallaxSelection for M5 with 9arcmins
# par_center = 0.14009025675968623

bin_n = 1000

for i, sub_label in enumerate(data_sub_sets_label): 
  
  mask = masks[i] 
  # label = label_gaia_cluster_members+sub_label
  label = 'N: '+str(np.size(data_x[mask]))+ sub_label
  print('\n'+label_gaia_cluster_members+sub_label+':')
  args =  [ data_x[mask] ]
  kwargs =  {
              'bins' : bin_n,\
              'text':label,\
              'show_text': lorentz_text_on,\
              'show_plot': lorentz_plot_on 
              }
  
  # Derive lorentz-fit
  [popt_lorentz,perr_lorentz,lorentz_curve,hist_data] = \
      p_e.lorentz_profile_hist_fit(*args,**kwargs)
      
  
  dist_eval_array[i] = (label,
                       color_selection[i],
                       popt_lorentz[0],
                       popt_lorentz[1],
                       popt_lorentz[2]
                       )
  
  if 'distance' in data_x_label:

    d_plot = dist_eval_array[i]['lor_cen']
    
    if search_radius_on is True:
      # Calc search radius in kPc using d in kPc-----------
      SR_plot = calcSearchRadius.calcSearchRadius(d_plot, SR_angle)
      print('search radius in kPc: '+str(SR_plot))
      # ---------------------------------------------------
    
    SR_min = d_plot-SR_plot
    SR_max = d_plot+SR_plot
    
    SR_min_label='SR'+str(SR_plot)
    SR_max_label='SR'+str(SR_plot)
    
    # search radius in line of sight converted to distance 
    dist_bound_mask = ((gaia_cluster_members['distance'] > SR_min)&(dist_x <  SR_max))
    dist_bound = gaia_cluster_members[dist_bound_mask]
    dist_bound_name = str(dist_bound.size)+' stars in radial search radius: '+str(SR_plot)+' kPc'
    
    print(' '+dist_bound_name)

    
    
  elif 'parallax' in data_x_label:  
    par_plot = dist_eval_array[i]['lor_cen']
    
    d_plot = 1 / par_plot
    
    if search_radius_on is True:
      # Calc search radius in kPc using d in kPc-----------
      SR_plot = calcSearchRadius.calcSearchRadius(d_plot, SR_angle)
      print(' search radius in kPc: '+str(SR_plot))
      # ---------------------------------------------------
      
    par_plot_min = 1 / (d_plot+SR_plot)
    par_plot_max = 1 / (d_plot-SR_plot)
    
    par_delta_min =   par_plot - par_plot_min
    par_delta_max = - par_plot + par_plot_max
    
    SR_min = par_plot_min
    SR_max = par_plot_max
    
    SR_min_label='SR min '+'{0:1.3g}'.format(par_delta_min)
    SR_max_label='SR max '+'{0:1.3g}'.format(par_delta_max)
    
    # search radius in line of sight converted to parallax
    par_bound_mask = ((gaia_cluster_members['parallax'] > SR_min)&(gaia_cluster_members['parallax'] <  SR_max))
    par_bound_array = gaia_cluster_members[par_bound_mask]
    par_bound_name = str(par_bound_array.size)+' stars in radial search radius [min,max]: '+'['+str(SR_min)+','+str(SR_max)+'] mas'
    
    print(' '+par_bound_name)

  
  SR_data = SR_min,SR_min_label,SR_max,SR_max_label
  
  # Multiple dictionaries in a list
  # Needs to be created each loop, otherwise copy necessary!
  keys = ['data_set','lor_curve','hist_data','SR_data']
  dist_eval_plot_dict = dict.fromkeys(keys)
  
  dist_eval_plot_dict['data_set']  = label
  dist_eval_plot_dict['lor_curve'] = lorentz_curve
  dist_eval_plot_dict['hist_data'] = hist_data
  dist_eval_plot_dict['SR_data']   = SR_data
  dist_eval_plot_dict_list.append(dist_eval_plot_dict)



  


#%% ### Structure of structured numpy arrays 

# dist_eval_array
# dist_eval_bound_array
  
# np_dtype = [                         
#             ('data_set',  'S10'      ), 
#             ('color',     'S10'      ), 
#             ('lor_amp',   '<f8'         ),
#             ('lor_cen',   '<f8'         ),
#             ('lor_width', '<f8'         )            
#            ]

# gaia_cluster_members
# data_sub_sets_label = ['full','10% error','0.01kPc error']
# masks               = [mask_dist_valid,
#                        mask_dist_valid & mask_dist_error_01kpc,
#                        mask_dist_valid & mask_dist_error_01frac]
#%% ### Visualization of chosen data sets distance and adjoint error distribution--------------


# data Lorentz selection & error plot -----------------------------------
if par_Lor_sel_plot_on      is True:
  
  data_ax_label       =  data_x_label+' '+data_x_unit
  data_error_ax_label =  data_x_label+' '+data_x_unit
  
  hist_ax_label = 'number [N]'
  error_ax_label = data_x_label+' error '+data_x_unit
  
  
  # Create data evaluation plot
  dist_fig, dist_axs = plt.subplots(1,2,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[2,3]},)
  
  hist_ax = dist_axs[1]
  error_ax = dist_axs[0]
  
  
  # for i, data_set in enumerate(data_sub_sets_label):
  if True:
    i = 0
    # data_set = dist_eval_array[i]
    data = data_x[masks[i]]
    # data_error = data_error_x[masks[i]]
       
    # Multiple dictionaries
    N, bins, patches = \
      hist_ax.hist( \
                          data,\
                            bins = (dist_eval_plot_dict_list[i])['hist_data'][1],
                            label = dist_eval_array[i]['data_set'],
                              color=   dist_eval_array[i]['color'],\
                                # alpha = 1 / (i+1),\
                                alpha = 0.3
                                )
        
    
        
    
    # Plot Search Radius boundaries into histogram
    if search_radius_on is True:
      
      SR_min,SR_min_label,SR_max,SR_max_label = dist_eval_plot_dict_list[i]['SR_data']
         
      
      hist_ax.axvline(SR_min, color=t_fit, label=SR_min_label, \
                                linewidth= 0.8, linestyle = '--')
      hist_ax.axvline(SR_max, color=t_fit, label=SR_max_label, \
                                linewidth= 0.8, linestyle = '-.') 
      
      
    # Plot Lorentz fit into histogram  
    args   =  (
              # With multiple dictionaries
              (dist_eval_plot_dict_list[i])['lor_curve'][0],
              (dist_eval_plot_dict_list[i])['lor_curve'][1]
              )
    
    kwargs = {
              'label' : 'Lorentz Fit ' ,\
              # 'color' : dist_eval_array[i]['color']\
              'color' :t_fit\
              }
    
    hist_ax.plot(*args,**kwargs)  
    
          
    # Plot Lorentz width as fwhm/2 to compare to Search Radius
    lor_delta = dist_eval_array[i]['lor_width']
    lor_center = dist_eval_array[i]['lor_cen']
    lor_amp = dist_eval_array[i]['lor_amp']
    
    hist_ax.scatter(lor_center, lor_amp, color= t_fit, label='lor center'+'{0:10.3g}'.format(lor_center), marker = '+' ,s=100 )
    # hist_ax.axvline(lor_center, ymax=lor_amp,color=t_fit, label='lor center'+'{0:10.3g}'.format(lor_center), \
    #                           linewidth= 0.6, linestyle = '-')
      
    # dummy_1 = '{0:10.3g}'.format(lor_delta)
    
    hist_ax.axvline(lor_center-lor_delta , color=t_fit, label='lor width'+'{0:10.3g}'.format(lor_delta), \
                              linewidth= 0.4, linestyle = '-')
    hist_ax.axvline(lor_center+lor_delta , color=t_fit,  \
                              linewidth= 0.4, linestyle = '-')
  if parallax_eval_on is True:
    # for parallax
    hist_ax.set_xlim(-.5,.5)
  
  hist_ax.set_xlabel(data_ax_label)
  hist_ax.set_ylabel(hist_ax_label)
  hist_ax.legend()
  
  
  
  
  for i, data_set in enumerate(data_sub_sets_label):
    markers = ['.','.','+']    
    marker_size = [.5,.5,35]
    data = data_x[masks[i]]
    data_error = data_error_x[masks[i]]
    color_list = [dist_eval_array[0]['color'],'black','black']               
    label_list = dist_eval_array[1]['data_set']
    
    error_ax.scatter(\
                        data,\
                          data_error,\
                            label=dist_eval_array[i]['data_set'],\
                              marker= markers[i],\
                              s = marker_size[i],
                              # s=.5,\
                                # alpha = 1 / (i+1),\
                                alpha = 0.3,\
                                  # color= dist_eval_array[i]['color'],
                                  color= color_list[i] )  
  
  
  if error_selection is True:
    error_ax.axhline(error_upper_bound, color='black',  \
                   # label = dist_eval_array[1]['data_set'],
                            linewidth= 0.4, linestyle = '-')
  
  if parallax_eval_on is True:
                               
    # parallax 
    error_ax.set_ylim( 0,.6)


  error_ax.set_xlabel(data_error_ax_label)
  error_ax.set_ylabel(error_ax_label)
  error_ax.legend()
    
  if savefig is True:
    dist_fig.savefig(fig_sub_path+file_name+'_distEval.jpg')
  


  


  
  
  
  
  
  







  