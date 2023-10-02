#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 19:40:59 2022

@author: tom

Improvements:
  - Clean unnecessary comments and lines of code
  - Write documentation and improve in code documentation.
  
  - Include calcSearchRadius() and convArcminToRad() functions from
    gaiaDistEval.py
  - Bonus: Call functions in gaiaParallaxSelection.py
"""

import csv
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from zero_point import zpt

from scipy.stats import kde

import scipy as scipy
from scipy import optimize

### ------------ Set default plot style----------------------------------------
#https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure


import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style()




#%% Lorentz profile fit function
# https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.optimize.curve_fit.html
# https://github.com/emilyripka/BlogRepo/blob/master/181119_PeakFitting.ipynb




def Lorentzian(x, amp1, cen1, wid1):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2))
  
def lorentz_profile_hist_fit(data,bins = 5000,plot=[],text='', show_text = False, show_plot = True):
  
  data_hist, data_hist_edges = np.histogram(data,bins = bins)
  
  lorentz_x = (data_hist_edges[:-1] + data_hist_edges[1:]) /2
  popt_lorentz, pcov_lorentz = \
   scipy.optimize.curve_fit( Lorentzian, lorentz_x, data_hist)
  
  resolution = bins*10
  lorentz_x = np.linspace(lorentz_x.min(), lorentz_x.max(),resolution)
  
  lorentz_y = Lorentzian(lorentz_x, *popt_lorentz)
  
  # Source: scipy.optimize.curve.fit
  # pcov_lorentz is the estimated covariance of popt. 
  # The diagonals provide the variance of the parameter estimate. 
  # To compute one standard deviation errors on the parameters use:
  # perr = np.sqrt(np.diag(pcov))
  
  # perr_lorentz = one standard deviation /1 Sigma errors
  perr_lorentz = np.sqrt(np.diag(pcov_lorentz))
  
  
  # Note: lor_width* 2 = full widht_half_max
  lor_ampl   = popt_lorentz[0]
  lor_center = popt_lorentz[1]
  lor_width  = popt_lorentz[2]
  
  lor_err_ampl   = perr_lorentz[0]
  lor_err_center = perr_lorentz[1]
  lor_err_width  = perr_lorentz[2]
  
  if show_text is True:
    print('')
    print("-------------Peak -------------")
    print(text)
    print ("amplitude = %0.3f (+/-) %0.3f" % (popt_lorentz[0], perr_lorentz[0]))
    print ("center    = %0.3f (+/-) %0.3f" % (popt_lorentz[1], perr_lorentz[1]))
    print ("width     = %0.3f (+/-) %0.3f" % (popt_lorentz[2], perr_lorentz[2]))
    print ("--------------------------------")

  if show_plot is True: 
    l_fig,l_ax = plt.subplots(1,1)
    # l_ax.hist(data, bins = bins)
    l_ax.hist(data, bins = data_hist_edges)
    l_ax.plot(lorentz_x,lorentz_y, color = 'orange')

  return [popt_lorentz,perr_lorentz,[lorentz_x,lorentz_y],[data_hist, data_hist_edges]]

#%% Lorentz profile peak evaluation using average
def lorentz_profile_average_in_width(data, popt_lorentz,perr_lorentz, search_factor = False, text='',show=False):
  
  # Note: lor_width* 2 = full widht_half_max
  lor_center = popt_lorentz[1]
  lor_width  = popt_lorentz[2]
  
  lor_center_err = perr_lorentz[1]
  lor_width_err  = perr_lorentz[2]
  
  
  if search_factor is False:
    search_factor = lor_width
    
  mask = (data <= lor_center+search_factor) & \
         (data >= lor_center-search_factor) 

  peak_average = np.average(data[mask])
  
    
  if show is True:
    print('')
    print('--------------Average--------------')
    print(text)
    print ("average center = %0.3f            " % (peak_average))
    print ("lorentz with   = (+/-) %0.3f" % (lor_width))
    print ("search factor =  (+/-) %0.3f" % (search_factor))
    print ("search area = %0.3f to  %0.3f" % (peak_average-search_factor, peak_average+search_factor))
    print ("--------------------------------")
  
  return [peak_average,search_factor,mask]




#%% calculate distance 
def distance_calc(parallax, parallax_error = 0,show=False,text=''):

  # parallax in milli arc sec per year (mas)
  # distance in parsec (kpc)
  distance = 1 / parallax
  
  # Source for error propagation: 
  # https://www.statology.org/error-propagation/
  # Note: Get another source ?
  
  # Right way to calculate error propagation
  # error_array = (parallax_error / parallax )**2
  # distance_err =  np.abs( distance ) * np.sqrt( np.sum(error_array) )
  
  # Simplified error propagation
  distance_err = np.abs( distance ) * (parallax_error / parallax )
  
  
  if show is True:
    
    par_unit = 'mas/year'
    dist_unit = 'kpc'
    
    label_par = '{0:s} at:{1:10.3f} (+/-) {3:0.3f} {2:s}'.format(r' parallax',parallax,par_unit,parallax_error)
    label_dist ='{0:s} at:{1:10.3f} (+/-) {3:0.3f} {2:s}'.format(r' distance',distance,dist_unit, distance_err)

    print('')
    print('--------------Distance--------------')
    print(text)
    print(label_par)
    print(label_dist)
    print ("--------------------------------")
    
  return [distance, distance_err]
    

#%% Test Script
# Testing functionality of parallax_evaluation.py
# Disable to decrease plot and print


#Load data
# [parallax_data,parallax_error_data] = np.load('selected_parallax_array.npy')


# par_bins = 5000
# label_parallax_data = 'test'

# [popt_lorentz,perr_lorentz,lorentz_curve] = \
#   lorentz_profile_hist_fit(\
#                                parallax_data,\
#                                  bins = par_bins,\
#                                    text=label_parallax_data,\
#                                      show_text=True,\
#                                        show_plot = False \
#                                          )
    
# [peak_average, search_factor, mask_parallax_data_parallax_range] = \
#   lorentz_profile_average_in_width(\
#                                        parallax_data, \
#                                          popt_lorentz,\
#                                            perr_lorentz,\
#                                              text=label_parallax_data,\
#                                                show= True)
    
# [distance, distance_err] = \
#   distance_calc(\
#                     peak_average, \
#                       parallax_error = search_factor,\
#                         text= label_parallax_data)

#%%



