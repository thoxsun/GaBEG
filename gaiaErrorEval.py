#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 21:29:04 2023

@author: tom

import as: 
  from gaiaErrorEval import trgb_conv_int
"""

import numpy as np
from scipy import stats

def trgb_conv_int(data, sigma = np.nan , cl=0.95):
  
  mean  = np.mean(data)
     
  if np.isnan(sigma):
    sigma = np.std(data)
  
    
  conf_interval_margin = stats.norm.interval(cl, loc=mean, scale=sigma ) 
  conf_interval_margin_delta = conf_interval_margin[1]- conf_interval_margin[0]
  delta  = conf_interval_margin_delta /2

  n_in_cl = ((data >= conf_interval_margin[0]) & (data <= conf_interval_margin[1])).sum() / float(data.size)
  print('{:0.2%} of the single draws are in conf_int'.format(n_in_cl))
  print('Mean and cl delta: {0:0.3} $\pm$ {1:0.3}'.format(mean,delta))
      
  return [mean, delta, conf_interval_margin]