#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Oct 28 14:39:13 2022

@author: Thomas Suntinger

Import recommendation:
  import plotStyle  

Functionality:
  - plotStyle can set global figure and axis configurations for plots 
    created by matplotlib. 
  - It can be executed as script
  - or by calling its function plotStyle.set_plot_style() 
  
Improvements:
  
  
"""


### ------------ Set default plot style----------------------------------------
#https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure
# 

fig_width = 5.47
fig_height = 5
  
def set_plot_style(fig_width = 5.47,fig_height = 5):
  
  import matplotlib.pyplot as plt
  
  SMALL_SIZE = 9
  # SMALL_SIZE = 13 
  MEDIUM_SIZE = 11 
  # MEDIUM_SIZE = 18 
  BIGGER_SIZE = 13 
  # BIGGER_SIZE = 20 
  
  plt.rc('font', size=SMALL_SIZE) 
  # controls default text sizes 
  
  plt.rc('axes', titlesize=BIGGER_SIZE) 
  # fontsize of the axes title 
  
  plt.rc('axes', labelsize=MEDIUM_SIZE) 
  # fontsize of the x and y labels 
  
  plt.rc('xtick', labelsize=MEDIUM_SIZE) 
  # fontsize of the tick labels 
  
  plt.rc('ytick', labelsize=MEDIUM_SIZE) 
  # fontsize of the tick labels 
  
  plt.rc('legend', fontsize=SMALL_SIZE) 
  # legend fontsize 
  
  plt.rc('legend', framealpha=1.) 
  # legend non transparent
  
  plt.rc('figure', titlesize=BIGGER_SIZE) 
  # fontsize of the figure title
  
  plt.rc('lines', linewidth=0.8)
  

  
  # Workaround for minus sign not included in "Computer Modern" font used by LaTex
  # Computer Modern = Serif.cmr10
  # https://github.com/matplotlib/matplotlib/pull/18397
  
  import matplotlib
  matplotlib.rcParams.update({
      'font.family': 'cmr10',
      'axes.formatter.use_mathtext': True})
  
  plt.rc('font', family='Serif',serif="cmr10") 
    
  # dots per inch
  plt.rc('figure', dpi=300 ) 
  
  # plt.rc('figure', autolayout = True ) 
  plt.rcParams['figure.constrained_layout.use'] = True
  
  plt.rcParams.update({ "text.usetex": True })
  
  #-----------------fig sizes -------------------------------------------------
  
  # single column width (21 picas, 8.9 cm, 3.5 inches)
  # plt.rc('figure', figsize=[3.5,fig_height]) 
  
  # a 1.5-column width (30 picas, 12.7 cm, 5 inches)
  # plt.rc('figure', figsize=[5,fig_height]) 
  
  # a full two-column width (43 picas, 18.2 cm, 7.2 inches)
  # plt.rc('figure', figsize=[7.2,fig_height]) 
  
  # Figurewith = 395pt = 5.47inch for MasterThesis Template
  # plt.rc('figure', figsize=[5.47,fig_height]) 
  
  #----------------------------------------------------------------------------
  plt.rc('figure', figsize=[fig_width,fig_height]) 
  
  
  ready_to_publish = True
  
  return [ready_to_publish, fig_height, fig_width]


import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style()