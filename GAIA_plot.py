#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:31:33 2022

@author: Thomas Suntinger

Import recommendation:
  import GAIA_plot as gp

Functionality:
  
  - singleClusterMemberSelectionPlot()
      Creates a plot with spacial and CMD visualization of a cluster member
      selection.
  Includes two functions not ready for puplication standards:
    - cmd_plot() plots:
        - all data sets in one color magnitude diagram  
        - one color magnitude diagram for each data set
    - cmd_plot_full() takes a different input list format, restructures it
      and calls cmd_plot() to plot the data
    
Improvement:
  - Too messy, needs restructuring!
  - Only used for quick plots
  - For publication purposes cmds are created individually, 
    not with this functions!

"""

import matplotlib.pyplot as plt
import numpy as np

#%% ### ------------ Set default plot style------------------------------------
#https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure

# Note: fig_width is fixed in plotStyle.py for MasterThesis specifications
#       fig_height can and should be adjusted to special needs of individual plots
#       fig_height of 5inch is approximatelly half of height on A4 MasterThesis      
fig_height = 5
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)


def singleClusterMemberSelectionPlot(
                                     memberSelection,
                                     label_memberSelection,
                                     memberSelectionColor,
                                     
                                     noMemberSelection,
                                     label_noMemberSelection,
                                     noMemberSelectionColor,
                                     
                                     cluster_coord_center,
                                     cluster_search_radius,
                                     
                                     brightness_threshold_on = False,
                                     brightness_threshold = 0,
                                     
                                     savefig = True,
                                     fig_sub_path= 'gfx/',
                                     fig_name = 'cluster_selection'
                                     ):
  
  marker_size = 0.4
  marker_alpha = 0.4
  
  marker_size_lost = 30
  marker_alpha_lost = 1
  marker = '+'
  
  # import matplotlib.gridspec as gridspec
  
  pos_fig, pos_ax = plt.subplots(1,2,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[2,1]},sharex=False,sharey=False)
  
  # gs = gridspec.GridSpec(1, 2,width_ratios=[2,1])
  
  pos_ax[0].scatter(memberSelection['ra'],memberSelection['dec'],label=label_memberSelection,s=marker_size,alpha=marker_alpha,color=memberSelectionColor)
  pos_ax[0].scatter(noMemberSelection['ra'],noMemberSelection['dec'],label=label_noMemberSelection,marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=noMemberSelectionColor)
  
  circle_cluster = plt.Circle(cluster_coord_center, cluster_search_radius, color='black', fill=False)
  pos_ax[0].add_patch(circle_cluster)
  
  pos_ax[0].set_ylabel(r'dec')
  pos_ax[0].set_xlabel(r'ra')
  
  pos_ax[0].legend()
  pos_ax[0].set_aspect('equal', adjustable='box')
  
  
  
  
  pos_ax[1].scatter(memberSelection['bp_rp'],memberSelection['phot_g_mean_mag'],label=label_memberSelection, s=marker_size,alpha=marker_alpha,color=memberSelectionColor)  
  pos_ax[1].scatter(noMemberSelection['bp_rp'],noMemberSelection['phot_g_mean_mag'],label=label_noMemberSelection,marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=noMemberSelectionColor)  
  
  if brightness_threshold_on is True:
    pos_ax[1].axhline(brightness_threshold, linestyle = 'dashed', color = 'k',linewidth=0.5, zorder = 1.5)  
    
  
  pos_ax[1].set_xlabel('BP-RP [mag]')
  pos_ax[1].set_ylabel('phot_g_mean_mag [mag]')
  pos_ax[1].invert_yaxis()
  # pos_ax[1].legend(loc = 'lower center')
  
    
  # pos_ax[0,i].set_aspect('equal', adjustable='box')
  # pos_ax[1].set_aspect('equal')
  
    
  if savefig is True:
    pos_fig.savefig(fig_sub_path+fig_name+'_pos_cmd')

  return



def cmd_plot_full(full_data, \
             labels,\
               colors=[], \
                 data_background = [], \
                  label_data_background = [],\
                  color_data_background ='black',\
                   super_title = 'CMD', \
                    all_in_one = True, \
                     ready_to_publish = False,\
                       y_labels = [],\
                         markers = [],\
                           savefig = False,\
                             fig_sub_path = ''
                       ):
      
  data_list = []
  for i, el in enumerate(full_data):
    data_list.append( [full_data[i]['bp_rp'],full_data[i]['phot_g_mean_mag']] )
  
  if len(data_background) > 0:
    data_background = [data_background['bp_rp'],data_background['phot_g_mean_mag']]
    
  
  return cmd_plot(data_list, \
              labels=labels, 
               colors=colors, \
                data_background = data_background, \
                 label_data_background = label_data_background,\
                 color_data_background =color_data_background,\
                  super_title=super_title,
                    all_in_one=all_in_one,
                      ready_to_publish=ready_to_publish,\
                        y_labels = y_labels,\
                          markers = markers,\
                            savefig = savefig, \
                              fig_sub_path = fig_sub_path
                          )

  

def cmd_plot(data_list, \
             labels,\
               colors=[], \
                 data_background = [], \
                  label_data_background = [],\
                  color_data_background ='black',\
                   super_title = 'CMD', \
                    all_in_one = True, \
                     ready_to_publish = False,\
                       y_labels = [],\
                         markers =[],\
                           savefig = False,\
                             fig_sub_path = ''
                       ):
          
  if not any(markers):
    markers = np.full( shape=len(data_list), fill_value='.')
    marker_size = 1
  else:
    marker_size = 10
      
  if all_in_one is True:
    cmd_fig, cmd_ax = plt.subplots(1,1,sharex=True,sharey=True)
    if ready_to_publish is False:
      cmd_fig.suptitle(super_title)
    
    # if list is not empty
    if data_background: 
      if not label_data_background:
        label_data_background = 'full data set'
      cmd_ax.scatter(data_background[0],data_background[1],\
                     label= label_data_background ,\
                       s=1,color=color_data_background,\
                         alpha=0.3, zorder = 0.5,)
        
    for i, data in enumerate(data_list):
      if labels is float:
        label = '{0:s}: {1:1.2e}'.format(r'level',labels[i])
      else:
        label = labels[i]
      cmd_ax.scatter(data[0],data[1],label= label ,s=marker_size,color=colors[i],alpha=0.7, marker = markers[i])
    
    if y_labels:
      for i, el in enumerate(y_labels):
        cmd_ax[i].set_ylabel(el)  
    else:
      cmd_ax.set_ylabel('G-band mean magnitude [mag]')  
      
    cmd_ax.set_xlabel('BP-RP [mag]')
    cmd_ax.invert_yaxis()
    cmd_ax.legend(loc='lower right')

    
  else:
    n_plots = len(data_list)
    cmd_fig, cmd_ax = plt.subplots(1,n_plots,sharex=True,sharey=True)
    if ready_to_publish is False:
      cmd_fig.suptitle(super_title) 
    cmd_ax[0].invert_yaxis()
    
    if y_labels:
      for i, el in enumerate(y_labels):
        cmd_ax[i].set_ylabel(el)  
    else:
      cmd_ax[0].set_ylabel('G-band mean magnitude [mag]')  
    
            
    for i, data in enumerate(data_list):
      # if list is not empty
      if data_background: 
        if not label_data_background:
          label_data_background = 'full data set'
        cmd_ax[i].scatter(data_background[0],data_background[1],\
                       label= 'full data set' ,\
                         s=1,color=color_data_background,\
                           alpha=0.3, zorder = 0.5)
          
      if labels is float: 
        label = '{0:s}: {1:1.2e}'.format(r'Isopleth',labels[i])
      else:
        label = labels[i]
      cmd_ax[i].scatter(data[0],data[1],label= label ,s=marker_size,color=colors[i],alpha=0.7, marker = markers[i])
      cmd_ax[i].set_xlabel('BP-RP [mag]')
      cmd_ax[i].legend(loc='lower right')
  
  if savefig is True:
    cmd_fig.savefig(fig_sub_path+super_title)
  
  return [cmd_fig,cmd_ax,all_in_one]




