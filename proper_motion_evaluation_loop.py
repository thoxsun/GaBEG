#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 21:35:28 2022

@author: Thomas Suntinger

Import recommendation:
  import proper_motion_evaluation as pme
  
Input:
  - load: Gaia query data set as csv file (e.g. from GAIA_query_circle.py)

Dependencies:
  - csv
  - proper_motion_evaluation.py
  - GAIA_load.py
  
Dependencies for visualization:
  - plotStyle.py

Functionality:
  - Scirpt estimates probability density distribution in 2D phase space 
    (proper motion right ascension (pmra) and declination (pmdec)) 
    with a Gaussian Kernel Density Estimator (KDE) by calling 
    propability_density_estimation() function repeatedly.
  - Each new iteration uses pre selected data from iteration before to improve
    the results. 
  - A variety of plots can be made and saved, if adjoined parameters 
    are set to True (see proper_motion_evaluation.py for more details). 
    

Structure ( & parameters):
  - Init 
  - Plot settings 
      (figure height)
  - Load /Store data settings 
      (path-, file- & data set names)
  - toggle proper motion evaluation plots 
      (see proper_motion_evaluation.py)
  - proper motion evaluation settings 
      (scaling, iterations, labels of isopleths of interest)
  - round 0
  - additional iterations

Output:
  - store: csv files for each iteration with pm center values ([ra,dec] in mas)
  - print: pm center values ([ra,dec] in mas)
  - store: csv files for each iteration and isopleth of interest

Output plots:
  - curve of growth and KDE as default (see Structure:Plot settings)
  
  
Improvements:
  - Clean unnecessary comments and lines of code
  - Write documentation and improve in code documentation.
  
"""

import csv
import proper_motion_evaluation as pme


#%% ### Plot settings -------------------------------------------


fig_height = 4
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)




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

file_sub_path = 'data/pmData/'
from GAIA_save import save_gaia as save_gaia
save_data = True


    
#%% ### toggle proper motion evaluation plots ---------------------------------

proper_motion_evaluation_on = True
fig_name = file_name
fig_sub_path = 'gfx/pm_eval/'


# ------------------------ proper motion plot on/off---------------------------
if proper_motion_evaluation_on is True:
  pm_eval_plot_on = True
  cmd_from_pm_eval_plot_on = True

  levels_spacing_plot   = False
  contour_plot          = False
  proper_motion_surface_plot_on = False
  
  savefig = True

else:
  pm_eval_plot_on = False
  cmd_from_pm_eval_plot_on = False  
  
  levels_spacing_plot   = False
  contour_plot          = False
  proper_motion_surface_plot_on = False
  
  savefig = False
  
  
  
  
  
  
  
  
  
#%% ### proper motion evaluation settings ----------------------------  


# levels_scaling = 'lin' OR 'log' OR 'loginv'

# M5 needs 'log'
levels_scaling = 'log'  
levels_of_interesst = 0    


# Number of additional iterations to first one
iterations = 3

# Labels of isopleths in each iteration
iteration_list_labels = ['pm_iso_05','pm_iso_1']
iteration_list = []



#%% round 0
i = 0
fig_name = file_name +'_round'+str(i)

x = gaia_dr3_array['pmra']
y = gaia_dr3_array['pmdec']

[levels, density_sums, normed_density_sums, \
 mask_sample_in_levels_of_interesst, \
 levels_of_interesst, contour_labels, \
   [pm_ra_kde_center,pm_dec_kde_center], \
    [pm_ra_gauss_center,pm_dec_gauss_center] ] = \
pme.propability_density_estimation(x,y, \
                     levels_scaling             = levels_scaling, \
                     levels_of_interesst        = levels_of_interesst, \
                     curve_of_growth_plot       = pm_eval_plot_on ,\
                     contour_plot_of_interesst  = pm_eval_plot_on,\
                    
                     levels_spacing_plot        = levels_spacing_plot ,\
                     contour_plot               = contour_plot,\
                     proper_motion_surface_plot_on = proper_motion_surface_plot_on ,\
                     
                     ready_to_publish           = ready_to_publish,
                     savefig                    = savefig,
                     fig_name = fig_name,
                     fig_sub_path = fig_sub_path
                     )

pm_eval_name = file_name +'_round'+str(i) + '_pm-center.csv'
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
  
print ("pm kde center [ra,dec]:    [%0.3f , %0.3f] mas" % (pm_ra_kde_center, pm_dec_kde_center))
print ("pm gauss center [ra,dec]:  [%0.3f , %0.3f] mas" % (pm_ra_gauss_center, pm_dec_gauss_center))


# data sample with pm enclosed by isopleth with 1.0 prob.max
gaia_dr3_array_pm_iso_05 =  gaia_dr3_array[mask_sample_in_levels_of_interesst[0]]  
gaia_dr3_array_pm_iso_1 =  gaia_dr3_array[mask_sample_in_levels_of_interesst[1]]  
label_gaia_dr3_array_pm_iso_05 = contour_labels[0]

isos = [gaia_dr3_array_pm_iso_05,gaia_dr3_array_pm_iso_1]
iteration_list.append( isos )

for j, iso in enumerate(isos):
  iso_name = file_name +'_round'+str(i) + '_'+iteration_list_labels[j]
  save_gaia(iso,iso_name,file_sub_path=file_sub_path)




#%% additional rounds
while i < iterations:
  
  x = gaia_dr3_array_pm_iso_05['pmra']
  y = gaia_dr3_array_pm_iso_05['pmdec']
  
  # Needs to be reseted
  levels_of_interesst=0    
  
  i += 1
  fig_name = file_name +'_round'+str(i)

  
  import proper_motion_evaluation as pme
  [levels, density_sums, normed_density_sums, \
   mask_sample_in_levels_of_interesst, \
     levels_of_interesst, contour_labels, \
       [pm_ra_kde_center,pm_dec_kde_center], \
        [pm_ra_gauss_center,pm_dec_gauss_center] ] = \
    pme.propability_density_estimation(x,y, \
                         levels_scaling             = levels_scaling, \
                         levels_of_interesst        = levels_of_interesst, \
                         curve_of_growth_plot       = pm_eval_plot_on ,\
                         contour_plot_of_interesst  = pm_eval_plot_on,\
                        
                         levels_spacing_plot        = levels_spacing_plot ,\
                         contour_plot               = contour_plot,\
                         proper_motion_surface_plot_on = proper_motion_surface_plot_on ,\
                         
                         ready_to_publish           = ready_to_publish,
                         savefig                    = savefig,
                         fig_name = fig_name,
                         fig_sub_path = fig_sub_path
                         )
  
  print ("pm kde center [ra,dec]:    [%0.3f , %0.3f] mas" % (pm_ra_kde_center, pm_dec_kde_center))
  print ("pm gauss center [ra,dec]:  [%0.3f , %0.3f] mas" % (pm_ra_gauss_center, pm_dec_gauss_center))
  
  
  gaia_dr3_array_pm_iso_1 =  gaia_dr3_array_pm_iso_05[mask_sample_in_levels_of_interesst[1]] 
  gaia_dr3_array_pm_iso_05 =  gaia_dr3_array_pm_iso_05[mask_sample_in_levels_of_interesst[0]]  
  
  isos = [gaia_dr3_array_pm_iso_05,gaia_dr3_array_pm_iso_1] 
  iteration_list.append( isos )
  
  for j, iso in enumerate(isos):
    iso_name = file_name +'_round'+str(i) + '_'+iteration_list_labels[j]
    save_gaia(iso,iso_name,file_sub_path=file_sub_path)
 
  pm_eval_name = file_name +'_round'+str(i) + '_pm-center.csv'
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
 

# for k, isos in enumerate(iteration_list):
#   for j, iso in enumerate(isos):
#     iso_name = file_name +'_round'+str(k) + '_'+iteration_list_labels[j]
#     save_gaia(iso,iso_name,file_sub_path=file_sub_path)

  




