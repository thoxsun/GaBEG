#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:34:11 2023

@author: tom
"""




#%% Init

import csv
import numpy as np
# import os
# import glob
import matplotlib.pyplot as plt
# from zero_point import zpt

# from scipy.stats import kde

import GAIA_plot as gp

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

# dataDR3 - viaux13 convolution
t07 = 'orangered'

# Av extinction models
t08 = 'lightgrey'

# RGB
t09 = 'teal'
t_RGB = t07
# AGB
t10 = 'purple'
t_AGB = t09

# AGB RGB
t_ARGB = t10


#%% ### ------------ Load /Store data settings ---------------------------------------

folder_name = 'data/clusterMemberSelection/'

## -------------- GAIA data ---------------------------------------------------

data_set_selected = 'Gaia_M5_r-9arcmin_zpc_pm_clusterMemberSelection_distEval_absG-band_mbol_bcLinAprx'
file_name = data_set_selected

## ------------ load data into variables --------------------------------------

from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

gaia_cluster_members = file_data
gaia_cluster_members.sort(order = 'mbol_mag') 
label_gaia_cluster_members = 'cluster members'
label_gaia_cluster_members_n = 'cluster members ('+str(gaia_cluster_members.size)+' stars)'

print()
print(label_gaia_cluster_members_n)


#---------Load TRGB_M5_viaux13.csv --------------------------------
viaux13_set_selected = 'TRGB_M5_Viaux13.csv'

viaux13_data = load_query('data/trgbEvalData/' + viaux13_set_selected)
label_viaux13_data = 'viaux13'
label_viaux13_n = label_viaux13_data+' ('+str(viaux13_data.size)+' stars)'

print()
print(label_viaux13_n)



#---------Store GAIA data --------------------------------

file_sub_path = 'data/trgbEvalData/'
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### ----------- Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/trgb_eval/'


# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)

viaux13_full_data_load_and_print  = False
empirical_fit_cmd_plot            = False
hist_cmd_rgb_sel_plot             = False
pos_comparison_plot               = False
hist_cmd_pos_rgb_sel_plot         = True




#%% Viaux13 GC M5 members 

# Note: Used, to verify coordinates of four brightest stars in Viaux13 TRGB selection.

if viaux13_full_data_load_and_print is True:

  #---------Load photm5_Viaux13_catalog.csv --------------------------------
  viaux13_all_stars_name = 'photm5_Viaux13_catalog.csv'
  
  viaux13_all_stars = load_query('data/trgbEvalData/' + viaux13_all_stars_name)
  label_viaux13_all_stars = 'viaux13_all'
  label_viaux13_all_stars_n = label_viaux13_all_stars+' ('+str(viaux13_all_stars.size)+' stars)'
  
  print()
  print(label_viaux13_all_stars_n)
  
  viaux13_all_stars['Imag'].sort()
  
  mask = np.logical_not((viaux13_all_stars['Vmag'] == 99.999) | (viaux13_all_stars['Imag'] == 99.999))
  data_list = []
  data_list.append( 
    [
     (viaux13_all_stars['Vmag']-viaux13_all_stars['Imag'])[mask],
      viaux13_all_stars['Imag'][mask]
    ] 
    )
  # --------Plot CMD ------------------------------------------------------
  gp.cmd_plot(data_list, \
                labels='Viaux13',\
                  colors='k', \
                 #   data_background = [], \
                 #    label_data_background = [],\
                 #    color_data_background ='black',\
                      super_title = 'CMD Viaux13', \
                 #      all_in_one = True, \
                 #       ready_to_publish = False,\
                          y_labels = 'I',\
                            x_labels = 'V-I',\
                            # markers =['.'],\
                 #             savefig = False,\
                 #               fig_sub_path = ''
                         )


#%% viaux13 cross referencation



# Transfomr Viaux13 coordinates to Gaia coordinates ---------------------------

viaux13_pos_trans = np.zeros(np.size(viaux13_data))

# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn
viaux13_data = rfn.append_fields(viaux13_data, 'ra_deg_J2016', viaux13_pos_trans, dtypes='<f8')
viaux13_data = rfn.append_fields(viaux13_data, 'dec_deg_J2016', viaux13_pos_trans, dtypes='<f8')

#------------------------------------------------------------------
#convert coordinates into diverent equinox (=epoch)
# https://docs.astropy.org/en/stable/coordinates/transforming.html

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic, FK5, ICRS
from astropy.time import Time

for i, el in enumerate(viaux13_data):
  fk5c = SkyCoord(el['ra'], el['dec'],frame=FK5(equinox=Time('J2000')))
  fk5_2016 = FK5(equinox=Time(2016, format='jyear'))
  z = fk5c.transform_to(fk5_2016)
  
  z = fk5c.transform_to(ICRS)
  viaux13_data['ra_deg_J2016'][i]  = z.ra.deg
  viaux13_data['dec_deg_J2016'][i] = z.dec.deg

# '02h31m49.09s', '+89d15m50.8s'  
# fk5c = SkyCoord('02h31m49.09s', '+89d15m50.8s',
#                 frame=FK5(equinox=Time('J2000')))



# cross-match functions---------------------------- ---------------------------

def points_in_circle(x, y, center_x, center_y, radius):
    distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    return [distance, distance <= radius]


def gaiaCrossRef(data_set, star_ra,star_dec, radius = 1/3600):
  
  [distance, mask] = points_in_circle(
                                      data_set['ra'], 
                                      data_set['dec'],
                                      star_ra, 
                                      star_dec, 
                                      radius
                                      )
  
  return [data_set[mask],distance[mask],mask,distance.min()]

# Cross-match Viaux13 coordinates to Gaia coordinates -------------------------

viaux13_crossmatch = []
viaux13_crossmatch_details = []
viaux13_crossmatch_offset = np.zeros([np.size(viaux13_data),4])

offset_calibration = False
if offset_calibration is True:
  viaux13_crossmatch_offset = np.zeros([np.size(viaux13_data)+1,4])
  offset_ra_calibration = 0.0
  offset_dec_calibration = 0.0

# cross match radius in degree:
# crossmatch_radius = 1/3600 # 1 arcsec
crossmatch_radius = 0.07 # 0.07*3600
# crossmatch_radius = 1/60 #10 arcsec


viaux_I = np.zeros_like(gaia_cluster_members['ra'])
viaux_I[:] = np.nan
viaux_nr = np.zeros_like(gaia_cluster_members['ra'],dtype=np.uint)
# viaux_nr [:] = np.nan

for i, el in enumerate(viaux13_data):
  [data_set,distance,mask,distance_min] = gaiaCrossRef(
                                                      gaia_cluster_members, 
                                                      el['ra_deg_J2016'],
                                                      el['dec_deg_J2016'], 
                                                      radius=crossmatch_radius
                                                      )
  
  
  gaia_star_sel = data_set[(distance == distance_min)]
  
  mask = (gaia_star_sel['source_id'] == gaia_cluster_members['source_id'])
  (viaux_I)[mask] = el['ap_I_mag']
  (viaux_nr)[mask] = el['Nr']
  
  offset_ra  = gaia_star_sel['ra'] - el['ra_deg_J2016']
  offset_dec = gaia_star_sel['ra'] - el['ra_deg_J2016']
  viaux13_crossmatch_offset[i,:] = el['Nr'],offset_ra, offset_dec, distance_min
  
  # viaux13_crossmatch_details.append([el,data_set,distance,mask,distance_min])
  viaux13_crossmatch_details.append([el,gaia_star_sel,distance_min, offset_ra, offset_dec] )
  viaux13_crossmatch.append(gaia_star_sel )
  
  
  
  if offset_calibration is True:
    if i < 2:
      offset_ra_calibration += offset_ra
      offset_dec_calibration += offset_dec
    elif i == 2 :
      offset_ra_calibration /= 2
      offset_dec_calibration /= 2
      
      [data_set,distance,mask,distance_min] = gaiaCrossRef(
                                                          gaia_cluster_members, 
                                                          el['ra_deg_J2016'] -offset_ra_calibration,
                                                          el['dec_deg_J2016']-offset_dec_calibration, 
                                                          radius=crossmatch_radius
                                                          )
      offset_ra  = gaia_star_sel['ra'] - el['ra_deg_J2016'] -offset_ra_calibration
      offset_dec = gaia_star_sel['ra'] - el['ra_deg_J2016'] -offset_dec_calibration
      viaux13_crossmatch_offset[i+1,:] = el['Nr'],offset_ra, offset_dec, distance_min
      
      # viaux13_crossmatch_details.append([el,data_set,distance,mask,distance_min])
      viaux13_crossmatch_details.append([el,gaia_star_sel,distance_min, offset_ra, offset_dec] )
      viaux13_crossmatch.append(gaia_star_sel )
    
  
# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'Viaux13_Nr', viaux_nr, dtypes='<f8')
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'Viaux13_I', viaux_I, dtypes='<f8')

  
  
#%% RGB AGB selection preparation ---------------------------------------------  



# preparation of data and bounds ----------------------------------------

# bolometric magnitude Mbol
# see gaiaBCg_approximation.py
mbol_trgb_lower_bound = -2.4030841435619883
# mbol_trgb_lower_bound     = -2.40
mask_bc_aprx = gaia_cluster_members['bcg_aprx_flag'] == True
mbol_trgb_lower_bound_unc = gaia_cluster_members['mbol_mag_error'][mask_bc_aprx][0]


# Manualy chosen to exclude HB in Mbol CMD
upper_bound_hor_branch  = -0.1 
  

# Gaia DR3 with BCg derivation only (no BC approximation included)
mask_bcg = np.logical_not(gaia_cluster_members['bcg_aprx_flag'])
gaia_cluster_members_bc = gaia_cluster_members[mask_bcg]
gaia_cluster_members_bc_rgb = gaia_cluster_members_bc[(gaia_cluster_members_bc['mbol_mag'] < mbol_trgb_lower_bound)]
print()
print('gaia_cluster_members_bc: '+str(gaia_cluster_members_bc.size))
print('gaia_cluster_members_bc_rgb: '+str(gaia_cluster_members_bc_rgb.size))
print()



# Gaia DR3 selection for separation evaluation ----------------------------

# Chose lower Mbol thresholds 
mbol_rgb_agb_lower_bound = upper_bound_hor_branch
# mbol_rgb_agb_lower_bound = mbol_trgb_lower_bound

# Create a reference (dr3_above_lower_bound_selection) to the chosen data set
dr3_above_lower_bound_selection  = gaia_cluster_members[(   gaia_cluster_members['mbol_mag'] < mbol_rgb_agb_lower_bound)]
# dr3_above_lower_bound_selection  = gaia_cluster_members_bc[(gaia_cluster_members_bc['mbol_mag'] < mbol_rgb_agb_lower_bound)]

# --------------------------------------------------------------------------

# sort array with brightest stars on top
dr3_above_lower_bound_selection.sort(order = 'mbol_mag')
dr3_selection_brightest = dr3_above_lower_bound_selection['mbol_mag'][0]
dr3_selection_lowest    = dr3_above_lower_bound_selection['mbol_mag'][-1]




#%% ### RGB AGB selection ----------------------------------------------------
# Note: Make RGB AGB fit and separation function 
#          from this section with data_sample as input

# def rgb_agb_separation(data_sample):
if True:
  data_sample = dr3_above_lower_bound_selection
  
  # sort array with brightest stars on top
  data_sample.sort(order = 'mbol_mag')
  data_sample_brightest = data_sample['mbol_mag'][0]
  data_sample_dimmest    = data_sample['mbol_mag'][-1]
# -------------------------------------------------------------------  


# RGB star color t09
# AGB star color t10

# x = np.linspace(0.8, 2.1,1000)
# x = np.linspace(0.8, 2.1,1000)

b = 1.95+0.16+0.01
# b = 1.95
x = np.linspace(1.1, b,1000)


# brightest_rgb_star = 11.415 from cmd 
y_1 = data_sample_brightest
# y_1 = dr3_selection_brightest

# TRGB estimate
# y_0 = y_1 - 0.04
y_0 = y_1 - 0.1

# RGB empirical fit from Viaux et.al 2013 (DOI: 10.1051/0004-6361/201322004)
# I = I0 + 3.83 * (1.95 - col_temp)**2.5
a = 3.83-0.5+0.7
a = 3.83-0.2
b = b
c = 2.5-0.2
# c = 2.5

y = y_0 + a * (b - x)**c
# y = y_0 + (3.83-0.5) * (1.95+0.16 - x)**2.5
# y = y_0 + (3.83-0.0) * (1.95 +0.0 - x)**2.5

# diff equation
y_d = - c * a * (b - x)**(c-1)
# y_d = - 2.5 * (3.83-0.5) * (1.95+0.16 - x)**1.5
# y_d = - 2.5 * (3.83-0.0) * (1.95+0.00 - x)**1.5

# normal line of 
# d_r = 0.1
# d_r = 0.3
d_r = 0.35

# d_r_agb = 0.166
# d_r_agb = 0.17
# d_r_rgb = 0.035
# d_r_agb = 0.095
d_r_agb = d_r

# d_r_rgb = 0.035
# d_r_rgb = 0.0656

d_r_left = 0.35
d_r_right = 0.2
d_r_right = 10.2

d_r_rgb_left  = 0.03
d_r_rgb_right_true = 0.1
d_r_rgb_right = 10.1


# sx = data_sample['bp_rp']
# sy = data_sample['phot_g_mean_mag']

sx = data_sample['bp_rp']
sy = data_sample['mbol_mag']

s_r = np.zeros_like(sx)


data_sample_rgb_agb_selection_mask = np.full(sx.size, False)
data_sample_rgb_selection_mask = np.full(sx.size, False)
data_sample_agb_selection_mask = np.full(sx.size, False)

for i, el in enumerate(data_sample):
  d_x = (y - sy[i])
  d_y = (x - sx[i])
  
  l = np.sqrt(d_x**2+d_y**2)
    
  mask = (np.equal(l,l.min()))
  if np.greater_equal(d_y[mask], 0):
    # star right from fit-> negative (y-axis inverted in cmd)
    s_r[i] = -l.min()
  else:
    # star left from fit -> positive (y-axis inverted in cmd)
    s_r[i] = l.min()
    
  

  
  # if np.less_equal(l.min(), d_r):
  if np.less_equal(s_r[i], d_r_right ) & np.greater_equal(s_r[i], -d_r_left ):
    data_sample_rgb_agb_selection_mask |= np.isin(data_sample['source_id'],el['source_id'] )
    
    # Print distance to empirical fit 
    # print(str(el['source_id'])+': '+str(s_r[i]))
    
    # print absolute magnitude 
    # print(str(el['source_id'])+': '+str(data_sample['mbol_mag'][i]))
  
  # if np.less_equal(l.min(), d_r_rgb ):
  #   data_sample_rgb_selection_mask |= np.isin(data_sample['source_id'],el['source_id'] )
  if np.less_equal(s_r[i], d_r_rgb_right ) & np.greater_equal(s_r[i], -d_r_rgb_left ):
    data_sample_rgb_selection_mask |= np.isin(data_sample['source_id'],el['source_id'] )
  elif np.less_equal(s_r[i], -d_r_rgb_left ) & np.greater_equal(s_r[i], -d_r_agb ):
    
    data_sample_agb_selection_mask |= np.isin(data_sample['source_id'],el['source_id'] )
    # Print distance to empirical fit 
    # print(str(el['source_id'])+': '+str(s_r[i]))
    
    # print absolute magnitude 
    # print(str(el['source_id'])+': '+str(data_sample['mbol_mag'][i]))




# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn
data_sample = rfn.append_fields(data_sample, 'fit_delta', s_r)

# Note: data_sample is a reference on one of following arrays 
#       created fields are therefore added to the original array
# data_sample = dr3_above_lower_bound_selection
# data_sample = gaia_cluster_members
# data_sample = gaia_cluster_members_bc


# Store selections ------------------------------------------------------------

    
# RGB selection
label_rgb_selection     = 'delta left from empirical fit : '+ str(-d_r_rgb_left)
rgb_selection           = data_sample[data_sample_rgb_selection_mask]
rgb_selection_distance  = s_r[data_sample_rgb_selection_mask]
save_gaia(rgb_selection, file_name +'_rgb_selection' ,file_sub_path=file_sub_path)

# AGB selection
label_agb_selection     = str(-d_r_agb)+ '< delta fit < '+ str(-d_r_rgb_left)
agb_selection           = data_sample[data_sample_agb_selection_mask]
agb_selection_distance  = s_r[data_sample_agb_selection_mask]
save_gaia(rgb_selection, file_name +'_agb_selection' ,file_sub_path=file_sub_path)

# AGB RGB selection d_r < 0.35
label_agb_rgb_selection = 'delta to empirical fit: '+ str(d_r)
agb_rgb_selection = data_sample[data_sample_rgb_agb_selection_mask]
agb_rgb_selection_distance = s_r[data_sample_rgb_agb_selection_mask]
save_gaia(rgb_selection, file_name +'_agb_rgb_selection' ,file_sub_path=file_sub_path)

agb_rgb_selection_brightest = agb_rgb_selection['mbol_mag'][0]
agb_rgb_selection_lowest    = agb_rgb_selection['mbol_mag'][-1]

trgb_mbol = agb_rgb_selection_brightest
trgb_star = agb_rgb_selection[0]
print('Brightest star in the mix: '+str(agb_rgb_selection['mbol_mag'].min()))
print('Brightest star in the mix: '+str(trgb_mbol))



#%% RGB eval fit plot 

if empirical_fit_cmd_plot is True:
    list_comparison = []
    comparison_labels = [label_gaia_cluster_members,'empirical RGB fit']
    comparison_colors = [t02,t09]
    comparison_markers = ['']
    
    data_background = [gaia_cluster_members['bp_rp'],gaia_cluster_members['mbol_mag']]
    list_comparison.append([agb_rgb_selection['bp_rp'],agb_rgb_selection['mbol_mag']])
    list_comparison.append([x,y])
    # list_comparison.append([x,y_d])

    
    import GAIA_plot as gp
    
    # Note: cmd plot savefig in function!!!
    [cmd_fig,cmd_ax,all_in_one] = gp.cmd_plot(list_comparison, \
                labels=comparison_labels, 
                 colors=comparison_colors, \
                    data_background = data_background,\
                    color_data_background = t00,\
                    label_data_background = label_gaia_cluster_members,\
                    super_title=file_name+'_CMD_clust_sel_fit',
                     markers = comparison_markers,
                      all_in_one=True,
                        ready_to_publish=True,
                          savefig = savefig,
                            fig_sub_path = fig_sub_path
                          )
      
#%% Store selection as csv table for LaTex

star_list = []
list_header = [
               r'Viaux13 Nr', 
               r'$ M_{I} $', 
               r'GAIA DR3 ID', 
               r'G3 Nr', 
               r'$ m_{G} $',
               r'$ M_{bol} $ ',
               r'$\Updelta M_{bol}$ ',
               r'$ \Updelta fit $'
               ]
star_list.append(list_header)
for i, el in enumerate(agb_rgb_selection[0:4]):
  star_list.append([
                    str(el['Viaux13_Nr']),\
                    '{0:5.3f}'.format(el['Viaux13_I']),\
                    el['source_id'] ,\
                    i+1 ,\
                    '{0:5.3f}'.format(el['phot_g_mean_mag']),\
                    '{0:5.2f}'.format(el['mbol_mag']),\
                    '{0:5.2f}'.format(el['mbol_mag_error']),\
                    '{0:5.3f}'.format(el['fit_delta']),\
                      ]
                    )
# create new csv file
with open(file_sub_path + data_set_selected+'_trgb_candidates_selection.csv', 'w', newline='') as csvfile:
  # create new csv object
  csv_write = csv.writer(csvfile, delimiter=',')
  # store data in new csv file
  csv_write.writerows(star_list)
  
#%% Store slection as LaTex formatted text file 
#Note: paste text from file to latex document
list_header = [
               r'Viaux13 Nr', 
               r'$ M_{I} $', 
               r'GAIA DR3 ID', 
               r'G3 Nr',
               r'$ m_{G} $',
               r'$ M_{bol} $ ',
               r'$\Updelta M_{bol}$ ',
               r'$ \Updelta fit $'
               ]

list_header_laTex =''
for i, el in enumerate(list_header):
  list_header_laTex += el 
  list_header_laTex += ' & '
list_header_laTex += ' \\\ '

star_list_laTex = []
star_list_laTex.append([list_header_laTex])

for i, el in enumerate(agb_rgb_selection[0:4]):
  star_list_laTex.append([
      str(el['Viaux13_Nr'])                         + ' & ' + 
      '{0:5.3f}'.format(el['Viaux13_I'])            + ' & ' + 
      str(el['source_id'])                          + ' & ' + 
      str(i+1)                                      + ' & ' + 
      '{0:5.3f}'.format(el['phot_g_mean_mag'])      + ' & ' + 
      '{0:5.2f}'.format(el['mbol_mag'])             + ' & ' +
      '{0:5.2f}'.format(el['mbol_mag_error'])       + ' & ' +
      '{0:5.3f}'.format(el['fit_delta'])            + ' \\\ '
      ])

# create new txt file
with open(file_sub_path + data_set_selected+'_trgb_candidates_selection.txt', 'w', newline='') as txtfile:
  # create new csv object
  csv_write = csv.writer(txtfile, delimiter=',')
  # store data in new csv file
  csv_write.writerows(star_list_laTex)



#%% Histogram + CMD comparison

    
    # data_background = [gaia_cluster_members['bp_rp'],gaia_cluster_members['phot_g_mean_mag']]
    # list_comparison.append([agb_rgb_selection['bp_rp'],agb_rgb_selection['phot_g_mean_mag']])
if hist_cmd_rgb_sel_plot is True:

  marker_size = 0.4
  marker_alpha = 0.4
  
  marker_size_lost = 30
  marker_alpha_lost = 1
  marker = '+'
  
  n_bin = 28*2
  
  pos_fig, pos_ax = plt.subplots(1,2,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[1,1]},sharex=False,sharey=False)
  
  pos_ax_0 = pos_ax[0]  
  pos_ax_1 = pos_ax[1]  
  # pos_ax_2 = pos_ax[2]  
  
  # pos_ax_0.set_ylabel('N')
  
  
  # pos_ax_0.set_ylim(0,par_prob_fit.max()+0.001)
  
  
  # pos_ax_0 = pos_ax_0.twinx()
  label= 'selected RGB AGB stars'
  n,bins,patches = pos_ax_0.hist(agb_rgb_selection_distance,color = t07, bins = n_bin, label=label,alpha=0.6)
  
  print('bin range: '+'{0:0.2f} mag with {1:1.0f} bins'.format((bins[1]-bins[0]),n_bin))
  
  # pos_ax_0.hist(par_kde_1D,color = t06, bins = par_bins,label=label_par_kde_1D)
  
  # pos_ax_0.hist(par_bound_array['parallax'],color = t07, bins = par_bins,\
                        # label=par_bound_name
                        # )
  
  
  pos_ax_0.axvline(0,
                        label = 'empirical fit',\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '-' \
                            )
  
  # pos_ax_0.axvline(par_max,
  #                       # label = par_bound_name,\
  #                       color = 'k',\
  #                         linewidth=0.6, \
  #                           zorder = 4,\
  #                         linestyle = '--' \
  #                           )
  
  # pos_ax_0.axvline(par_min,
  #                       # label = par_bound_name,\
  #                       color = 'k',\
  #                         linewidth=0.6, \
  #                           zorder = 4,\
  #                         linestyle = '--' \
  #                           )
  
    
    
  # delta = 1.5
  # pos_ax_0.set_xlim(par_x_center-delta,par_x_center+delta)
  pos_ax_0.set_xlabel('distance [mag]')
  pos_ax_0.set_ylabel('count of stars [N]')
  # pos_ax_0.legend(loc='upper left')
  # pos_ax_0.legend()
  # pos_ax_0.legend()
  # pos_ax_0.set_aspect('equal', adjustable='box')
  
  
  
  # pos_ax[1].plot(x,y-brightness_modul,label = 'empirical fit', color = t_fit)
  pos_ax[1].plot(x,y,label = 'empirical fit', color = t_fit)
  
  # All cluster members
  label = 'cluster members'
  pos_ax[1].scatter(gaia_cluster_members['bp_rp'],gaia_cluster_members['mbol_mag'],label=label,s=1,color=t04)
  
  # cluster members above threshold
  # label = '{0:s}: {1:2.1e}'.format(r'cluster members above: ',mbol_rgb_agb_lower_bound)
  # pos_ax[1].scatter(dr3_above_lower_bound_selection['bp_rp'],dr3_above_lower_bound_selection['mbol_mag'],label=label,s=1,color=t04)
  
  # rgb selection derived by fit
  label= 'selected RGB AGB stars'
  # pos_ax[1].scatter(agb_rgb_selection['bp_rp'],agb_rgb_selection['mbol_mag'],color=t09,label= label,s=1)   
  pos_ax[1].scatter(agb_rgb_selection['bp_rp'],agb_rgb_selection['mbol_mag'],label = label, marker = '+',s=10,alpha=1,color=t07)
  
  # # rgb agb selection with BC
  # label= 'BC RGB AGB stars'
  # pos_ax[1].scatter(gaia_cluster_members_bc['bp_rp'],gaia_cluster_members_bc['mbol_mag'],label = label, marker = 'o',facecolor = 'none',s=10,alpha=1,color='lightseagreen')
  
  # Viaux13 crossmatch
  label= 'Viaux13 stars'
  for el in  viaux13_crossmatch: 
    pos_ax[1].scatter(el['bp_rp'],el['mbol_mag'],label = label, marker = 'o',facecolor = 'none',s=10,alpha=1,color='lightseagreen')
  
  
  # Threshold as line
  pos_ax[1].axhline(mbol_rgb_agb_lower_bound, linestyle = 'dashed', label= 'bound: '+str(mbol_rgb_agb_lower_bound)+' mag', color = 'k',linewidth=0.5, zorder = 1.5)  
  
  # pos_ax[1].set_ylim(y_0-0.2,-1.5)
  # pos_ax[1].set_xlim(0.6,2.2)
  
  pos_ax[1].set_xlabel('BP-RP [mag]')
  pos_ax[1].set_ylabel('mbol_mag [mag]')
  pos_ax[1].invert_yaxis()
  # pos_ax[1].legend(loc = 'lower center')
  
    
  # pos_ax[0,i].set_aspect('equal', adjustable='box')
  # pos_ax[1].set_aspect('equal')
  pos_ax[1].set_aspect(.5)
  
  # fig_sub_path = 'gfx/clus_sel/'
  fig_name = data_set_selected
  
  if savefig is True:
    pos_fig.savefig(fig_sub_path+fig_name+'_rgb_sel.jpg')
  
  
  
  
  

#%% ### TRGB selection spatial distribution and CMD plot ------------

if pos_comparison_plot is True:
   
  

  # Cluster center coordinates
  coord_center_M5 = 229.638417,2.081028
  # Search cone radius in unit degrees
  radius_M5 = 0.15
  
  marker_size = 0.4
  marker_alpha = 0.4
  
  marker_size_lost = 30
  marker_alpha_lost = 1
  marker = '+'
  
  
  # import matplotlib.gridspec as gridspec
  
  pos_fig, pos_ax = plt.subplots(1,2,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[2,1]},sharex=False,sharey=False)
  
  # gs = gridspec.GridSpec(1, 2,width_ratios=[2,1])
  
  pos_ax[0].scatter(agb_rgb_selection['ra'],agb_rgb_selection['dec'],label=label_agb_rgb_selection,s=marker_size,alpha=marker_alpha,color=t_ARGB)
  # pos_ax[0].scatter(gaia_dr3_array_radVel_excluded['ra'],gaia_dr3_array_radVel_excluded['dec'],label=label_gaia_dr3_array_radVel_excluded,marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t02)
  pos_ax[0].scatter(agb_rgb_selection['ra'][0:4],agb_rgb_selection['dec'][0:4],label='',marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t02)
  
  for i, el in enumerate(viaux13_crossmatch):
    
    if offset_calibration is True:
      label_offset_cal = str(viaux13_crossmatch_offset[i,0])
    else:
      label_offset_cal = viaux13_data['Nr'][i]
      
    pos_ax[0].scatter(el['ra'],el['dec'],label='',marker = 'o',facecolor = 'none',s=10,alpha=1,color='lightseagreen')
    # add text to markers with plt.annotate()  
    pos_ax[0].annotate( label_offset_cal, \
                     xy=(el['ra'],el['dec']), \
                     textcoords='offset pixels', \
                     xytext=(15-15*i,30),
                     size= 14
                    )
      
  circle_cluster = plt.Circle(coord_center_M5, radius_M5, color='black', fill=False)
  pos_ax[0].add_patch(circle_cluster)
  
  pos_ax[0].set_ylabel(r'dec')
  pos_ax[0].set_xlabel(r'ra')
  
  # pos_ax[0].legend()
  pos_ax[0].set_aspect('equal', adjustable='box')
  
  
 
  
  pos_ax[1].scatter(agb_rgb_selection['bp_rp'],agb_rgb_selection['mbol_mag'],label=label_agb_rgb_selection, s=marker_size,alpha=marker_alpha,color=t_ARGB)  
  # pos_ax[1].scatter(gaia_dr3_array_radVel_excluded['bp_rp'],gaia_dr3_array_radVel_excluded['phot_g_mean_mag'],label=label_gaia_dr3_array_radVel_excluded,marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t02)  
  pos_ax[1].scatter(agb_rgb_selection['bp_rp'][0:4],agb_rgb_selection['mbol_mag'][0:4],label='',marker = marker,s=marker_size_lost,alpha=marker_alpha_lost,color=t02)  
  
  for i, el in enumerate(viaux13_crossmatch):
    pos_ax[1].scatter(el['bp_rp'],el['mbol_mag'],label='',marker = 'o',facecolor = 'none',s=10,alpha=1,color='lightseagreen')  
    
    if offset_calibration is True:
      label_offset_cal = str(viaux13_crossmatch_offset[i,0])
    else:
      label_offset_cal = viaux13_data['Nr'][i]
      
    pos_ax[1].annotate( 
                      label_offset_cal, \
                     xy=(el['bp_rp'],el['mbol_mag']), \
                     textcoords='offset pixels', \
                     xytext=(-15,5),
                     size= 14
                    )
  pos_ax[1].axhline(mbol_rgb_agb_lower_bound, linestyle = 'dashed', color = 'k',linewidth=0.5, zorder = 1.5)  
    
  
  pos_ax[1].set_xlabel('BP-RP [mag]')
  pos_ax[1].set_ylabel('mbol_mag [mag]')
  pos_ax[1].invert_yaxis()
  # pos_ax[1].legend(loc = 'lower center')
  
    
  # pos_ax[0,i].set_aspect('equal', adjustable='box')
  pos_ax[1].set_aspect('equal')
 
  
  
  fig_name = data_set_selected
  
  if savefig is True:
    pos_fig.savefig(fig_sub_path+fig_name+'_TRGB_pos_cmd')





#%% Pos + Histogram + CMD comparison

    
    # data_background = [gaia_cluster_members['bp_rp'],gaia_cluster_members['phot_g_mean_mag']]
    # list_comparison.append([agb_rgb_selection['bp_rp'],agb_rgb_selection['phot_g_mean_mag']])
if hist_cmd_pos_rgb_sel_plot is True:

  marker_size = 0.4
  marker_alpha = 0.4
  
  marker_size_lost = 30
  marker_alpha_lost = 1
  marker = '+'
  
  n_bin = 28 *2
  
  # Cluster center coordinates
  coord_center_M5 = np.array([229.638417,2.081028])
  # Search cone radius in unit degrees
  radius_M5 = 0.15
  
  
  # multi_fig, multi_ax = plt.subplots(2,2,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[1,1]},sharex=False,sharey=False)
  # pos_ax = multi_ax[0]  
  # hist_ax = multi_ax[1]  
  # cmd_ax = multi_ax[2]  
  
  # multi_fig, multi_ax = plt.subplots(figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[1,1]})
  multi_fig = plt.figure(figsize=(fig_width,fig_height*.6))
  
  from matplotlib.gridspec import GridSpec
  gs = GridSpec(2, 2, figure=multi_fig,width_ratios=[1, 2], height_ratios=[1,2])
  
  # First subplot: top-left
  pos_ax  = multi_fig.add_subplot(gs[1,0])
  hist_ax = multi_fig.add_subplot(gs[0,0])
  cmd_ax  = multi_fig.add_subplot(gs[:,1])
  
  
  #pos_ax-----------------------------------------------------------
  # pos_ax.scatter(agb_rgb_selection['ra'],agb_rgb_selection['dec'],
  #                label=label_agb_rgb_selection,
  #                # s=marker_size,
  #                 s=1,
  #                 alpha=marker_alpha,
  #                color=t_ARGB)
  
  # pos_ax.scatter(rgb_selection['ra'],rgb_selection['dec'],
  #                label=label_rgb_selection,
  #                # s=marker_size,
  #                s=1,
  #                 alpha=marker_alpha,
  #                color=t_RGB)
  
  # pos_ax.scatter(agb_selection['ra'],agb_selection['dec'],
  #                label=label_agb_selection,
  #                # s=marker_size,
  #                s=1,
  #                 alpha=marker_alpha,
  #                color=t_AGB)
  
  pos_ax.scatter(dr3_above_lower_bound_selection['ra'] -coord_center_M5[0],
                 dr3_above_lower_bound_selection['dec']-coord_center_M5[1],
                  label='stars above threshold',
                 # s=marker_size,
                 s=1,
                 alpha=marker_alpha,
                 color=t00)
  
  pos_ax.scatter(
                  agb_rgb_selection['ra'][0:4]-coord_center_M5[0],
                  agb_rgb_selection['dec'][0:4]-coord_center_M5[1],
                 label='',
                 marker = marker,
                 s=marker_size_lost,
                 alpha=marker_alpha_lost,
                 color=t02
                 )
  
  annotate_spacing= [(-35,10),(-10,20),(10,10),(5,10)]
  for i, el in  enumerate(annotate_spacing): 
    # add text to markers with plt.annotate()  
    pos_ax.annotate( str(i+1), \
                      xy=(
                        agb_rgb_selection['ra'][i] -coord_center_M5[0],
                        agb_rgb_selection['dec'][i]-coord_center_M5[1]
                        )
                      ,
                      textcoords='offset pixels', \
                      xytext= annotate_spacing[i],
                      size= 14
                    )
  
  
  for i, el in enumerate(viaux13_crossmatch):
    
    if offset_calibration is True:
      label_offset_cal = str(viaux13_crossmatch_offset[i,0])
    else:
      label_offset_cal = viaux13_data['Nr'][i]
      
    pos_ax.scatter(el['ra']-coord_center_M5[0],
                   el['dec']-coord_center_M5[1],
                   label='',marker = 'o',facecolor = 'none',
                   s=10,alpha=1,
                   color='lightseagreen',
                   zorder = 0.5
                   )
    # add text to markers with plt.annotate()  
    # pos_ax.annotate( label_offset_cal, \
    #                  xy=(el['ra']-coord_center_M5[0],
    #                      el['dec']-coord_center_M5[1]), \
    #                  textcoords='offset pixels', \
    #                  xytext=((0,0,-25)[i],20), \
    #                    size= 14
    #                 )
      
  circle_cluster = plt.Circle(coord_center_M5[:]-coord_center_M5[:], 
                              radius_M5, 
                              color='black', 
                              fill=False)
  pos_ax.add_patch(circle_cluster)
  
  rx = 0.0,0.02,0.04
  for  el in rx:
    circle_cluster = plt.Circle(coord_center_M5-coord_center_M5, 
                              # radius_M5/5, 
                              el, 
                              color='black', 
                              fill=False,
                              linestyle = '-',
                              linewidth = 0.2,
                              zorder = 0.5
                              )
    pos_ax.add_patch(circle_cluster)
  
  pos_ax.axhline(y=coord_center_M5[1]-coord_center_M5[1], linestyle = '-', color = 'k',linewidth=0.2, zorder = 0.5)  
  pos_ax.axvline(x=coord_center_M5[0]-coord_center_M5[0], linestyle = '-', color = 'k',linewidth=0.2, zorder = 0.5)  
  
  pos_ax.set_ylabel(r'dec')
  pos_ax.set_xlabel(r'ra')
  
  r_circle = 3
  pos_ax.set_xlim(coord_center_M5[0]-coord_center_M5[0]-radius_M5/r_circle,
                  coord_center_M5[0]-coord_center_M5[0]+radius_M5/r_circle
                  )
  pos_ax.set_ylim(coord_center_M5[1]-coord_center_M5[1]-radius_M5/r_circle,
                  coord_center_M5[1]-coord_center_M5[1]+radius_M5/r_circle
                  )
  
  pos_ax.set_xticks([-0.05,0.0,0.05])
  
  
  
  # pos_ax.legend()
  # pos_ax.set_aspect('equal', adjustable='box')
  pos_ax.set_aspect('equal')
  
 
  
 
  
  # hist_ax -----------------------------------------------------------
  
  # hist_ax.set_ylabel('N')
  
  
  # hist_ax.set_ylim(0,par_prob_fit.max()+0.001)
  
  
  # hist_ax = hist_ax.twinx()
  label= 'selected RGB AGB stars'
  n,bins,patches = hist_ax.hist(agb_rgb_selection_distance[agb_rgb_selection_distance< d_r_rgb_right_true],
                                # color = t07, 
                                color = t_ARGB, 
                                bins = n_bin, 
                                label=label,
                                alpha=0.
                                )
  n,bins,patches = hist_ax.hist(rgb_selection_distance,
                                color = t_RGB, 
                                bins = bins, 
                                label=label,
                                alpha=1.0
                                )
  
  n,bins,patches = hist_ax.hist(agb_selection_distance,
                                color = t_AGB, 
                                bins = bins, 
                                label=label,
                                alpha=1.
                                )
  
  print('bin range: '+'{0:0.2f} mag with {1:1.0f} bins'.format((bins[1]-bins[0]),n_bin))  
  
  hist_ax.axvline(0,
                        label = 'empirical fit',\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '-' \
                            )
  hist_ax.axvline(-d_r_rgb_left,
                        label = 'empirical fit',\
                        color = 'k',\
                          linewidth=0.6, \
                            zorder = 4,\
                          linestyle = '--' \
                            )
  # hist_ax.axvline(d_r_rgb_right,
  #                       label = 'empirical fit',\
  #                       color = 'k',\
  #                         linewidth=0.6, \
  #                           zorder = 4,\
  #                         linestyle = '--' \
  #                           )  
  # hist_ax.axvline(-d_r_agb,
  #                       label = 'empirical fit',\
  #                       color = 'k',\
  #                         linewidth=0.6, \
  #                           zorder = 4,\
  #                         linestyle = '-.' \
  #                           )  
    
  # delta = 1.5
  # hist_ax.set_xlim(par_x_center-delta,par_x_center+delta)
  hist_ax.set_xlabel(r' $\Delta$ [mag]')
  hist_ax.set_ylabel('N')
  # hist_ax.legend(loc='upper left')
  # hist_ax.legend()
  # hist_ax.legend()
  # hist_ax.set_aspect('equal', adjustable='box')
  
  # cmd_ax-----------------------------------------------------------
  
  
  
  x1, x2, y1, y2 = 1.7, 2.12, -3.52, -3.02
  ratio = np.abs(x1-x2) / np.abs(y1-y2)
  # ratio *= 1/2
  
  height = 0.38
  # width  = 0.60
  width  = height*ratio
  offset = 0.2
  
  # inset axes....
  from mpl_toolkits.axes_grid1.inset_locator import inset_axes
  
  axins = cmd_ax.inset_axes([0.18, 0.6, width, height])
  # axins.set_aspect(.5)
  # subregion of the original image
  axins.set_xlim(x1, x2)
  axins.set_ylim(y1, y2)
  axins.invert_yaxis()
  
  
  # G3 selection from bright (1) to dim (4)
  annotate_spacing= [(-30,-10),(15,-40),(-30,10),(5,-40)]
  for i, el in  enumerate(annotate_spacing): 
    # add text to markers with plt.annotate()  
    axins.annotate( str(i+1), \
                      xy=(agb_rgb_selection['bp_rp'][i],agb_rgb_selection['mbol_mag'][i]), \
                      textcoords='offset pixels', \
                      xytext= annotate_spacing[i],
                      size= 14
                    )
  
  # Viaux13 crossmatch
  # label= 'Viaux13 stars'
  # annotate_spacing= [(10,-10),(-40,-15),(15,-30)]
  # for i, el in  enumerate(viaux13_crossmatch): 
    
  #   if offset_calibration is True:
  #     label_offset_cal = str(viaux13_crossmatch_offset[i,0])
  #   else:
  #     label_offset_cal = viaux13_data['Nr'][i]
          
  #   # add text to markers with plt.annotate()  
  #   axins.annotate( label_offset_cal, \
  #                     xy=(el['bp_rp'],el['mbol_mag']), \
  #                     textcoords='offset pixels', \
  #                     xytext= annotate_spacing[i],
  #                     size= 14
  #                   )
    
  
  for i, ax_el in enumerate([cmd_ax,axins]):
    # ax_el.plot(x,y-brightness_modul,label = 'empirical fit', color = t_fit)
    ax_el.plot(x,y,label = 'empirical fit', color = t_fit, linewidth=0.5, alpha = 1,zorder=1.5)
    
    # d_r = 0.35
    # # d_r_rgb = 0.035
    # d_r_rgb = 0.0656
    # d_r_agb = 0.166
    # -0.023 = x**2 + y**2
    # ax_el.plot(x,y+np.sqrt(d_r**2-x**2),label = 'empirical fit', color = t_fit, linewidth=0.5, alpha = 1,zorder=1.5,linestyle='--')
    # ax_el.plot(x,y+0.077,label = 'empirical fit', color = t_fit, linewidth=0.5, alpha = 1,zorder=1.5,linestyle='--')
    
    # All cluster members
    label = 'cluster members'
    ax_el.scatter(gaia_cluster_members['bp_rp'],
                  gaia_cluster_members['mbol_mag'],
                  label=label,
                  s=marker_size,
                  color=t00
                  )
    
    
    # cluster members above threshold
    # label = '{0:s}: {1:2.1e}'.format(r'cluster members above: ',mbol_rgb_agb_lower_bound)
    # ax_el.scatter(dr3_above_lower_bound_selection['bp_rp'],
    #               dr3_above_lower_bound_selection['mbol_mag'],
    #               label=label,
    #               s=marker_size,
    #               color=t00
    #               )
    
    
    
    # rgb selection derived by fit
    # label= 'selected RGB AGB stars'
    # # ax_el.scatter(agb_rgb_selection['bp_rp'],agb_rgb_selection['mbol_mag'],color=t09,label= label,s=1)   
    # ax_el.scatter(agb_rgb_selection['bp_rp'],agb_rgb_selection['mbol_mag'],
    #               label = label, 
    #               s=marker_size,
    #               alpha=1,
    #               # color=t07
    #               color=t_ARGB
    #               )
    
    label= 'selected RGB stars'
    ax_el.scatter(rgb_selection['bp_rp'],rgb_selection['mbol_mag'],
                  label = label, 
                  s=1,
                  alpha=1,
                  color=t_RGB
                  )
    
    label= 'selected AGB stars'
    ax_el.scatter(agb_selection['bp_rp'],agb_selection['mbol_mag'],
                  label = label, 
                  s=4,
                  marker = 'x',
                  alpha=1,
                  color=t_AGB
                  )
    
    # G3 selection
    ax_el.scatter(agb_rgb_selection['bp_rp'][0:4],
                  agb_rgb_selection['mbol_mag'][0:4],
                  label = 'brightest in the mix', 
                  marker = marker,
                  s = marker_size_lost,
                  alpha = marker_alpha_lost,
                  color='black',
                  zorder = 0.2
                  )
    
    # # rgb agb selection with BC
    # label= 'BC RGB AGB stars'
    # ax_el.scatter(gaia_cluster_members_bc['bp_rp'],
    #               gaia_cluster_members_bc['mbol_mag'],
    #               label = label, 
    #               marker = 'o',
    #               facecolor = 'none',
    #               s=10,
    #               alpha=1,
    #               color='lightseagreen'
    #               )
    
    
    

    
    # Viaux13 crossmatch
    label= 'Viaux13 stars'
    annotate_spacing= [(10,-10),(-40,-15),(15,-30)]
    
    for i, el in  enumerate(viaux13_crossmatch): 
      
      if offset_calibration is True:
        label_offset_cal = str(viaux13_crossmatch_offset[i,0])
      else:
        label_offset_cal = viaux13_data['Nr'][i]
        
      ax_el.scatter(el['bp_rp'],el['mbol_mag'],
                    label = label, marker = 'o',
                    facecolor = 'none',
                    s=10,
                    alpha=1,
                    color='lightseagreen',
                    zorder = 0.1
                    )
      # add text to markers with plt.annotate()  
      # ax_el.annotate( label_offset_cal, \
      #                   xy=(el['bp_rp'],el['mbol_mag']), \
      #                   textcoords='offset pixels', \
      #                   xytext= annotate_spacing[i],
      #                   size= 14
      #                 )
    
    
    # Threshold as line
    ax_el.axhline(mbol_rgb_agb_lower_bound, 
                  linestyle = 'dashed', 
                  label= 'bound: '+str(mbol_rgb_agb_lower_bound)+' mag', 
                  color = 'k',
                  linewidth=0.5, 
                  zorder = 1.5
                  )  

  

  # cmd_ax.set_ylim(y_0-0.1,1.4)  
  cmd_ax.set_ylim(ymax=1.4)
  # cmd_ax.set_ylim(y_0-0.2,-1.5)
  # cmd_ax.set_xlim(0.6,2.2)
  # cmd_ax.set_xlim(-1.0,2.2)
  
  cmd_ax.set_xlabel('BP-RP [mag]')
  cmd_ax.set_ylabel(r'$M_{bol} $[mag]')
  cmd_ax.invert_yaxis()
  # cmd_ax.legend(loc = 'lower center')
  
    
  # cmd_ax[0,i].set_aspect('equal', adjustable='box')
  # cmd_ax.set_aspect('equal')
  cmd_ax.set_aspect(.5)
  
  # fig_sub_path = 'gfx/clus_sel/'
  fig_name = data_set_selected
  
  if savefig is True:
    multi_fig.savefig(fig_sub_path+fig_name+'_agb_rgb_sep.jpg')
  
  
  
  


