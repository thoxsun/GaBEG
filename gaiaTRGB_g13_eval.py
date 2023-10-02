#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 01:34:11 2023

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

folder_name = 'data/trgbEvalData/'

## -------------- GAIA data ---------------------------------------------------

data_set_selected = 'Gaia_M5_r-9arcmin_zpc_pm_clusterMemberSelection_distEval_absG-band_mbol_bcLinAprx_rgb_selection'
file_name = data_set_selected

## ------------ load data into variables --------------------------------------

from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

rgb_selection = file_data
rgb_selection.sort(order = 'mbol_mag') 
label_rgb_selection= 'rgb_selection'
label_rgb_selection_n = 'rgb_selection ('+str(rgb_selection.size)+' stars)'

print()
print(label_rgb_selection_n)


#---------Load TRGB_candidates_selection.csv --------------------------------

trgb_star      = rgb_selection[0].copy()
trgb_star_cl95 = rgb_selection[0].copy()

# Test deepcopy
# trgb_star_cl95['designation']=';)'

set_selected = 'Gaia_M5_r-9arcmin_zpc_pm_clusterMemberSelection_distEval_absG-band_mbol_bcLinAprx_trgb_candidates_selection.csv'

trgb_selection = load_query('data/trgbEvalData/' + set_selected)
label_trgb_selection = 'trgb selection'
label_trgb_selection_n = label_trgb_selection +' ('+str(trgb_selection.size)+' stars)'

print()
print(label_trgb_selection_n)


#---------Store GAIA data --------------------------------

file_sub_path = 'data/trgbEvalData/'
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### ----------- Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/trgb_eval/'


# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)



g13_interpolation_function_plot     = True
g13_eval_plot                       = True
pretty_g13_eval_plot                = True


#%% BCg upper bound derivation based on metallicity deviation for trgb_star


# Metallicity ([Fe/H],[alpha/Fe]) deviation comparison for BCg of trgb_star ---


# Using gdr3bcg to get bc_g
# source: https://gitlab.oca.eu/ordenovic/gaiadr3_bcg
import gdr3bcg.bcg as bcg
table=bcg.BolometryTable()
# Conversion from IAU B2 2015 convention to Straniero20
offset = 0.01

teff_max = float(trgb_star['teff_gspphot'])
logg_max= float(trgb_star['logg_gspphot'])



# [fe/H] ----------------------------------------------------
# Harris 1996 (2010 edition)
# source: https://physics.mcmaster.ca/~harris/mwgc.dat
label_harris10 = 'Harris2010'
feh_harris10 = -1.29

# Carretta2010
# Doi: https://doi.org/10.1051/0004-6361/200913451
label = 'Carretta2010'
feh_carretta = -1.340 

feh_array = np.array([feh_harris10,feh_carretta])

from gaiaErrorEval import trgb_conv_int
[feh_mean, feh_mean_95cl_delta, conf_interval_margin] = trgb_conv_int(feh_array, cl=0.95)

# feh_array = np.array([feh_harris10,feh_mean,feh_carretta])  
feh_array = np.array([conf_interval_margin[0],feh_mean,conf_interval_margin[1]])  

# [alpha/Fe] ----------------------------------------------------
# For information on afe = 0.2 \pm 0.2 see BCg approximation
afe_array = np.array([0.0,0.2,0.4])  

# Calc BC ----------------------------------------------------
bc_g_array = np.full((feh_array.size,afe_array.size),np.nan)
for i, feh in enumerate(feh_array):
  for j, afe in enumerate(afe_array):
  
    stars_bc = [teff_max, logg_max, feh, afe]
    bc_g_array[i,j] = table.computeBc(stars_bc,offset)  

bc_g_array_dev = bc_g_array - bc_g_array[1,1]
print('Largest deviation for all [feh,afe]-combinationsis from mean is {0:0.3f}'.format(np.max(np.abs(bc_g_array_dev))))  

  
# Calculate BC TRGB and unc----------------------------------------------------

# BC 
# feh_mean with 95% CL from Harris2010 & Carretta2010
# afe full range from BCg table 0.0 to 0.4
afe_mean  = 0.2
star_bc   = [teff_max, logg_max, feh_mean, afe_mean]
trgb_bc_g = table.computeBc(stars_bc,offset)

# BC unc with 95cl from feh and full range from afe
# afe full range from BCg table 0.0 to 0.4
afe_unc   = 0.2
stars_bc = [teff_max, logg_max, feh_mean+feh_mean_95cl_delta, afe_mean+ afe_unc]
trgb_bc_g_max = table.computeBc(stars_bc,offset)
trgb_bc_g_unc = trgb_bc_g_max - trgb_bc_g

print()
print('BC with feh from Harris2010 and full range from afe:')
print('BC = ({0:0.3f} \pm {1:0.3f}) mag'.format(trgb_star['bcg_mag'],trgb_star['bcg_mag_error']))
print('Updated BC with mean from literature and unc 95cl from feh and full range from afe:')
print('BC = ({0:0.3f} \pm {1:0.3f}) mag'.format(trgb_bc_g,trgb_bc_g_unc))
print()
print()

# Updated BC is slightly increased in value and unc in third digit after comma
# Note: tested metallicity does not influence BCg significantly
trgb_star['bcg_mag']       = trgb_bc_g
trgb_star['bcg_mag_error'] = trgb_bc_g_unc

trgb_star_cl95['bcg_mag']       = trgb_bc_g
trgb_star_cl95['bcg_mag_error'] = trgb_bc_g_unc

#%% Mbol and Mbol_unc derivation with 95% CL and unc maximization for upper bound derivation

trgb_star_mg_unc = 0.0
# Source: 
# https://www.cosmos.esa.int/web/gaia/science-performance
trgb_star_mg_unc = 0.2

# error, chosen in gaiaExt.py module
distance_modulus_error = trgb_star['distance_modulus_starhorse_error_cl95']

# conservative addition: 
distance_modulus_error = trgb_star['ag_starhorse_mean_error_cl95']    + trgb_star['true_distance_modulus_error']

# squared error propagation
distance_modulus_error = np.sqrt(trgb_star['ag_starhorse_mean_error_cl95']**2 + trgb_star['true_distance_modulus_error']**2 )


# Mg = mg - distance_modulo 
# Note: Starhorse 95% CL as unc 
trgb_star_Mg     = trgb_star['phot_g_mean_mag'] - trgb_star['distance_modulus_starhorse'] 
trgb_star_Mg_unc = np.sqrt(distance_modulus_error**2 +trgb_star_mg_unc**2) 

print('Mg delta {0:0.4f}'.format(trgb_star_Mg - trgb_star['abs_phot_g_mean_mag']))
print('Mg_unc delta {0:0.4f}'.format(trgb_star_Mg_unc - trgb_star['abs_phot_g_mean_mag_error']))

# Mbol = Mg + BC
# Note: bcg_mag is slightly increased due to metallicity mean plus 95%CL
trgb_star_Mbol     = trgb_star_Mg     + trgb_star['bcg_mag']
trgb_star_Mbol_unc = np.sqrt(trgb_star_Mg_unc**2 + trgb_star['bcg_mag_error']**2)

print('Mbol delta {0:0.4f}'.format(trgb_star_Mbol - trgb_star['mbol_mag']))
print('Mbol delta_unc {0:0.4f}'.format(trgb_star_Mbol_unc - trgb_star['mbol_mag_error']))
print()
print()

trgb_star_cl95['mbol_mag']       = trgb_star_Mbol
trgb_star_cl95['mbol_mag_error'] = trgb_star_Mbol_unc


#%% (not needed) preparation of data and bounds ----------------------------------------

# bolometric magnitude Mbol
# see gaiaBCg_approximation.py
# mbol_trgb_lower_bound = -2.4030841435619883
# mbol_trgb_lower_bound     = -2.40
# mask_bc_aprx = rgb_selection['bcg_aprx_flag'] == True
# mbol_trgb_lower_bound_unc = rgb_selection['mbol_mag_error'][mask_bc_aprx][0]

# Manualy chosen to exclude HB in Mbol CMD
# upper_bound_hor_branch  = -0.1 

#%% Mbol TRGB estimation approximation Straniero20

range_Straniero20 = 2.5

limit_Straniero20 = trgb_star['mbol_mag'] +range_Straniero20
mask = rgb_selection['mbol_mag'] < limit_Straniero20
N_trgb_smaller_2_5 = rgb_selection[mask].size

limit_Straniero20_upper_bound = trgb_star['mbol_mag'] -trgb_star['mbol_mag_error'] +range_Straniero20
mask = rgb_selection['mbol_mag'] < limit_Straniero20_upper_bound
N_trgb_smaller_2_5_upper_bound = rgb_selection[mask].size
  
print('Number of RGB stars within 2.5mag from brightest: ')
print(N_trgb_smaller_2_5)
print(r'Number of RGB stars within 2.5mag - Delta Mbol from brightest: ')
print(N_trgb_smaller_2_5_upper_bound)
print()

def trgb_approximation_straniero20(N):
  # Straniero20 TRGB approximation
  # Note N must be inside [50,300]
  x_straniero = np.log10(N)
  
  m_bol_delta_tip = 0.03824 * x_straniero**3 - 0.1627 * x_straniero**2 + 0.09774 * x_straniero + 0.179 
  
  mdm_N = -0.07664 * x_straniero**3 + 0.56869*x_straniero**2 - 1.43798*x_straniero + 1.24946
  
  return [m_bol_delta_tip, mdm_N]


N = np.arange(50,300,1,dtype = float)
# function characteristic
[m_bol_delta_tip_func, mdm_N_func] = trgb_approximation_straniero20(N)
# TRGB brightness increase from Mbol_max with unc
[m_bol_delta_tip, mdm_N] = trgb_approximation_straniero20(N_trgb_smaller_2_5)
# TRGB brightness increase from Mbol_max + Mbol_max unc to calc upper bound
[m_bol_delta_tip_upper_bound, mdm_N_upper_bound] = trgb_approximation_straniero20(N_trgb_smaller_2_5_upper_bound)

print('m_bol_delta_tip - m_bol_delta_tip_upper_bound = {0:0.4f}'.format(m_bol_delta_tip - m_bol_delta_tip_upper_bound))
print('m_bol_delta_tip unc - m_bol_delta_tip_upper_bound unc = {0:0.4f}'.format(mdm_N - mdm_N_upper_bound))
print()
print()


#%% Mbol TRGB estimation approximation Straniero20 plot

fig,ax = plt.subplots(1)
ax.set_title('Straniero20 Approximation')
ax.scatter(N, m_bol_delta_tip_func,s=1,label='delta_tip')
ax.scatter(N, mdm_N_func,s=1,label='mdm')


ax.scatter(N_trgb_smaller_2_5 , m_bol_delta_tip,s=1,label='GC M5 delta_tip')
ax.scatter(N_trgb_smaller_2_5, mdm_N,s=1,label='GC M5 mdm')

ax.scatter(N_trgb_smaller_2_5_upper_bound, m_bol_delta_tip_upper_bound,s=1,label='GC M5 delta_tip upper bound')
ax.scatter(N_trgb_smaller_2_5_upper_bound, mdm_N_upper_bound,s=1,label='GC M5 mdm upper bound')

fig.legend()




#%% Mbol TRGB & uncertainty derivation
# Note:  
# print(trgb_star['mbol_mag'] == trgb_star_cl95['mbol_mag'])
# print(trgb_star['mbol_mag_error'] != trgb_star_cl95['mbol_mag_error'])

# Mbol TRGB 
m_bol_trgb = trgb_star['mbol_mag'] - m_bol_delta_tip
m_bol_trgb_unc = np.sqrt(trgb_star['mbol_mag_error']**2+mdm_N**2)
# Mbol TRGB 95CL unc
m_bol_trgb_unc_cl95 = np.sqrt(trgb_star_cl95['mbol_mag_error']**2+mdm_N**2)

print()
print('Estimated bolometric TRGB brightness for M5:')
print('{0:0.3f} \pm {1:0.3f}'.format(m_bol_trgb,m_bol_trgb_unc))
print()
print('Estimated bolometric TRGB brightness for M5 with CL95 as unc:')
print('{0:0.3f} \pm {1:0.3f}'.format(m_bol_trgb,m_bol_trgb_unc_cl95))
print()
print('m_bol_delta_tip')
print('{0:0.3f} \pm {1:0.3f}'.format(m_bol_delta_tip,mdm_N))


# Mbol TRGB upper bound
# $ (M_{bol})_{TRGB_{80}} $ = -3.784 mag with $ ((M_{bol})_{max_{80}} $ = -3.74 mag and $ < \updelta >_{N=80} $ = 0.039 mag
m_bol_trgb_upper_bound_cl95 = trgb_star['mbol_mag']- m_bol_trgb_unc_cl95 - m_bol_delta_tip_upper_bound 
# $ \Updelta (M_{bol})_{TRGB_{80}} $ = 0.044 mag with $ \Updelta ((M_{bol})_{max_{80}} $ = 0.00 mag and $ MDM_{N=80} $ = 0.044 mag 
m_bol_trgb_upper_bound_cl95_unc = 0.0 + mdm_N_upper_bound


# $ (M_{bol})_{TRGB_{up}} $ = $ (M_{bol})_{TRGB_{80}} $ - $ \Updelta (M_{bol})_{TRGB_{80}} $ = -3.784 -0.044 = -3.825 mag 
m_bol_trgb_upper_bound_cl95_final = m_bol_trgb_upper_bound_cl95 - m_bol_trgb_upper_bound_cl95_unc 
# $ \Updelta (M_{bol})_{TRGB_{up}} $ = $ (M_{bol})_{TRGB_{up}} $ - $ (M_{bol})_{TRGB} $ = -3.825 - (-3.498)  =   0.327 mag 
m_bol_trgb_upper_bound_cl95_final_unc = m_bol_trgb_upper_bound_cl95_final - m_bol_trgb 

print()
print()
print('Estimated upper bound bolometric TRGB brightness for M5 with unc CL95:')
print('{0:0.3f}'.format(m_bol_trgb_upper_bound_cl95_final))
print()
print('Estimated bolometric TRGB brightness for M5 with upper bound as unc:')
print('{0:0.3f} \pm {1:0.3f}'.format(m_bol_trgb,m_bol_trgb_upper_bound_cl95_final_unc))
print()
print()

# trgb_star['mbol_mag'] = m_bol_trgb
# trgb_star['mbol_mag_error'] = m_bol_trgb_upper_bound_cl95_final_unc

# m_bol_trgb     = trgb_star['mbol_mag'] 
# m_bol_trgb_unc = trgb_star['mbol_mag_error']
# m_bol_trgb_unc_cl95 = trgb_star['mbol_mag_error']

m_bol_trgb_upper_bound_cl95 = m_bol_trgb_upper_bound_cl95_final
m_bol_trgb_upper_bound_cl95_unc = m_bol_trgb_upper_bound_cl95_final_unc

# # Mbol TRGB upper bound - faulty, delete in future
# m_bol_trgb_upper_bound = trgb_star['mbol_mag'] - m_bol_delta_tip_upper_bound

# m_bol_trgb_upper_bound_unc_cl95 = np.sqrt(trgb_star_cl95['mbol_mag_error']**2+mdm_N_upper_bound**2)
# m_bol_trgb_upper_bound_cl95 = trgb_star['mbol_mag'] - m_bol_delta_tip_upper_bound - m_bol_trgb_upper_bound_unc_cl95

# m_bol_trgb_upper_bound_cl95 = trgb_star['mbol_mag'] - m_bol_delta_tip_upper_bound - m_bol_trgb_upper_bound_unc_cl95 

# m_bol_trgb_upper_bound_cl95_unc = m_bol_trgb - m_bol_trgb_upper_bound_cl95

# print()
# print()
# print('Estimated upper bound bolometric TRGB brightness for M5 with unc CL95:')
# print('{0:0.3f}'.format(m_bol_trgb_upper_bound_cl95))
# print()
# print('m_bol_delta_tip:')
# print('{0:0.3f}'.format(m_bol_delta_tip_upper_bound + m_bol_trgb_upper_bound_unc_cl95))
# print()
# print()






# Mbol TRGB & Mbol_TRGB upper bound dependencies:
# # Mbol TRGB and dependencies
# m_bol_trgb:         Mbol_TRGB_M5
# m_bol_trgb_unc:     Mbol_TRGB_M5 unc 1Sigma Starhorse 
# m_bol_trgb_unc_cl95: 
#             - Mbol TRGB M5 unc 95CL Starhorse 
#             - Updated BC with: 
#                 -mean feh \pm 95cl unc from Harris2010 & Carretta2010
#                 -afe full range 
# # Mbol TRGB upper bound and dependencies:
# m_bol_trgb_upper_bound_cl95 depends on:
#             - Mbol_TRGB_M5  
#             - m_bol_delta_tip_upper_bound (dependent on N RGB stars in RGB selection)
#                   -N_trgb_smaller_2_5_upper_bound
#             - m_bol_trgb_upper_bound_unc_cl95 (dependent on N RGB stars in RGB selection)
#                   -m_bol_trgb_unc_cl95 
#                   -mdm_N_upper_bound



#%% Straniero20 results for M5

# VI + JK with ZAHB distance scale
m_bol_trgb_Str20 = np.mean([-3.60,-3.63])
m_bol_trgb_Str20_unc = np.mean([0.25,0.26])




#%% Start of g13eval:
  # Note: Move to own file?

#%% Choose TRGB Mbol

# Note: necessary, if moved to new file
# m_bol_trgb     = trgb_star['mbol_mag'] 
# m_bol_trgb_unc = trgb_star['mbol_mag_error']


#%% Metallicity deviation comparison [M/H]

# [M/H] evaluation from Straniero20--------------------------
def mh_eval_Straniero20(feh, afe):
  
  # mh = feh + np.log10(b+a*10**(afe))
  a = 0.6695 
  b = 0.3305
  
  mh = feh + np.log10(b+a*10**(afe))

  return mh






# [fe/H] ----------------------------------------------------
# Harris 1996 (2010 edition)
# source: https://physics.mcmaster.ca/~harris/mwgc.dat
label_harris10 = 'Harris2010'
feh_harris10 = -1.29

# Carretta2010
# Doi: https://doi.org/10.1051/0004-6361/200913451
label = 'Carretta2010'
feh_carretta = -1.340 

feh_array = np.array([feh_harris10,feh_carretta])

from gaiaErrorEval import trgb_conv_int
[feh_mean, feh_mean_95cl_delta, conf_interval_margin] = trgb_conv_int(feh_array, cl=0.95)

# feh_mean with 95% CL from Harris2010 & Carretta2010
# feh_mean
# feh_mean_95cl_delta
# feh_array = np.array([feh_harris10,feh_mean,feh_carretta])  
feh_array = np.array([conf_interval_margin[0],feh_mean,conf_interval_margin[1]])  

# [alpha/Fe] ----------------------------------------------------
# For information on afe = 0.2 \pm 0.2 see BCg approximation
afe_array = np.array([0.0,0.2,0.4])  
# afe full range from BCg table 0.0 to 0.4
afe_mean  = 0.2
afe_unc   = 0.2
afe_upper_limit = afe_mean+ afe_unc



# Calc [M/H] ----------------------------------------------------
mh_array = np.zeros((feh_array.size,afe_array.size))

for i, feh in enumerate(feh_array):
  
  print()
  print('-------------------------')
  print('feh: {0:0.3f}'.format(feh))
  
  for j, afe in enumerate(afe_array):
      
    mh_array[i,j] = mh_eval_Straniero20(feh, afe)
    print() 
    print('afe: {0:0.3f}'.format(afe))
    print(' mh: {0:0.3f}'.format(mh_array[i,j]))

print('-------------------------')
print('-------------------------')
mh_array_dev = mh_array - mh_array[1,1]
print('mh range:')
print(mh_array)
print()
print('deviation from mh center in mh range:')
print(mh_array_dev)


print('-------------------------')
print('-------------------------')

symmetry_testing = True
if symmetry_testing == True:
  
  mh_feh_array_afe_unc = (mh_array[0,:] - mh_array[2,:])/2
  mh_afe_array_feh_unc = (mh_array[:,0] - mh_array[:,2])/2
  
  print('mh unc with fixed feh symmetric: ')
  print('row1: max-center; row2: (max-min)/2')
  print(str(mh_array[1,:]-mh_array[2,:]))
  print(str(mh_feh_array_afe_unc))
  print()
  # print('mh unc with fixed afe symmetric: {0:0.3f} == {1:0.3f}'.format(mh_array[:,1]-mh_array[:,2],mh_afe_array_feh_unc))
  print('mh unc with fixed afe symmetric: ')
  print('row1: max-center; row2: (max-min)/2')
  print(str(mh_array[:,1]-mh_array[:,2]))
  print(str(mh_afe_array_feh_unc))

print('-------------------------')
print('-------------------------')



#%% [M/H] estimation with [Fe/H] & [alpha/Fe] using Straniero20

# mh_G3 = mh_G3_carretta
# mh_G3 = mh_G3_afe04

afe_G3 = afe_mean
afe_G3_unc = afe_unc

feh_G3 = feh_mean
feh_G3_unc = feh_mean_95cl_delta

mh_G3             =  mh_array[1,1]
mh_G3_margin      = [mh_array[0,0],mh_array[2,2]]
mh_G3_margin_dev  = [mh_array_dev[0,0],mh_array_dev[2,2]]
mh_G3_unc         =  np.abs(mh_G3 - mh_G3_margin).max()


print(' [M/H] derived with Straniero2020 Eq.2 in dex: ')
print('mh_G3: '+str(mh_G3))
print('mh_G3_margin_dev: '+str(mh_G3_margin_dev))
print('mh_G3_unc: '+str(mh_G3_unc))
print('-------------------------')
print('-------------------------')

#%% [M/H] estimation with [Fe/H] & [alpha/Fe] using Straniero20 plot

feh = feh_mean
feh_unc = feh_mean_95cl_delta

afe = afe_mean
afe_unc = afe_unc

mh = mh_eval_Straniero20(feh, afe)


afe_func = np.linspace(afe-afe_unc,afe+afe_unc,100)
feh_func = np.linspace(feh-feh_unc,feh+feh_unc,100)

mh_func_feh_fix = mh_eval_Straniero20(feh, afe_func)
mh_func_afe_fix = mh_eval_Straniero20(feh_func, afe)

mh_func_feh_fix_dev = mh_eval_Straniero20(feh, afe_func) - mh
mh_func_afe_fix_dev = mh_eval_Straniero20(feh_func, afe) - mh

afe_func_dev = afe_func - afe
afe_range = 2*afe_unc
afe_scaled = afe_func_dev / afe_range

feh_func_dev = feh_func - feh
feh_range = 2*feh_unc
feh_scaled = feh_func_dev / feh_range


fig,axs = plt.subplots(2,1)
ax = axs[0]
ax_1 = axs[1]

ax.scatter(afe_scaled, mh_func_feh_fix,s=1,label='mh(afe) with fixed feh:'+str(feh))
ax.scatter(feh_scaled, mh_func_afe_fix,s=1,label='mh(feh) with fixed afe:'+str(afe))

ax_1.scatter(afe_scaled, mh_func_feh_fix_dev,s=1,label='fixed feh:'+str(feh))
ax_1.scatter(feh_scaled, mh_func_afe_fix_dev,s=1,label='fixed afe:'+str(afe))


ax.set_ylabel('[M/H]')  
ax.set_xlabel('scale')
# ax.invert_yaxis()
ax_1.set_ylabel('[M/H]')  
ax_1.set_xlabel('scale')
ax.legend()



#%% Prepare Straniero20 data

# Load reference table from Straniero20 
# Note: add source!
file_sub_path = 'data/trgbEvalData/'
file_name = 'Straniero20Table1.csv'

file_data = load_query(file_sub_path+file_name)

# Sorted by metallicity 
file_data.sort(order='MH')
j = 9
i = 0
z = []

while i < file_data.size:
# for i, el in enumerate(file_data):
  dummy = np.copy(file_data[i:i+j])
  i += j
  z.append(dummy)
  
# Sorted by metallicity   
z_1 = z[0]
z_2 = z[1]
z_3 = z[2]


# Sorted by g_13
file_data.sort(order='g_13')
j = 3
i = 0
g_13 = []
  
while i < file_data.size:
# for i, el in enumerate(file_data):
  dummy = np.copy(file_data[i:i+j])
  i += j
  g_13.append(dummy)

#%% Mbol([M/H]) interpolation function for all 9 g13 values provided by Straniero20

#  Note: For each g13 value, one polynomial fit Mbol([M/H]) of degree 2 is derived
#         mh_mbol_poly2_fit_list includes all 9 fits for 9 g13 values
#         adjoined g13 values are stored in g_13_poly_mh_G3

# [M/H]    array from Straniero20
# Mbol     array from Straniero20
# g13      array from Straniero20


def mh_mbol_poly2_fit(g_13):

  # Calc g13([M/H]=[M/H]_M5,Mbol) poly fit
  g_13_poly_mh_G3 = np.zeros(len(g_13))
  mh_mbol_poly2_fit_list = []
  
  # for each g13 value, one polynomial fit of degree 2 is derived for Mbol([M/H])
  for i, el in enumerate(g_13):
    
    mh_poly = el['MH'][:]
    mbol_poly = el['Mbol'][:]
    
    # polynomial fit with order = 2 for sets with g13 = const
    poly_mh_mbol = np.poly1d(np.polyfit(mh_poly, mbol_poly, 2))
    
    mh_mbol_poly2_fit_list.append(poly_mh_mbol)
    g_13_poly_mh_G3[i] = el['g_13'][0]

  
  return [g_13_poly_mh_G3,mh_mbol_poly2_fit_list]

[g_13_poly_mh_G3,mh_mbol_poly2_fit_list] = mh_mbol_poly2_fit(g_13)

#%% Mbol([M/H]) interpolation of GC M5 with [M/H]_M5

# [M/H]_M5 scalar derived from GC M5

# mh_G3             =  mh_array[1,1]
# mh_G3_margin      = [mh_array[0,0],mh_array[2,2]]
# mh_G3_margin_dev  = [mh_array_dev[0,0],mh_array_dev[2,2]]
# mh_G3_unc         =  np.abs(mh_G3 - mh_G3_margin).max()

# g_13_poly_mh_G3 = Includes adjoined g13 values

mbol_poly_mh_G3 = np.zeros(len(g_13))
mbol_poly_mh_G3_margin_max = np.zeros(len(g_13))
mbol_poly_mh_G3_margin_min = np.zeros(len(g_13))

mh_G3_margin_min = np.min(mh_G3_margin)
mh_G3_margin_max = np.max(mh_G3_margin)

for i, mh_mbol_poly2 in enumerate(mh_mbol_poly2_fit_list):
  mbol_poly_mh_G3[i]            = mh_mbol_poly2(mh_G3)
  mbol_poly_mh_G3_margin_max[i] = mh_mbol_poly2(mh_G3_margin_max)
  mbol_poly_mh_G3_margin_min[i] = mh_mbol_poly2(mh_G3_margin_min)

# Uncertainty estimates of poly2 fits:
if True:
  # mbol_poly_mh_G3_margin_max
  # mbol_poly_mh_G3_margin_min
  # mbol_poly_mh_G3
  
  min_dev = mbol_poly_mh_G3 - mbol_poly_mh_G3_margin_min
  max_dev = mbol_poly_mh_G3 - mbol_poly_mh_G3_margin_max 
  
  min_dev_mean = min_dev.mean()
  min_dev_std = min_dev.std()
  
  max_dev_mean = max_dev.mean()
  max_dev_std = max_dev.std()
  
  mbol_poly_mh_G3_unc = np.max(np.array(max_dev_mean,min_dev_mean))


#%% Mbol([M/H]_M5) metallicity comparison plot

   

if g13_interpolation_function_plot is True:
  
  mbol_fig, mbol_ax = plt.subplots(1,1,
                                   figsize = (fig_width,fig_height*0.4),
                                   sharex=(True),
                                   sharey=False)
  
  if ready_to_publish is False:
    mbol_fig.suptitle('Comparison of TRGB M_bol to Straniero20')
          
  
  mbol_ax.set_ylabel(r'$ M_{bol} [mag] $')  
  mbol_ax.set_xlabel(r'$ [\frac{M}{H}] $')
  mbol_ax.invert_yaxis()
  
  
  # Mbol([M/H]_M5)
  mh_G3_scatter = np.full(g_13_poly_mh_G3.size,mh_G3)
  label_z = "[M/H] %0.3f  dex" % (mh_G3)
  mbol_ax.scatter(mh_G3_scatter,mbol_poly_mh_G3,
                label = label_z,
                color=t_fit,linestyle ='-',marker='x')
  # Mbol([M/H]_M5_max) 
  mh_G3_scatter = np.full(g_13_poly_mh_G3.size,mh_G3_margin_max)
  label_z = "[M/H]_max %0.3f  dex" % (mh_G3_margin_max)
  mbol_ax.scatter(mh_G3_scatter,mbol_poly_mh_G3_margin_max,
                label = label_z,
                color=t_fit,linestyle ='-',marker='x')
  # Mbol([M/H]_M5_min)
  mh_G3_scatter = np.full(g_13_poly_mh_G3.size,mh_G3_margin_min)
  label_z = "[M/H]_min %0.3f  dex" % (mh_G3_margin_min)
  mbol_ax.scatter(mh_G3_scatter,mbol_poly_mh_G3_margin_min,
                label = label_z,
                color=t_fit,linestyle ='-',marker='x')
    
    
  for i, el in enumerate(z_1):
    label_z = "g_13 %0.3f  dex" % (z_3['g_13'][0])
    mbol_ax.plot([z_3['MH'][i],z_2['MH'][i],z_1['MH'][i]], [z_3['Mbol'][i],z_2['Mbol'][i],z_1['Mbol'][i]],
                 # label = label_z,
                 color=t_fit,linestyle ='',marker='^')
    
    z_i = np.array([z_3['MH'][i],z_2['MH'][i],z_1['MH'][i]])
    mh_i = np.linspace(z_i.min(), z_i.max(),100)
    
    mbol_i = mh_mbol_poly2_fit_list[i](mh_i)
    
    mbol_ax.plot(mh_i,mbol_i,
                 # label = label_z,
                 color=t_fit,linestyle ='-',marker='')
    
  
  mbol_ax.legend()
  
   
    
  if savefig is True:
    mbol_fig.savefig(fig_sub_path+'M5_mbol_mh_trgb_poly2_fit_based_on_Straniero20')



#%% testing poly fit func g13(Mbol_TRGB_M5,[M/H]_M5) with poly2 fit Mbol(g13, [M/H]_M5 ) as input

# Note: Not in use because
#       Poly fit with orders 2 to 9 and higher do not resemble data set! 
poly_fit_true = False

if poly_fit_true is True:
  # Calc Mbol_MH_G3(g13) poly fit
  poly_g13_mbol = np.poly1d(np.polyfit(mbol_poly_mh_G3,g_13_poly_mh_G3 , 3))
      
  # Calc g13 for (Mbol)TRGB and upper bound
  g_13_poly_mbol_mh_G3_TRGB_M5       = poly_g13_mbol(m_bol_trgb)
  g_13_poly_mbol_mh_G3_TRGB_M5_upper = poly_g13_mbol(m_bol_trgb - m_bol_trgb_unc)
  g_13_poly_mbol_mh_G3_TRGB_M5_lower = poly_g13_mbol(m_bol_trgb + m_bol_trgb_unc)
  
  mbol_poly_g13_mbol_mh_G3 = np.linspace(-5.0,-3.0,100)
  mbol_poly_g13_mbol_mh_G3_limited_edition = g_13_poly_mh_G3
  g_13_poly_mbol_mh_G3       = poly_g13_mbol(mbol_poly_g13_mbol_mh_G3 )





#%% g13(Mbol_TRGB,[M/H]_const) linear interpolation function of poly2 fit Mbol(g13, [M/H]_const )

# g13(Mbol_TRGB,[M/H]_const) linear interpolation 
def g13_mbol_lin_interpol(m_bol_fit,g_13_poly_mh_G3,mbol_poly_mh_G3,label=''):
  # Calc g13 for (Mbol)TRGB 
  #
  # m_bol_fit is (Mbol)TRGB with metallicity mh from observed GC ------
  # Note: mh must be used in Mbol([M/H]) poly fit to derive 
  #        g_13_poly_mh_G3 & mbol_poly_mh_G3
  #
  # g_13_poly_mh_G3,mbol_poly_mh_G3 ---------------------------
  # From Mbol([M/H]) poly fit with metallicity mh from observed GC 
  
  
  mask_brighter =  (m_bol_fit <= mbol_poly_mh_G3)
  mask_dimmer = (m_bol_fit >= mbol_poly_mh_G3)
  
  mbol_poly_larger  = mbol_poly_mh_G3[mask_brighter]
  mbol_poly_smaller = mbol_poly_mh_G3[mask_dimmer]
  if np.any(mbol_poly_larger):
    mbol_poly_larger  = mbol_poly_larger[-1]
    g_13_poly_mh_G3_larger  = (g_13_poly_mh_G3[mask_brighter]  )[-1]
  
    if np.any(mbol_poly_smaller):
      mbol_poly_smaller = mbol_poly_smaller[0]
      g_13_poly_mh_G3_smaller = (g_13_poly_mh_G3[mask_dimmer])[0]
      
      k = (mbol_poly_smaller - mbol_poly_larger ) / (g_13_poly_mh_G3_smaller - g_13_poly_mh_G3_larger) 
      
      # m_bol_fit = g_13_lin_interpol_trgb * k + mbol_poly_larger
      g_13_lin_interpol_trgb = (m_bol_fit-mbol_poly_smaller) / k + g_13_poly_mh_G3_smaller
      # g_13_lin_interpol_trgb = (m_bol_fit-mbol_poly_smaller) / k + g_13_poly_mh_G3_smaller
      print('-------------------------------------------------------')
      print()
      print('Evaluate g13 for: '+label)
      print('(Mbol)TRGB to fit: '+str(m_bol_fit))
      print('g13 upper bound: '+str(g_13_lin_interpol_trgb))
      print()
      print('-------------------------------------------------------')
    else:
      print('-------------------------------------------------------')
      print('Evaluate g13 for: '+label)
      print('(Mbol)TRGB to fit: '+str(m_bol_fit))
      print('Mbol so high, that it leads to g13 > 6!')
      print('No data for approximation available!')
      print('TRGB g13 set to = 6 with: '+str(mbol_poly_mh_G3.min()))
      g_13_lin_interpol_trgb = 6.0
      print()
      print('-------------------------------------------------------')
  else:
    print('-------------------------------------------------------')
    print()
    print('Evaluate g13 for: '+label)
    print('(Mbol)TRGB to fit: '+str(m_bol_fit))
    print('Mbol so low, that g13 < 0!')
    print('g13 cannot be negative per definition!')
    print('TRGB g13 set to 0 with: '+str(mbol_poly_mh_G3.max()))
    g_13_lin_interpol_trgb = 0.0
    print()
    print('-------------------------------------------------------')
  return g_13_lin_interpol_trgb



# -----------------------------------------------------------------------------
# Testing function:
  
# m_bol_fit = -3.9613
m_bol_fit = -3.62-0.26
m_bol_fit = m_bol_trgb
m_bol_fit = m_bol_trgb + m_bol_trgb_unc

g_13_lin_interpol_trgb = g13_mbol_lin_interpol(m_bol_fit,g_13_poly_mh_G3,mbol_poly_mh_G3,label ='Test')  





#%% g13(Mbol_TRGB_M5,[M/H]_M5) linear interpolation function of poly2 fit Mbol(g13, [M/H]_M5 ) Calculation

# -----------------------------------------------------------------------------  
# Calculation of g13(Mbol_TRGB_M5,[M/H]_M5) linear interpolation 

# Note: To get upper g13 bound, mh_margin_min needs to be used (see plots)
# Note: To get lower g13 bound, mh_margin_max needs to be used (see plots)


# g13(Mbol_TRGB_M5,[M/H]_M5) linear interpolation 
g_13_M5_TRGB_G3 = g13_mbol_lin_interpol(
                                        m_bol_trgb,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3,
                                        label='GC M5'
                                        )  

#-------------------------------------------------------------------------------------------------
# g13(Mbol_TRGB_M5,[M/H]_M5) linear interpolation with m_bol_trgb_unc and mbol_poly_mh_G3_unc
g_13_M5_TRGB_G3_limit_dim = g13_mbol_lin_interpol(
                                        m_bol_trgb+m_bol_trgb_unc,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_max,
                                        label='GC M5 dimmer limit'
                                        )

g_13_M5_TRGB_G3_limit_bright = g13_mbol_lin_interpol(
                                        m_bol_trgb-m_bol_trgb_unc,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min,
                                        label='GC M5 brighter limit'
                                        )
g_13_M5_TRGB_G3_limit_bright_upper = g13_mbol_lin_interpol(
                                        m_bol_trgb-m_bol_trgb_unc,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min+mbol_poly_mh_G3_unc,
                                        label='GC M5 brighter limit with upper mbol poly fit unc based on mh_G3_unc'
                                        )
g_13_M5_TRGB_G3_limit_bright_lower = g13_mbol_lin_interpol(
                                        m_bol_trgb-m_bol_trgb_unc,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min-mbol_poly_mh_G3_unc,
                                        label='GC M5 brighter limit lower mbol poly fit unc based on mh_G3_unc'
                                        )


#-------------------------------------------------------------------------------------------------
# g13(Mbol_TRGB_M5,[M/H]_M5) linear interpolation with m_bol_trgb_unc_cl95 and mbol_poly_mh_G3_unc

g_13_M5_TRGB_G3_cl95_limit_dim = g13_mbol_lin_interpol(
                                        m_bol_trgb+m_bol_trgb_unc_cl95,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_max,
                                        label='GC M5 dimmer limit'
                                        )

g_13_M5_TRGB_G3_cl95_limit_bright = g13_mbol_lin_interpol(
                                        m_bol_trgb-m_bol_trgb_unc_cl95,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min,
                                        label='GC M5 cl95 brighter limit'
                                        )
g_13_M5_TRGB_G3_cl95_limit_bright_upper = g13_mbol_lin_interpol(
                                        m_bol_trgb-m_bol_trgb_unc_cl95,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min+mbol_poly_mh_G3_unc,
                                        label='GC M5 cl95 brighter limit with upper mbol poly fit unc based on mh_G3_unc'
                                        )
g_13_M5_TRGB_G3_cl95_limit_bright_lower = g13_mbol_lin_interpol(
                                        m_bol_trgb-m_bol_trgb_unc_cl95,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min-mbol_poly_mh_G3_unc,
                                        label='GC M5 cl95 brighter limit with lower mbol poly fit unc based on mh_G3_unc'
                                        )


# Uncertainty estimates of linear interpolation:
if True:

  print('For M5 TRGB G3:')
  print(g_13_M5_TRGB_G3)
  
  print()
  print('g13 for M5 TRGB G3 dimmer limit:')
  print(g_13_M5_TRGB_G3)
  print(g_13_M5_TRGB_G3_cl95_limit_dim)
  
  print()
  print('g13 for M5 TRGB G3 brighter limit:')
  print(g_13_M5_TRGB_G3)
  print(g_13_M5_TRGB_G3_cl95_limit_bright)
  print()
  print('..with mbol poly fit unc based on mh G3 unc =')
  print(mbol_poly_mh_G3_unc)
  print('yielding brighter lower limit:')
  print(g_13_M5_TRGB_G3_cl95_limit_bright_lower)
  print('...and yielding brighter upper limit:')
  print(g_13_M5_TRGB_G3_cl95_limit_bright_upper)
  
  
  min_dev = g_13_M5_TRGB_G3_cl95_limit_bright_lower -g_13_M5_TRGB_G3_cl95_limit_bright
  max_dev = g_13_M5_TRGB_G3_cl95_limit_bright_upper -g_13_M5_TRGB_G3_cl95_limit_bright  
    
  g_13_M5_TRGB_G3_cl95_limit_bright_unc = np.max(np.abs(np.array(max_dev,min_dev)))
  
  print()
  print('g13 min for g13 brighter limit:')
  print(min_dev)
  print()
  print('g13 max for g13 brighter limit:')
  print(max_dev)
  print()
  


#-------------------------------------------------------------------------------------------------
# g13(Mbol_TRGB_M5 upper bound,[M/H]_M5) linear interpolation with mbol_poly_mh_G3_unc and mbol_poly_mh_G3_unc
g_13_M5_TRGB_G3_upper_bound = g13_mbol_lin_interpol(
                                        m_bol_trgb_upper_bound_cl95,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min,
                                        label='GC M5 upper bound'
                                        )  

g_13_M5_TRGB_G3_upper_bound_upper = g13_mbol_lin_interpol(
                                        m_bol_trgb_upper_bound_cl95,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min+mbol_poly_mh_G3_unc,
                                        label='GC M5 upper bound with upper mbol poly fit unc based on mh_G3_unc'
                                        )  
g_13_M5_TRGB_G3_upper_bound_lower = g13_mbol_lin_interpol(
                                        m_bol_trgb_upper_bound_cl95,
                                        g_13_poly_mh_G3,
                                        mbol_poly_mh_G3_margin_min-mbol_poly_mh_G3_unc,
                                        label='GC M5 upper bound with lower mbol poly fit unc based on mh_G3_unc'
                                        )  

# Uncertainty estimates of linear interpolation:
if True:
  print('For M5 upper bound:')
  print(m_bol_trgb_upper_bound_cl95)
  
  print()
  print('g13 for M5 upper bound:')
  print(g_13_M5_TRGB_G3_upper_bound)
  
  print()
  print('..with mbol poly fit unc based on mh G3 unc =')
  print(mbol_poly_mh_G3_unc)
  print('yielding lower limit:')
  print(g_13_M5_TRGB_G3_upper_bound_lower)
  print('...and yielding brighter upper limit:')
  print(g_13_M5_TRGB_G3_upper_bound_upper)

  g_13_M5_TRGB_G3_cl95_limit_bright_max_dev = g_13_M5_TRGB_G3_cl95_limit_bright-g_13_M5_TRGB_G3_cl95_limit_bright_lower
  g_13_M5_TRGB_G3_cl95_limit_bright_min_dev = g_13_M5_TRGB_G3_cl95_limit_bright  - g_13_M5_TRGB_G3_cl95_limit_bright_upper
    
  g_13_M5_TRGB_G3_cl95_limit_bright_unc = np.max(np.abs(np.array(g_13_M5_TRGB_G3_cl95_limit_bright_max_dev,g_13_M5_TRGB_G3_cl95_limit_bright_min_dev)))
  
  print()
  print('g13 min for g13 brighter limit:')
  print(min_dev)
  print()
  print('g13 max for g13 brighter limit:')
  print(max_dev)
  print()


  



# Note: All lower bounds yield negative g13 
#       g13_lower := 0



#%% Definitions and Dependencies of input and results for g13 eval

# --------------------- Definitions and dependencies --------------------------
# [M/H]_M5:
# mh_G3             =  mh_array[1,1]
# mh_G3_margin      = [mh_array[0,0],mh_array[2,2]]
# mh_G3_margin_dev  = [mh_array_dev[0,0],mh_array_dev[2,2]]
# mh_G3_unc         =  np.abs(mh_G3 - mh_G3_margin).max()
# mh_G3_margin_min = np.min(mh_G3_margin)
# mh_G3_margin_max = np.max(mh_G3_margin)

# Mbol(g13, [M/H]_M5 ) poly2 fit:
# g_13_poly_mh_G3             = Includes adjoined g13 values
# mbol_poly_mh_G3             = Includes adjoined Mbol values 
# mbol_poly_mh_G3_margin_max  = Includes adjoined Mbol values
# mbol_poly_mh_G3_margin_min  = Includes adjoined Mbol values




# Mbol TRGB & Mbol_TRGB upper bound dependencies:
# # Mbol TRGB and dependencies
# m_bol_trgb:         Mbol_TRGB_M5
# m_bol_trgb_unc:     Mbol_TRGB_M5 unc 1Sigma Starhorse 
# m_bol_trgb_unc_cl95: 
#             - Mbol TRGB M5 unc 95CL Starhorse 
#             - Updated BC with: 
#                 -mean feh \pm 95cl unc from Harris2010 & Carretta2010
#                 -afe full range 
# # Mbol TRGB upper bound and dependencies:
# m_bol_trgb_upper_bound_cl95:
#             - Mbol_TRGB_M5  
#             - m_bol_delta_tip_upper_bound (dependent on N RGB stars in RGB selection)
#                   -N_trgb_smaller_2_5_upper_bound
#             - m_bol_trgb_upper_bound_unc_cl95 (dependent on N RGB stars in RGB selection)
#                   -m_bol_trgb_unc_cl95 
#                   -mdm_N_upper_bound


# --------------------- Results --------------------------------------

# Input: mbol_poly_mh_G3
# g_13_M5_TRGB_G3

# Input: mbol_poly_mh_G3_margin_max
# g_13_M5_TRGB_G3_limit_dim
# g_13_M5_TRGB_G3_cl95_limit_dim


# Input: mbol_poly_mh_G3_margin_min
#   g_13_M5_TRGB_G3_limit_bright
#   g_13_M5_TRGB_G3_cl95_limit_bright

#   g_13_M5_TRGB_G3_upper_bound






#%% g13 interpolate Straniero2020 comparison plot

if g13_eval_plot is True:
  
  mbol_fig, mbol_ax = plt.subplots(1,1,
                                   figsize = (fig_width,fig_height*0.4),
                                   sharex=(True),
                                   sharey=False)
  
  if ready_to_publish is False:
    mbol_fig.suptitle('Comparison of TRGB M_bol to Straniero20')
          
  
  mbol_ax.set_ylabel('$ M_{bol} [mag] $')  
  mbol_ax.set_xlabel('$ g_{13} $')
  mbol_ax.invert_yaxis()
  
  # mbol_ax.set_xticks()
  # mbol_ax.set_xlim()  
  # mbol_ax.set_ylim(top=15,bottom=11)  
  
  label_z = "[M/H] %0.3f  dex" % (mh_G3)
  
  if poly_fit_true is True:
    mbol_ax.plot(g_13_poly_mbol_mh_G3,mbol_poly_g13_mbol_mh_G3,
                  label = label_z,
                  color='green',linestyle ='-',marker='')
  
  label_z = "[M/H] %0.3f  dex" % (z_3['MH'][0])
  mbol_ax.plot(g_13_poly_mh_G3,mbol_poly_mh_G3,
                # label = label_z,
               color='darkgreen',linestyle ='-',marker='*')
  
  
  mbol_ax.scatter(g_13_lin_interpol_trgb,m_bol_fit,
                  label = ' ', 
                  color = 'darkgreen',
                  zorder=4, 
                  marker='+',
                  s = 50)
  
  label_z = "[M/H] %0.3f  dex" % (z_3['MH'][0])
  mbol_ax.plot(z_3['g_13'],z_3['Mbol'],
                label = label_z,
               color=t_fit,linestyle ='-',marker='^')
  
  label_z = "[M/H] %0.3f  dex" % (z_2['MH'][0])
  mbol_ax.plot(z_2['g_13'],z_2['Mbol'],
                label = label_z,
               color=t_fit,linestyle ='-',marker='x')
  
  label_z = "[M/H] %0.3f  dex" % (z_1['MH'][0])
  mbol_ax.plot(z_1['g_13'],z_1['Mbol'],
                label = label_z,
               color=t_fit,linestyle ='-',marker='o')
  
  label_m_bol_trgb = r"$[\frac{M}{H}]$ = (%0.2f $\pm$ %0.2f), $M_{bol}$=(%0.2f $\pm$ %0.2f) mag" % (-0.99, 0.20, m_bol_trgb, m_bol_trgb_unc)
  
  
  mbol_ax.axhline(m_bol_trgb+m_bol_trgb_unc, 
                  color=t_fit,linewidth=0.6,linestyle='--')
  mbol_ax.axhline(m_bol_trgb-m_bol_trgb_unc, 
                  color=t_fit,linewidth=0.6,linestyle='--')
  mbol_ax.axhline(m_bol_trgb, label = label_m_bol_trgb, 
                  color=t_fit,linewidth=0.6,linestyle='-')
  
  mbol_ax.axhspan(m_bol_trgb-m_bol_trgb_unc, m_bol_trgb+m_bol_trgb_unc, color=t00,alpha=0.4)  #horizontal shading
  
  # VI + JK with ZAHB distance scale
  m_bol_trgb_Str20 = np.mean([-3.60,-3.63])
  m_bol_trgb_Str20_unc = np.mean([0.25,0.26])
  
  label_m_bol_trgb_Str20 = r"$[\frac{M}{H}]$ = (%0.2f $\pm$ %0.2f), $M_{bol}$=(%0.2f $\pm$ %0.2f) mag" % (-1.11 ,0.20 , m_bol_trgb_Str20, m_bol_trgb_Str20_unc)
  # mbol_ax.axhline(m_bol_trgb_Str20, label = label_m_bol_trgb_Str20, 
  #                 color=t_fit,linewidth=0.6,linestyle='-')
  mbol_ax.scatter(0,m_bol_trgb_Str20,label = label_m_bol_trgb_Str20, color = t_lost,zorder=4, marker='.',s = 1)
  mbol_ax.errorbar(0,m_bol_trgb_Str20,yerr=m_bol_trgb_Str20_unc,xerr=[[0],[2.30]], color = t_lost,zorder=4,fmt='-.', capsize=5)
  
  mbol_ax.axhline(m_bol_trgb_Str20 - m_bol_trgb_Str20_unc, label = ' ', 
                  color=t_lost,linewidth=0.6,linestyle='-')
  
  mbol_ax.axvline(2.30, label = ' ', 
                  color=t_lost,linewidth=0.6,linestyle='-')
  
  # mbol_ax.set_xlim(right = 2.5, left = -0.1)
  # mbol_ax.set_ylim(top = -3.89,bottom=-3.2)
  mbol_ax.legend()
  
   
    
  if savefig is True:
    mbol_fig.savefig(fig_sub_path+data_set_selected+'_mbol_trgb_Straniero20_comp_g13_full')

#%% Pretty Mbol comparison plot Straniero2020


if pretty_g13_eval_plot is True:
  
  mbol_fig, mbol_ax = plt.subplots(1,1,
                                   figsize = (fig_width,fig_height*0.4),
                                   sharex=(True),
                                   sharey=False)
  
  if ready_to_publish is False:
    mbol_fig.suptitle('Comparison of TRGB M_bol to Straniero20')
          
  
  mbol_ax.set_ylabel('$ M_{bol} [mag] $')  
  mbol_ax.set_xlabel('$ g_{13} $')
  mbol_ax.invert_yaxis()
  
  # mbol_ax.set_xticks()
  # mbol_ax.set_xlim()  
  # mbol_ax.set_ylim(top=15,bottom=11)  
  
  label_z = "[M/H] %0.3f  dex" % (z_3['MH'][0])
  mbol_ax.plot(z_3['g_13'],z_3['Mbol'],
               # label = label_z,
               color=t_fit,linestyle ='-',marker='^')
  
  label_z = "[M/H] %0.3f  dex" % (z_2['MH'][0])
  mbol_ax.plot(z_2['g_13'],z_2['Mbol'],
               # label = label_z,
               color=t_fit,linestyle ='-',marker='x')
  
  label_z = "[M/H] %0.3f  dex" % (z_1['MH'][0])
  mbol_ax.plot(z_1['g_13'],z_1['Mbol'],
               # label = label_z,
               color=t_fit,linestyle ='-',marker='o')
  
  
  
  label_m_bol_trgb = r"$[\frac{M}{H}]$ = (%0.2f $\pm$ %0.2f), $M_{bol}$=(%0.2f $\pm$ %0.2f) mag" % (-0.99, 0.20, m_bol_trgb, m_bol_trgb_unc)
    
  mbol_ax.axhline(m_bol_trgb+m_bol_trgb_unc, 
                  color=t_fit,linewidth=0.6,linestyle='--')
  mbol_ax.axhline(m_bol_trgb-m_bol_trgb_unc, 
                  color=t_fit,linewidth=0.6,linestyle='--')
  mbol_ax.axhline(m_bol_trgb, label = label_m_bol_trgb, 
                  color=t_fit,linewidth=0.6,linestyle='-')
  
  mbol_ax.axhspan(m_bol_trgb-m_bol_trgb_unc, m_bol_trgb+m_bol_trgb_unc, color=t00,alpha=0.4)  #horizontal shading
  
  # VI + JK with ZAHB distance scale
  m_bol_trgb_Str20 = np.mean([-3.60,-3.63])
  m_bol_trgb_Str20_unc = np.mean([0.25,0.26])
  
  label_m_bol_trgb_Str20 = r"$[\frac{M}{H}]$ = (%0.2f $\pm$ %0.2f), $M_{bol}$=(%0.2f $\pm$ %0.2f) mag" % (-1.11 ,0.20 , m_bol_trgb_Str20, m_bol_trgb_Str20_unc)
  # mbol_ax.axhline(m_bol_trgb_Str20, label = label_m_bol_trgb_Str20, 
  #                 color=t_fit,linewidth=0.6,linestyle='-')
  mbol_ax.scatter(0,m_bol_trgb_Str20,label = label_m_bol_trgb_Str20, color = t_lost,zorder=4, marker='.',s = 1)
  mbol_ax.errorbar(0,m_bol_trgb_Str20,yerr=m_bol_trgb_Str20_unc,xerr=[[0],[2.30]], color = t_lost,zorder=4,fmt='-.', capsize=5)
  
  mbol_ax.set_xlim(right = 2.5, left = -0.1)
  mbol_ax.set_ylim(top = -3.89,bottom=-3.2)
  mbol_ax.legend()
  
   
    
  if savefig is True:
    mbol_fig.savefig(fig_sub_path+data_set_selected+'_mbol_trgb_Straniero20_comp_g13_pretty')

#%% Mbol([M/H]_M5) metallicity comparison plot

   

if g13_interpolation_function_plot is True:
  
  mbol_fig, mbol_ax = plt.subplots(1,1,
                                   figsize = (fig_width,fig_height*0.4),
                                   sharex=(True),
                                   sharey=False)
  
  if ready_to_publish is False:
    mbol_fig.suptitle('Comparison of TRGB M_bol to Straniero20')
          
  
  mbol_ax.set_ylabel(r'$ M_{bol} [mag] $')  
  mbol_ax.set_xlabel(r'$ [\frac{M}{H}] $')
  mbol_ax.invert_yaxis()
  
  
  # mask = (z_3['g_13'] < 2.6)  
  u,v,w = 0,-1,4
  # zp1 = z_1[0:6:2]
  
  zp1 = z_1[u::w]
  zp2 = z_2[u::w]
  zp3 = z_3[u::w]
  
  label_z = r"$[\frac{M}{H}] \ \ \ \ \  \  \ \ = %0.2f$ " % (zp3['MH'][0])
  mbol_ax.axvline(zp3['MH'][0], 
                color=t_fit,
                label = label_z,              
                alpha = 0.6,
                zorder=1,
                linewidth=0.8,
                linestyle='-.')
  
  label_z = r"$[\frac{M}{H}] \ \ \ \ \  \  \ \ = %0.2f$ " % (zp2['MH'][0])                
  mbol_ax.axvline(zp2['MH'][0], 
                color=t_fit,
                label = label_z,              
                alpha = 0.6,
                zorder=1,
                linewidth=0.6,
                linestyle=':')
  
  label_z = r"$[\frac{M}{H}] \  \ \ \ \  \  \ \ = %0.2f$ " % (zp1['MH'][0])
  mbol_ax.axvline(zp1['MH'][0], 
                  label = label_z,
                  alpha = 0.6,
                  zorder=1,
                  color=t_fit,
                  linewidth=0.6,
                  linestyle='--')
  
  
  
  label_z = r"$ M_{bol}([\frac{M}{H}],g_{13_{const}} ) $ "
  mbol_ax.plot([],[],
                label = label_z,
               color=t_fit,
               linestyle ='-',
               marker='',
               linewidth=0.4
               # alpha = 0.6
               )

  
  for i, el in enumerate(zp1):
    
          
    mbol_ax.plot([zp3['MH'][i],zp2['MH'][i],zp1['MH'][i]], [zp3['Mbol'][i],zp2['Mbol'][i],zp1['Mbol'][i]],
                  # label = label_z,
                 color=t_fit,
                 linestyle ='',
                  marker='',
                  markersize = 2
                 )
    
    z_i = np.array([zp3['MH'][i],zp2['MH'][i],zp1['MH'][i]])
    mh_i = np.linspace(z_i.min()-0.3, z_i.max()+0.3,100)
    
    mbol_i = (mh_mbol_poly2_fit_list[u::w])[i](mh_i)
    
    mbol_ax.plot(mh_i,mbol_i,
                 # label = label_z,
                 color=t_fit,
                 linestyle ='-',
                 marker='',
                 linewidth=0.4
                 # alpha = 0.6
                 )
    
    label_z = r"$g_{13}$ = %0.1f " % (zp3['g_13'][i])
    # label_z = r" %0.1f " % (zp3['g_13'][i])
    mbol_ax.annotate( label_z, \
                      xy=(mh_i[25],mbol_i[25]), \
                     textcoords='offset pixels', \
                     xytext=(-60,18),
                     size= 8
                    )
  
  marker_size = 52
  # Mbol([M/H]_M5_max) 
  mh_G3_scatter = np.full(g_13_poly_mh_G3.size,mh_G3_margin_max)
  label_z = r"$[\frac{M}{H}]_{G3_{max}} \  = %0.2f $" % (mh_G3_margin_max)
  mbol_ax.scatter(mh_G3_scatter[u::w],mbol_poly_mh_G3_margin_max[u::w],
                label = label_z,
                color=t_fit,
                linestyle ='-',
                marker='*',
                s = marker_size
                )
  # Mbol([M/H]_M5)
  mh_G3_scatter = np.full(g_13_poly_mh_G3.size,mh_G3)
  label_z = r"$[\frac{M}{H}]_{G3} \ \ \ \  \  = %0.2f $" % (mh_G3)
  mbol_ax.scatter(mh_G3_scatter[u::w],mbol_poly_mh_G3[u::w],
                label = label_z,
                color=t_fit,
                linestyle ='-',
                marker='.',
                s = marker_size
                )
  
  # Mbol([M/H]_M5_min)
  mh_G3_scatter = np.full(g_13_poly_mh_G3.size,mh_G3_margin_min)
  label_z = r"$[\frac{M}{H}]_{G3_{min}} \ = %0.2f$ " % (mh_G3_margin_min)
  mbol_ax.scatter(mh_G3_scatter[u::w],mbol_poly_mh_G3_margin_min[u::w],
                label = label_z,
                color=t_fit,
                linestyle ='-',
                marker='+',
                s = marker_size
                )
  
  
  
  # mbol_ax.set_ylim(bottom= -3.3)
  # mbol_ax.set_xlim(left= -2.3, right = -0.2)  
  
  mbol_ax.legend(ncol = 1,  bbox_to_anchor=(1, 1.045))
  
   
    
  if savefig is True:
    mbol_fig.savefig(fig_sub_path+'M5_mbol_mh_trgb_poly2_fit_based_on_Straniero20_pretty')




#%% g13 interpolate Straniero2020 comparison plot

if g13_eval_plot is True:
  
  mbol_fig, mbol_ax = plt.subplots(1,1,
                                   figsize = (fig_width,fig_height*0.5),
                                   sharex=(True),
                                   sharey=False)
  
  if poly_fit_true is True:
    label_z = r"$[\frac{M}{H}]$ %0.3f  dex" % (mh_G3)
    mbol_ax.plot(
                  g_13_poly_mh_G3,
                  mbol_poly_mh_G3,
                  label = label_z,
                  color=t_fit,
                  linestyle ='-',
                  marker='o'
                  )        
  
  if ready_to_publish is False:
    mbol_fig.suptitle('Comparison of TRGB M_bol to Straniero20')
    
  mbol_ax.set_ylabel('$ M_{bol} [mag] $')  
  mbol_ax.set_xlabel('$ g_{13} $')
  mbol_ax.invert_yaxis()
  
  
  
  marker = '.'
  marker_size = 6
  
  
  # [M/H]_M5: 
  label_z = r"$[\frac{M}{H}]_{G3_{max}} \ = %0.2f$ " % (mh_G3_margin_max)
  mbol_ax.plot(
                g_13_poly_mh_G3,
                mbol_poly_mh_G3_margin_max,
                label = label_z,
                color=t_fit,
                linestyle ='-',
                linewidth=0.6,
                zorder = 3,
                # marker=marker,
                marker='*',
                markersize = marker_size
                )
    
  label_z_sum = r"$ [\frac{M}{H}]_{G3} \ $ = $%0.2f < %0.2f < %0.2f$ " \
                      % (mh_G3_margin_min, \
                         mh_G3, \
                         mh_G3_margin_max)
  label_z = r"$[\frac{M}{H}]_{G3} \ \ \ \  \  = %0.2f $" % (mh_G3)
  mbol_ax.plot(
                g_13_poly_mh_G3,
                mbol_poly_mh_G3,
                label = label_z,
                # label = label_z_sum,
                color=t_fit,
                linestyle ='-',
                linewidth=0.6,
                zorder = 3,
                marker=marker,
                markersize = marker_size
                )
  
    
  label_z = r"$[\frac{M}{H}]_{G3_{min}} \  = %0.2f $" % (mh_G3_margin_min)
  mbol_ax.plot(
                g_13_poly_mh_G3,
                mbol_poly_mh_G3_margin_min,
                label = label_z,
                color=t_fit,
                linestyle ='-',
                linewidth=0.6,
                zorder = 3,
                # marker=marker,
                marker='+',
                markersize = marker_size
                )
  
  
  # Straniero20 Mbol_TRGB for mh based on feh_Carreta2010 and afe = 0.4
  # VI + JK with ZAHB distance scale
  m_bol_trgb_Str20 = np.mean([-3.60,-3.63])
  m_bol_trgb_Str20_unc = np.mean([0.25,0.26])
  
  label_m_bol_trgb_Str20 = r"$ (M_{bol})_{S20}  = (%0.2f\  \pm %0.2f)$ mag" % (m_bol_trgb_Str20, m_bol_trgb_Str20_unc)  
  label_mh_Str20 = r"$[\frac{M}{H}]_{S20} = (%0.2f \pm %0.2f)$ dex" % (-1.11, 0.20)  
  # label_m_bol_trgb_Str20 = r"$ (M_{bol})_{S20} $ " 
  # mbol_ax.axhline(m_bol_trgb_Str20+m_bol_trgb_Str20_unc, 
  #                 color=t_fit,
  #                 alpha = 0.6,
  #                 zorder=1,
  #                 linewidth=0.6,
  #                 linestyle='--')
  # mbol_ax.axhline(m_bol_trgb_Str20-m_bol_trgb_Str20_unc, 
  #                 color=t_fit,
  #                 alpha = 0.6,
  #                 zorder=1,
  #                 linewidth=0.8,
  #                 linestyle='--')
  mbol_ax.axhline(m_bol_trgb_Str20, 
                  label = label_m_bol_trgb_Str20, 
                  alpha = 0.6,
                  zorder=1,
                  color=t_fit,
                  linewidth=0.6,
                  linestyle='-')
  #horizontal shading
  mbol_ax.axhspan(
                  m_bol_trgb_Str20-m_bol_trgb_Str20_unc, 
                  m_bol_trgb_Str20+m_bol_trgb_Str20_unc, 
                  # label = label_mh_Str20,
                  color=t00,
                  alpha=0.4)  
  
  
  
  
  
  
  
  # Mbol TRGB:
    # unc range
  label_m_bol_trgb = r"$ (M_{bol})_{G3} = (%0.2f < %0.2f < %0.2f$) mag" \
                      % (m_bol_trgb+m_bol_trgb_unc_cl95, \
                         m_bol_trgb, \
                         m_bol_trgb_upper_bound_cl95)
    # unc upper & lower bound
  label_m_bol_trgb = r"$ (M_{bol})_{G3} = (%0.2f^{-%0.2f}_{+%0.2f})$ mag" \
                      % (m_bol_trgb, \
                         m_bol_trgb-m_bol_trgb_upper_bound_cl95, \
                         m_bol_trgb_unc_cl95)
    # unc set to upper bound                      
  label_m_bol_trgb = r"$ (M_{bol} )_{G3} \ = (%0.2f\  \pm {%0.2f})$ mag" \
                      % (m_bol_trgb, \
                         m_bol_trgb_unc_cl95, \
                         )                        
                        
  # label_m_bol_trgb = r"$ (M_{bol})_{G3}$ "
  
  # Mbol TRGB & unc
  yerr = [[(m_bol_trgb_unc_cl95)],[m_bol_trgb_unc_cl95]] 
  mbol_ax.errorbar(
                    g_13_M5_TRGB_G3,
                    m_bol_trgb,
                    yerr=yerr,
                    # xerr=[[0],[2.30]], 
                    label = label_m_bol_trgb, 
                    color = t_lost,
                    zorder=4, 
                    fmt='.',
                    markersize = 10,
                    capsize=5
                    )
  # mbol_ax.hlines( 
  #                 y=m_bol_trgb-m_bol_trgb_unc_cl95, 
  #                 xmin=g_13_M5_TRGB_G3, 
  #                 xmax=g_13_M5_TRGB_G3_upper_bound, 
  #                 # label=label_g_13_M5_TRGB_G3_upper_bound,
  #                 color=t_lost, 
  #                 linestyle='--', 
  #                 zorder=4,
  #                 linewidth=0.6
  #                 )
  
  # mbol_ax.scatter( 
  #                 g_13_M5_TRGB_G3_upper_bound, 
  #                 m_bol_trgb_upper_bound_cl95, 
  #                 # label=label_g_13_M5_TRGB_G3_upper_bound,
  #                 color = t_lost,
  #                 zorder=5, 
  #                 marker=r'$\mathbf{x}$',
  #                 s = 30
  #                 )             
  
  # g13 for Mbol TRGB upper bound:
  label_g_13_M5_TRGB_G3_upper_bound = r"$ (g_{13})_{G3} \ $ = ($%0.2f <  (g_{13})_{G3} < %0.2f $)" \
                            % (g_13_M5_TRGB_G3_cl95_limit_dim, \
                               g_13_M5_TRGB_G3_upper_bound)
  label_g_13_M5_TRGB_G3_upper_bound = r"$ (g_{13})_{G3} \ \ \ \leq \ \ \ \ %0.2f $" \
                            % (g_13_M5_TRGB_G3_upper_bound)
                            
                            
  label_g_13_M5_TRGB_G3_upper_bound = r"$ (g_{13} )_{G3} \ \ \ = \ \ %0.2f^{+%0.2f}_{-%0.2f}$ " \
                            % (g_13_M5_TRGB_G3_upper_bound, 
                               np.abs(g_13_M5_TRGB_G3_cl95_limit_bright_max_dev),
                               np.abs(g_13_M5_TRGB_G3_cl95_limit_bright_min_dev)
                               )
  # label_g_13_M5_TRGB_G3_upper_bound = r"$ (g_{13})_{G3} \leq %0.2f$, $(M_{up} )_{G3} \ = %0.2f$ mag " \
  #                           % (g_13_M5_TRGB_G3_upper_bound, m_bol_trgb_upper_bound_cl95 )
  
  # label_g_13_M5_TRGB_G3_upper_bound = r"$ (g_{13})_{G3} $"  
  
  label_M5_TRGB_G3_upper_bound = r"$ (M_{up} )_{G3} \ = %0.2f $ mag " \
                            % (m_bol_trgb_upper_bound_cl95
                               )
  mbol_ax.axhline( 
                  y=m_bol_trgb_upper_bound_cl95, 
                  label=label_M5_TRGB_G3_upper_bound,
                  color=t_lost, 
                  linestyle='-', 
                  zorder=4,
                  linewidth=0.6
                  )
  
  mbol_ax.axvline( 
                  x=g_13_M5_TRGB_G3_upper_bound, 
                  label=label_g_13_M5_TRGB_G3_upper_bound,
                  color=t_lost, 
                  linestyle='--', 
                  zorder=4,
                  linewidth=0.6
                  )
  
  #vertical shading
  mbol_ax.axvspan(
                  g_13_M5_TRGB_G3_upper_bound + g_13_M5_TRGB_G3_cl95_limit_bright_max_dev, 
                  g_13_M5_TRGB_G3_upper_bound + g_13_M5_TRGB_G3_cl95_limit_bright_min_dev, 
                  # label = label_mh_Str20,
                  color=t_lost,
                  alpha=0.2)  
  
  
  # mbol_ax.scatter( 
  #                 g_13_M5_TRGB_G3_upper_bound, 
  #                 m_bol_trgb_upper_bound_cl95, 
  #                 label=label_g_13_M5_TRGB_G3_upper_bound,
  #                 color = t_lost,
  #                 zorder=5, 
  #                 marker=r'$\mathbf{x}$',
  #                 s = 30
  #                 )        
  
     
    
  mbol_ax.set_xlim(right = 2.5, left = -0.1)
  mbol_ax.set_ylim(top = -3.89,bottom=-3.2)
  mbol_ax.set_xlim(right = 3, left = -0.1)
  mbol_ax.set_ylim(top = -4.0,bottom=-3.15)
  
  
  
  # Get the handles and labels for the legend
  handles, labels = mbol_ax.get_legend_handles_labels()
  
  # Reorder the handles and labels
  new_order = np.arange(0,len(handles))
  # new_order[-1] = new_order[-1]-1 
  # new_order[-2] = new_order[-2]+1 
  
  new_order_1 = new_order[0:3]
  handles_1 = [handles[i] for i in new_order_1]
  labels_1 = [labels[i] for i in new_order_1]
  legend1 = mbol_ax.legend(handles_1, labels_1,ncol=1, loc= 'lower left', bbox_to_anchor=(0.14, 0))
  
  new_order_2 = new_order[3:]
  new_order_2[:] = 3,6,4,5
  handles_2 = [handles[i] for i in new_order_2]
  labels_2 = [labels[i] for i in new_order_2]
  legend2 =mbol_ax.legend(handles_2, labels_2,ncol=1, loc= 'lower left', bbox_to_anchor=(0.48, 0))
  
  
  
  
  # Create the legend with the new order
  # mbol_ax.legend(handles, labels,ncol=2)
  # mbol_ax.legend(ncol=2)
  
  # Create two separate legends
  # legend1 = mbol_ax.legend(loc='upper left')
  # legend2 = mbol_ax.legend(loc='upper right')
  
  # Add the legends to the figure
  mbol_ax.add_artist(legend1)
  
  # Set labels for the two legends
  # legend1.get_texts()[0].set_text('Legend 1')
  # legend2.get_texts()[0].set_text('Legend 2')
  
   
    
  if savefig is True:
    mbol_fig.savefig(fig_sub_path+data_set_selected+'_mbol_trgb_g13')