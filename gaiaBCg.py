#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 12:38:43 2023

@author: tom
"""


import csv
import numpy as np
import os
import glob
import matplotlib.pyplot as plt

import GAIA_plot as gp

#%% ### ------------ Set default plot style------------------------------------
#https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure

# Note: fig_width is fixed in plotStyle.py for MasterThesis specifications
#       fig_height can and should be adjusted to special needs of individual plots
#       fig_height of 5inch is approximatelly half of height on A4 MasterThesis      
fig_height = 5
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)



### ------------ Set default plot colors --------------------------------------

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

#%% ### ------------ Load /Store data settings ---------------------------------------

folder_name = 'data/clusterMemberSelection/'

## -------------- GAIA data ---------------------------------------------------
# data_set_selected = 'Gaia_M3_r-18arcmin.csv'
data_set_selected = 'Gaia_M5_r-9arcmin_zpc_pm_clusterMemberSelection_distEval_absG-band'
# data_set_selected = 'Gaia_M67_r-25arcmin.csv'

file_name = data_set_selected

## ------------ load data into variables --------------------------------------
import GAIA_load as gl
from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

gaia_cluster_members = file_data
label_gaia_cluster_members = 'cluster members'
label_gaia_cluster_members_n = 'cluster members ('+str(gaia_cluster_members.size)+' stars)'

#---------Load BCg table --------------------------------
bc_table_path = 'data/bcEvalData/bc_dr3_feh_all'

file_data = gl.load_query(bc_table_path)

bcg_table = file_data
label_bcg_table = 'BC G-filter'
label_bcg_table_n = 'BC G-filter ('+str(bcg_table.size)+' stars)'

#%%

#---------Store GAIA data --------------------------------

file_sub_path = folder_name
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### ----------- Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/clus_sel/'


# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)

bcg_table_plot_on = True



#%% Definitions



# point is a list of 4 elements : [teff, logg, metallicity ([Fe/H]), alpha/Fe]
# offset is an optional floating value (0 by default), see below for details.


 
#%% ### RGB AGB star selection -------------------------------------


rgb_agb_lower_bound = 12.3

#StarHorse extinction
# brightness_modul = 14.55912732455515
# abs_rgb_agb_lower_bound = rgb_agb_lower_bound-brightness_modul 

mask_rgb_agb = (gaia_cluster_members['phot_g_mean_mag'] < rgb_agb_lower_bound)

#%% M5 NGC 5904 [Fe/H] literature values and uncertainty estimation with 95cl

# [fe/H] ----------------------------------------------------
# Harris 1996 (2010 edition)  (used by Viaux13)
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



#%% # cluster members with gspphot evaluation--------------------------------------

print()
# cluster members with Teff
mask_teff = np.logical_not( np.isnan(gaia_cluster_members['teff_gspphot']) )
gaia_cluster_members_teff = gaia_cluster_members[mask_teff]
print('gaia_cluster_members_teff: '+str(gaia_cluster_members_teff.size))
gaia_cluster_members_teff_rgb = gaia_cluster_members_teff[(gaia_cluster_members_teff['phot_g_mean_mag'] < rgb_agb_lower_bound)]
print('gaia_cluster_members_teff_rgb: '+str(gaia_cluster_members_teff_rgb.size))

# cluster members with metallicity [M/H]
mask_mh = np.logical_not( np.isnan(gaia_cluster_members['mh_gspphot']) )
gaia_cluster_members_mh = gaia_cluster_members[mask_mh]
print('gaia_cluster_members_mh: '+str(gaia_cluster_members_mh.size))
gaia_cluster_members_mh_rgb = gaia_cluster_members_teff[(gaia_cluster_members_mh['phot_g_mean_mag'] < rgb_agb_lower_bound)]
print('gaia_cluster_members_mh_rgb: '+str(gaia_cluster_members_mh_rgb.size))

# cluster members with surface gravity log(g)
mask_logg = np.logical_not(np.isnan(gaia_cluster_members['logg_gspphot']))
gaia_cluster_members_logg = gaia_cluster_members[mask_logg]
print('gaia_cluster_members_logg: '+str(gaia_cluster_members_logg.size))
gaia_cluster_members_logg_rgb = gaia_cluster_members_logg[(gaia_cluster_members_logg['phot_g_mean_mag'] < rgb_agb_lower_bound)]
print('gaia_cluster_members_logg_rgb: '+str(gaia_cluster_members_logg_rgb.size))


#%% # cluster members with gspspec evaluation--------------------------------------

print()

# cluster members with ionized iron - alpha-elements [alpha/Fe]
mask_alphafe = np.logical_not(np.isnan(gaia_cluster_members['alphafe_gspspec']))
gaia_cluster_members_alphafe = gaia_cluster_members[mask_alphafe]
print('gaia_cluster_members_alphafe: '+str(gaia_cluster_members_alphafe.size))
gaia_cluster_members_alphafe_rgb = gaia_cluster_members_alphafe[(gaia_cluster_members_alphafe['phot_g_mean_mag'] < rgb_agb_lower_bound)]
print('gaia_cluster_members_alphafe_rgb: '+str(gaia_cluster_members_alphafe_rgb.size))

# cluster members with neutral iron [Fe/M]
mask_fem = np.logical_not( np.isnan(gaia_cluster_members['fem_gspspec']) )
gaia_cluster_members_fem_gspspec = gaia_cluster_members[mask_fem]
print('gaia_cluster_members_fem_gspspec: '+str(gaia_cluster_members_fem_gspspec.size))
gaia_cluster_members_fem_gspspec_rgb = gaia_cluster_members_fem_gspspec[(gaia_cluster_members_fem_gspspec['phot_g_mean_mag'] < rgb_agb_lower_bound)]
print('gaia_cluster_members_fem_gspspec_rgb: '+str(gaia_cluster_members_fem_gspspec_rgb.size))

# -----------------------------------------------------------------------------
#%% BC data selection 5951

mask_rgb_agb = (gaia_cluster_members['phot_g_mean_mag'] < rgb_agb_lower_bound)

print()
print('gaia_cluster_members: '+str(gaia_cluster_members.size))
# gaia_cluster_members_rgb = gaia_cluster_members[(gaia_cluster_members['phot_g_mean_mag'] < rgb_agb_lower_bound)]
gaia_cluster_members_rgb = gaia_cluster_members[mask_rgb_agb]
print('gaia_cluster_members_rgb: '+str(gaia_cluster_members_rgb.size))

# gsp_spec population very low!!! Change with DR4/5?
# mask_bc = mask_teff & mask_mh & mask_logg & mask_alphafe & mask_fem

# gsp_phot population adequat?
mask_bc = mask_teff & mask_mh & mask_logg
# Necessary literature values
# feh[:]= feh_harris10
# alphafe = ?


print()
# gaia_cluster_members_bc = gaia_cluster_members
gaia_cluster_members_bc = gaia_cluster_members[mask_bc]
print('gaia_cluster_members_bc: '+str(gaia_cluster_members_bc.size))

gaia_cluster_members_bc_rgb = gaia_cluster_members_bc[(gaia_cluster_members_bc['phot_g_mean_mag'] < rgb_agb_lower_bound)]
print('gaia_cluster_members_bc: '+str(gaia_cluster_members_bc_rgb.size))

# Note: calc bc for all valid cluster members
gaia_cluster_members_bc = gaia_cluster_members_bc
mask_bc = mask_bc 

# Note: calc bc for valid cluster members above rgb agb threshold (Saves computing time! )
# gaia_cluster_members_bc = gaia_cluster_members_bc_rgb
# mask_bc = mask_bc & mask_rgb_agb


#%% Prepare data selection

# Using Gaia data (all [val] have unit dex,)

teff = gaia_cluster_members_bc['teff_gspphot']
logg = gaia_cluster_members_bc['logg_gspphot']

# Global metallicity [M/H] from GSP-Spec MatisseGauguin using RVS spectra and Monte Carlo realisations (float, Abundances[dex])
# Median global metallicity (assuming source is a single star) inferred by GSP-Spec MatisseGauguin (Recio-Blanco et al. 2022) from RVS spectra. 
mh = gaia_cluster_members_bc['mh_gspphot'] 


# alphafe_gspspec : Abundance of alpha-elements [alpha/Fe] with respect to iron from GSP-Spec MatisseGauguin using RVS 
alpha_Fe = gaia_cluster_members_bc['alphafe_gspspec'] 
# Abundance of neutral iron [Fe/M] from GSP-Spec MatisseGauguin using RVS spectra and Monte Carlo realisations, applied to the individual N lines of the element, given in fem_gspspec_nlines (float, Abundances[dex])
# Median abundance of neutral iron (assuming source is a single star) from RVS spectra and Monte Carlo realisations derived using MatisseGauguin (Recio-Blanco et al. 2022) atmospheric parameters and the Gauguin algorithm, applied to the individual N lines of the element, where the number of lines is given in fem_gspspec_nlines. 
fem = gaia_cluster_members_bc['fem_gspspec'] 


# The neutral iron abundance [Fe/H] is obtained by [Fe/H]=[Fe/M]+[M/H]. (float, Abundances[dex])
feh = fem + mh



# %%  BC Table

bcg_table_matrix = [0,0,0,0,0]
bcg_table_matrix = [  
                                  bcg_table['teff'], 
                                  bcg_table['logg'],
                                  bcg_table['feh'],
                                  bcg_table['alphafe'],
                                  bcg_table['bcg']
                                  ]

teff_m, logg_m, feh_m, alphaFe_tab, bcg_m = bcg_table_matrix 

# Fe/H error estimation for missing fem data in Gaia query
feh_tab_min  = bcg_table['feh'].min()
feh_tab_max  = bcg_table['feh'].max()
feh_tab_mean = np.abs(feh_tab_min-feh_tab_max) / 2

fem_nan = feh_tab_mean - mh
fem_nan_min = feh_tab_min - mh
fem_nan_max = feh_tab_max - mh



#alpha/Fe error estimation
afe_tab_min  = bcg_table['alphafe'].min()
afe_tab_max  = bcg_table['alphafe'].max()
afe_tab_mean = np.abs(afe_tab_min-afe_tab_max) / 2



fig_name = label_bcg_table

if bcg_table_plot_on is True:
  
  fig, ax = plt.subplots(2,2)
  s = 0.5
  ax[0,0].scatter(teff_m,   bcg_m,s=s,label='teff')
  ax[0,0].set_xlabel('teff')
  ax[0,1].scatter(logg_m,   bcg_m,s=s,label='logg')
  ax[0,1].set_xlabel('logg')
  ax[1,0].scatter(feh_m,    bcg_m,s=s,label='Fe/H')
  ax[1,0].set_xlabel('Fe/H')
  ax[1,1].scatter(alphaFe_tab,bcg_m,s=s,label='Alph/Fe')
  ax[1,1].set_xlabel('Alph/Fe')
  # fig.legend()



#%% Set literature values for BC eval

# Necessary literature values
# feh_array = np.array([conf_interval_margin[0],feh_mean,conf_interval_margin[1]])  
# alpha_Fe =...

# afe_tab_max
# afe_tab_min
# afe_tab_mean

#%% Calculate BCg

# Using gdr3bcg to get bc_g
# source: https://gitlab.oca.eu/ordenovic/gaiadr3_bcg
import gdr3bcg.bcg as bcg
table=bcg.BolometryTable()

# prepare BC arras with NaN values
bc_all = np.copy(gaia_cluster_members['teff_gspphot'])
bc_all[:] = np.nan
bc_all_unc = np.copy(bc_all)
m_bol_all     = np.copy(bc_all)
m_bol_all_unc = np.copy(bc_all)

# prepare BC arrays without NaN values
bc_g = np.copy(gaia_cluster_members_bc['teff_gspphot'])
bc_g[:] = np.nan
bc_g_unc = np.copy(bc_g)

bc_g_mean = np.copy(bc_g)
bc_g_max = np.copy(bc_g)
bc_g_min = np.copy(bc_g)

# prepare M_bol arrays without NaN values
# Note: Not necessary, reassigned below!
m_bol = np.copy(bc_g)
m_bol_unc = np.copy(bc_g)

# Conversion from IAU B2 2015 convention to Straniero20
offset = 0.01

# Calc BC ------------------------------------------------------

gsp_phot_spec_eval_true = False
gsp_phot_eval_true = True

# Evaluation with gsp_spec values (sample size too small)
if gsp_phot_spec_eval_true is True:
  for i, el in enumerate(bc_g):
    # stars_bc = [teff[0], logg[0], feh[0], alpha_Fe[0]]
    stars_bc = [teff[i], logg[i], feh[i], alpha_Fe[i]]
    bc_g[i] = table.computeBc(stars_bc,offset)
    # print(str(gaia_cluster_members_bc['designation'][i])+' BC in mag: '+str(bc_g[i]))
    bc_g_unc = 0

# Evaluation only with gsp_phot, literature values, and alpha_Fe min/max
if gsp_phot_eval_true is True:
  for i, el in enumerate(bc_g):
    # stars_bc = [teff[0], logg[0], feh[0], alpha_Fe[0]]
    stars_bc = [teff[i], logg[i], feh_array[1], afe_tab_mean]
    bc_g_mean[i] = table.computeBc(stars_bc,offset)
    print(str(gaia_cluster_members_bc['designation'][i])+' BC in mag: '+str(bc_g_mean[i]))
  
  
    # stars_bc = [teff[0], logg[0], feh[0], alpha_Fe[0]]
    stars_bc = [teff[i], logg[i], feh_array[0], afe_tab_min]
    bc_g_min[i] = table.computeBc(stars_bc,offset)
    # print(str(gaia_cluster_members_bc['designation'][i])+' BC min: '+str(bc_g_min[i]))
  
  
    # stars_bc = [teff[0], logg[0], feh[0], alpha_Fe[0]]
    stars_bc = [teff[i], logg[i], feh_array[2], afe_tab_max]
    bc_g_max[i] = table.computeBc(stars_bc,offset)
    # print(str(gaia_cluster_members_bc['designation'][i])+' BC max: '+str(bc_g_max[i]))

  # Test symmetric error distribution
  error_ratio = np.abs(bc_g_mean -bc_g_max) / np.abs(bc_g_mean -bc_g_min)
  error_ratio_mean = error_ratio.mean()
  error_ratio_std  = error_ratio.std()
  
  # Note: If error distribution symmetric
  bc_g_unc = np.abs(bc_g_mean -bc_g_max)
  bc_g     = bc_g_mean




#%% # Calc M_bol -------------------------------------------------------

# BC = M_bol - M_G
# M_bol = BC + M_G
m_bol = bc_g + (gaia_cluster_members_bc['abs_phot_g_mean_mag'])

# quadratic error?!
# M_bol_unc = np.squrt(BC_unc² + M_G_unc²)
m_bol_unc = np.sqrt(bc_g_unc**2 + ((gaia_cluster_members_bc['abs_phot_g_mean_mag_error']))**2 )

# print(m_bol)

#%% # Store in gaia_cluster_members_bc ------------------------------------

bc_all[mask_bc] = bc_g
bc_all_unc[mask_bc] = bc_g_unc

m_bol_all[mask_bc] = m_bol
m_bol_all_unc[mask_bc] = m_bol_unc

# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'bcg_mag', bc_all, dtypes='<f8')
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'bcg_mag_error', bc_all_unc , dtypes='<f8')
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'mbol_mag', m_bol_all, dtypes='<f8')
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'mbol_mag_error', m_bol_all_unc , dtypes='<f8')

save_gaia(gaia_cluster_members, file_name +'_mbol' ,file_sub_path=file_sub_path)




#%%  BCg visualization plot 
# Note: Includes an approximation of BC but gaiaBCg_approximation.py should
#       be used for that! 
# Note: Clean up in future!

import numpy as np
from sklearn.linear_model import LinearRegression
coord_in_plot = -3.2
slope = 1
delta = 0.1

# mask = mask_rgb_agb & mask_bc
mask = mask_bc
x, y = (gaia_cluster_members['mbol_mag'])[mask],(gaia_cluster_members['abs_phot_g_mean_mag'])[mask]

model = LinearRegression().fit(x.reshape(-1, 1), y)
r_sq = model.score(x.reshape(-1, 1), y)
k = model.coef_[0]
intercept = model.intercept_

filter_outlier = (x<1.1) & (y>-0.1) & (y > x*k+intercept)

filter_m = np.logical_not(filter_outlier)

x_f,y_f = x[filter_m],y[filter_m]
model_f = LinearRegression().fit(x_f.reshape(-1, 1), y_f)
r_sq_f = model.score(x_f.reshape(-1, 1), y_f)
k_f = model_f.coef_[0]
intercept_f = model_f.intercept_


model_1 = LinearRegression().fit(x[(y<-0.1)].reshape(-1, 1), y[(y<-0.1)])
# r_sq_1 = model.score(x.reshape(-1, 1), y)
k_1 = model_1.coef_[0]
intercept_1 = model_1.intercept_

model_2 = LinearRegression().fit(x[(x>1.1)].reshape(-1, 1), y[(x>1.1)])
# r_sq_2 = model.score(x.reshape(-1, 1), y)
k_2 = model_2.coef_[0]
intercept_2 = model_2.intercept_

x_bcg, y_Mg = (gaia_cluster_members['bcg_mag'])[mask],(gaia_cluster_members['abs_phot_g_mean_mag'])[mask]


# fig_bc, axs_bc = plt.subplots(2,2,sharey=True)
# ax_Mg       = axs_bc[0,1]
# ax_bc       = axs_bc[0,0]
# ax_cmd      = axs_bc[1,0]
# ax_cmd_bol  = axs_bc[1,1]

fig_bc, axs_bc = plt.subplots(1,4,figsize=(fig_width,fig_height*0.6),gridspec_kw={'width_ratios':[1,1,3,3]},sharey=True)
ax_Mg       = axs_bc[1]
ax_bc       = axs_bc[0]
ax_cmd      = axs_bc[2]
ax_cmd_bol  = axs_bc[3]

ax_Mg.scatter(x[filter_m],y[filter_m],s=0.5,color='black')
ax_Mg.scatter(x[filter_outlier],y[filter_outlier],color='red',s=0.5,alpha=0.4)

# ax_Mg.axline([coord_in_plot,coord_in_plot],slope=slope,color ='black',linestyle='-')

# ax_Mg.axline([0+coord_in_plot,intercept_1+coord_in_plot*k_1],slope=k_1,color ='green',linestyle='--')
# ax_Mg.axline([0+coord_in_plot,intercept_2+coord_in_plot*k_2],slope=k_2,color ='blue',linestyle='-.')
ax_Mg.axline([0+coord_in_plot,intercept  +coord_in_plot*k  ],slope=k,  color = 'red')

ax_Mg.axline([0+coord_in_plot,intercept_f  +coord_in_plot*k_f  ],slope=k_f,  color = 'blue')
# ax_Mg.invert_yaxis()
ax_Mg.set_xlabel('$M_{bol}$')
ax_Mg.set_ylabel('$M_{G}$')



ax_bc.scatter(x_bcg[filter_m],y_Mg[filter_m],s=0.5,color='black')
ax_bc.scatter(x_bcg[filter_outlier],y_Mg[filter_outlier],color='red',s=0.5,alpha=0.4)

# ax_bc.invert_yaxis()
ax_bc.set_xlabel('$BC_{G}$')
ax_bc.set_ylabel('$M_{G}$')



M_g  = (gaia_cluster_members['abs_phot_g_mean_mag'])[mask]


# m_bol  = (gaia_cluster_members['mbol_mag'])[mask]
# m_bol_unc  = (gaia_cluster_members['mbol_mag_error'])[mask]
bp_rp = (gaia_cluster_members['bp_rp'])[mask]

ax_cmd.scatter(gaia_cluster_members['bp_rp'],gaia_cluster_members['abs_phot_g_mean_mag'],color='green',s=0.5)
ax_cmd.scatter(bp_rp[filter_m],M_g[filter_m],color='black',s=0.5)
ax_cmd.scatter(bp_rp[filter_outlier],M_g[filter_outlier],color='red',s=0.5)

# ax_cmd.invert_yaxis()
ax_cmd.set_xlabel('bp-rp')
ax_cmd.set_ylabel('$M_{G}$')



mbol  = (gaia_cluster_members['mbol_mag'])[mask]
bp_rp = (gaia_cluster_members['bp_rp'])[mask]

ax_cmd_bol.scatter(bp_rp[filter_m],mbol[filter_m],color='black',s=0.5)
ax_cmd_bol.scatter(bp_rp[filter_outlier],mbol[filter_outlier],color='red',s=0.5)


ax_cmd_bol.invert_yaxis()
ax_cmd_bol.set_xlabel('bp-rp')
ax_cmd_bol.set_ylabel('$M_{bol}$')



fig_test,ax_test = plt.subplots(1,1)
ax_test.scatter(bc_g_mean,bc_g_unc,s=0.6,label='BCg')
ax_test.scatter(m_bol,bc_g_unc,s=0.6, label='Mbol')
ax_test.set_ylabel('BCg uncertainty')
ax_test.legend()