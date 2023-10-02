 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 22:08:24 2023

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

data_set_selected = 'Gaia_M5_r-9arcmin_zpc_pm_clusterMemberSelection_distEval_absG-band_mbol'


file_name = data_set_selected

## ------------ load data into variables --------------------------------------
import GAIA_load as gl
from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

gaia_cluster_members = file_data
label_gaia_cluster_members = 'cluster members'
label_gaia_cluster_members_n = 'cluster members ('+str(gaia_cluster_members.size)+' stars)'



#%%

#---------Store GAIA data --------------------------------

file_sub_path = folder_name
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### ----------- Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/clus_sel/'
fig_name = data_set_selected

# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)

int_stel_ext_plot_on = True
spat_ext_plot = False
CMD_abs_plot_on = True

lin_reg_plots_on = True
 
#%% ### RGB AGB star selection -------------------------------------


rgb_agb_lower_bound = 12.3

# brightness module from StarHorse extinction with cl95 unc 
brightness_modul = gaia_cluster_members['distance_modulus_starhorse'][0]
brightness_modul_unc = gaia_cluster_members['distance_modulus_starhorse_error_cl95'][0]

abs_rgb_agb_lower_bound = rgb_agb_lower_bound-brightness_modul 

mask_rgb_agb = (gaia_cluster_members['phot_g_mean_mag'] < rgb_agb_lower_bound)

#%% Find linear regression for BC_G(M_G)

import numpy as np
from sklearn.linear_model import LinearRegression

mask_bc_aprx = np.isnan(gaia_cluster_members['bcg_mag'])
mask_bc = np.logical_not( mask_bc_aprx )


# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'bcg_aprx_flag', mask_bc_aprx, dtypes='bool')

#%%
coord_in_plot = -3.2
slope = 1
delta = 0.1

# mask = mask_rgb_agb & mask_bc
mask = mask_bc
mbol = (gaia_cluster_members['mbol_mag'])[mask]  
Mg = (gaia_cluster_members['abs_phot_g_mean_mag'])[mask]
bcg = (gaia_cluster_members['bcg_mag'])[mask]

# mask = mask_rgb_agb & mask_bc_aprx
mask = mask_bc_aprx
mbol_aprx = (gaia_cluster_members['mbol_mag'])[mask]  
Mg_aprx = (gaia_cluster_members['abs_phot_g_mean_mag'])[mask]
Mg_aprx_unc = (gaia_cluster_members['abs_phot_g_mean_mag_error'])[mask]
bcg_aprx = (gaia_cluster_members['bcg_mag'])[mask]


#%% # First attempt to find Linear Regression Mbol

model_mbol_1 = LinearRegression().fit(mbol.reshape(-1, 1), Mg)
r_sq = model_mbol_1.score(mbol.reshape(-1, 1), Mg)
k_mbol_1 = model_mbol_1.coef_[0]
i_mbol_1 = model_mbol_1.intercept_

# First attempt to find Linear Regression Mbol
model_bcg_1 = LinearRegression().fit(bcg.reshape(-1, 1), Mg)
r_sq = model_bcg_1.score(bcg.reshape(-1, 1), Mg)
k_bcg_1 = model_bcg_1.coef_[0]
i_bcg_1 = model_bcg_1.intercept_


if lin_reg_plots_on is True:
  fig_linReg_1, ax_linReg_1 = plt.subplots(1,2,sharey=True)
  
  
  ax_Mg = ax_linReg_1[0]
  ax_bc = ax_linReg_1[1]
  
    
  ax_Mg.scatter(mbol,Mg,s=0.5,color='black')
  
  ax_Mg.axline([coord_in_plot,coord_in_plot],slope=slope,color ='black',linewidth = 0.6,linestyle='-.')
  ax_Mg.axline([0+coord_in_plot,i_mbol_1 + coord_in_plot*k_mbol_1  ],slope=k_mbol_1,  color = 'red')
  ax_Mg.axhline(abs_rgb_agb_lower_bound,color = 'black',linewidth = 0.6,linestyle ='--')
  ax_Mg.set_xlabel('$M_{bol}$')
  ax_Mg.set_ylabel('$M_{G}$')

  ax_bc.scatter(bcg,Mg,s=0.5,color='black')
  ax_bc.axvline(0, color = 'black',linewidth = 0.6,linestyle ='-.')
  ax_bc.axline([0+coord_in_plot,i_bcg_1 + coord_in_plot*k_bcg_1  ],slope=k_bcg_1,  color = 'red')
  ax_bc.axhline(abs_rgb_agb_lower_bound,color = 'black',linewidth = 0.6,linestyle ='--')
  
  ax_bc.invert_yaxis()
  ax_bc.set_xlabel('$BC_{G}$')
  ax_bc.set_ylabel('$M_{G}$')
  
  if savefig is True:
    fig_linReg_1.savefig(fig_sub_path+fig_name+'_bcLinAprx_1.jpg')


#%%
# 

m_bol_outlier_limit =  1.1
Mg_outlier_limit    = -0.1

filter_outlier = (mbol< m_bol_outlier_limit) & (Mg> Mg_outlier_limit) \
                                              & (Mg > mbol*k_mbol_1+i_mbol_1)
filter_m = np.logical_not(filter_outlier)

mbol_sel = mbol[filter_m]
Mg_sel = Mg[filter_m]
bcg_sel_poly = bcg[filter_m & (Mg<-0.5)] 
Mg_sel_poly = bcg[filter_m & (Mg<-0.5)] 

model_f = LinearRegression().fit(mbol_sel.reshape(-1, 1), Mg_sel)
r_sq_f = model_f.score(mbol_sel.reshape(-1, 1), Mg_sel)
k_mbol_2 = model_f.coef_[0]
i_mbol_2 = model_f.intercept_


  

poly_bc_Mg = np.poly1d(np.polyfit(bcg_sel_poly, Mg_sel_poly, 2))
poly_bc = np.linspace(np.min(bcg_sel_poly),np.max(bcg_sel_poly),100)
poly_Mg = poly_bc_Mg(poly_bc)

if lin_reg_plots_on is True:
  
  fig_linReg_2, ax_linReg_2 = plt.subplots(2,2,sharey=True,sharex=False)
  
  
  ax_Mg_1 = ax_linReg_2[0,0]
  ax_bc_1 = ax_linReg_2[0,1]
  
  
  ax_Mg_1.scatter(mbol,Mg,s=0.5,color='black')
  ax_Mg_1.axline([coord_in_plot,coord_in_plot],slope=slope,color ='black',linewidth = 0.6,linestyle='-.')
  ax_Mg_1.axline([0+coord_in_plot,i_mbol_1 + coord_in_plot*k_mbol_1  ],slope=k_mbol_1,  color = 'red')
  ax_Mg_1.axhline(abs_rgb_agb_lower_bound,color = 'black',linewidth = 0.6,linestyle ='--')
  

  ax_bc_1.scatter(bcg,Mg,s=0.5,color='black')
  ax_bc_1.axvline(0, color = 'black',linewidth = 0.6,linestyle ='-.')
  ax_bc_1.axline([0+coord_in_plot,i_bcg_1 + coord_in_plot*k_bcg_1  ],slope=k_bcg_1,  color = 'red')
  ax_bc_1.axhline(abs_rgb_agb_lower_bound,color = 'black',linewidth = 0.6,linestyle ='--')
  
    
  
  
  ax_Mg_2 = ax_linReg_2[1,0]
  ax_bc_2 = ax_linReg_2[1,1]
  
  ax_Mg_2.scatter(mbol[filter_m],Mg[filter_m],s=0.5,color='black')
  ax_Mg_2.scatter(mbol[filter_outlier],Mg[filter_outlier],color='red',s=0.5,alpha=0.4)
  
  ax_Mg_2.axline([coord_in_plot,coord_in_plot],slope=slope,color ='black',linewidth = 0.6,linestyle='-.')
  ax_Mg_2.axline([0+coord_in_plot,i_mbol_2 + coord_in_plot*k_mbol_2  ],slope=k_mbol_2,  color = 'green')
  ax_Mg_2.axhline(abs_rgb_agb_lower_bound,color = 'black',linewidth = 0.6,linestyle ='--')
  
  ax_bc_2.scatter(bcg[filter_m],Mg[filter_m],s=0.5,color='black')
  ax_bc_2.scatter(bcg[filter_outlier],Mg[filter_outlier],color='red',s=0.5,alpha=0.4)
  
  ax_bc_2.plot(poly_bc,poly_Mg,  color = 'green')
  ax_bc_2.axhline(abs_rgb_agb_lower_bound,color = 'black',linewidth = 0.6,linestyle ='--')
  
  
  ax_bc_2.invert_yaxis()
    
  ax_Mg_1.set_ylabel('$M_{G}$')
  ax_Mg_2.set_ylabel('$M_{G}$')
  ax_Mg_2.set_xlabel('$M_{bol}$')
  ax_bc_2.set_xlabel('$BC_{G}$')
  
  if savefig is True:
    fig_linReg_2.savefig(fig_sub_path+fig_name+'_bcLinAprx_2.jpg')
  
  
#%% # Approximation of stars without BC and Mbol values

# Calculate approximated Mbol for stars without BC using lin fit -----
k = k_mbol_1
intercept = i_mbol_1
# lin_fit_color = t_fit
lin_fit_color = t_lost

# Mbol approximation derived from following formula
# Mg = k*Mbol + intercept 
mbol_aprx = (Mg_aprx - intercept)/k


# Calculate approximated Mbol uncertainty---------------------------------

# propagation of uncertainty ----------------------------
mbol_aprx_fit_unc = brightness_modul_unc / np.abs(k)
print('mbol_aprx_fit_unc: '+str(mbol_aprx_fit_unc))


# Statistical uncertainty evaluation from scattering ----
# Calculate approximated Mbol for stars with BC
Mbol_aprx_sel = (Mg_sel-intercept)/k
# Calculate deviation to derived Mbol by the gaiaBCg algorithm
mbol_aprx_deviation = np.abs(mbol_sel - Mbol_aprx_sel)

# Derive standard deviation and mean from mbol_aprx_deviation
# Use ( mean + 3 * std ) as uncertainty estimate
# 3* std includes 99.7% of deviations in normal distribution
mbol_aprx_fit_deviation_unc = np.mean(mbol_aprx_deviation) + 3*np.std(mbol_aprx_deviation)
print('mbol_aprx_fit_deviation_unc: '+str(mbol_aprx_fit_deviation_unc))

# Sum of uncertainties for Mbol ----------------------------
mbol_aprx_unc = mbol_aprx_fit_deviation_unc + mbol_aprx_fit_unc
# mbol_aprx_unc = np.sqrt(mbol_aprx_fit_deviation_unc**2 + mbol_aprx_fit_unc**2)
print('mbol_aprx_unc: '+str(mbol_aprx_unc))


# Calculate approximated BC and uncertainty ----------------
bcg_aprx     = mbol_aprx - Mg_aprx
# d(Mbol)**2 = np.sqrt(d(BC)**2 + d(MG)**2)
bcg_aprx_unc = np.sqrt(mbol_aprx_unc**2 - brightness_modul_unc**2)
print('bcg_aprx_unc: '+str(bcg_aprx_unc))



# Calc rgb agb mbol lower bound and uncertainty ------------
mbol_rgb_agb_lower_bound     = (abs_rgb_agb_lower_bound-intercept)/k
mbol_rgb_agb_lower_bound_unc = mbol_aprx_unc
print('mbol lower bound: '+str(mbol_rgb_agb_lower_bound))

# Print rgb agb mbol lower bound and uncertainty ------------
label_mg_lower_bound   = "$ M_{g}   = (%0.2f \pm %0.2f) $ mag" % (abs_rgb_agb_lower_bound,brightness_modul_unc)
label_mbol_fit         = "$ M_{bol} = (M_{g} - %0.2f) / %0.2f)$ mag" % (intercept,k)
label_mbol_fit_func    = "$ M_{bol} = (M_{g} - i) / k)$ mag"
label_mbol_lower_bound = "$ M_{bol} = (%0.2f \pm %0.2f) $ mag" % (mbol_rgb_agb_lower_bound,mbol_aprx_unc)

print()
print('AGB RGB Threshold:')
print(label_mg_lower_bound)
print(label_mbol_fit_func)
print(label_mbol_fit       )
print(label_mbol_lower_bound)



# Store approximated values in file ---------------------------
gaia_cluster_members['bcg_mag'][mask_bc_aprx] = bcg_aprx  
gaia_cluster_members['bcg_mag_error'][mask_bc_aprx] = bcg_aprx_unc
gaia_cluster_members['mbol_mag'][mask_bc_aprx] = mbol_aprx
gaia_cluster_members['mbol_mag_error'][mask_bc_aprx] = mbol_aprx_unc

save_gaia(gaia_cluster_members, file_name +'_bcLinAprx' ,file_sub_path=file_sub_path)

print_bcg_eval = False
if print_bcg_eval is True:
  print()
  bcg_unc = gaia_cluster_members['bcg_mag_error'][(gaia_cluster_members['bcg_aprx_flag']==False)]
  print(str(bcg_unc.mean()))
  print(str(bcg_unc.std()))
  print(str(bcg_unc.min()))
  print(str(bcg_unc.max()))
  bcg_eval_fig, bcg_eval_ax = plt.subplots(1,1)
  bcg_eval_ax.hist(bcg_unc)
  
  # Note: bcg is not normaly distributed -> mean & std meaningless!
  print()
  bcg_eval = gaia_cluster_members['bcg_mag'][(gaia_cluster_members['bcg_aprx_flag']==False)]
  print(str(bcg_eval.mean()))
  print(str(bcg_eval.std()))
  print(str(bcg_eval.min()))
  print(str(bcg_eval.max()))
  bcg_eval_fig, bcg_eval_ax = plt.subplots(1,1)
  bcg_eval_ax.hist(bcg_eval)

#%% BC CMD plot


fig_bc, axs_bc = plt.subplots(1,3,
                              figsize=(fig_width,fig_height*0.6),
                              gridspec_kw={'width_ratios':[1,1,1]},
                              sharey=True
                              )
# Linear Fit plot with (x,y) = (M_bol, Mg)
ax_Mg       = axs_bc[0]
# CMD without BC (y = Mg)
ax_cmd      = axs_bc[1]
# CMD with BC (y = M_bol)
ax_cmd_bol  = axs_bc[2]

# ax_Mg.set_aspect('equal', adjustable='box')
# ax_cmd.set_aspect('equal', adjustable='box')
# ax_cmd_bol.set_aspect('equal', adjustable='box')

# ax_Mg.set_aspect('equal')
# ax_cmd.set_aspect('equal')
# ax_cmd_bol.set_aspect('equal')

#-------------------------- Data visualization preparation --------------------

alpha = 1

# BC with gaiaBCg.py (_sel)
# Note: color pm_selection (t04='purple', dot)
bp_rp_sel       = ((gaia_cluster_members['bp_rp'])[mask_bc])[filter_m]

markersize_sel = 2.
marker_sel = '.'
color_sel = t02
# color_sel = t_fit

# BC with gaiaBCg.py (_sel) & linear aproximation
# Note: color pm_selection (t04='purple', dot)
m_bol_full_bc   = gaia_cluster_members['mbol_mag']
M_g_full_bc     = gaia_cluster_members['abs_phot_g_mean_mag']
bp_rp_full_bc   = gaia_cluster_members['bp_rp']

markersize_full_bc = 0.5
marker_full_bc = '.'
color_full_bc = t02
# color_full_bc = t_fit

# BC with gaiaBCg.py, but not nearby linear fit (_sel_outlier)
# Note: color = (t04, t_lost, t_fit), marker='+'
bp_rp_outlier   = ((gaia_cluster_members['bp_rp'])[mask_bc])[filter_outlier]

markersize_sel_outlier = 20
marker_sel_outlier = '+'
color_sel_outlier = t_lost
# color_sel_outlier = t02

# BC linear fit approximation (_aprx)
# Note: color t04, dot with t_fit ring
bp_rp_aprx      = (gaia_cluster_members['bp_rp'])[mask_bc_aprx]

markersize_aprx = 2.
marker_aprx = '.'
color_aprx = lin_fit_color
color_aprx = t01
# facecolor_aprx = 'none'
facecolor_aprx = color_aprx






# ----------------------- plot data -------------------------------------------
# Linear Fit plot with (x,y) = (M_bol, Mg) -------------------

x1, x2, y1, y2 = 0.5, 1.2, 0.5, 1.3
ratio = np.abs(x1-x2) / np.abs(y1-y2)

height = 0.40
width  = height*ratio
offset = 0.06

# inset axes....

axins = ax_Mg.inset_axes([1-(offset*ratio/2+width), 1-offset/2-height, width, height])
# subregion of the original image
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.invert_yaxis()

# x1, x2, y1, y2 = 0.0, 1.0, 0.0, 2.0
# ratio = np.abs(x1-x2) / np.abs(y1-y2)
# height = 0.30
# offset = 0.06
# # inset axes....
# width = height*ratio
# axins1 = ax_Mg.inset_axes([offset, 0.05+offset*ratio, width, height])
# # subregion of the original image
# axins1.set_xlim(x1, x2)
# axins1.set_ylim(y1, y2)
# axins1.invert_yaxis()


for i, ax_el in enumerate([ax_Mg,axins]):

  # # BC with gaiaBCg.py, but not nearby linear fit (_sel_outlier)
  # ax_el.scatter(mbol_aprx,Mg_aprx,
  #               s=markersize_sel_outlier,
  #               color=color_sel_outlier,
  #               alpha= alpha,
  #               marker=marker_sel_outlier
  #               )
  
  # BC with gaiaBCg.py (_sel)
  ax_el.scatter(mbol_sel, Mg_sel,
                color=color_sel,
                alpha=alpha,
                s=markersize_sel,
                # marker = marker_sel
                )
  
  # # full BC with gaiaBCg.py (_sel) & linear aproximation
  # ax_el.scatter(m_bol_full_bc, M_g_full_bc,
  #               color=color_full_bc,
  #               alpha=alpha,
  #               s=markersize_full_bc,
  #               marker = marker_full_bc
  #               )
  
  # BC with gaiaBCg.py, but not nearby linear fit (_sel_outlier)
  ax_el.scatter(mbol[filter_outlier],Mg[filter_outlier],
                color=color_sel_outlier,
                s=markersize_sel_outlier,
                alpha=alpha,
                marker=marker_sel_outlier
                )
  
  
  # Chosen linear fit for BC approximation
  ax_el.axline([0+coord_in_plot,intercept  +coord_in_plot*k  ],slope=k,  
               color = lin_fit_color 
               )
  
  
  # Chosen horizontal border for _sel_outlier
  x0,x1 = ax_el.get_xlim()
  x_min = M_g_full_bc.min()*20
  x_max = (Mg_outlier_limit - intercept) / k
  
  ax_el.hlines(
               Mg_outlier_limit,
               xmin = x_min, xmax=x_max ,  
               color = color_sel_outlier ,
               linewidth = 0.5,
               linestyle = '--'
               )
  ax_el.set_xlim(x0,x1)
  
  # Chosen vertical border for _sel_outlier
  y0,y1=ax_el.get_ylim()
  y_min = m_bol_full_bc.max()*20
  y_max = m_bol_outlier_limit* k + intercept
  
  ax_el.vlines(m_bol_outlier_limit,
               ymin = y_min, ymax=y_max ,  
               color = color_sel_outlier,
               linewidth = 0.5,
               linestyle = '--'
               )
  ax_el.set_ylim(y0,y1)
  
  ax_el.axhline(abs_rgb_agb_lower_bound,color = 'black',linewidth = 0.6,linestyle ='--')




# CMD without BC (y = Mg)--------------------------------


# BC linear fit approximation (_aprx)
ax_cmd.scatter(bp_rp_aprx, Mg_aprx,
               s = markersize_aprx ,
               # marker = marker_aprx,
               # facecolor = facecolor_aprx,
               color = color_aprx 
               )

# BC with gaiaBCg.py (_sel)
ax_cmd.scatter(bp_rp_sel,Mg_sel,
                color=color_sel,
                alpha=alpha,
                s=markersize_sel,
                # marker = marker_sel
                )

# # full BC with gaiaBCg.py (_sel) & linear aproximation
# ax_cmd.scatter(bp_rp_full_bc, M_g_full_bc,
#               color=color_full_bc,
#               alpha=alpha,
#               s=markersize_full_bc,
#               marker = marker_full_bc
#               )


# BC with gaiaBCg.py, but not nearby linear fit (_sel_outlier)
ax_cmd.scatter(bp_rp_outlier,Mg[filter_outlier],
               s=markersize_sel_outlier,
               color=color_sel_outlier,
               alpha= alpha,
               marker=marker_sel_outlier
               )
ax_cmd.axhline(abs_rgb_agb_lower_bound,
               color = 'black',
               linewidth = 0.6,
               linestyle ='--'
               )

# CMD with BC (y = M_bol) ----------------------------------

# BC linear fit approximation (_aprx)
ax_cmd_bol.scatter(bp_rp_aprx,mbol_aprx,
                   s = markersize_aprx ,
                   # marker = marker_aprx,
                   # facecolor = facecolor_aprx,
                   color = color_aprx                    
                   )

# BC with gaiaBCg.py (_sel)
ax_cmd_bol.scatter(bp_rp_sel,mbol_sel,
                    color=color_sel,
                    alpha=alpha,
                    s=markersize_sel,
                    # marker = marker_sel
                    )

# # full BC with gaiaBCg.py (_sel) & linear aproximation
# ax_cmd_bol.scatter(bp_rp_full_bc, m_bol_full_bc,
#               color=color_full_bc,
#               alpha=alpha,
#               s=markersize_full_bc,
#               marker = marker_full_bc
#               )


# BC with gaiaBCg.py, but not nearby linear fit (_sel_outlier)
ax_cmd_bol.scatter(bp_rp_outlier,mbol[filter_outlier],
                   s=markersize_sel_outlier,
                   color=color_sel_outlier,
                   alpha= alpha,
                   marker=marker_sel_outlier
                   )

ax_cmd_bol.axhline(abs_rgb_agb_lower_bound,
               color = 'black',
               linewidth = 0.6,
               linestyle ='--'
               )

ax_cmd_bol.axhline(mbol_rgb_agb_lower_bound,
               color = 'black',
               linewidth = 0.6,
               linestyle ='-.'
               )

# Ax labels -------------------------------------------------

ax_Mg.set_xlabel('$M_{bol}$')
ax_Mg.set_ylabel('$M_{G}$')

# ax_bc.invert_yaxis()
ax_bc.set_xlabel('$BC_{G}$')
ax_bc.set_ylabel('$M_{G}$')

# ax_cmd.invert_yaxis()
ax_cmd.set_xlabel('bp-rp')
ax_cmd.set_ylabel('$M_{G}$')

ax_cmd_bol.invert_yaxis()
ax_cmd_bol.set_xlabel('bp-rp')
ax_cmd_bol.set_ylabel('$M_{bol}$')

ax_cmd_bol.set_ylim(ymin=2.2)


# fig_test,ax_test = plt.subplots(1,1)
# ax_test.scatter(bc_g_mean,bc_g_unc,s=0.6)
# ax_test.scatter(m_bol,bc_g_unc,s=0.6)


  
if savefig is True:
  fig_bc.savefig(fig_sub_path+fig_name+'_bcLinAprx_cmd.jpg')
  
  
