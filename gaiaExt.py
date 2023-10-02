#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 21:13:20 2023

@author: Thomas Suntinger

Import recommendation:

  
Input:
  - load: data set of a cluster as csv file (e.g. from gaiaClusterMemberSelection.py)
  
Dependencies:
  - numpy
  - astropy
  - GAIA_load.py
  - GAIA_save.py
  
Dependencies for visualization:
  - matplotlib
  - plotStyle.py
  
Functionality:
  - Script loads a cluster member data set and evaluates absolute magnitude of
    stars, using one of multiple interstellar extinction models. 
  - Alternatively, one can also use a convolution of extinction models
    to derive a distance modulus.
  - The derived distance modulus of the interstellar extinction models are 
    compared in a plot. 
  - Description of the different extinction models can be found in adjoined
    'Structure:*'. 
  
Structure ( & parameters):
  - Init 
  - Set default plot style
  - Set default plot colors
  - Load /Store data settings
      (distance & error of selected cluster)
  - Store figures settings 
  - Toggle plots on/off
  - Interstellar Extinction 
      - Distance modulus without extinction correction 
      - Data preparation
  - 3.0) M5 distance modulus determination using RR-Lyr stars 
  - Sort out in future (comments, to be deleted if not needed!)
  - 3.a) Derivation, using a reference value for abs.magn of RCs (Chan&Bovy2022)
      (search box for horizontal branch, reference value for abs.magn of RCs)
  - 3.b) Derivation using GALEX
  - Amores&Lepine(2005) - Model A
  - Amores&Lepine(2005) - Model A2
  - Bayestar(2019)
  - Chen et al.(2019) - check out for converting to A_G!!!
  - 3.c) Derivation using GAIA DR3 extinction A_g
  - 3.c.2) Calculate abs.magn for every star individually
  - 3.c.2) Calculate abs.magn for every star individually plots
  - 3.d) Include interstellar extinction from starhors2021 catalouge in comparison
  - 3.d) Include interstellar extinction from starhors2021 catalouge in comparison plots
  - sort switcharoo in future (to be deleted if implementation does not help anyone)
  - Interstellar Extinction comparison
      (selected extinction model)
  - Interstellar Extinction comparison plot
  - CMD plot with Absolute Magn
    
Output:
  - store: csv file cluster member selection with absolute brightness & adjoined error 
  - store: csv file cluster member selection with Starhorse extinction & distance modulus
              (Ag,Ag error std & 95% Confidence level, mu, mu error std % 95% CL)
  - print: Distance modulus without extinction
  - print: Value and name of chosen distance modulus
  
Output plots:
  - Interstellar extinction comparison plot
  - CMD with absolute magnitudes after extinction correction
  - 3D Spatial interstellar extinction distribution of Gaia and starhors data sets
  
Improvements:
  - Clean unnecessary comments and lines of code
  - Write documentation and improve in code documentation.

  - Bonus: script set to globular cluster M5, changing cluster involves changing
           nearly the whole script. One might generate a copy and adjust it for
           another cluster
       - Make script executable for different clusters  
           - e.g. use solely starhors cataloge to derive extinction
           - Wait for DR4/DR5 and redo evaluation with Gaia extinction 
       - Rewritte 'Structure:sort switcharoo in future'

"""

#%% Init

import numpy as np
import matplotlib.pyplot as plt

# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn

#%% ### Set default plot style------------------------------------
#https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure

# Note: fig_width is fixed in plotStyle.py for MasterThesis specifications
#       fig_height can and should be adjusted to special needs of individual plots
#       fig_height of 5inch is approximatelly half of height on A4 MasterThesis      
fig_height = 5
import plotStyle
[ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height=fig_height)



#%% ### Set default plot colors --------------------------------------

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


#%% ### Load /Store data settings ---------------------------------------

folder_name = 'data/clusterMemberSelection/'

## -------------- GAIA data ---------------------------------------------------
# data_set_selected = 'Gaia_M3_r-18arcmin.csv'
data_set_selected = 'Gaia_M5_r-9arcmin_zpc_pm_clusterMemberSelection_distEval'
# data_set_selected = 'Gaia_M67_r-25arcmin.csv'

file_name = data_set_selected

## ------------ load data into variables --------------------------------------

from GAIA_load import load_query as load_query

#---------Load GAIA query --------------------------------

file_data = load_query(folder_name + data_set_selected)

gaia_cluster_members = file_data
label_gaia_cluster_members = 'cluster members'
label_gaia_cluster_members_n = 'cluster members ('+str(gaia_cluster_members.size)+' stars)'

#--------- Global distance & error --------------------------------

# Baumgardt et al. 2018 mean
d       = 7.57
d_error = 0.13


#---------Load Starhorse GAIA query --------------------------------

# Needs to be queried from the website:
# https://gaia.aip.de/query/
starhorse_set_selected = 'sth_query_M5_radius_32arcmin.csv'

starhorse_data = load_query('data/queryData/' + starhorse_set_selected)
label_starhorse_data = 'StarHorse query'
label_starhorse_n = label_starhorse_data+' ('+str(starhorse_data.size)+' stars)'



#---------Store GAIA data --------------------------------

file_sub_path = folder_name
from GAIA_save import save_gaia as save_gaia
save_data = True


#%% ### Store figures settings ------------------------------------------

savefig = True
fig_sub_path = 'gfx/clus_sel/'


# if savefig is True:
#   cmd_fig.savefig(file_path+file_name)

#%% ### Toggle plots on/off
int_stel_ext_plot_on = True
spat_ext_plot = False
CMD_abs_plot_on = False


#%% ### Interstellar Extinction ------------------------------------

# Uncertainty method toggle
# conservative_unc = True
conservative_unc = False

# Starhorse catalog uncertainty
# If confidence level as uncertainty is True:
#   - int.st.ext.plot will include cl as unc
#   - abs.mag unc will be set to cl as unc
cl_as_unc = False


# ---------------- Distance modulus without extinction correction -------

# distance modulus using only parallax, not corrected for extinction
# dist_mod_DR3_dist = - 5 * np.log10(1000/parallax) + 5
dist_mod_DR3_dist = - 5 * np.log10(1000 * d) + 5
dist_mod_DR3_dist = - dist_mod_DR3_dist

# Uncertainty in parallax evaluation
# dist_mod_DR3_dist_unc = parallax_error

label_dist_mod_DR3_dist = 'distance modulus'

mu_min = - 5 * np.log10(1000 * (d+d_error) ) + 5
mu_max = - 5 * np.log10(1000 * (d-d_error) ) + 5

mu_min_error = dist_mod_DR3_dist - mu_min 
mu_max_error = mu_max - dist_mod_DR3_dist 


# dist_mod_DR3_dist_unc = [mu_min_error,mu_max_error]
dist_mod_DR3_dist_unc = np.mean([mu_min_error,mu_max_error])

unc = np.full(gaia_cluster_members.size, dist_mod_DR3_dist)
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'true_distance_modulus', unc , dtypes='<f8')

unc = np.full(gaia_cluster_members.size, dist_mod_DR3_dist_unc )
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'true_distance_modulus_error', unc, dtypes='<f8')

# ---------------- Data preparation ------------------------------

# gaia_cluster_members = gaia_dr3_array_special
gaia_cluster_members = gaia_cluster_members
label_data_DR3 = label_gaia_cluster_members

ext_array_length = 6
int_stel_ext_array = np.zeros(ext_array_length ,dtype=float)
int_stel_ext_array_unc = np.zeros(ext_array_length ,dtype=float)

dist_mod_array     = np.copy(int_stel_ext_array)
dist_mod_unc_array = np.copy(int_stel_ext_array)
dist_mod_model     = ['','','','','','','',]
dist_mod_model_note= dist_mod_model.copy()

#%% # 3.0) M5 distance modulus determination using RR-Lyr stars 
#Note: first attempt, more research necessary!

# source: https://ui.adsabs.harvard.edu/abs/2011MNRAS.416.1056C/abstract
# distance_modulus = 14.44
# distance_modulus = 14.50


#%% # ------------------- Sort out in future --------------------------------------
# # Absolute Magnitude using distance from GAIA EDR3 parallaxes
# # m-M=5*log(distance in pc)-5 (+ A_lambda)

# M = m - 5 * np.log10(distance) + 5 
# M = m - 5 * np.log10(1 / parallax in arcsec) + 5 
# For parallax in mas
# M = m - 5 * np.log10(1000 / parallax in marcsec) + 5


# currently gaia_dr3_array = file_data
# phot_g_mean_mag   = gaia_dr3_array['phot_g_mean_mag']
# par               = gaia_dr3_array['parallax']
# bp_rp             = gaia_dr3_array['bp_rp']
                            
# phot_g_mean_mag   = gaia_dr3_array_special['phot_g_mean_mag']
# par               = gaia_dr3_array_special['parallax']
# bp_rp             = gaia_dr3_array_special['bp_rp']

# phot_g_mean_mag   = gaia_dr3_array_pm_selection['phot_g_mean_mag']
# par               = gaia_dr3_array_pm_selection['parallax']
# bp_rp             = gaia_dr3_array_pm_selection['bp_rp']

# abs_phot_g_mean_mag = phot_g_mean_mag   - 5 * np.log10(1000/parallax) + 5

# # distance_modulus = m-M = 5*log(distance in pc)-5
# # M = m - distance_modulus
# # source: https://ui.adsabs.harvard.edu/abs/2011MNRAS.416.1056C/abstract
# distance_modulus = 14.44
# abs_phot_g_mean_mag_RR =  phot_g_mean_mag - distance_modulus

# offset = abs_phot_g_mean_mag - abs_phot_g_mean_mag_RR

# ------------------- End of Sort out in future -------------------------------



#%% # 3.a) Derivation, using a reference value for abs.magn of RCs (Chan&Bovy2022)

# Ref. value from source
m_abs_g_ref     = 0.435
m_abs_g_ref_unc = 0.004


# ---------- Filter ap.mag from RCs in GAIA data sample ------------
# Note: Select rectangle arround horizontal branch, get stars and calc mean,std...

# For M5
m_ap_max = 15.2
m_ap_min = 14.85
b_r_min = -0.1
b_r_max = 0.93

# For M5
# m_ap_max = 15.28
# m_ap_min = 14.938
# b_r_min = 0.26
# b_r_max = 0.85

# Filter area of high probability RC population
 
m_ap_mask = (np.less(gaia_cluster_members['phot_g_mean_mag'], m_ap_max) & \
              np.greater(gaia_cluster_members['phot_g_mean_mag'], m_ap_min) )

b_r_mask = (np.less(gaia_cluster_members['bp_rp'], b_r_max) & \
              np.greater(gaia_cluster_members['bp_rp'], b_r_min) )

# Calculate mean apparent magn for data area
m_ap_mean_rc_3a = np.mean(gaia_cluster_members['phot_g_mean_mag'][(m_ap_mask & b_r_mask)])
# Use standarddeviation as uncertainty
m_ap_mean_rc_3a_unc = np.std(gaia_cluster_members['phot_g_mean_mag'][(m_ap_mask & b_r_mask)])

# Calc distance modulus and use highest uncertainty as new uncertainty
dist_mod_rc_3a      = m_ap_mean_rc_3a - m_abs_g_ref

# dist_mod_rc_3a_unc  = np.max([m_ap_mean_rc_3a_unc , m_abs_g_ref_unc] )
# dist_mod_rc_3a_unc  = m_ap_mean_rc_3a_unc + m_abs_g_ref_unc
# dist_mod_rc_3a_unc  = np.sqrt(m_ap_mean_rc_3a_unc**2 + m_abs_g_ref_unc**2)

if conservative_unc is True:
  dist_mod_rc_3a_unc  = m_ap_mean_rc_3a_unc + m_abs_g_ref_unc
  a_g_rc_3a_unc =  dist_mod_rc_3a_unc + dist_mod_DR3_dist_unc  
else:
  dist_mod_rc_3a_unc  = np.sqrt(m_ap_mean_rc_3a_unc**2 + m_abs_g_ref_unc**2)
  a_g_rc_3a_unc =  np.sqrt(dist_mod_DR3_dist_unc**2 + dist_mod_rc_3a_unc**2)

a_g_rc_3a =  dist_mod_DR3_dist - dist_mod_rc_3a



i = 5
int_stel_ext_array[i] = a_g_rc_3a
int_stel_ext_array_unc[i] = a_g_rc_3a_unc

dist_mod_array[i]     = dist_mod_rc_3a
dist_mod_unc_array[i] = dist_mod_rc_3a_unc
# dist_mod_model[i]     = 'HB fitting with abs.M_RC by Chan \& Bovy2020'
dist_mod_model[i]     = 'Horizontal branch fitting'
dist_mod_model_note[i]= '$A_{V}$ Error bar is stand.dev. of stars in HB'

# ----------------------------------------------------------------------



#%% # 3.b) Derivation using GALEX


#------------------------------------------------------------------
#convert coordinates into diverent equinox (=epoch)
# https://docs.astropy.org/en/stable/coordinates/transforming.html
#
# Note: No difference between 2016 and 2000 for galactic coordinates

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic, FK5, ICRS
from astropy.time import Time

J2016_coord_data = SkyCoord(gaia_cluster_members['ra'], gaia_cluster_members['dec'],
                        unit=(u.degree, u.degree), frame = ICRS)

J2016_coord_data.transform_to(FK5(equinox=Time('J2016')))
J2000_coord_data = J2016_coord_data.copy()
J2000_coord_data.transform_to(FK5(equinox=Time('J2000')))

# ICRS is close to mas to FK5 Epoch J2000
J2000_coord_data.galactic
J2016_coord_data.galactic
#------------------------------------------------------------------

ra_M5 = 229.638417
dec_M5 = 2.081028
J2016_coord_M5 = SkyCoord(ra_M5, dec_M5, unit=(u.degree, u.degree), frame = ICRS)

J2016_coord_M5.transform_to(FK5(equinox=Time('J2016')))
J2000_coord_M5 = J2016_coord_M5.copy()
J2000_coord_M5.transform_to(FK5(equinox=Time('J2000')))
# ICRS is close to mas to FK5 Epoch J2000
J2000_coord_M5.galactic
J2016_coord_M5.galactic

# Use this coordinates for GALEX querries:
l_M5 = 3.85872974  #longitude
b_M5 = 46.79636532 #latitude

# 7.366075209090368
# 7.008540961250695




#%% # Amores&Lepine(2005) - Model A
# Valid for distances:
  
# distance = 7.341343903517526 kpc
# distance = 7.366075209090368 kpc
# distance = 7.008540961250695 kpc
# a_V_GALEX = 0.1198

# distance = 7.57 kPc
a_V_GALEX = 0.120
a_V_GALEX_std = 0

# Figure out, how to convert or how to use E(B-V) to adjust for G-band
a_g_GALEX     = a_V_GALEX
a_g_GALEX_unc = 0

dist_mod_GALEX_3b_Amores      =  a_g_GALEX + dist_mod_DR3_dist

if conservative_unc is True:  
  dist_mod_GALEX_3b_Amores_unc  = a_g_GALEX_unc + dist_mod_DR3_dist_unc
else:
  dist_mod_GALEX_3b_Amores_unc  = np.sqrt(a_g_GALEX_unc**2 + dist_mod_DR3_dist_unc**2)



i = 0
int_stel_ext_array[i] = a_g_GALEX
int_stel_ext_array_unc[i] = a_g_GALEX_unc

dist_mod_array[i]     = dist_mod_GALEX_3b_Amores
dist_mod_unc_array[i] = dist_mod_GALEX_3b_Amores_unc
dist_mod_model[i]     = 'Amores \& Lepine(2005) - A'
dist_mod_model_note[i]= '$A_{V}$ A_V not converted'


#%% # Amores&Lepine(2005) - Model A2
# Valid for distances:
# distance = 7.341343903517526 kpc
# distance = 7.366075209090368 kpc
# distance = 7.008540961250695 kpc

# a_V_GALEX = 0.121270 
# a_V_GALEX_std = 0

# distance = 7.57 kPc
a_V_GALEX = 0.121
a_V_GALEX_std = 0

# Figure out, how to convert or how to use E(B-V) to adjust for G-band
a_g_GALEX     = a_V_GALEX
a_g_GALEX_unc = a_V_GALEX_std

if conservative_unc is True:
  dist_mod_GALEX_3b_Amores_A2_unc  = a_g_GALEX_unc+ dist_mod_DR3_dist_unc
else:
  dist_mod_GALEX_3b_Amores_A2_unc  = np.sqrt(a_g_GALEX_unc**2 + dist_mod_DR3_dist_unc**2)

# Calc distance modulus and use highest uncertainty as new uncertainty
dist_mod_GALEX_3b_Amores_A2      =  a_g_GALEX + dist_mod_DR3_dist



i = 1
int_stel_ext_array[i] = a_g_GALEX
int_stel_ext_array[i] = a_g_GALEX_unc

dist_mod_array[i]     = dist_mod_GALEX_3b_Amores_A2
dist_mod_unc_array[i] = dist_mod_GALEX_3b_Amores_A2_unc
dist_mod_model[i]     = 'Amores \& Lepine(2005) - A2'
dist_mod_model_note[i]= '$A_{V}$ not converted'



#%% # Bayestar(2019)
# Valid for distances:
# distance = 7.341343903517526 kpc
# distance = 7.366075209090368 kpc
# distance = 7.008540961250695 kpc
# distance = 7.57 kpc

# To convert Bayestar2019 (=a_V_GALEX) to extinction, multiply by the coefficients given in Table 1 of Green et al. (2019). 
# source: https://dustmaps.readthedocs.io/en/latest/maps.html#bayestar 
# source: http://argonaut.skymaps.info/usage#units 
a_V_GALEX = 0.1100 
a_V_GALEX_std = 0.0036 

# Figure out, how to convert or how to use E(B-V) to adjust for G-band
a_g_GALEX     = a_V_GALEX
a_g_GALEX_unc = a_V_GALEX_std

# Calc distance modulus and use highest uncertainty as new uncertainty
dist_mod_GALEX_3b_Bayestar      =  a_g_GALEX + dist_mod_DR3_dist

if conservative_unc is True:
  dist_mod_GALEX_3b_Bayestar_unc  = a_g_GALEX_unc + dist_mod_DR3_dist_unc
else:
  dist_mod_GALEX_3b_Bayestar_unc  = np.sqrt(a_g_GALEX_unc**2 + dist_mod_DR3_dist_unc**2)

i = 2
int_stel_ext_array[i] = a_g_GALEX
int_stel_ext_array_unc[i] = a_g_GALEX_unc

dist_mod_array[i]     = dist_mod_GALEX_3b_Bayestar
dist_mod_unc_array[i] = dist_mod_GALEX_3b_Bayestar_unc
dist_mod_model[i]     = 'Bayestar19'
dist_mod_model_note[i]= '$A_{V}$ not converted'










#%% # Chen et al.(2019) - check out for converting to A_G!!!
# Extinction in A_G using the relation provided by the authors: A_G = E(G-K) + 1.987*E(H-K) (Chen et al. 2019)
# nan means no valid values (coordinates outside the coverage of the map).
# Link for the paper: https://academic.oup.com/mnras/article/483/4/4277/5237718
# The values with subscript -- in -- means the input data and with - out -- as obtained in the map

# l_in, l_out, b_in, b_out, dist. in, dist. out, A_G, A_G_error Max. distance available
# ------------------------------------------------------------------
# 3.8587 nan 46.7964 nan 7.341 nan nan nan nan 

# Chen does not include M5 position extinction values







#%% # 3.c) Derivation using GAIA DR3 extinction A_g

# Only for test cases, disable if done! 
# data_DR3 = gaia_dr3_array

a_g_DR3_mask = np.logical_not(np.isnan(gaia_cluster_members['ag_gspphot'] ))
a_g_DR3_data = (gaia_cluster_members['ag_gspphot'][(a_g_DR3_mask)])


# 3.c.1) Calculate mean interstellar extinction for data area
a_g_DR3 = np.mean(gaia_cluster_members['ag_gspphot'][(a_g_DR3_mask)])
# Use standarddeviation as uncertainty
a_g_DR3_unc = np.std(gaia_cluster_members['ag_gspphot'][(a_g_DR3_mask)])

# Calc distance modulus and use highest uncertainty as new uncertainty
dist_mod_DR3_3c1      =  a_g_DR3 + dist_mod_DR3_dist

if conservative_unc is True:
  dist_mod_DR3_3c1_unc  = a_g_DR3_unc + dist_mod_DR3_dist_unc
else:
  dist_mod_DR3_3c1_unc  = np.sqrt(a_g_DR3_unc**2 + dist_mod_DR3_dist_unc**2)

i = 3
int_stel_ext_array[i] = a_g_DR3
int_stel_ext_array_unc[i] = a_g_DR3_unc

dist_mod_array[i]     = dist_mod_DR3_3c1
dist_mod_unc_array[i] = dist_mod_DR3_3c1_unc
dist_mod_model[i]     = 'GAIA DR3 extinction'
dist_mod_model_note[i]= '$A_{G}$ Error bar is stand.dev.'




#%% # 3.c.2) Calculate abs.magn for every star individually
# Note: Appropriate, if parallax is also mean parallax? If time, test it ;)


# distance modulus of every star individually with constant parallax
dist_mod_DR3_3c2 = a_g_DR3_data - dist_mod_DR3_dist
 
# spat_ext_plot =True
if spat_ext_plot is True:
 
 hist_gaia_ext_fig, hist_gaia_ext_ax = plt.subplots()

 N, bins, patches = \
    hist_gaia_ext_ax.hist( \
                          a_g_DR3_data,\
                            bins = 100,\
                             color=t04,\
                                label=dist_mod_model[4])    


 hist_gaia_ext_ax.set_xlabel('ag_gspphot in mag')
 hist_gaia_ext_ax.set_ylabel('count of stars')
 hist_gaia_ext_ax.legend(loc='upper left')



N_bins = 100
x = gaia_cluster_members['ra'][(a_g_DR3_mask)]
y = gaia_cluster_members['dec'][(a_g_DR3_mask)]
z = dist_mod_DR3_3c2
z = a_g_DR3_data

# x_delta = (x.max()-x.min()) / (N_bins)
# y_delta = (y.max()-y.min()) / (N_bins)
# xx, yy = np.mgrid[x.min()-x_delta:x.max()+x_delta:N_bins*1j, y.min()-y_delta:y.max()+y_delta:N_bins*1j]
xx, yy = np.mgrid[x.min():x.max():N_bins*1j, y.min():y.max():N_bins*1j]
zz = np.zeros_like(xx)
# zz = zz[:-1,:-1]

x_delta = (xx[1,1] - xx[0,1])/2
y_delta = (yy[1,1] - yy[1,0])/2

for i, el in enumerate(xx[:,1]):
  for j, el in enumerate(yy[1,:]):    
    

    mask_x_box = ( (xx[i  ,1]-x_delta) <= x) & ( (xx[i,1] +x_delta) >  x) 
    mask_y_box = ( (yy[1  ,j]-y_delta) <= y) & ( (yy[1,j] +y_delta) >  y)
        
    mask_in_box = mask_x_box & mask_y_box
    zz[i,j] = np.average(z[mask_in_box])
    
# zz[(np.isnan(zz))] = 0.0 
zz -= dist_mod_DR3_dist


# ------------ Interpolation Test of extinction data---------------------
# Note: Not fruitfull! 

# from scipy import interpolate as inpo

# Note: Gets messy with negative extinction values
# f = inpo.interp2d(xx[:,1], yy[1,:],zz, kind = 'cubic')
# zz_inpo = f(xx[:,1], yy[1,:])

# Note: Does only show same characteristic as in averaging boxes
# f = inpo.interp2d(x, y,z, kind = 'cubic')
# z_inpo = f(xx[:,1], yy[1,:])
# zz_inpo = z_inpo.reshape(xx.shape)
  
# ------------------------------------------------------------------------

#%% # 3.c.2) Calculate abs.magn for every star individually plots
if spat_ext_plot is True:
  fig_3d = plt.figure()
  ax_3d = fig_3d.add_subplot(projection='3d')
  surf = ax_3d.plot_surface(xx, \
                            yy, \
                            zz, \
                              rstride=1, cstride=1, linewidth=0, antialiased=False)
  ax_3d.set_xlabel('ra')
  ax_3d.set_ylabel('dec')
  ax_3d.set_zlabel('distance modulus [mag]')
  ax_3d.set_title(dist_mod_model[4])
  
  # plot two  
  colmesh_gaia_ext_fig,colmesh_gaia_ext_axe =    plt.subplots()
  im_mesh = \
    colmesh_gaia_ext_axe.pcolormesh(xx, yy, zz, \
                                          shading='gouraud', \
                                            cmap='winter_r',\
                                              zorder = 0)
  col_bar = colmesh_gaia_ext_fig.colorbar(im_mesh, ax = colmesh_gaia_ext_axe,\
                             # ticks = col_ticks,\
                             #   fraction = 0.046,pad = 0.04
                             label = 'distance modulus [mag]'
                               )
  # # colmesh_gaia_ext_axe.scatter(x[(z<=.10)],y[(z<=.1)], label = r'<0.10:'+str(y[(z<.10)].size),s=1, color='orange')
  # colmesh_gaia_ext_axe.scatter(x[(z>1.0)],y[(z>1.0)], label = r'>1.00:'+str(y[(z>1.0)].size),s=1, color='white')
  # mask = (z>0.25)&(z<=1.0)
  # colmesh_gaia_ext_axe.scatter(x[mask],y[mask], label = r'> 0.25: '+str((y[mask]).size),s=1,color = 'grey')
  # mask = (z>0.1)&(z<=0.25)
  # colmesh_gaia_ext_axe.scatter(x[mask],y[mask], label = r'> 0.10:'+str((y[mask]).size),s=1,color = 'black')
  colmesh_gaia_ext_axe.set_xlabel('ra')
  colmesh_gaia_ext_axe.set_ylabel('dec')
  colmesh_gaia_ext_axe.legend()
  colmesh_gaia_ext_axe.set_title(dist_mod_model[4])





#%%  3.d) Include interstellar extinction from starhors2021 catalouge in comparison
# https://data.aip.de/projects/starhorse2021.html

# sth_mask = np.where(starhorse_data['g_id']==gaia_cluster_members['source_id'],starhorse_data['g_id'],gaia_cluster_members['source_id'])
sth_mask = np.isin(starhorse_data['g_id'],gaia_cluster_members['source_id'])
sth_selected = starhorse_data[sth_mask]

dr3_sth_mask = np.isin(gaia_cluster_members['source_id'],sth_selected['g_id'])
dr3_sth_selected = gaia_cluster_members[dr3_sth_mask]



# No Nan values in data set
a_g_sth_mask = np.logical_not(np.isnan(sth_selected['ag50'] ))
a_g_sth_data = (sth_selected['ag50'][(a_g_sth_mask)])


# 3.d.1) Calculate mean interstellar extinction for data area
a_g_sth = np.mean(a_g_sth_data)
# Use standarddeviation as uncertainty
a_g_sth_unc_std = np.std(a_g_sth_data)

dist_mod_sth_3d1      =  a_g_sth + dist_mod_DR3_dist



from gaiaErrorEval import trgb_conv_int
a_g_sth, a_g_sth_unc_cl95, a_g_sth_margin_cl95 = \
        trgb_conv_int(a_g_sth, 
                      sigma = a_g_sth_unc_std, 
                      cl = 0.95
                      )

  

unc = np.full(gaia_cluster_members.size, a_g_sth)
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'ag_starhorse_mean', unc , dtypes='<f8')

unc = np.full(gaia_cluster_members.size, a_g_sth_unc_std)
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'ag_starhorse_mean_error_std', unc , dtypes='<f8')

unc = np.full(gaia_cluster_members.size, a_g_sth_unc_cl95)
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'ag_starhorse_mean_error_cl95', unc , dtypes='<f8')

if conservative_unc is True:
  dist_mod_sth_3d1_unc_std   = a_g_sth_unc_std      + dist_mod_DR3_dist_unc
  dist_mod_sth_3d1_unc_cl95  = a_g_sth_unc_cl95 + dist_mod_DR3_dist_unc
else:
  dist_mod_sth_3d1_unc_std   = np.sqrt(a_g_sth_unc_std**2 + dist_mod_DR3_dist_unc**2)
  dist_mod_sth_3d1_unc_cl95  = np.sqrt(a_g_sth_unc_cl95**2 + dist_mod_DR3_dist_unc**2)


if cl_as_unc is True:
  a_g_sth_unc  = a_g_sth_unc_cl95
  dist_mod_sth_3d1_unc = dist_mod_sth_3d1_unc_cl95
  print('Starhorse extinction uncertainty set to 95% confidence level')
else:
  a_g_sth_unc = a_g_sth_unc_std
  dist_mod_sth_3d1_unc = dist_mod_sth_3d1_unc_std
  print('Starhorse extinction uncertainty set to one standard deviation')

 
# Global interstellar extinction and distance correction for absolute magnitude
unc = np.full(gaia_cluster_members.size, dist_mod_sth_3d1)
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'distance_modulus_starhorse', unc , dtypes='<f8')

unc = np.full(gaia_cluster_members.size, dist_mod_sth_3d1_unc_std)
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'distance_modulus_starhorse_error_std', unc , dtypes='<f8')

unc = np.full(gaia_cluster_members.size, dist_mod_sth_3d1_unc_cl95)
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'distance_modulus_starhorse_error_cl95', unc , dtypes='<f8')



i = 4
int_stel_ext_array[i] = a_g_sth
int_stel_ext_array_unc[i] = a_g_sth_unc

dist_mod_array[i]     = dist_mod_sth_3d1
dist_mod_unc_array[i] = dist_mod_sth_3d1_unc
dist_mod_model[i]     = 'StarHorse extinction'
dist_mod_model_note[i]= '$A_{G}$ Error bar is stand.dev.'












#--------------------------------------------------------

# Source with same dimensions:
  # sth_selected           from starhorse catalouge
  # dr3_sth_selected       from gaia catalouge
#--------------------------------------------------------



# distance modulus of every star individually with constant parallax
dist_mod_sth_3d2 = a_g_sth_data - dist_mod_DR3_dist

if spat_ext_plot is True:
 
 hist_gaia_ext_fig, hist_gaia_ext_ax = plt.subplots()

 N, bins, patches = \
    hist_gaia_ext_ax.hist( \
                          a_g_sth_data,\
                            bins = 50,\
                             color=t07,\
                                label=dist_mod_model[5])    
 # N, bins, patches = \
 #         hist_gaia_ext_ax.hist( \
 #                               a_g_DR3_data,\
 #                                 bins = 100,\
 #                                  color=t04,\
 #                                     label=dist_mod_model[4])    


 
 hist_gaia_ext_ax.set_xlabel('ag_50 in mag')
 hist_gaia_ext_ax.set_ylabel('count of stars')
 hist_gaia_ext_ax.legend(loc='upper right')



N_bins = 100
x = dr3_sth_selected['ra']
y = dr3_sth_selected['dec']
z = dist_mod_sth_3d2
z = a_g_sth_data

xx, yy = np.mgrid[x.min():x.max():N_bins*1j, y.min():y.max():N_bins*1j]
zz = np.zeros_like(xx)
# zz = zz[:-1,:-1]

x_delta = (xx[1,1] - xx[0,1])/2
y_delta = (yy[1,1] - yy[1,0])/2

for i, el in enumerate(xx[:,1]):
  for j, el in enumerate(yy[1,:]):    
    mask_x_box = ( (xx[i  ,1]-x_delta) <= x) & ( (xx[i,1] +x_delta) >  x) 
    mask_y_box = ( (yy[1  ,j]-y_delta) <= y) & ( (yy[1,j] +y_delta) >  y)
        
    mask_in_box = mask_x_box & mask_y_box
    zz[i,j] = np.average(z[mask_in_box])
    
# zz[(np.isnan(zz))] = 0.0 
zz -= dist_mod_DR3_dist

# ------------ Interpolation Test of extinction data---------------------
# Note: Not fruitfull! 

# from scipy import interpolate as inpo

# Note: Gets messy with negative extinction values
# f = inpo.interp2d(xx[:,1], yy[1,:],zz, kind = 'cubic')
# zz_inpo = f(xx[:,1], yy[1,:])

# Note: Does only show same characteristic as in averaging boxes
# f = inpo.interp2d(x, y,z, kind = 'cubic')
# z_inpo = f(xx[:,1], yy[1,:])
# zz_inpo = z_inpo.reshape(xx.shape)
  
# ------------------------------------------------------------------------

#%%  3.d) Include interstellar extinction from starhors2021 catalouge in comparison plots

if spat_ext_plot is True:
  fig_3d = plt.figure()
  ax_3d = fig_3d.add_subplot(projection='3d')
  surf = ax_3d.plot_surface(xx, \
                            yy, \
                            zz, \
                              rstride=1, cstride=1, linewidth=0, antialiased=False,\
                                label = dist_mod_model[5])
  ax_3d.set_xlabel('ra')
  ax_3d.set_ylabel('dec')
  ax_3d.set_zlabel('distance modulus [mag]')
  ax_3d.set_title(dist_mod_model[5])
  
  # second plot
  colmesh_gaia_ext_fig,colmesh_gaia_ext_axe =    plt.subplots()
  im_mesh = \
    colmesh_gaia_ext_axe.pcolormesh(xx, yy, zz, \
                                          shading='gouraud', \
                                            # cmap=plt.cm.BuGn_r,\
                                            cmap=  'winter_r',\
                                              zorder = 0)
  col_bar = colmesh_gaia_ext_fig.colorbar(im_mesh, ax = colmesh_gaia_ext_axe,\
                             # ticks = col_ticks,\
                             #   fraction = 0.046,pad = 0.04
                             label = 'distance modulus [mag]'
                               )
  # # colmesh_gaia_ext_axe.scatter(x[(z<=.10)],y[(z<=.1)], label = r'<0.10:'+str(y[(z<.10)].size),s=1, color='orange')
  # colmesh_gaia_ext_axe.scatter(x[(z>1.0)],y[(z>1.0)], label = r'> 1.00:'+str(y[(z>1.0)].size),s=1, color='white')
  # mask = (z>0.25)&(z<=1.0)
  # colmesh_gaia_ext_axe.scatter(x[mask],y[mask], label = r'> 0.25: '+str((y[mask]).size),s=1,color = 'grey')
  # mask = (z>0.1)&(z<=0.25)
  # colmesh_gaia_ext_axe.scatter(x[mask],y[mask], label = r'> 0.10: '+str((y[mask]).size),s=1,color = 'black')
  colmesh_gaia_ext_axe.set_xlabel('ra')
  colmesh_gaia_ext_axe.set_ylabel('dec')
  # colmesh_gaia_ext_axe.legend()
  colmesh_gaia_ext_axe.set_title(dist_mod_model[5])


#%% # sort switcharoo in future
# if False:
#   i=0
#   int_stel_ext_xd           = int_stel_ext_array[i] 
#   dist_mod_xd                = dist_mod_array[i]
#   mod_name_xd                = dist_mod_model[i]
#   dist_mod_unc_xd            = dist_mod_unc_array[i]
  
#   i=5
#   int_stel_ext_xd5           = int_stel_ext_array[i] 
#   dist_mod_xd5                = dist_mod_array[i]
#   mod_name_xd5                = dist_mod_model[i]
#   dist_mod_unc_xd5            = dist_mod_unc_array[i]
  
#   i=5
#   int_stel_ext_array[i]   = int_stel_ext_xd           
#   dist_mod_array[i]       = dist_mod_xd             
#   dist_mod_model[i]       = mod_name_xd                
#   dist_mod_unc_array[i]   = dist_mod_unc_xd           
  
#   i=0
#   int_stel_ext_array[i]   = int_stel_ext_xd5           
#   dist_mod_array[i]       = dist_mod_xd5             
#   dist_mod_model[i]       = mod_name_xd5                
#   dist_mod_unc_array[i]   = dist_mod_unc_xd5           
  
#   color_list = [t02,t02,t02,t04,t07,t_lost]
# else:
#   color_list = [t_lost,t02,t02,t02,t04,t07]






#%% # Interstellar Extinction comparison -------------------------

color_list = [t02,t02,t02,t04,t07,t_lost]

# int_stel_ext_array[i] = a_g_DR3
# dist_mod_array[i]     = dist_mod_DR3_3c1
# dist_mod_unc_array[i] = dist_mod_DR3_3c1_unc
# dist_mod_model[i]     = 'GAIA DR3 extinction A_g'
# dist_mod_model_note[i]= ''

# Note: One aproach to evaluate int.st.ext. is averaging all models
#       HB-fitting excluded from average - method prone to strong variation of 
#       result due to search-window definition

# ------------------Int.st.ext. selection / evaluation ---------------

# # Mean extinction of models (1-5) -----------------------------
# int_stel_ext            = np.mean(int_stel_ext_array[1:]) 
# dist_mod                = np.mean(dist_mod_array[1:]     )

# # dist_mod_unc_max        = np.max(dist_mod_unc_array[1:])
# # dist_mod_unc_std        = np.std(dist_mod_unc_array[1:])
# # dist_mod_unc            = np.max([dist_mod_unc_max,dist_mod_unc_std])

# n = dist_mod_unc_array[1:].size
# dist_mod_unc        = np.sum(dist_mod_unc_array[1:]) / n

# int_stel_ext_unc = np.sum(int_stel_ext_array_unc[1:]) / n

# mod_name = '(mean of models 1 to 5)'
# # ------------------------------------------------------------

# Model i is taken as extinction -------------------------------
i = 4
int_stel_ext            = int_stel_ext_array[i] 
dist_mod                = dist_mod_array[i]
mod_name                = dist_mod_model[i]
dist_mod_unc            = dist_mod_unc_array[i]

if cl_as_unc is True: 
  unc = '95% CL error' 
else: 
  unc = 'standard deviation error'
# --------------------------------------------------------------

print('')
print('-------------Interstellar Extinction Evaluation ---------------')
print('')
print('Distance modulus without extinction: ')
print ("(m-M) = %0.3f (+/-) %0.3f mag" % (dist_mod_DR3_dist, dist_mod_DR3_dist_unc))

print('')
print('Chosen distance modulus: '+mod_name + ' with '+ unc)
print ("(m-M) = %0.3f (+/-) %0.3f mag" % (dist_mod, dist_mod_unc))
print('')


# Global interstellar extinction and distance correction for absolute magnitude
abs_mag = gaia_cluster_members['phot_g_mean_mag'] - dist_mod
abs_mag_unc = np.zeros_like(abs_mag)
abs_mag_unc[:] = dist_mod_unc

# Extend existing structured numpy nd array with fields
# Source: https://stackoverflow.com/questions/25427197/numpy-how-to-add-a-column-to-an-existing-structured-array
import numpy.lib.recfunctions as rfn
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'abs_phot_g_mean_mag', abs_mag, dtypes='<f8')
gaia_cluster_members = rfn.append_fields(gaia_cluster_members, 'abs_phot_g_mean_mag_error', abs_mag_unc , dtypes='<f8')





save_gaia(gaia_cluster_members, file_name +'_absG-band' ,file_sub_path=file_sub_path)





#%% # Interstellar Extinction comparison plot -------------------------

ext_x = np.arange(ext_array_length)

# int_stel_ext_plot_on = True
if int_stel_ext_plot_on                     is True:
  
  
  ext_fig, ext_ax = plt.subplots(1,1,figsize = (fig_width,fig_height*0.7),sharex=(True),sharey=False)
  
  if ready_to_publish is False:
    ext_fig.suptitle('Comparison of interstellar extinction models')
          
  
  ext_ax.set_ylabel('(m-M) [mag]')  
  ext_ax.set_xlabel('extinction model')
  ext_ax.set_xticks(ext_x)
  ext_ax.set_xlim(-0.5,5.5)  
  # ext_ax.set_ylim(top=15,bottom=11)  
  # ext_ax.set_xlabel('extinction model')
   
  # ext_ax.axhline(dist_mod+dist_mod_unc, color=t_fit,linewidth=0.6,linestyle='--')
  # ext_ax.axhline(dist_mod-dist_mod_unc, color=t_fit,linewidth=0.6,linestyle='--')
  
  ext_ax2 = ext_ax.twinx()
  
  # Distance modulus 
  # renamed to $ (5 log(d)-5)$
  label_dist_mod_DR3_dist = r'$ \mu_{0} $'
  label_dist_mod_val =" %0.2f $\pm$ %0.2f mag" % (dist_mod_DR3_dist, dist_mod_DR3_dist_unc)
  ext_ax2.plot([],label = label_dist_mod_val+' ('+label_dist_mod_DR3_dist+')',color=t06,)
  ext_ax.axhline(dist_mod_DR3_dist, color=t06, linewidth=0.6,linestyle='-')
  ext_ax.axhspan(-mu_min, -mu_max, color=t06,alpha=0.4)  #horizontal shading
  
  # mean (m-M)
  if True:
    dist_mod_mean           = np.mean(dist_mod_array[0:-1])
    dist_mod_mean_unc       = 3 * np.std (dist_mod_array[0:-1])
                                       
    label_mean_dist_mod_val = " %0.2f $\pm$ %0.2f mag" % (dist_mod_mean, dist_mod_mean_unc)
    ext_ax2.plot([],label = label_mean_dist_mod_val+' '+'($ \overline{\mu}_{[0,4]} $)',color=t_fit,)
    ext_ax.axhline(dist_mod_mean, color=t_fit,linewidth=0.6)
    ext_ax.axhspan(dist_mod_mean-dist_mod_mean_unc, dist_mod_mean+dist_mod_mean_unc, color=t_fit,alpha=0.1)
  
  ext_ax2.axis('off')
  
  # color_list = [t_lost,t02,t02,t02,t04,t07]
  
  # label_ext = '{0:10.3f}'.format(dist_mod)
  # ext_ax.annotate( label_ext, \
  #                   xy=(0.25,dist_mod), \
  #                   textcoords='offset pixels', \
  #                   arrowprops=dict(arrowstyle='->', facecolor='black'), \
  #                   xytext=(-50,70) \
  #                 )
  
  for i, entry in enumerate(dist_mod_array):
    # ext_ax.scatter(ext_x[i],entry,color=color_list[i],label=str(i) +' '+ dist_mod_model[i]+' - '+dist_mod_model_note[i],s=3 )   
    ext_ax.scatter(ext_x[i],entry,color=color_list[i],label=str(i) +' '+ dist_mod_model[i],s=3 )   
    ext_ax.errorbar(ext_x[i],entry,yerr= dist_mod_unc_array[i], color = color_list[i], fmt='o')   
    
    # label_ext = '{0:10.3f}'.format(entry)
    label_ext = '({0:1.2f}$\pm${1:1.2f})'.format(entry,dist_mod_unc_array[i])
    
    ext_ax.annotate( label_ext, \
                      xy=(ext_x[i],entry), \
                      textcoords='offset pixels', \
                      # xytext=(16,46), \
                        xytext=(-100,120 + 40*np.mod(i,2)), \
                    )
  # ext_ax.legend()
  leg1 = ext_ax.legend(loc='lower left')
  leg1.remove()
  ext_ax2.legend(loc='upper left')
  ext_ax2.add_artist(leg1)
  
  if savefig is True:
    ext_fig.savefig(fig_sub_path+data_set_selected+'_extComp')
  



#%% CMD plot with Absolute Magn-------------------------------------------

      
if CMD_abs_plot_on                     is True:
  cmd_abs_fig, cmd_abs_ax = plt.subplots(1,2,figsize = (10,5),sharex=(False),sharey=False)
  cmd_abs_fig.suptitle('CMD corrected using interstellar extinction')
  
  cmd_abs_ax[0].set_ylabel('absolute photometric magnitude G')
  label = '{0:s}: {1:5.3e}'.format(r'Using distance modulus in mag: ',dist_mod)
  cmd_abs_ax[0].scatter(gaia_cluster_members['bp_rp'],gaia_cluster_members['abs_phot_g_mean_mag'],label=label,s=1,color=t04)
  
  
  cmd_abs_ax[1].set_ylabel('absolute photometric magnitude G')  
  cmd_abs_ax[1].scatter(gaia_cluster_members['bp_rp'][(gaia_cluster_members['phot_g_mean_mag'] < 12.3)],gaia_cluster_members['abs_phot_g_mean_mag'][(gaia_cluster_members['phot_g_mean_mag'] < 12.3)],label=label,s=1,color=t04)
  
  
  
  for i in [0,1]:
    cmd_abs_ax[i].set_xlabel('(blue-red) in mag')
    cmd_abs_ax[i].invert_yaxis()
    # cmd_abs_ax[i].legend(loc = 'lower center')
  cmd_abs_ax[0].legend(loc = 'lower center')

