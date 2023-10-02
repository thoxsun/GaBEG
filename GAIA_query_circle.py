#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 11:25:15 2022

@author: Thomas Suntinger

Input:
  - Coordinates of 
  - and search cone radius for cluster of interesst 

Dependencies:
  - astroquery.gaia

Dependencies for visualization:
  - plotStyle.py
  - GAIA_plot.py
  - GAIA_load.py

Functionality:
  - Download and store Gaia data from archive in csv-file.
  - Query uses ICRS Coordinates and radius for circular search window.
  - Plots a CMD and saves it, if visualization is turned on.

Structure (& parameters):
  - Init (including visualization on/off)
  - Gaia query
  - CMD visualization 

Output:
  - store: csv file with queried data
  - print: csv path and file name 
  - print: number of stars in sample

Output plots:
  - CMD of sample 
  

Improvements:
  - Include convolution with astrophysical_parameters
  - and use it to download necessary parameters for the Gaia BC algorithm 
  
"""

#%%
from astroquery.gaia import Gaia

# Toggle plot and save a CMD of retrieved data
cmd_plot_on = True

if cmd_plot_on == True:
  fig_height = 4
  import plotStyle as plotStyle
  [ready_to_publish, fig_height, fig_width] = plotStyle.set_plot_style(fig_height = fig_height)



#%% Gaia query 

query_all_joint= """SELECT 
                    gaia_source.source_id,
                    gaia_source.designation,
                    
                    gaia_source.ra, 
                    gaia_source.dec,
                    gaia_source.l, 
                    gaia_source.b,
                    gaia_source.ref_epoch,
                    
                    gaia_source.parallax,
                    gaia_source.parallax_error,
                    gaia_source.pm,
                    gaia_source.pmra,
                    gaia_source.pmdec,
                    gaia_source.phot_g_mean_mag,
                    gaia_source.bp_rp,
                    gaia_source.radial_velocity,
                    
                    gaia_source.ruwe,
                    gaia_source.phot_variable_flag,
                    gaia_source.non_single_star,
                    
                    gaia_source.teff_gspphot,
                    gaia_source.logg_gspphot,
                    
                    gaia_source.distance_gspphot,
                    gaia_source.mh_gspphot,
                    gaia_source.azero_gspphot,
                    gaia_source.ag_gspphot,
                    gaia_source.ebpminrp_gspphot,

                    gaia_source.astrometric_params_solved,
                    gaia_source.nu_eff_used_in_astrometry,
                    gaia_source.pseudocolour,
                    gaia_source.ecl_lat,
                    
                    gaiadr3.astrophysical_parameters.alphafe_gspspec,
                    gaiadr3.astrophysical_parameters.mh_gspphot,
                    gaiadr3.astrophysical_parameters.fem_gspspec
                    
FROM gaiadr3.gaia_source 
LEFT OUTER JOIN gaiadr3.astrophysical_parameters USING (source_id)
WHERE 
CONTAINS(
	POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),
 	CIRCLE(
		'ICRS',
		COORD1(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
		COORD2(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
        0.15)
)=1
            AND gaia_source.pmra IS NOT NULL AND abs(gaia_source.pmra)>0 \
            AND gaia_source.pmdec IS NOT NULL AND abs(gaia_source.pmdec)>0 \
            AND abs(gaia_source.pmra_error/gaia_source.pmra)<0.10 \
            AND abs(gaia_source.pmdec_error/gaia_source.pmdec)<0.10 ;",


"""

# https://de.wikipedia.org/wiki/Join_(SQL)#Innerer_Verbund

# astrophysical_parameters.alphafe_gspspec,
# astrophysical_parameters.mh_gspphot,
# astrophysical_parameters.fem_gspspec,
# INNER JOIN gaiadr3.astrophysical_parameters AS astrophysical_parameters USING (source_id)

#---------------------------------M3-------------------------------------------
# Querry for M3
# Radius is 18arcmin 
# 	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(205.548417,28.377278,.1100,-.1400,-2.6400,-147.2000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(205.548417,28.377278,.1100,-.1400,-2.6400,-147.2000,2000,2016.0)),
# 		0.3)
# output_file = "Gaia_M3_r-18arcmin"

# Querry for M3
# Radius is 8.5arcmin 
#  	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(205.548417,28.377278,.1100,-.1400,-2.6400,-147.2000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(205.548417,28.377278,.1100,-.1400,-2.6400,-147.2000,2000,2016.0)),
# 		0.141666667)
# output_file = "Gaia_M3_r-8_5arcmin"


# Querry for M3
# Radius is 32arcmin 
# 	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(205.548417,28.377278,.1100,-.1400,-2.6400,-147.2000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(205.548417,28.377278,.1100,-.1400,-2.6400,-147.2000,2000,2016.0)),
# 		0.5333333333333333)
# output_file = "Gaia_M3_r-32arcmin"

#---------------------------------M5-------------------------------------------

# Coordinates of M5 resolved by GAIA
# Radius is 8arcmin 
#  	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
#         0.13333333333333333)
# output_file = "Gaia_M5_r-8arcmin"

# Coordinates of M5 resolved by GAIA
# Radius is 9arcmin EPOCH_PROP_POS
#  	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
#         0.15)
output_file = "Gaia_M5_r-9arcmin"



# Slight offset from M5 center
# Radius is 8arcmin 
#  	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(229.638957,2.081528,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(229.638957,2.081528,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
#         0.13333333333333333)
# original offset
# CIRCLE('ICRS',229.638957,2.081528,0.13333333333333333)
# output_file = "Gaia_M5_r-8arcmin_slight_offset"

# Coordinates of M5 resolved by GAIA
# Radius is 32arcmin 
#  	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
#         0.5333333333333333)
# output_file = "Gaia_M5_r-32arcmin"

# Coordinates of M5 resolved by GAIA
# Radius is 300arcmin
# 	CIRCLE('ICRS',229.638417,2.081028,5.)
# 	CIRCLE(
# 		'ICRS',
# 		COORD1(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(229.638417,2.081028,.1410,4.0600,-9.8900,53.8000,2000,2016.0)),
#         5)
# output_file = "Gaia_M5_r-300arcmin"

#-----------------------------------M67----------------------------------------

# Coordinates of M67 resolved by GAIA

# Radius is 13arcmin
#  	CIRCLE('ICRS',
#         COORD1(EPOCH_PROP_POS(132.846,11.814,1.1325,-10.9737,-2.9396,33.9200,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(132.846,11.814,1.1325,-10.9737,-2.9396,33.9200,2000,2016.0)),
#         0.21666666666666667)
# output_file = "Gaia_M67_r-13arcmin"

# Radius is 25arcmin
#  	CIRCLE('ICRS',
#         COORD1(EPOCH_PROP_POS(132.846,11.814,1.1325,-10.9737,-2.9396,33.9200,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(132.846,11.814,1.1325,-10.9737,-2.9396,33.9200,2000,2016.0)),
#         0.4166666666666667)
# output_file = "Gaia_M67_r-25arcmin"

# Coordinates of M67 resolved by GAIA
# Radius is 32arcmin
#  	CIRCLE('ICRS',
#         COORD1(EPOCH_PROP_POS(132.846,11.814,1.1325,-10.9737,-2.9396,33.9200,2000,2016.0)),
# 		COORD2(EPOCH_PROP_POS(132.846,11.814,1.1325,-10.9737,-2.9396,33.9200,2000,2016.0)),
#         0.5333333333333333)
# output_file = "Gaia_M67_r-32arcmin"



# Random spot for testing
# CIRCLE('ICRS',228.638417,2.181028,0.13333333333333333)
# output_file = "Gaia_random_test_area"


output_file_path = 'data/queryData/'
save_path_name = output_file_path+output_file+'.csv'
print('file-path: ' + save_path_name)
job_gaiadr3 = Gaia.launch_job_async(query=query_all_joint, dump_to_file=True, output_format='csv', output_file = save_path_name)

#%% CMD
# ------------------ Visualize data -------------------------------------------


if cmd_plot_on == True:

  from GAIA_load import load_query as load_query
  
  file_data  = load_query(save_path_name)
  data_list = []
  data_list.append( [file_data['bp_rp'],file_data['phot_g_mean_mag']] )
  
  print(str(len(file_data))+' stars in sample')
  
  import GAIA_plot as gp
  fig_sub_path = 'gfx/gaiaQuery/'
  
  gp.cmd_plot(data_list, \
              labels=[output_file], 
                colors=['black'], \
                 # data_background = data_background,\
                 # color_data_background = t00,\
                  # label_data_background = label_data_background,\
                  super_title='Raw data from '+output_file,
                    all_in_one=True,
                      ready_to_publish=True,
                        savefig = True,\
                          fig_sub_path = fig_sub_path)
  
   