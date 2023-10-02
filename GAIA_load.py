#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 12:13:27 2023

@author: Thomas Suntinger

Import recommendation:
  import GAIA_load as gl
  from GAIA_load import load_query as load_query
  
Functionality:
  It loads and returns GAIA csv files as structured numpy arrays.

Improvement:
  
  
"""


import numpy as np
import os

# Load csv file as astropy.Table object
# Reading and writing tables
# https://docs.astropy.org/en/stable/table/io.html 

def load_query_astroTable(file_name, n=np.nan):
  file_name = '/'+ file_name
  file_path = os.getcwd() + file_name
  
  from astropy.table import Table
  file_data = Table.read(file_path, format='ascii.csv')
  
  return file_data


  
# Load csv file as structured numpy array
def load_query(file_name):
  
  if '.csv' != file_name[-4:]: 
    file_name = file_name+'.csv'
  file_name = '/'+ file_name
  file_path = os.getcwd() + file_name
    
  file_data = np.genfromtxt(file_path,skip_header = 0, encoding = None, delimiter = ',' ,dtype=None, names = True)
     
  return file_data