#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 10:46:00 2023

@author: Thomas Suntinger

Import recommendation:
  import GAIA_save as gs
  from GAIA_save import save_gaia as save_gaia

Functionality:
  - function, saving structured numpy arrays as GAIA csv files.
  
Improvement:
  


"""

# astropy.table
from astropy.table import Table

def save_gaia(file_data,file_name,file_sub_path=''):
  
  if '.csv' != file_name[-4:]: 
    file_name = file_name+'.csv'
  
  # Construct a table 
  # https://docs.astropy.org/en/stable/table/construct_table.html
  file_table = Table(file_data)  
  
  # Reading and writing tables
  # https://docs.astropy.org/en/stable/table/io.html 
  
  # Table write formats (format='ascii.csv', Note: .csv in file name necessary)
  # https://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers
  file_table.write(file_sub_path + file_name, format='ascii.csv',overwrite=True)

  return 