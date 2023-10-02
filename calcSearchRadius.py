#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 20:25:15 2023

@author: tom

Import recommentdation:
  import calcSearchRadius
"""

import numpy as np

#%% ### Functions for search radius calculation


def convertArcMinToRad(angle):
  # angle convert from arcmin to deg
  SR_angle = angle / 60
  # convert from deg to rad
  SR_angle = SR_angle *2*np.pi / 360
  
  return SR_angle

# d in any distance unit yields same distance unit for the search radius
# SR_angle in unit rad
def calcSearchRadius(d,SR_angle):
  SR = (d) * np.tan(SR_angle) 
  return SR



