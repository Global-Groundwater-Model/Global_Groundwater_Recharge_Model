#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 13:16:21 2022

@author: sara.nazari
"""


## ----------------------------------------------------------------------------- ##
#     Consider only negative values of TEva that show the evapotranspiration      #
## ----------------------------------------------------------------------------- ##

import xarray as xr
import pandas as pd


    
print('Start: Deriving Teva from Raw Data')

fn_TEva_Dims = snakemake.input[0]
fn_TEva = snakemake.output[0]
        
ds_TEva = xr.open_dataset(fn_TEva_Dims)   
ds_TEva_new = xr.where(ds_TEva.Teva > 0, 0, -1*ds_TEva.Teva)
ds_TEva_new.to_netcdf(fn_TEva, mode='w', format = 'NETCDF4',engine='netcdf4')

ds_TEva.close
ds_TEva_new.close

print('Finish: Deriving Teva from Raw Data')


