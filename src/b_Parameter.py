#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 11:33:34 2022

@author: sara.nazari
"""

## ------------------------------------------------------------------------------------- ##
#  Calculate b parameter describing the sub-grid distribution of soil moisture capacities #
## ------------------------------------------------------------------------------------- ##


import numpy as np
import xarray as xr

print('Start: Soil Moisture Shape Parameter Calculation') 

ds_topo_std = xr.open_dataset(snakemake.input[0], decode_times=False)
ds_beta = xr.open_dataset(snakemake.input[1], decode_times=False)
ds_TopSoil_Depth = xr.open_dataset(snakemake.input[2], decode_times=False)
     
topo_std = ds_topo_std.__xarray_dataarray_variable__
beta = ds_beta.__xarray_dataarray_variable__
topo_min = 100
topo_max = 1000
b_oro = (topo_std[0,:,:] - topo_min)/(topo_std[0,:,:] + topo_max)
b = np.empty(shape=(1,ds_topo_std.y.shape[0],ds_topo_std.x.shape[0]))
    
for i in range(0,ds_topo_std.y.shape[0]):
    for j in range(0,ds_topo_std.x.shape[0]):
        if b_oro [i,j]> 0.01:
            b[0,i,j] = beta[0,i,j] + b_oro[i,j]
    
        else:
            b[0,i,j] = beta[0,i,j]
        j += 1
    i += 1
    print(i)
        
ar_b = xr.DataArray(b, dims = ['band', 'latitude', 'longitude'],
                    coords = [ds_topo_std.coords['band'],ds_TopSoil_Depth.coords['y'], ds_TopSoil_Depth.coords['x']],
                    attrs=dict(description="shape parameter b describing the sub-grid distribution of soil moisture capacities"))
    
ds_b = ar_b.to_dataset(name='b_Parameter')
ds_b.attrs['title']='shape parameter b'
ds_b.to_netcdf(snakemake.output[0], mode='w', format = 'NETCDF4',engine='netcdf4')
   
ds_beta.close
ds_topo_std.close
ds_TopSoil_Depth.close
ds_b.close

print('Finish: Soil Moisture Shape Parameter Calculation') 
