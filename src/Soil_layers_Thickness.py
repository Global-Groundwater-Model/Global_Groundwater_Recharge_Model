#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 10:00:01 2022

@author: sara.nazari
"""
## ------------------------------------------------------------------------------------- ##
#                 Calculate 2 soil layers and aquifer thickness                           #
## ------------------------------------------------------------------------------------- ##


import xarray as xr
import numpy as np

print('Start: Soil layers thickness calculation') 

fn_Total_Depth = snakemake.input[0]
fn_Topsoil_Depth = snakemake.output[0]
fn_Subsoil_Thickness = snakemake.output[1]
fn_Aquifer_Thickness = snakemake.output[2]


ds_Total_Depth = xr.open_rasterio(fn_Total_Depth)
ds_new = ds_Total_Depth[0,:,:].where(ds_Total_Depth[0,:,:] != -99999.0, np.nan)

Z1 = ds_Total_Depth[0,:,:].where(ds_Total_Depth[0,:,:] <= 30 , 30)
Z1 = Z1.where(Z1 != -99999.0, np.nan)
ar_depth_topsoil = xr.DataArray(Z1, dims = ['y','x'], coords =[ds_Total_Depth.coords['y'],ds_Total_Depth.coords['x']], 
                                attrs=dict(description="Depth and Thickness of the TopSoil layer",
                                           units='cm'))
ds_depth_topsoil = ar_depth_topsoil.to_dataset(name='Z1')
ds_depth_topsoil.to_netcdf(fn_Topsoil_Depth, mode='w', format = 'NETCDF4',engine='netcdf4')

Z2 = xr.where(Z1 != np.nan, ds_new[:,:] - Z1[:,:],np.nan)
Z2= xr.where (Z2 > 170, 170, Z2)
ar_Z2 = xr.DataArray(Z2, dims = ['y','x'], coords =[ds_Total_Depth.coords['y'],ds_Total_Depth.coords['x']],
                     attrs=dict(description="Thickness of the SubSoil layer",
                                units='cm'))
ds_Z2 = ar_Z2.to_dataset(name='Z2')
ds_Z2.to_netcdf(fn_Subsoil_Thickness, mode='w', format = 'NETCDF4',engine='netcdf4')

Depth_Subsoil = xr.where(Z2 != 0 , Z2 + 30, 0)

Aquifer_Thickness = xr.where(Depth_Subsoil != np.nan,ds_new - Depth_Subsoil ,np.nan)
ar_aquifer_thickness = xr.DataArray(Aquifer_Thickness, dims = ['y','x'], coords =[ds_Total_Depth.coords['y'],ds_Total_Depth.coords['x']],
                                attrs=dict(description="Aquifer Thickness",
                                           units='cm'))
ds_aquifer_thickness = ar_aquifer_thickness.to_dataset(name='Z3')
ds_aquifer_thickness.to_netcdf(fn_Aquifer_Thickness, mode='w', format = 'NETCDF4',engine='netcdf4')


print('Finish: Soil layers thickness calculation') 

