#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 20:55:41 2022

@author: sara.nazari
"""

## ------------------------------------------------------------------------------------------------- ##
#   In this code the following procedures are done for saturated and residual water content data     #         
# 1-The raw tif data of the residual and saturated water content of the subsoil are modified and     # 
#   saved in netcdf format with the thickness dimension to match the dimensions and                  #
#   avoid unreal zero values mainly in the boundaries.                                               #
# 2-There are some grids between ( 79.75 43.65, 81.55 43.65, 80.45 45, 82.35 45)                     # 
#   coordinates that have high values for saturated water content,                                   #
#   So here also this value are changed to the maximum value available for other grids.              #
# 3-Calculate maximum and minimum SubSoil water storage                                              #
#   Capacities with output from first step                                                           #                 
## ------------------------------------------------------------------------------------------------- ##

import xarray as xr
import numpy as np

print('Start: SubSoil Properties Calculation') 

fn_SubSoil_tetaSat = snakemake.input[0]
fn_SubSoil_tetaRes =snakemake.input[1]
fn_Thickness_SubSoil =snakemake.input[2]

fn_WC_sat_SubSoil = snakemake.output [0]
fn_WC_res_SubSoil = snakemake.output [1]
fn_SWC_Max_SubSoil = snakemake.output[2]
fn_SWC_Min_SubSoil = snakemake.output[3]

ds_SubSoil_tetaSat = xr.open_rasterio(fn_SubSoil_tetaSat)
ds_SubSoil_tetaRes = xr.open_rasterio(fn_SubSoil_tetaRes)
ds_Z2 = xr.open_dataset(fn_Thickness_SubSoil)

WC_sat = ds_SubSoil_tetaSat [0,:,:]*0.0001 # Multiply with 0.0001 because of the raw data necessary convention
WC_sat = xr.where(WC_sat<0, np.nan, 
                  xr.where( WC_sat > 1, 0.72, WC_sat))

WC_res = ds_SubSoil_tetaRes [0,:,:]*0.0001 # Multiply with 0.0001 because of the raw data necessary convention
WC_res = xr.where(WC_res<0, np.nan, WC_res)

ar_WC_sat = xr.DataArray(WC_sat, dims = ['y','x'], coords =[ds_Z2.coords['y'],ds_Z2.coords['x']], 
                                attrs=dict(description="Saturated SubSoil Water Content",
                                            units='m3/m3'))
ds_WC_sat = ar_WC_sat.to_dataset(name='WC_Sat_SubSoil')
ds_WC_sat.to_netcdf(fn_WC_sat_SubSoil, mode='w', format = 'NETCDF4',engine='netcdf4')

ar_WC_res = xr.DataArray(WC_res, dims = ['y','x'], coords =[ds_Z2.coords['y'],ds_Z2.coords['x']], 
                                attrs=dict(description="Residual SubSoil Water Content",
                                            units='m3/m3'))
ds_WC_Res = ar_WC_res.to_dataset(name='WC_Res_SubSoil')
ds_WC_Res.to_netcdf(fn_WC_res_SubSoil, mode='w', format = 'NETCDF4',engine='netcdf4')

Z2_m = np.empty(shape=(ds_Z2.y.shape[0],ds_Z2.x.shape[0]))
Z2_m [:,:] = ds_Z2.Z2*0.01

WC_dif = ds_WC_sat.WC_Sat_SubSoil - ds_WC_Res.WC_Res_SubSoil
SWC_Max = Z2_m*WC_dif
SWC_Min = ds_WC_Res.WC_Res_SubSoil*Z2_m

ar_SWC_Max = xr.DataArray(SWC_Max, dims = ['y','x'], coords =[ds_Z2.coords['y'],ds_Z2.coords['x']], 
                                attrs=dict(description="Maximum SubSoil Water Storage Capacity",
                                            units='m'))
ds_SWC_Max = ar_SWC_Max.to_dataset(name='SWC_Max_SubSoil')
ds_SWC_Max.to_netcdf(fn_SWC_Max_SubSoil, mode='w', format = 'NETCDF4',engine='netcdf4')

ar_SWC_Min = xr.DataArray(SWC_Min, dims = ['y','x'], coords =[ds_Z2.coords['y'],ds_Z2.coords['x']], 
                                attrs=dict(description="Minimum SubSoil Water Storage Capacity",
                                            units='m'))
ds_SWC_Min = ar_SWC_Min.to_dataset(name='SWC_Min_SubSoil')
ds_SWC_Min.to_netcdf(fn_SWC_Min_SubSoil, mode='w', format = 'NETCDF4',engine='netcdf4')


print('Finish: SubSoil Properties Calculation') 
