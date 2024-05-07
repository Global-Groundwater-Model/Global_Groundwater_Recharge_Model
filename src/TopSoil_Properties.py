#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:47:31 2022

@author: sara.nazari
"""
## --------------------------------------------------------------------------------- ##
#   This code performs two operations on saturated and residual water content data   #        
#   1-The raw tif topsoil residual and saturated water content data is modified and  #
#     saved in netcdf format with thickness dimension to match the dimensions and    #
#     avoid unreal zero values mainly in the boundaries.                             #
#   2- Calculate maximum and minimum topsoil water storage                           #
#      capacities with first step output                                             #                 
## --------------------------------------------------------------------------------- ##

import xarray as xr
import numpy as np

print('Start: Topsoil Properties Calculation') 

fn_TopSoil_tetaSat = snakemake.input[0]
fn_TopSoil_tetaRes =snakemake.input[1]
fn_Thickness_TopSoil =snakemake.input[2]

fn_WC_sat_TopSoil = snakemake.output [0]
fn_WC_res_TopSoil = snakemake.output [1]
fn_SWC_Max_TopSoil = snakemake.output[2]
fn_SWC_Min_TopSoil = snakemake.output[3]


ds_TopSoil_tetaSat = xr.open_rasterio(fn_TopSoil_tetaSat)
ds_TopSoil_tetaRes = xr.open_rasterio(fn_TopSoil_tetaRes)
ds_Z1 = xr.open_dataset(fn_Thickness_TopSoil)

WC_sat = ds_TopSoil_tetaSat [0,:,:]*0.0001 # Multiply with 0.0001 because of the raw data necessary convention
WC_sat = xr.where( WC_sat<0, np.nan, WC_sat)

WC_res = ds_TopSoil_tetaRes [0,:,:]*0.0001 # Multiply with 0.0001 because of the raw data necessary convention
WC_res = xr.where(WC_res<0, np.nan, WC_res)

ar_WC_sat = xr.DataArray(WC_sat, dims = ['y','x'], coords =[ds_Z1.coords['y'],ds_Z1.coords['x']], 
                                attrs=dict(description="Saturated TopSoil Water Content",
                                            units='m3/m3'))
ds_WC_sat = ar_WC_sat.to_dataset(name='WC_Sat_TopSoil')
ds_WC_sat.to_netcdf(fn_WC_sat_TopSoil, mode='w', format = 'NETCDF4',engine='netcdf4')

ar_WC_res = xr.DataArray(WC_res, dims = ['y','x'], coords =[ds_Z1.coords['y'],ds_Z1.coords['x']], 
                                attrs=dict(description="Residual TopSoil Water Content",
                                            units='m3/m3'))
ds_WC_Res = ar_WC_res.to_dataset(name='WC_Res_TopSoil')
ds_WC_Res.to_netcdf(fn_WC_res_TopSoil, mode='w', format = 'NETCDF4',engine='netcdf4')

Z1_m = np.empty(shape=(ds_Z1.y.shape[0],ds_Z1.x.shape[0]))
Z1_m [:,:] = ds_Z1.Z1*0.01

WC_dif = ds_WC_sat.WC_Sat_TopSoil - ds_WC_Res.WC_Res_TopSoil
SWC_Max = Z1_m*WC_dif
SWC_Min = ds_WC_Res.WC_Res_TopSoil*Z1_m

ar_SWC_Max = xr.DataArray(SWC_Max, dims = ['y','x'], coords =[ds_Z1.coords['y'],ds_Z1.coords['x']], 
                                attrs=dict(description="Maximum TopSoil Water Storage Capacity",
                                            units='m'))
ds_SWC_Max = ar_SWC_Max.to_dataset(name='SWC_Max_TopSoil')
ds_SWC_Max.to_netcdf(fn_SWC_Max_TopSoil, mode='w', format = 'NETCDF4',engine='netcdf4')

ar_SWC_Min = xr.DataArray(SWC_Min, dims = ['y','x'], coords =[ds_Z1.coords['y'],ds_Z1.coords['x']], 
                                attrs=dict(description="Minimum TopSoil Water Storage Capacity",
                                            units='m'))
ds_SWC_Min = ar_SWC_Min.to_dataset(name='SWC_Min_TopSoil')
ds_SWC_Min.to_netcdf(fn_SWC_Min_TopSoil, mode='w', format = 'NETCDF4',engine='netcdf4')



print('Finish: Topsoil Properties Calculation') 
