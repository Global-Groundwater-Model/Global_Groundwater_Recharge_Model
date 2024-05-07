#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 22:03:12 2022

@author: sara.nazari
"""

## ----------------------------------------------------------------------- ##
#    Calculate Sum of Daily IMERG Rainfall (Pr) and ERA5 Snowmelt (Smlt)   #
## ----------------------------------------------------------------------- ##



import xarray as xr
import numpy as np
import geopandas as gpd
from shapely.geometry import mapping
   
print('Start: Calculate daily sum rainfall Pr and snowmelt Smlt') 

fn_Pr =  snakemake.input[0]
fn_Smlt = snakemake.input[1]
fn_Sum_Pr_Smlt = snakemake.output[0]
    
dsPr = xr.open_dataset(fn_Pr, decode_times=False)
dsSmlt = xr.open_dataset(fn_Smlt)

time = dsPr.time[:]
Latitude = dsPr.latitude[:]
Longitude = dsPr.longitude[:]
Num_time = time.shape[0]
Num_Lat = Latitude.shape[0]
Num_Lon = Longitude.shape[0] 
Pr = np.empty(shape=(Num_time,Num_Lat,Num_Lon))
Smlt = np.empty(shape=(Num_time,Num_Lat,Num_Lon))
    
Pr [:,:,:] = dsPr.Pr
Smlt_reindex = dsSmlt.reindex(latitude=list(reversed(dsSmlt.latitude)))
Smlt[:,:,:] = Smlt_reindex.Smlt

print(1)

Pr_Smlt =  Pr [:,:,:] + (Smlt[:,:,:]*1000)
ar_Pr_Smlt = xr.DataArray(Pr_Smlt,
                            dims = ['time', 'latitude', 'longitude'],
                            coords = [dsSmlt.coords['time'],dsPr.coords['latitude'], dsPr.coords['longitude']],
                            attrs=dict(description="Sum of daily total rainfall Pr (IMERG) and snowmelt Smlt (ERA5)",
                                        units="mm/day"))
ds_Pr_Smlt = ar_Pr_Smlt.to_dataset(name='Pr_Smlt')

print(2)

ds_Pr_Smlt_int = ds_Pr_Smlt.interpolate_na(dim=("time"), method="linear")
ds_Pr_Smlt_int = ds_Pr_Smlt_int.interpolate_na(dim=("latitude"), method="linear")
ds_Pr_Smlt_int = ds_Pr_Smlt_int.interpolate_na(dim=("longitude"), method="linear")

ds_Pr_Smlt_int.attrs['title']='Daily Rainfall and Snowmelt'
ds_Pr_Smlt_int = ds_Pr_Smlt_int.rio.write_crs("epsg:4326")

World_Basin_Shape = snakemake.input[2]
World_Shape = gpd.read_file(World_Basin_Shape, crs="epsg:4326")
clip_ds_Pr_Smlt_int = ds_Pr_Smlt_int.rio.clip(World_Shape.geometry.apply(mapping), World_Shape.crs, drop=False)
clip_ds_Pr_Smlt_int.to_netcdf(fn_Sum_Pr_Smlt, mode='w', format = 'NETCDF4',engine='netcdf4')

dsPr.close
dsSmlt.close
ds_Pr_Smlt.close
ds_Pr_Smlt_int.close
clip_ds_Pr_Smlt_int.close

print('Finish: Calculate daily sum rainfall Pr and snowmelt Smlt') 



