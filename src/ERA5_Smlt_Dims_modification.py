#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 16:06:21 2022

@author: sara.nazari
"""

## ------------------------------------------------------------------ ##
#            Change the dimension of Snowmelt from ERA5               #
#           to the same dimension as IMERG rainfall dataset           #
## ------------------------------------------------------------------ ##

import xarray as xr
import numpy as np

print('Start: Smlt Dimension Modification') 

fnPr = snakemake.input[0]
fnSmlt = snakemake.input[1]

dsPr = xr.open_dataset(fnPr)
dsSmlt = xr.open_dataset(fnSmlt)


Time_Steps = dsPr.time
Latitude = dsPr.latitude
Longitude = dsPr.longitude

Num_Time_Steps = Time_Steps.shape[0]
Num_Lat = Latitude.shape[0]
Num_Lon = Longitude.shape[0] 

nsmlt = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon)) 
nnSmlt = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon)) 

dsSmlt = dsSmlt.smlt[:]
i = 0

# Average of two neighbour grids in lat dim

while i < Num_Lat:
    
    nsmlt[:,i,:] = (dsSmlt[:,i,:] + dsSmlt[:,i+1,:])/2
    i +=1
    
# Average of two neighbour grids in lon dim   
j = 0    
while j < Num_Lon-1:
    
    nnSmlt[:,:,j] = (nsmlt[:,:,j] + nsmlt[:,:,j+1])/2
    j +=1

dsnSmlt = xr.DataArray(nnSmlt,dims = ['time', 'latitude', 'longitude'], 
                          coords = [dsSmlt.coords['time'],-1*dsPr.latitude, 
                                    dsPr.longitude],
                          attrs=dict(description="Matching Dimensions of Daily ERA5 Snowmelt and IMERG Total Preceipitation",
                                     units="m/day"))                          
ds = dsnSmlt.to_dataset(name='Smlt')
ds.to_netcdf(snakemake.output[0], mode='w', format = 'NETCDF4',engine='netcdf4')
print('Finish: Smlt Dimension Modification') 
