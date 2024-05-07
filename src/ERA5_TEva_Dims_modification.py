#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 15:22:48 2022

@author: sara.nazari
"""

## -------------------------------------------------------------------- ##
#  Change the dimension of Total Evapotranspiration (Teva) from ERA5    #
#           to the same dimension as IMERG rainfall dataset             #
## -------------------------------------------------------------------- ##

import xarray as xr
import numpy as np

print('Start: Teva Dimension Modification') 

fnPre = snakemake.input[0]
fnTEva = snakemake.input[1]

dsPr = xr.open_dataset(fnPr)
dsEva = xr.open_dataset(fnTEva)


Time_Steps = dsPr.time
Latitude = dsPr.latitude
Longitude = dsPr.longitude

Num_Time_Steps = Time_Steps.shape[0]
Num_Lat = Latitude.shape[0]
Num_Lon = Longitude.shape[0] 

nEva = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon)) 
nnEva = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon)) 

dsEva = dsEva.e[:]
i = 0

# Average of two neighbour grids in lat dim

while i < Num_Lat:
    
    nEva[:,i,:] = (dsEva[:,i,:] + dsEva[:,i+1,:])/2
    i +=1
    
# Average of two neighbour grids in lon dim   
j = 0    
while j < Num_Lon-1:
    
    nnEva[:,:,j] = (nEva[:,:,j] + nEva[:,:,j+1])/2
    j +=1

dsnEva = xr.DataArray(nnEva,dims = ['time', 'latitude', 'longitude'], 
                          coords = [dsEva.coords['time'],-1*dsPr.latitude, 
                                    dsPr.longitude],
                          attrs=dict(description="Matching Dimensions of Daily ERA5 Total Evapotranspiration and IMERG Total Preceipitation",
                                     units="m/day"))                          
ds = dsnEva.to_dataset(name='Teva')
ds.to_netcdf(snakemake.output[0], mode='w', format = 'NETCDF4',engine='netcdf4')

print('Finish: Teva Dimension Modification') 
