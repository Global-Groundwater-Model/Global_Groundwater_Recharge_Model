#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:53:59 2023

@author: sara.nazari
"""

## ------------------------------------------------------------------ ##
#  Merg the daily GPM IMERG Final Precipitation L3 1 day 0.1 degree x #
#  0.1 degree (GPM_3IMERGDF) data to create each year daily rainfall, #
#                        variable name: Pr                            #
## ------------------------------------------------------------------ ##

import xarray as xr
import pandas as pd
import numpy as np


print('Start: Merg Daily Precipitation Based on IMERG Data to Create Each Year from '+str(snakemake.params[0])+' to '+str(snakemake.params[0]) + ' Daily Precipitation')

# Setting Coordinate and time steps based on first dataset of IMERG
fn_first_Day = snakemake.input[0]
ds_first_Day = xr.open_dataset(fn_first_Day)
ds_first_Day = xr.open_dataset(fn_first_Day)
Latitude = ds_first_Day.lat
Longitude = ds_first_Day.lon


Num_Lat = Latitude.shape[0]
Num_Lon = Longitude.shape[0] 
time = 0

year = snakemake.params[0]
df_time = pd.DataFrame({'Date':pd.date_range(str(year)+'-01-01', str(year)+'-12-31', freq='D')})
Num_Time = df_time.shape[0]
Daily_Pr = np.empty(shape=(Num_Time,Num_Lat,Num_Lon))

print('start',year)
        
for timesteps in range (0,Num_Time):
    fn = snakemake.params[1] + str(df_time.Date.dt.strftime('%Y%m%d').astype(int)[timesteps])+'.nc'                 
    ds_next_day = xr.open_dataset(fn)
    Daily_Pr[time,:,:] = ds_next_day.precipitationCal[:].transpose('time','lat','lon')
    time = time +1
    print(time)
    
ar_Annual_Pr = xr.DataArray(Daily_Pr, dims = ['time','latitude','longitude'],coords =[df_time.Date[:], Latitude, Longitude], 
                                attrs=dict(description="IMERG_Daily_Precipitation", units='mm/day'))                    
ds_Annual_Pr = ar_Annual_Pr.to_dataset(name='Pr')
ds_Annual_Pr.attrs['title'] = 'GPM IMERG Final Precipitation L3 1 day 0.1 degree x 0.1 degree (GPM_3IMERGDF)'   
ds_Annual_Pr.to_netcdf(snakemake.output[0], mode='w', format = 'NETCDF4',engine='netcdf4')

print('finish',year)
    
print('Finish: Merg Daily Precipitation Based on IMERG Data to Create Each Year from '+str(year)+' to '+str(year) + ' Daily Precipitation')
