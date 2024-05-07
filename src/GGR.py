#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 17:02:32 2023

@author: sara.nazari
"""


import xarray as xr
import numpy as np

print('Start: Global Groundwater Recharge Model') 

def GGR(fn_S_TopSoil_Pre,fn_S_SubSoil_Pre,fn_Sum_Pr_Smlt,fn_Teva,year,fn_S_top,fn_R_run,fn_R_sub,fn_CR,fn_S_Sub,fn_R_gw):
             
    # Soil Properties Data
    ## TopSoil Properties
    fn_SWC_Max_TopSoil = snakemake.input[4]
    fn_SWC_Min_TopSoil = snakemake.input[5]
    fn_b = snakemake.input[6]
    fn_TopSoil_N = snakemake.input[7]
    fn_TopSoil_Ksat = snakemake.input[8]
    fn_TopSoil_tetaSat = snakemake.input[9]
    fn_TopSoil_tetaRes = snakemake.input[10]
    fn_Topsoil_Thickness = snakemake.input[11]
        
    ## SubSoil Properties
    fn_SWC_Max_SubSoil = snakemake.input[12]
    fn_SWC_Min_SubSoil = snakemake.input[13]
    fn_SubSoil_N = snakemake.input[14]
    fn_SubSoil_Ksat = snakemake.input[15]
    fn_SubSoil_tetaSat = snakemake.input[16]
    fn_SubSoil_tetaRes = snakemake.input[17]
    fn_Subsoil_Thickness = snakemake.input[18]
             
    # Read input data
    ds_Pr_Smlt = xr.open_dataset(fn_Sum_Pr_Smlt)
    ds_Teva = xr.open_dataset(fn_Teva)
    
    ds_SWC_Max_TopSoil = xr.open_dataset(fn_SWC_Max_TopSoil)
    ds_SWC_Min_TopSoil = xr.open_dataset(fn_SWC_Min_TopSoil)
    ds_b = xr.open_dataset(fn_b)
    ds_N_TopSoil = xr.open_rasterio(fn_TopSoil_N)
    ds_Ksat_TopSoil = xr.open_rasterio(fn_TopSoil_Ksat)
    ds_tetaSat_TopSoil = xr.open_rasterio(fn_TopSoil_tetaSat)
    ds_tetares_TopSoil = xr.open_rasterio(fn_TopSoil_tetaRes)
    ds_Z1 = xr.open_dataset(fn_Topsoil_Thickness)
    ds_SWC_Min_SubSoil = xr.open_dataset(fn_SWC_Min_SubSoil)
    ds_SWC_Max_SubSoil = xr.open_dataset(fn_SWC_Max_SubSoil)
    ds_N_SubSoil = xr.open_rasterio(fn_SubSoil_N)
    ds_Ksat_SubSoil = xr.open_rasterio(fn_SubSoil_Ksat)
    ds_tetaSat_SubSoil = xr.open_rasterio(fn_SubSoil_tetaSat)
    ds_tetares_SubSoil = xr.open_rasterio(fn_SubSoil_tetaRes)
    ds_Z2 = xr.open_dataset(fn_Subsoil_Thickness)

    # For Meteo Input Data Set an Empty Array Equals to the Input Meteo
    # Reindex the Latitude of the Sum Pr and Smlt to match other array's latitudes
    ds_Pr_Smlt_reindex = ds_Pr_Smlt.reindex(latitude=list(reversed(ds_Pr_Smlt.latitude)))
        
    # Setting Coordinate and time steps, considering devision  
    Time_Steps = ds_Pr_Smlt.time
    Latitude = ds_Pr_Smlt_reindex.latitude[:]
    Longitude = ds_Pr_Smlt.longitude[:]
    Num_Time_Steps = Time_Steps.shape[0]
    Num_Lat = Latitude.shape[0]
    Num_Lon = Longitude.shape[0] 
    Time_Steps_S  = np.empty(shape=(Num_Time_Steps+1), dtype='datetime64[ns]')
    Time_Steps_S[1:Num_Time_Steps+1] = Time_Steps
        
    lat1 = ds_Pr_Smlt_reindex.latitude[0]
    lat2 = ds_Pr_Smlt_reindex.latitude[1799]
    lon1 = ds_Pr_Smlt_reindex.longitude[0]
    lon2 = ds_Pr_Smlt_reindex.longitude[3599]
               
    # Get Values from input datasets and assign them to defined empty arrays 
    
    Teva = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon))
    Teva[:,:,:] = ds_Teva.Teva[:,:,:]*1000 # Change the unit to mm

    Pr_smlt = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon))
    Pr_smlt[:,:,:] = ds_Pr_Smlt_reindex.Pr_Smlt[:,:,:]
    
    ## TopSoil                      
    SWC_Max =np.empty(shape=(Num_Lat,Num_Lon))
    SWC_Max[:,:] = ds_SWC_Max_TopSoil.SWC_Max_TopSoil[:,:]*1000 # Change the unit to mm
    
    SWC_Min = np.empty(shape=(Num_Lat,Num_Lon))
    SWC_Min[:,:] = ds_SWC_Min_TopSoil.SWC_Min_TopSoil[:,:]*1000 # Change the unit to mm
    
    b = np.empty(shape=(Num_Lat,Num_Lon))
    b[:,:] = ds_b.b_Parameter[0,:,:]
    
    tetaSat_TopSoil = np.empty(shape=(Num_Lat,Num_Lon))
    tetaSat_TopSoil = ds_tetaSat_TopSoil[0,:,:]*0.0001 # Multiple by 0.0001 because of dataset unit correction 
    tetares_TopSoil = np.empty(shape=(Num_Lat,Num_Lon))
    tetares_TopSoil = ds_tetares_TopSoil[0,:,:]*0.0001 
    N_TopSoil = np.empty(shape=(Num_Lat,Num_Lon))
    N_TopSoil[:,:] = ds_N_TopSoil[0,:,:]*0.0001 # Multiple by 0.0001 because of dataset unit correction
    Ksat_TopSoil = np.empty(shape=(Num_Lat,Num_Lon))
    Ksat_TopSoil[:,:] = ds_Ksat_TopSoil[0,:,:]*0.001 # Multiple by 0.0001 because of dataset unit correction and by 10 change cm to mm
    
    Z1 = np.empty(shape=(Num_Lat,Num_Lon))
    Z1[:,:] = ds_Z1.Z1[:,:]*10 # Change cm to mm
    
    ## SubSoil
    SWC_Min_SubSoil = np.empty(shape=(Num_Lat,Num_Lon))
    SWC_Min_SubSoil[:,:] = ds_SWC_Min_SubSoil.SWC_Min_SubSoil[:,:]*1000  # Change the unit to mm

    SWC_Max_SubSoil = np.empty(shape=(Num_Lat,Num_Lon))
    SWC_Max_SubSoil[:,:] = ds_SWC_Max_SubSoil.SWC_Max_SubSoil[:,:]*1000  # Change the unit to mm

    tetaSat_SubSoil = np.empty(shape=(Num_Lat,Num_Lon))
    tetaSat_SubSoil = ds_tetaSat_SubSoil[0,:,:]*0.0001 # Multiple by 0.0001 because of dataset unit correction 
    tetares_SubSoil = np.empty(shape=(Num_Lat,Num_Lon))
    tetares_SubSoil = ds_tetares_SubSoil[0,:,:]*0.0001       
    N_SubSoil = np.empty(shape=(Num_Lat,Num_Lon))
    N_SubSoil[:,:] = ds_N_SubSoil[0,:,:]*0.0001
    Ksat_SubSoil = np.empty(shape=(Num_Lat,Num_Lon))
    Ksat_SubSoil[:,:] = ds_Ksat_SubSoil[0,:,:]*0.001
    
    Z2 = np.empty(shape=(Num_Lat,Num_Lon))
    Z2[:,:] = ds_Z2.Z2[:,:]*10
    
    
    # Calculate Constant values
    
    ## TopSoil
    SWC_Max_Min = SWC_Max - SWC_Min
    b1 = b+1
    b2 = 1/b1
    b_dif_W = b1*SWC_Max_Min
    M_TopSoil  = np.empty(shape=(Num_Lat,Num_Lon))
    N_1_TopSoil = 1 / N_TopSoil [:,:]
    M_TopSoil [:,:] = 1-(N_1_TopSoil)
    M_1_TopSoil = -1*(1/M_TopSoil)
    dif_tetaSatRes_TopSoil = np.empty(shape=(Num_Lat,Num_Lon))
    dif_tetaSatRes_TopSoil[:,:] = tetaSat_TopSoil[:,:] - tetares_TopSoil[:,:]  
    
    ## SubSoil
    M_SubSoil  = np.empty(shape=(Num_Lat,Num_Lon))
    N_1_SubSoil = 1 / N_SubSoil [:,:]
    M_SubSoil [:,:] = 1-(N_1_SubSoil)
    M_1_SubSoil = -1*(1/M_SubSoil)
    dif_tetaSatRes_SubSoil = np.empty(shape=(Num_Lat,Num_Lon))
    dif_tetaSatRes_SubSoil[:,:] = tetaSat_SubSoil[:,:] - tetares_SubSoil[:,:] 
    
    # Assign empty and initial values for the ouputs
    
    ## TopSoil
    
    S_TopSoil = np.empty(shape=(Num_Time_Steps+1,Num_Lat,Num_Lon))
    ds_TopSoil_Pre = xr.open_dataset(fn_S_TopSoil_Pre)
    Stop = ds_TopSoil_Pre.S_top[-1,:,:]
    S_TopSoil[0,:,:] = Stop.loc[lat1:lat2,lon1:lon2]
    S_TopSoil_check = np.empty(shape=(Num_Time_Steps+1,Num_Lat,Num_Lon))
    S_TopSoil_check [0,:,:] = Stop.loc[lat1:lat2,lon1:lon2]
    
    Runoff = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon))
    R_1 = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon))
    
    ## SubSoil
    S_SubSoil = np.empty(shape=(Num_Time_Steps+1,Num_Lat,Num_Lon))
    ds_SubSoil_Pre = xr.open_dataset(fn_S_SubSoil_Pre)
    Ssub = ds_SubSoil_Pre.S_sub[-1,:,:]
    S_SubSoil[0,:,:] = Ssub.loc[lat1:lat2,lon1:lon2]
    S_SubSoil_check = np.empty(shape=(Num_Time_Steps+1,Num_Lat,Num_Lon))
    S_SubSoil_check [0,:,:] = Ssub.loc[lat1:lat2,lon1:lon2]
    
    R_2 = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon))
    CR = np.empty(shape=(Num_Time_Steps,Num_Lat,Num_Lon))
    
    for time in range(1,Num_Time_Steps+1):
        
        # Calculate Runoff
        Cond_1 = S_TopSoil[time-1,:,:] - SWC_Min[:,:]
        Coef_W_1 = (SWC_Max[:,:] - S_TopSoil[time-1,:,:])
        Coef_W_2 = Coef_W_1[:,:]/SWC_Max_Min[:,:]    
        Coef_W_3 = np.power(Coef_W_2[:,:],b2[:,:])
        Cond_2 = (b1[:,:])*(SWC_Max_Min[:,:])*(Coef_W_3[:,:])
        Runoff[time-1,:,:] = xr.where(Pr_smlt[time-1,:,:] < Cond_1[:,:], 0,
                                    (xr.where(Pr_smlt[time-1,:,:] >= Cond_2[:,:], Pr_smlt[time-1,:,:] - Coef_W_1[:,:],
                                            (Pr_smlt[time-1,:,:] - Coef_W_1[:,:] + SWC_Max_Min[:,:] * np.power((Coef_W_3[:,:]-(Pr_smlt[time-1,:,:] /b_dif_W[:,:])),b1[:,:])))))

        # Calculate R_1
        teta_TopSoil = S_TopSoil[time-1,:,:]/Z1[:,:]  # When Z1 is zero, it will be nan
        effective_saturation_TopSoil = np.where(teta_TopSoil < tetares_TopSoil , 0,
                                                ((teta_TopSoil - tetares_TopSoil[:,:])/dif_tetaSatRes_TopSoil[:,:]))
        elm_1  = np.power(effective_saturation_TopSoil, M_1_TopSoil[:,:])
        Alpha_h_TopSoil  = np.power((elm_1 -1),N_1_TopSoil[:,:])
        Coef_TopSoil_1 = np.power(Alpha_h_TopSoil,(N_TopSoil[:,:]-1))
        Coef_TopSoil_2 = np.power(Alpha_h_TopSoil,N_TopSoil[:,:])
        Coef_TopSoil_3 = np.power(1+Coef_TopSoil_2,(-1*M_TopSoil[:,:]))
        numer_TopSoil  = np.power((1-(Coef_TopSoil_1*Coef_TopSoil_3)),2)
        denom_TopSoil  = np.power((1 + Coef_TopSoil_2),(M_TopSoil[:,:]/2))
        R_1_Value = Ksat_TopSoil[:,:]*(numer_TopSoil/denom_TopSoil)
        R_1[time-1,:,:] = xr.where(effective_saturation_TopSoil != 0, R_1_Value, 0) # add this condition to consider when the effective saturation of topsoil is zero, then R1 is zero                              
        
        # Calculate R_2 
        teta_SubSoil = S_SubSoil[time-1,:,:]/Z2[:,:] # When Z2 is zero, it will be nan
        effective_saturation_SubSoil = np.where(teta_SubSoil < tetares_SubSoil , 0,
                                                ((teta_SubSoil - tetares_SubSoil[:,:])/dif_tetaSatRes_SubSoil[:,:]))
        elm_2  = np.power(effective_saturation_SubSoil, M_1_SubSoil[:,:])
        Alpha_h_SubSoil = np.power((elm_2 -1),N_1_SubSoil[:,:])
        Coef_SubSoil_1 = np.power(Alpha_h_SubSoil,(N_SubSoil[:,:]-1))
        Coef_SubSoil_2 = np.power(Alpha_h_SubSoil,N_SubSoil[:,:])
        Coef_SubSoil_3 = np.power(1+Coef_SubSoil_2,(-1*M_SubSoil[:,:]))
        numer_SubSoil = np.power((1-(Coef_SubSoil_1*Coef_SubSoil_3)),2)
        denom_SubSoil = np.power((1 + Coef_SubSoil_2),(M_SubSoil[:,:]/2))
        R_2[time-1,:,:] = xr.where(effective_saturation_SubSoil != 0, Ksat_SubSoil[:,:]*(numer_SubSoil/denom_SubSoil),0) # add this condition to consider when the effective saturation of subsoil is zero, then R2 is zero 
        
        # Calculate CR
        #CR[time-1,:,:] = R_2[time-1,:,:] * (1-(S_TopSoil[time-1,:,:]/SWC_Max[:,:])) # When Z1, or SWC_Max is zero, it will be nan
        CR[time-1,:,:] =xr.where(effective_saturation_TopSoil<effective_saturation_SubSoil, (R_2[time-1,:,:] *effective_saturation_SubSoil*(1-effective_saturation_TopSoil)),
                                 0)     
        
        # Calculate S of TopSoil, SubSoil, and Aquifer
        S_TopSoil_check [time,:,:] = S_TopSoil[time-1,:,:] + Pr_smlt[time-1,:,:]  - Teva[time-1,:,:]- Runoff[time-1,:,:] - R_1[time-1,:,:] + CR[time-1,:,:]
         
        ## Check S_TopSoil to be in between maximum and minimum possible water storage
        S_TopSoil[time,:,:] = xr.where(S_TopSoil_check [time,:,:] <= SWC_Min[:,:], SWC_Min[:,:],
                                        xr.where(S_TopSoil_check [time,:,:] >= SWC_Max[:,:], SWC_Max[:,:],
                                                S_TopSoil_check [time,:,:]))
                                        
        # If S_TopSoil is less than SWC_Min decrease R_1, since the maximum possible R_1 is passed
        R_1[time-1,:,:] = xr.where((S_TopSoil_check [time,:,:] <= SWC_Min[:,:]) &  (R_1[time-1,:,:] > (SWC_Min[:,:]-S_TopSoil_check [time,:,:])),
                                    R_1[time-1,:,:] - (SWC_Min[:,:]-S_TopSoil_check [time,:,:]),
                                    R_1[time-1,:,:]) # possible to R-1 negative values
        # if S_topSoil is saturated so the rest should change to runoff 
        Runoff[time-1,:,:] = xr.where(S_TopSoil_check [time,:,:] >= SWC_Max[:,:], Runoff[time-1,:,:] + (S_TopSoil_check [time,:,:] - SWC_Max[:,:]),
                                        Runoff[time-1,:,:])
        
        
        # Check S_SubSoil
        S_SubSoil_check[time,:,:]  = S_SubSoil[time-1,:,:] + R_1[time-1,:,:] - CR[time-1,:,:] - R_2[time-1,:,:]
                                                
        S_SubSoil[time,:,:] =  xr.where(S_SubSoil_check [time,:,:] <= SWC_Min_SubSoil[:,:], SWC_Min_SubSoil[:,:],
                                        xr.where(S_SubSoil_check [time,:,:] >= SWC_Max_SubSoil[:,:], SWC_Max_SubSoil[:,:],
                                                S_SubSoil_check [time,:,:]))
        # if S_Subsoil is less than Min then the R_2 value should decrease, so it should be updated and calculate again 
        R_2[time-1,:,:] =  xr.where(S_SubSoil_check [time,:,:] <= SWC_Min_SubSoil[:,:],
                                    ((S_SubSoil[time-1,:,:]- SWC_Min_SubSoil[:,:] + R_1[time-1,:,:])/(1+(S_TopSoil[time,:,:]/SWC_Max[:,:]))),
                                    R_2[time-1,:,:]) 
        R_2[time-1,:,:] = xr.where( R_2[time-1,:,:] <0, 0, R_2[time-1,:,:]) 
        
        # if S_Subsoil is less than Min, R_2 is updated so CR has to be updated too
        #CR[time-1,:,:] = xr.where(S_SubSoil_check [time,:,:] <= SWC_Min_SubSoil[:,:], (R_2[time-1,:,:]*(1-(S_TopSoil[time,:,:]/SWC_Max[:,:]))),
        #                           CR[time-1,:,:])
        CR[time-1,:,:] = xr.where(S_SubSoil_check [time,:,:] <= SWC_Min_SubSoil[:,:], (R_2[time-1,:,:] *effective_saturation_SubSoil*(1-effective_saturation_TopSoil)),
                                    CR[time-1,:,:])
        
        # When S_SubSoil is greater than Max it means R_1 has to be decrease based on diffrences between this greater value and maximum possible s_subsoil
        R_1[time-1,:,:] = np.around( xr.where((S_SubSoil_check [time,:,:] >= SWC_Max_SubSoil[:,:]) & (R_1[time-1,:,:] > (S_SubSoil_check [time,:,:]-SWC_Max_SubSoil[:,:])), 
                                    R_1[time-1,:,:] - (S_SubSoil_check [time,:,:]-SWC_Max_SubSoil[:,:]),
                                    R_1[time-1,:,:]), decimals= 10)
        
        # Now since R_1 has changed so again S_topSoil has to be calculated (where s_subsoil is greater than max, then use the updated R_1 to calculate S_topsoil) 
        S_TopSoil_check_2 = xr.where(S_SubSoil_check [time,:,:] >= SWC_Max_SubSoil[:,:], S_TopSoil[time-1,:,:] + Pr_smlt[time-1,:,:] - Teva[time-1,:,:]- Runoff[time-1,:,:] - R_1[time-1,:,:] + CR[time-1,:,:],
                                        S_TopSoil[time,:,:]) 
        
        
        # Since S_topSoil has increased now we have to check and see if this new S_Topsoil is greater than max S_topsoil or not
        S_TopSoil[time,:,:] = xr.where(S_TopSoil_check_2 >= SWC_Max[:,:], SWC_Max[:,:],
                                                S_TopSoil[time,:,:])
        
        # and if S-topsoil newly calculated is more than max then give the rest to runoff
        Runoff[time-1,:,:] = xr.where(S_TopSoil_check_2 >= SWC_Max[:,:], Runoff[time-1,:,:] + (S_TopSoil_check_2 - SWC_Max[:,:]),
                                        Runoff[time-1,:,:])
              
        print(time)
        
    ar_S_TopSoil = xr.DataArray(S_TopSoil[1:,:,:], dims = ['time','latitude','longitude'], coords =[Time_Steps,Latitude,Longitude],
                            attrs=dict(description="Cell-Averaged Topsoil Moisture",
                            units='mm/day'))                       
    ds_S_TopSoil = ar_S_TopSoil.to_dataset(name='S_top')                                                                              
    ds_S_TopSoil.to_netcdf(fn_S_top, mode='w', format = 'NETCDF4',engine='netcdf4')
        
    ar_Runoff= xr.DataArray(Runoff, dims = ['time','latitude','longitude'], coords =[Time_Steps,Latitude,Longitude],
                            attrs=dict(description="Runoff",
                            units='mm/day'))
    ds_Runoff = ar_Runoff.to_dataset(name='R_run')
    ds_Runoff.to_netcdf(fn_R_run, mode='w', format = 'NETCDF4',engine='netcdf4')
        
    ar_R_1 = xr.DataArray(R_1, dims = ['time','latitude','longitude'], coords =[Time_Steps,Latitude,Longitude],
                            attrs=dict(description="Subsoil_Recharge",
                            units='mm/day'))                       
    ds_R_1 = ar_R_1.to_dataset(name='R_sub')                                                                              
    ds_R_1.to_netcdf(fn_R_sub, mode='w', format = 'NETCDF4',engine='netcdf4')
        
    ## SubSoil
    ar_S_SubSoil = xr.DataArray(S_SubSoil[1:,:,:], dims = ['time','latitude','longitude'], coords =[Time_Steps,Latitude,Longitude],
                                attrs=dict(description="Cell-Averaged Subsoil Moisture",
                                units='mm/day'))                       
    ds_S_SubSoil = ar_S_SubSoil.to_dataset(name='S_sub')                                                                              
    ds_S_SubSoil.to_netcdf(fn_S_Sub, mode='w', format = 'NETCDF4',engine='netcdf4')
        
    ar_R_2 = xr.DataArray(R_2, dims = ['time','latitude','longitude'], coords =[Time_Steps,Latitude,Longitude],
                            attrs=dict(description="Groudwater_Recharge",
                            units='mm/day'))                       
    ds_R_2 = ar_R_2.to_dataset(name='R_gw')                                                                              
    ds_R_2.to_netcdf(fn_R_gw, mode='w', format = 'NETCDF4',engine='netcdf4')
        
    ar_CR = xr.DataArray(CR, dims = ['time','latitude','longitude'], coords =[Time_Steps,Latitude,Longitude],
                            attrs=dict(description="Capillary Rise from Subsoil to Topsoil",
                            units='mm/day'))                       
    ds_CR = ar_CR.to_dataset(name='CR')                                                                              
    ds_CR.to_netcdf(fn_CR, mode='w', format = 'NETCDF4',engine='netcdf4')
    
    
    ds_Runoff.close()
    ds_S_TopSoil.close()
    ds_S_SubSoil.close()
    ds_R_1.close()
    ds_R_2.close()
    ds_CR.close()
       
first_year = snakemake.params[0]
last_year = snakemake.params[1]
dir_Outputs = snakemake.params[2]

for year in range(first_year,last_year+1):
    if year == 2001:
            fn_S_TopSoil_Pre = snakemake.input[0] 
            fn_S_SubSoil_Pre = snakemake.input[1]
            fn_Sum_Pr_Smlt = snakemake.input[2]
            fn_Teva = snakemake.input[3]
    else:
        fn_S_TopSoil_Pre = dir_Outputs+'/GGR/S_top'+'/S_top_D_0101'+str(year-1)+'_3112'+str(year-1)+'_1.nc'
        fn_S_SubSoil_Pre = dir_Outputs+'/GGR/S_sub'+'/S_sub_D_0101'+str(year-1)+'_3112'+str(year-1)+'_1.nc'
        fn_Sum_Pr_Smlt = dir_Outputs+'/Meteo'+'/Sum_Pr_Smlt_mm_D_0101'+str(year)+'_3112'+str(year)+'.nc'
        fn_Teva = dir_Outputs+'/Meteo'+'/Teva_m_D_0101'+str(year)+'_3112'+str(year)+'.nc'
    
    fn_S_top = dir_Outputs+'/GGR/S_top'+'/S_top_D_0101'+str(year)+'_3112'+str(year)+'_1.nc'
    fn_R_run = dir_Outputs+'/GGR/R_run'+'/R_run_D_0101'+str(year)+'_3112'+str(year)+'_1.nc'
    fn_R_sub = dir_Outputs+'/GGR/R_sub'+'/R_sub_D_0101'+str(year)+'_3112'+str(year)+'_1.nc'
    fn_CR = dir_Outputs+'/GGR/CR'+'/CR_D_0101'+str(year)+'_3112'+str(year)+'_1.nc'
    fn_S_Sub = dir_Outputs+'/GGR/S_sub'+'/S_sub_D_0101'+str(year)+'_3112'+str(year)+'_1.nc'
    fn_R_gw = dir_Outputs+'/GGR/R_gw'+'/R_gw_D_0101'+str(year)+'_3112'+str(year)+'_1.nc'
    
    GGR(fn_S_TopSoil_Pre,fn_S_SubSoil_Pre,fn_Sum_Pr_Smlt,fn_Teva,year,fn_S_top,fn_R_run,fn_R_sub,fn_CR,fn_S_Sub,fn_R_gw)


print('Finish: Global Groundwater Recharge Model') 
