#Script que genera netcdf's con la producción de rocío por hora de WRF entre los años 1991 y 2020
# en la Península Ibérica.

import xarray as xr 
import numpy as np
from numba import njit
from scipy.interpolate import griddata
from datetime import datetime


#--------------------FUNCIÓN QUE CALCULA MEDIA ANUAL DE DEW YIELD P.U. TIEMPO-------------
@njit
def hpt(rh, t, h, n, V, sea, prep, dt):
    #rh en tanto por 1 (%/100)
    #t en ºC
    #h en km
    #n en octas
    #V en m/s
    #prep en mm
    #dt en horas

    #CÁLCULO DE T_d (Punto de Rocío) ---> Ec. Lawrence (https://www.mdpi.com/2073-4441/11/4/733)
    A=17.625
    B=243.04 #ºC
    num=B*(np.log(rh)+A*t/(B+t))
    den=A-np.log(rh)-A*t/(B+t)
    Td=num/den #PUNTO DE ROCÍO (ºC)

    #CÁLCULO DE RE (Energía disponible) (ec. 4 en Muselli et al.)
    RE1=1+0.204323*h-0.0238893*np.power(h,2)-(18.0132-1.04963*h+0.21891*np.power(h,2))*pow(10,-3)*Td
    RE2=(Td+273.15)/285
    RE3=1-n/8
    RE=0.37*RE1*np.power(RE2,4)*RE3

    #CÁLCULO DEL DEW YIELD POR UNIDAD DE TIEMPO
    #dt = periodo de medidas en horas
    h=dt/12*(0.06*(Td-t)+RE)

    #Cuando velocidad_viento>4.4m/s, Δh/Δt=0 (máscara 1)
    mask_v=np.where(V<4.4,V,0)
    #descartamos días con niebla, para ello consideramos que se cumplan dos de estas condiciones:
    #humedad relativa mayor o igual a 98%, diferencia entre Taire - T de rocío menores a 1 ºC y 
    #nubosidad de al menos 7 octas (máscara 2)
    mask_nubosidad=np.where(n>=7,1,0)
    mask_hum_rel=np.where(rh>=0.98,1,0)
    mask_temp=np.where((t-Td)<1,1,0)
    # Creamos una máscara que considera al menos dos de las tres condiciones
    mask_conditions = (mask_nubosidad + mask_hum_rel + mask_temp)
    mask_at_least_two = np.where(mask_conditions >= 2, np.nan, 1)

    h=mask_v*sea*mask_at_least_two*h

    #Suma diaria: hay 4 valores de h por día (cada 6h)
    for i in range(1, len(h), 4):
        h[i]=h[i]+h[i+1]+h[i+2]+h[i+3]

    #Cuando Δh/Δt<0 corresponde a evaporación, SE DESCARTA (máscara 3)
    mask_ev=np.where(h>0,h,0)
    #tenemos en cuenta solo días secos, precipitación<0.1mm (máscara 4)
    mask_pr=np.where(prep<0.1,1,0)
    #Solo tierra (máscara 5)
    sea=np.where(sea==0.0, np.nan, sea)

    h_def=mask_ev * mask_at_least_two * mask_pr * sea #DEW YIELD PER HOUR 

    return h_def * 24 * 365 #DEW YIELD PER YEAR (mm/año)
#-----------------------------------------------------------------------------------------------

#------------------FUNCIÓN PARA DETECTAR AÑOS BISIESTOS-------------------
# En los años bisiestos la dimensión de la variable temporal cambia 
# Si queremos optimizar extrayendo variables comunes antes del loop, 
# tenemos que distinguir años bisiestos
def bis_or_not(year):
    if (year % 400 == 0) and (year % 100 == 0):
        return True
    elif (year % 4 ==0) and (year % 100 != 0):
        return True
    else:
        return False
#--------------------------------------------------------------------------

#-------------------------------------MAIN--------------------------------------------- 

#------VARIABLES COMUNES------
fn=f'scripts\\data\\WRF\\Evaluation\\cloudfrac_d01_1991_time.nc'
wrff=xr.open_dataset(fn)
lons_W=wrff.coords['XLONG'].values #longitud
lats_W=wrff.coords['XLAT'].values #latitud
wrff.close()
fn=f'scripts\\data\\WRF\\hgt.nc'
wrff=xr.open_dataset(fn)
H=wrff.data_vars['HGT'].values/1000 #orografía (km)
wrff.close()
H_nob_W = np.repeat(H, 1460, axis=0)
H_b_W = np.repeat(H, 1464, axis=0)
fn=f'scripts\\data\\WRF\\landmask.nc'
wrff=xr.open_dataset(fn)
mask_sea=wrff.data_vars['LANDMASK'].values #máscara del mar
wrff.close()
mask_sea_nob_W = np.repeat(mask_sea, 1460, axis=0)
mask_sea_b_W = np.repeat(mask_sea, 1464, axis=0)
#-----------------------------

#Para la interpolación necesitamos:
fn='scripts\\data\\CERRA\\CERRA_IP_3h_1991.nc'
cerraf=xr.open_dataset(fn)
cerranb=cerraf.data_vars['lsm'].values
lons_C=cerraf.coords['longitude'].values
lats_C=cerraf.coords['latitude'].values
cerraf.close()
fn='scripts\\data\\CERRA\\CERRA_IP_3h_1992.nc'
cerraf=xr.open_dataset(fn)
cerrab=cerraf.data_vars['lsm'].values
cerraf.close()

for year in range(1991,2020+1):
    startTime = datetime.now() #para medir el tiempo de ejecución

    fn=f'scripts\\data\\WRF\\Evaluation\\cloudfrac_d01_{year}_time.nc'
    cloudfracf=xr.open_dataset(fn)
    N=cloudfracf.data_vars['CLDFRA'].values*8 #cobertura de nubes (octas)
    cloudfracf.close()
    fn=f'scripts\\data\\WRF\\Evaluation\\hur_d01_{year}_time.nc'
    hurf=xr.open_dataset(fn)
    RH=hurf.data_vars['hurs'].values #humedad relativa a 2m (tanto por 1)
    hurf.close()
    fn=f'scripts\\data\\WRF\\Evaluation\\t2_d01_{year}_time.nc'
    t2f=xr.open_dataset(fn)
    T=t2f.data_vars['T2'].values-(273.15) #temperatura aire a 2m (ºC)
    t2f.close()
    fn=f'scripts\\data\\WRF\\Evaluation\\windspeed_d01_{year}_time.nc'
    windspeedf=xr.open_dataset(fn)
    v=windspeedf.data_vars['W10'].values #velocidad del viento a 10m (m/s)
    windspeedf.close()
    fn=f'scripts\\precipitacion\\WRF\\pr_WRF_{year}.nc'
    prf=xr.open_dataset(fn)
    prep=prf.data_vars['precip'].values #precipitación (mm)
    prf.close()
    #variables comunes de años bisiestos
    bis=bis_or_not(year)
    if bis==True:
        H=H_b_W
        mask_sea=mask_sea_b_W
    else:
        H=H_nob_W
        mask_sea=mask_sea_nob_W
        
    h_WRF=hpt(RH, T, H, N, v, mask_sea, 6)

    #.................Interpolación de datos WRF....................
    #Creación de la matriz con dimensiones de CERRA
    if bis==True:
        [ntime, nlat, nlon]=cerrab.shape
        land_mask=cerrab
    else:
        [ntime, nlat, nlon]=cerranb.shape
        land_mask=cerranb
    h_WRF_regrid=np.full((ntime, nlat, nlon), np.nan)
    #De paso ya ajustamos también la máscara de mar
    land_mask=np.where(land_mask==0.0, np.nan, land_mask)

    #Máscara de valores válidos
    mask=~np.isnan(h_WRF[1,:,:])

    #Extraemos los datos válidos con la máscara
    lat = lats_W[mask]
    lon = lons_W[mask]
    h_mean = h_WRF[:,mask]

    #Interpolación de los datos WRF a la malla de CERRA
    for i in range(ntime):
        h_WRF_regrid[i,:,:] = griddata((lon,lat), h_mean[i,:], (lons_C, lats_C), method='linear')

    #Aplicamos la máscara de mar
    h_WRF = h_WRF_regrid*land_mask

    print(datetime.now() - startTime)
    #...............................................................

    #Creamos el archivo NETCDF
    #variable temporal:
    start=np.datetime64(f'{year}-01-01T00:00:00')
    end=np.datetime64(f'{year}-12-31T00:00:00')
    time = np.arange(start, end, np.timedelta64(24, 'h'))
    
    hours_since_start=(time-start).astype('timedelta64[h]').astype('float')

    nc_WRF= xr.Dataset(
        coords={
            'latitude': (['latitude', 'longitude'], lats_C),
            'longitude': (['latitude', 'longitude'], lons_C),
            'time': (['time'], hours_since_start)
        },
        data_vars={
            'dew_yield': (['time', 'latitude', 'longitude'], h_WRF)
        }
    )
    # Añadimos atributos a la coordenada temporal
    nc_WRF['time'].attrs['units'] = f'hours since {year}-01-01 00:00:00'
    nc_WRF['time'].attrs['calendar'] = 'gregorian'

    nc_WRF.to_netcdf(f'scripts\\dewyield_results\\WRF\\dew_yield_WRF_{year}.nc')

    print(f'Año {year} terminado')