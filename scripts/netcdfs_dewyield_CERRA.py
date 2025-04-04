#Script que genera netcdf's con la producción de rocío por hora de CERRA entre los años 1991 y 2020
# en la Península Ibérica.

import xarray as xr 
import numpy as np
from numba import njit
from scipy.interpolate import griddata

#--------------------FUNCIÓN QUE CALCULA MEDIA ANUAL DE DEW YIELD P.U. TIEMPO-------------
@njit
def hpt(rh, t, h, n, V, sea, dt):
    #CÁLCULO DE T_d (Punto de Rocío) ---> Ec. Lawrence (https://www.mdpi.com/2073-4441/11/4/733)
    A=17.625 #ºC
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

    #Cuando Δh/Δt<0 corresponde a evaporación, SE DESCARTA (máscara 1)
    mask_ev=np.where(h>0,h,0)
    #Cuando velocidad_viento>4.4m/s, Δh/Δt=0 (máscara 2)
    mask_v=np.where(V<4.4,V,0)
    #Solo tierra (máscara 3)
    sea=np.where(sea==0.0, np.nan, sea)

    h_def=mask_ev*mask_v*sea #DEW YIELD PER HOUR

    #MEDIA ANUAL
    return h_def
#-----------------------------------------------------------------------------------------------

#------------------FUNCIÓN PARA DETECTAR AÑOS BISIESTOS-------------------
# En los años bisiestos la dimensión de la variable temporal cambia 
# Si queremos optimizar extrayendo variables comunes antes del loop, 
# tenemos que distinguir años bisiestos
def bis_or_not(year):
    if year==2020: #por algún motivo a CERRA le falta 31-dic-2020
        return False
    else:
        if (year % 400 == 0) and (year % 100 == 0):
            return True
        elif (year % 4 ==0) and (year % 100 != 0):
            return True
        else:
            return False
#--------------------------------------------------------------------------

#--------------------------------------MAIN----------------------------------

#------VARIABLES COMUNES------
fn=f'scripts\\data\\CERRA\\CERRA_IP_3h_1991.nc'
cerraf=xr.open_dataset(fn)
H_nob_C=cerraf.data_vars['orog'].values/1000 #orografía (km)
mask_sea_nob_C=cerraf.data_vars['lsm'].values #máscara del mar
lons_C=cerraf.coords['longitude'].values #longitud
lats_C=cerraf.coords['latitude'].values #latitud
cerraf.close()
#años bisiestos
fn=f'scripts\\data\\CERRA\\CERRA_IP_3h_1992.nc'
cerraf=xr.open_dataset(fn)
H_b_C=cerraf.data_vars['orog'].values/1000 #orografía (km) (año bisiesto)
mask_sea_b_C=cerraf.data_vars['lsm'].values #máscara del mar (año bisiesto)
cerraf.close()
#-----------------------------

for year in range(1991,2020+1):
    fn=f'scripts\\data\\CERRA\\CERRA_IP_3h_{year}.nc'
    cerraf=xr.open_dataset(fn)
    RH=cerraf.data_vars['r2'].values/100 #humedad relativa a 2m (tanto por 1)
    T=cerraf.data_vars['t2m'].values-(273.15) #temperatura aire a 2m (ºC)
    N=cerraf.data_vars['tcc'].values*8/100 #cobertura de nubes (octas)
    v=cerraf.data_vars['si10'].values #velocidad del viento a 10m (m/s)
    cerraf.close()
    #variables comunes de años bisiestos
    bis=bis_or_not(year)
    if bis==True:
        H=H_b_C
        mask_sea=mask_sea_b_C
    else:
        H=H_nob_C
        mask_sea=mask_sea_nob_C

    h_mean_CERRA=hpt(RH, T, H, N, v, mask_sea, 6)

    #Creamos el archivo NETCDF
    #variable temporal:
    if year==2020:
        start=np.datetime64(f'{year}-01-01T00:00:00')
        end=np.datetime64(f'{year}-12-30T18:00:00') #por algún motivo a CERRA le falta 31-dic-2020
        time = np.arange(start, end+np.timedelta64(6,'h'), np.timedelta64(6, 'h'))
    else:
        start=np.datetime64(f'{year}-01-01T00:00:00')
        end=np.datetime64(f'{year}-12-31T18:00:00')
        time = np.arange(start, end+np.timedelta64(6,'h'), np.timedelta64(6, 'h'))

    
    hours_since_start=(time-start).astype('timedelta64[h]').astype('float')

    nc_CERRA= xr.Dataset(
        coords={
            'latitude': (['latitude', 'longitude'], lats_C),
            'longitude': (['latitude', 'longitude'], lons_C),
            'time': (['time'], hours_since_start)
        },
        data_vars={
            'dew_yield': (['time', 'latitude', 'longitude'], h_mean_CERRA)
        }
    )
    #Falta añadir los atributos de las coordenadas y la variable

    nc_CERRA.to_netcdf(f'scripts\\dewyield_results\\CERRA\\dew_yield_CERRA_{year}.nc')