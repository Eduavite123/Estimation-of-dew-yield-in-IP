#Script para calcular el "dew yield" por unidad de tiempo en PI cada 3h (1991-2020) DATOS CERRA
#(Muselli et al. -> https://www.mdpi.com/2073-4433/13/12/1974)

import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import pandas as pd

#medir tiempo ejecución (para optimización)
from datetime import datetime
startTime = datetime.now()

#--------------------FUNCIÓN QUE CALCULA MEDIA ANUAL DE DEW YIELD P.U. TIEMPO-------------
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
    return np.mean(h_def, axis=0)
#-----------------------------------------------------------------------------------------------

#--------------------------------------------PLOT--------------------------------------------------------------
def plot(lons, lats, h_mean, outroute):
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.LambertConformal(central_longitude=3)}, figsize=(10, 8))
    ax.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree()) #zoom a península ibérica
    ax.coastlines()

    c = ax.pcolormesh(lons, lats, h_mean, cmap='viridis', shading='auto', transform=ccrs.PlateCarree())
    cbar = plt.colorbar(c, ax=ax, orientation='horizontal', shrink=0.8, pad=0.05) 
    cbar.set_label('Producción de rocío por hora (mm/h)')

    ax.set_title(f'{year}')

    plt.savefig(outroute, dpi=300, bbox_inches='tight')

    #plt.show()
    plt.close(fig)
#-------------------------------------------------------------------------------------------------------------

#extraemos variables comunes
fn=f'scripts\\data\\CERRA\\CERRA_IP_3h_1991.nc'
cerraf=xr.open_dataset(fn)
H_nob=cerraf.data_vars['orog'].values/1000 #orografía (km)
mask_sea_nob=cerraf.data_vars['lsm'].values #máscara del mar
lons=cerraf.coords['longitude'].values #longitud
lats=cerraf.coords['latitude'].values #latitud
cerraf.close()
#en los años bisiestos la dimensión de la variable temporal es superior
#para optimizar extraemos variables comunes también de años bisiestos
fn=f'scripts\\data\\CERRA\\CERRA_IP_3h_1992.nc'
cerraf=xr.open_dataset(fn)
H_b=cerraf.data_vars['orog'].values/1000 #orografía (km) (año bisiesto)
mask_sea_b=cerraf.data_vars['lsm'].values #máscara del mar (año bisiesto)
cerraf.close()

#hacemos los plots de cada año
for year in range(1991,2020+1):
    #abrimos archivo netCDF y extraemos variables
    fn=f'scripts\\data\\CERRA\\CERRA_IP_3h_{year}.nc'
    cerraf=xr.open_dataset(fn)
    RH=cerraf.data_vars['r2'].values/100 #humedad relativa a 2m (tanto por 1)
    T=cerraf.data_vars['t2m'].values-(273.15) #temperatura aire a 2m (ºC)
    N=cerraf.data_vars['tcc'].values*8/100 #cobertura de nubes (octas)
    v=cerraf.data_vars['si10'].values #velocidad del viento a 10m (m/s)
    #cerramos archivo
    cerraf.close()

    #variables comunes de años bisiestos
    if year==2020: #por algún motivo a CERRA le falta 31-dic-2020
        H=H_nob ###########
        mask_sea=mask_sea_nob ##########
    else:
        if (year % 400 == 0) and (year % 100 == 0):
            H=H_b
            mask_sea=mask_sea_b
        elif (year % 4 ==0) and (year % 100 != 0):
            H=H_b
            mask_sea=mask_sea_b
        else:
            H=H_nob
            mask_sea=mask_sea_nob

    h_mean=hpt(RH, T, H, N, v, mask_sea, 6)
    route=f'figures\\CERRA\\dew_CERRA_{year}.png' #ruta de la gráfica
    plot(lons, lats, h_mean, route)

#tiempo ejecución (para optimización)
print(datetime.now() - startTime)