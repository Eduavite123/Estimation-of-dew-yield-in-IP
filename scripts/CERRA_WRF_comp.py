#Script que genera un plot para comparar la producción de rocío por hora en WRF y CERRA entre los años 1991 y 2020.

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
def plot_dual(lons_1, lats_1, lons_2, lats_2, h_mean1, h_mean2, title1, title2, general_title, outroute):
    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': ccrs.LambertConformal(central_longitude=3)}, figsize=(18, 8))
    fig.suptitle(general_title, fontsize=16, x=0.515)

    ax1.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree()) # zoom a península ibérica
    ax1.coastlines()
    c1 = ax1.pcolormesh(lons_1, lats_1, h_mean1, cmap='viridis', shading='auto', transform=ccrs.PlateCarree())
    ax1.set_title(title1)

    ax2.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree()) # zoom a península ibérica
    ax2.coastlines()
    c2 = ax2.pcolormesh(lons_2, lats_2, h_mean2, cmap='viridis', shading='auto', transform=ccrs.PlateCarree())
    ax2.set_title(title2)

    cbar = fig.colorbar(c1, ax=[ax1, ax2], orientation='horizontal', shrink=0.8, pad=0.05)
    cbar.set_label('Producción de rocío por hora (mm/h)')
    #fig.subplots_adjust(wspace=0.05, bottom=0.25, top=0.85)  

    plt.savefig(outroute, dpi=300, bbox_inches='tight')
    plt.close(fig)
#-------------------------------------------------------------------------------------------------------------

#------VARIABLES COMUNES------
#CERRA
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

#WRF
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

#hacemos los plots de cada año
for year in range(1991,2020+1):
    #CERRA
    fn=f'scripts\\data\\CERRA\\CERRA_IP_3h_{year}.nc'
    cerraf=xr.open_dataset(fn)
    RH=cerraf.data_vars['r2'].values/100 #humedad relativa a 2m (tanto por 1)
    T=cerraf.data_vars['t2m'].values-(273.15) #temperatura aire a 2m (ºC)
    N=cerraf.data_vars['tcc'].values*8/100 #cobertura de nubes (octas)
    v=cerraf.data_vars['si10'].values #velocidad del viento a 10m (m/s)
    cerraf.close()
    #variables comunes de años bisiestos
    if year==2020: #por algún motivo a CERRA le falta 31-dic-2020
        H=H_nob_C ###########
        mask_sea=mask_sea_nob_C ##########
    else:
        if (year % 400 == 0) and (year % 100 == 0):
            H=H_b_C
            mask_sea=mask_sea_b_C
        elif (year % 4 ==0) and (year % 100 != 0):
            H=H_b_C
            mask_sea=mask_sea_b_C
        else:
            H=H_nob_C
            mask_sea=mask_sea_nob_C
    h_mean_CERRA=hpt(RH, T, H, N, v, mask_sea, 6)

    #WRF
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
    #variables comunes de años bisiestos
    if (year % 400 == 0) and (year % 100 == 0):
        H=H_b_W
        mask_sea=mask_sea_b_W
    elif (year % 4 ==0) and (year % 100 != 0):
        H=H_b_W
        mask_sea=mask_sea_b_W
    else:
        H=H_nob_W
        mask_sea=mask_sea_nob_W
    h_mean_WRF=hpt(RH, T, H, N, v, mask_sea, 6)

    plot_dual(lons_C, lats_C, lons_W, lats_W, h_mean_CERRA, h_mean_WRF, 'CERRA', 'WRF', f'{year}', f'figures\\comp\\dew_CERRA_WRF_{year}.png')

#tiempo ejecución (para optimización)
print(datetime.now() - startTime)