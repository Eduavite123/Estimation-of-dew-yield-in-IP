#Script que genera un plot para comparar la media de producción de rocío por año en WRF y CERRA 
# entre los años 1991 y 2020.
#También genera archivos netcdf simples con las medias para un mejor manejo en sus análisis.
#Incluye la interpolación de los datos de WRF a la malla de CERRA.

import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import os
from scipy.interpolate import griddata

#--------------------------------------------PLOT--------------------------------------------------------------
def plot_dual(lons, lats, h_mean1, h_mean2, label, title1, title2, general_title, outroute):
    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': ccrs.LambertConformal(central_longitude=3)}, figsize=(18, 8))
    fig.suptitle(general_title, fontsize=16, x=0.515)

    ax1.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree()) # zoom a península ibérica
    ax1.coastlines()
    c1 = ax1.pcolormesh(lons, lats, h_mean1, cmap='viridis', shading='auto', transform=ccrs.PlateCarree())
    ax1.set_title(title1)

    ax2.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree()) # zoom a península ibérica
    ax2.coastlines()
    c2 = ax2.pcolormesh(lons, lats, h_mean2, cmap='viridis', shading='auto', transform=ccrs.PlateCarree())
    ax2.set_title(title2)

    cbar = fig.colorbar(c1, ax=[ax1, ax2], orientation='horizontal', shrink=0.8, pad=0.05)
    cbar.set_label(label)
    #fig.subplots_adjust(wspace=0.05, bottom=0.25, top=0.85)  

    plt.savefig(outroute, dpi=300, bbox_inches='tight')
    plt.close(fig)
#--------------------------------------------------------------------------------------------------------------

#----------------------------------------MEAN DEW YIELD--------------------------------------------------------
def mean_dewyield(input_folder):
    #Creamos lista con los archivos
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.nc')]
    #Cargamos y concatenamos los archivos
    mean_sum=0
    count=0
    for f in files:
        with xr.open_dataset(f) as ds:
            mean_sum += ds['dew_yield'].mean(dim='time')
            count += 1
    #Calculamos la media 
    mean_dewyield = mean_sum / count
    #Pasamos de mm/h a mm/año
    return mean_dewyield*24*365
#--------------------------------------------------------------------------------------------------------------

#------------------------------------------MAIN--------------------------------------------------
#Media del periodo CERRA
input_folder='scripts\\dewyield_results\\CERRA'
mean_dewyield_C=mean_dewyield(input_folder)
#Media del periodo WRF
input_folder='scripts\\dewyield_results\\WRF'
mean_dewyield_W=mean_dewyield(input_folder)

#extraemos latitudes y longitudes que serán necesarias
#WRF
fn='scripts\\dewyield_results\\WRF\\dew_yield_WRF_1991.nc'
wrff=xr.open_dataset(fn)
lons_WRF=wrff.coords['longitude'].values
lats_WRF=wrff.coords['latitude'].values
wrff.close()
#CERRA
fn='scripts\\dewyield_results\\CERRA\\dew_yield_CERRA_1991.nc'
cerraf=xr.open_dataset(fn)
lons_C=cerraf.coords['longitude'].values
lats_C=cerraf.coords['latitude'].values
cerraf.close()
#Usaremos también la máscara de mar en la interpolación
fn='scripts\\data\\CERRA\\CERRA_IP_3h_1991.nc'
cerraf=xr.open_dataset(fn)
land_mask=cerraf.data_vars['lsm'].values[0,:,:]
cerraf.close()
land_mask=np.where(land_mask==0.0, np.nan, land_mask)

#.................Interpolación de datos WRF....................
#Creación de la matriz con dimensiones de CERRA
[nlat, nlon]=mean_dewyield_C.shape
h_WRF_regrid=np.full((nlat, nlon), np.nan)

#Máscara de valores válidos
mask=~np.isnan(mean_dewyield_W.values)

#Extraemos los datos válidos con la máscara
lat = lats_WRF[mask]
lon = lons_WRF[mask]
h_mean = mean_dewyield_W.values[mask]

#Interpolación de los datos WRF a la malla de CERRA
h_WRF_regrid = griddata((lon,lat), h_mean, (lons_C, lats_C), method='linear')

#Aplicamos la máscara de mar
mean_dewyield_W = h_WRF_regrid*land_mask
#...............................................................

#Creamos el plot
leyenda='Producción de rocío (mm/año)'
title= 'Producción de rocío recolectable entre 1991 y 2020'
outroute='figures\\CERRA_WRF_dewyield.png'
plot_dual(lons_C, lats_C, mean_dewyield_C, mean_dewyield_W, leyenda, 'CERRA', 'WRF', title, outroute)

#........................ARCHIVOS NETCDF.........................
#CERRA
# nc_CERRA= xr.Dataset(
#     coords={
#         'latitude': (['latitude', 'longitude'], lats_C),
#         'longitude': (['latitude', 'longitude'], lons_C),
#     },
#     data_vars={
#         'mean_dew_yield': (['latitude', 'longitude'], mean_dewyield_C.values)
#     }
# )
# #Falta añadir los atributos de las coordenadas y la variable
# nc_CERRA.to_netcdf(f'scripts\\dewyield_results\\mean_dy_CERRA.nc')

# #WRF
# nc_WRF= xr.Dataset(
#     coords={
#         'latitude': (['latitude', 'longitude'], lats_C),
#         'longitude': (['latitude', 'longitude'], lons_C),
#     },
#     data_vars={
#         'mean_dew_yield': (['latitude', 'longitude'], mean_dewyield_W)
#     }
# )
# #Falta añadir los atributos de las coordenadas y la variable
# nc_WRF.to_netcdf(f'scripts\\dewyield_results\\mean_dy_WRF.nc')
#..................................................................