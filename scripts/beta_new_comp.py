import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import mannwhitneyu
from scipy.stats import f 
from datetime import datetime

#--------------------------------------------DUAL PLOT--------------------------------------------------------------
#Función que genera el gráfico (mapas de la PI) de comparación lado a lado
def plot_dual(lons, lats, h_mean1, h_mean2, label, title1, title2, outroute):
    # Calcular límites comunes para la barra de color
    vminimo = min(np.nanmin(h_mean1), np.nanmin(h_mean2))
    vmaximo = max(np.nanmax(h_mean1), np.nanmax(h_mean2))

    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': ccrs.LambertConformal(central_longitude=3)}, figsize=(18, 8))

    # Primer gráfico
    ax1.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree())  # Zoom a península ibérica
    ax1.coastlines()
    c1 = ax1.pcolormesh(lons, lats, h_mean1, cmap='viridis', shading='auto', transform=ccrs.PlateCarree(), vmin=vminimo, vmax=vmaximo)
    ax1.set_title(title1)

    # Segundo gráfico
    ax2.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree())  # Zoom a península ibérica
    ax2.coastlines()
    c2 = ax2.pcolormesh(lons, lats, h_mean2, cmap='viridis', shading='auto', transform=ccrs.PlateCarree(), vmin=vminimo, vmax=vmaximo)
    ax2.set_title(title2)

    # Barra de color compartida
    cbar = fig.colorbar(c1, ax=[ax1, ax2], orientation='horizontal', shrink=0.8, pad=0.05)
    cbar.set_label(label)

    plt.savefig(outroute, dpi=300, bbox_inches='tight')
    plt.close(fig)
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------PLOT--------------------------------------------------------------
#Función para crear mapas individuales
def plot(lons, lats, variable, label, color, outroute):
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.LambertConformal(central_longitude=3)}, figsize=(10, 8))
    ax.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree()) #zoom a península ibérica
    ax.coastlines()

    c = ax.pcolormesh(lons, lats, variable, cmap=color, shading='auto', transform=ccrs.PlateCarree())
    cbar = plt.colorbar(c, ax=ax, orientation='horizontal', shrink=0.8, pad=0.05) 
    cbar.set_label(label)

    # ax.set_title(f'{year}')

    plt.savefig(outroute, dpi=300, bbox_inches='tight')

    #plt.show()
    plt.close(fig)
#-------------------------------------------------------------------------------------------------------------

#----------------------------------------------MAIN-----------------------------------------------------------
#extraemos latitudes y longitudes para los plots
fn='scripts\\dewyield_results\\CERRA\\dew_yield_CERRA_1991.nc'
cerraf=xr.open_dataset(fn)
lons=cerraf.coords['longitude'].values
lats=cerraf.coords['latitude'].values
cerraf.close()
#por si no está aplicada la máscara de tierra, la cargamos también
fn='scripts\\data\\land_mask.nc'
cerraf=xr.open_dataset(fn)
land_mask=cerraf['lsm'].values
cerraf.close()
land_mask=np.where(land_mask==0, np.nan, land_mask) #convertimos a NaN los valores de mar

#Rutas y datos de los archivos a usar en el script
fn1='scripts\\dewyield_results\\dewyield_CERRA_daysum.nc'
fn2='scripts\\dewyield_results\\dewyield_WRF_daysum.nc'
selname1='dew_yield'
selname2='dew_yield'

land_mask_apply=False #True para aplicar, False para no

#nombres de los archivos de figuras y sus títulos y leyendas
label_comp='Producción de rocío (mm/año)' #leyenda comparación lado a lado
title_comp='Media de la producción de rocío entre 1991 y 2020' #título comparación lado a lado
out_name_comp='figures\\CERRA_WRF_dy.png' #nombre del archivo de comparación
title1='CERRA' #título del primer gráfico
title2='WRF' #título del segundo gráfico

#...........................PLOT DE COMPARACIÓN..........................
#Cargamos los sets con datos diarios
f1=xr.open_dataset(fn1)
mod1=f1[selname1].values
f1.close()
f2=xr.open_dataset(fn2)
mod2=f2[selname2].values
f2.close()

#Realizamos las medias 
if land_mask_apply==True:
    mean1=np.nanmean(mod1, axis=0) * land_mask
    mean2=np.nanmean(mod2, axis=0) * land_mask
else:
    mean1=np.nanmean(mod1, axis=0)
    mean2=np.nanmean(mod2, axis=0)

#Hacemos el plot
plot_dual(lons, lats, mean1, mean2, label_comp, title1, title2, out_name_comp)
#...........................................................................

