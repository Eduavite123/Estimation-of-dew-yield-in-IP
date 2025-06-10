import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import mannwhitneyu
from scipy.io import loadmat
import h5py

#--------------------------------------------DUAL PLOT--------------------------------------------------------------
#Función que genera el gráfico (mapas de la PI) de comparación lado a lado
def plot_dual(lons, lats, h_mean1, h_mean2, label, title1, title2, color, outroute):
    # Calcular límites comunes para la barra de color
    vminimo = min(np.nanmin(h_mean1), np.nanmin(h_mean2))
    vmaximo = max(np.nanmax(h_mean1), np.nanmax(h_mean2))

    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': ccrs.LambertConformal(central_longitude=3)}, figsize=(18, 8))

    # Primer gráfico
    ax1.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree())  # Zoom a península ibérica
    ax1.coastlines()
    c1 = ax1.pcolormesh(lons, lats, h_mean1, cmap=color, shading='auto', transform=ccrs.PlateCarree(), vmin=vminimo, vmax=vmaximo)
    ax1.set_title(title1)

    # Segundo gráfico
    ax2.set_extent([-10, 4.5, 35.97, 44.25], crs=ccrs.PlateCarree())  # Zoom a península ibérica
    ax2.coastlines()
    c2 = ax2.pcolormesh(lons, lats, h_mean2, cmap=color, shading='auto', transform=ccrs.PlateCarree(), vmin=vminimo, vmax=vmaximo)
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
fn='scripts\\data\\CERRA\\CERRA_daymean_total.nc'
cerraf=xr.open_dataset(fn)
lons=cerraf.coords['longitude'].values
lats=cerraf.coords['latitude'].values
cerraf.close()
lons=lons.transpose()
lats=lats.transpose()
#por si no está aplicada la máscara de tierra, la cargamos también
land_mask = np.load("scripts\\data\\PI_landmask.npy") #cargamos la máscara de tierra

land_mask_apply=True #True para aplicar, False para no

#nombres de los archivos de figuras y sus títulos y leyendas
label_comp='Producción de rocío (mm/año)' #leyenda comparación lado a lado
title_comp='Media de la producción de rocío entre 1991 y 2020' #título comparación lado a lado
out_name_comp='figures\\CERRA_WRF_dy.png' #nombre del archivo de comparación
title1='CERRA' #título del primer gráfico
title2='WRF' #título del segundo gráfico

out_name_dif='CERRA_WRF_dy_diff_cut.png' #nombre del archivo de diferencias relativas

#...........................PLOT DE COMPARACIÓN..........................
f = h5py.File('scripts\\dewyield_results\\dewyield_CERRA_data.mat','r')
data = f.get('daily_dewyield')
mod1 = np.array(data)
mod1=mod1*24*365
f = h5py.File('scripts\\dewyield_results\\dewyield_WRF_Evaluation_data_regridded.mat','r')
data = f.get('daily_dewyield')
mod2 = np.array(data)
mod2=mod2*24*365
print('Datos cargados')

#Realizamos las medias 
if land_mask_apply==True:
    mean1=np.nanmean(mod1, axis=2) * land_mask
    mean2=np.nanmean(mod2, axis=2) * land_mask
else:
    mean1=np.nanmean(mod1, axis=2)
    mean2=np.nanmean(mod2, axis=2)


#mean2=np.where(mean2>600, np.nan, mean2)
print('Medias calculadas')

#Hacemos el plot
plot_dual(lons, lats, mean1, mean2, label_comp, title1, title2, 'rainbow', out_name_comp)
print('Plot de comparación realizado')
#...........................................................................

#...........................TEST U DE MANN WHITNEY..........................
_, p_value = mannwhitneyu(mod1, mod2, alternative='two-sided', axis=2, nan_policy='omit')
print('Test de Mann Whitney realizado')

#Plot para revisar los valores de p-value en el mapa
plot(lons, lats, p_value, 'p-value', 'viridis', 'figures\\p_value_dy.png')
print('Plot de p-value realizado')

#Realizamos una máscara para celdas donde las diferencias son significativas
p_value_mask=np.where(p_value<0.05, 1, np.nan) #nivel de significancia 95%
#...........................................................................

#...........................DIFERENCIAS RELATIVAS...........................
#Cálculo de diferencias relativas
dif_rel=(mean2-mean1)/mean1*100
print('Diferencias relativas calculadas')

#Aplicamos la máscara de p-value
dif_rel=dif_rel*p_value_mask

dif_rel=np.where(dif_rel>200, np.nan, dif_rel)

#Plot de diferencias relativas
label_difrel='Diferencia relativa (%)'
outroute=f'figures\\{out_name_dif}'
plot(lons, lats, dif_rel, label_difrel, 'rainbow', outroute)
print('Plot de diferencias relativas realizado')
#...........................................................................