import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import mannwhitneyu
from scipy.io import loadmat
import h5py

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
fn='scripts\\data\\CERRA\\CERRA_IP_3h_1991.nc'
cerraf=xr.open_dataset(fn)
lons=cerraf.coords['longitude'].values
lats=cerraf.coords['latitude'].values
cerraf.close()
lats=lats.transpose() #los archivos mat tienen las coordenadas al revés
lons=lons.transpose()
#por si no está aplicada la máscara de tierra, la cargamos también
fn='scripts\\data\\land_mask.nc'
cerraf=xr.open_dataset(fn)
land_mask=cerraf['lsm'].values
cerraf.close()
land_mask=np.where(land_mask==0, np.nan, land_mask) #convertimos a NaN los valores de mar
land_mask=land_mask[0,:,:].transpose()

#Rutas y datos de los archivos a usar en el script
# fn1='scripts\\dewyield_results\\dewyield_CERRA_daysum.nc'
# fn2='scripts\\dewyield_results\\dewyield_WRF_daysum.nc'
# selname1='dew_yield'
# selname2='dew_yield'

land_mask_apply=True #True para aplicar, False para no

#nombres de los archivos de figuras y sus títulos y leyendas
label_comp='Producción de rocío (mm/año)' #leyenda comparación lado a lado
title_comp='Media de la producción de rocío entre 1991 y 2020' #título comparación lado a lado
out_name_comp='figures\\CERRA_WRF_dy.png' #nombre del archivo de comparación
title1='CERRA' #título del primer gráfico
title2='WRF' #título del segundo gráfico

out_name_dif='CERRA_WRF_dy_diff.png' #nombre del archivo de diferencias relativas

#...........................PLOT DE COMPARACIÓN..........................
#Cargamos los sets con datos diarios
# f1=xr.open_dataset(fn1)
# mod1=f1[selname1].values
# f1.close()
# f2=xr.open_dataset(fn2)
# mod2=f2[selname2].values
# f2.close()

f = h5py.File('scripts\\dewyield_results\\dewyield_CERRA_data.mat','r')
data = f.get('daily_dewyield')
mod1 = np.array(data)
f = h5py.File('scripts\\dewyield_results\\dewyield_WRF_Evaluation_data_regridded.mat','r')
data = f.get('daily_dewyield')
mod2 = np.array(data)
print('Datos cargados')

#Realizamos las medias 
if land_mask_apply==True:
    mean1=np.nanmean(mod1, axis=2) * land_mask
    mean2=np.nanmean(mod2, axis=2) * land_mask
else:
    mean1=np.nanmean(mod1, axis=2)
    mean2=np.nanmean(mod2, axis=2)

mean2=np.where(mean2>0.05, np.nan, mean2)
print('Medias calculadas')

#Hacemos el plot
plot_dual(lons, lats, mean1, mean2, label_comp, title1, title2, out_name_comp)
print('Plot de comparación realizado')
#...........................................................................

#...........................TEST U DE MANN WHITNEY..........................
_, p_value = mannwhitneyu(mod1, mod2, alternative='two-sided', axis=2, nan_policy='omit')
print(p_value.shape)
print('Test de Mann Whitney realizado')

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

dif_rel=np.where(dif_rel>600, np.nan, dif_rel)

#Plot de diferencias relativas
label_difrel='Diferencia relativa (%)'
outroute=f'figures\\{out_name_dif}'
plot(lons, lats, dif_rel, label_difrel, 'rainbow', outroute)
print('Plot de diferencias relativas realizado')
#...........................................................................