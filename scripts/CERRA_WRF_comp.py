# Script con varias funciones:
# - Genera un plot que compara la media de rocío recolectable al año en el periodo 1991-2020 de los datos CERRA 
#   y el modelo WRF.
# - Calcula y genera un plot de la diferencia relativa de las medias de los modelos. Realiza también un test U 
#   de Mann-Whitney.
# - Calcula las desviaciones típicas y genera un plot del cociente de las desviaciones de ambos modelos. Realiza
#   también un test F de Snedecor.

import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import mannwhitneyu
from scipy.stats import f 
import os
from datetime import datetime

#--------------------------------------------DUAL PLOT--------------------------------------------------------------
#Función que genera el gráfico de comparación entre los dos periodos
def plot_dual(lons, lats, h_mean1, h_mean2, label, title1, title2, general_title, outroute):
    # Calcular límites comunes para la barra de color
    vminimo = min(np.nanmin(h_mean1), np.nanmin(h_mean2))
    vmaximo = max(np.nanmax(h_mean1), np.nanmax(h_mean2))

    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': ccrs.LambertConformal(central_longitude=3)}, figsize=(18, 8))
    fig.suptitle(general_title, fontsize=16, x=0.515)

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

#----------------------------------------MEAN DEW YIELD--------------------------------------------------------
#Función que calcula la media de producción de rocío por año en todo el periodo
def mean_dewyield(input_folder):
    #Creamos lista con los archivos
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.nc')]
    #Cargamos y concatenamos los archivos
    datasets = [xr.open_dataset(f) for f in files]
    combined_ds = xr.concat(datasets, dim='time')
    # Calculamos la media ignorando valores NaN
    mean_dewyield = combined_ds['dew_yield'].mean(dim='time', skipna=True)
    # Cerramos los datasets para liberar memoria
    for ds in datasets:
        ds.close()    
    #Pasamos de mm/h a mm/año
    return mean_dewyield*24*365
#--------------------------------------------------------------------------------------------------------------

#----------------------------------------VAR DEW YIELD--------------------------------------------------------
#Función que calcula la varianza de producción de rocío por año en todo el periodo
def var_dewyield(input_folder):
    #Creamos lista con los archivos
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.nc')]
    #Cargamos y concatenamos los archivos
    datasets = [xr.open_dataset(f) for f in files]
    combined_ds = xr.concat(datasets, dim='time')
    #Necesitamos el número de elementos para los grados de libertad
    n_elements = combined_ds.sizes['time']
    # Calculamos la media ignorando valores NaN
    var_dewyield = combined_ds['dew_yield'].var(dim='time', skipna=True, ddof=1)
    # Cerramos los datasets para liberar memoria
    for ds in datasets:
        ds.close()    
    #Pasamos de mm/h a mm/año
    return var_dewyield*24*365, n_elements
#--------------------------------------------------------------------------------------------------------------

#------------------------------------------MAIN--------------------------------------------------
#....................PLOT COMPARACIÓN......................
#Media del periodo CERRA
input_folder='scripts\\dewyield_results\\CERRA'
mean_dewyield_C=mean_dewyield(input_folder)
#Media del periodo WRF
input_folder='scripts\\dewyield_results\\WRF'
mean_dewyield_W=mean_dewyield(input_folder)

#extraemos latitudes y longitudes para los plots
fn='scripts\\dewyield_results\\CERRA\\dew_yield_CERRA_1991.nc'
cerraf=xr.open_dataset(fn)
lons_C=cerraf.coords['longitude'].values
lats_C=cerraf.coords['latitude'].values
cerraf.close()

#Creamos el plot
leyenda='Producción de rocío (mm/año)'
title= 'Producción de rocío recolectable entre 1991 y 2020'
outroute='figures\\CERRA_WRF_dewyield.png'
plot_dual(lons_C, lats_C, mean_dewyield_C, mean_dewyield_W, leyenda, 'CERRA', 'WRF', title, outroute)
#.........................................................

#........DIFERENCIAS RELATIVAS..............
h_difrel=(mean_dewyield_W-mean_dewyield_C)/mean_dewyield_C*100

#Plot de diferencias relativas
label='Diferencia relativa (%)'
outroute='figures\\dif_relativas.png'
plot(lons_C, lats_C, h_difrel, label, 'plasma', outroute)
#..............................................

#..........TEST U DE MANN-WHITNEY..............
# Filtrar valores NaN de h_C y h_W para el test de Mann-Whitney
valid_mask = ~np.isnan(mean_dewyield_C) & ~np.isnan(mean_dewyield_W)
h_C_valid = mean_dewyield_C.values[valid_mask]
h_W_valid = mean_dewyield_W.values[valid_mask]

# Test U de Mann-Whitney
stat, p_value = mannwhitneyu(h_C_valid, h_W_valid, alternative='two-sided')
# Imprimir resultados
print(f'Estadístico U: {stat}')
print(f'Valor p: {p_value}')

# Interpretación
alpha = 0.05 #nivel de significancia 95%
if p_value < alpha:
    print("Las distribuciones son significativamente diferentes.")
else:
    print("No hay evidencia suficiente de que sean diferentes.")
#..................................................

#...........VARIANZAS..............
#Cálculo de varianzas
#CERRA
input_folder='scripts\\dewyield_results\\CERRA'
var_dewyield_C, dfC = var_dewyield(input_folder)
#WRF
input_folder='scripts\\dewyield_results\\WRF'
var_dewyield_W, dfW = var_dewyield(input_folder)

#Calculamos el cociente de varianzas
var_coc=var_dewyield_W/var_dewyield_C

#Plot de cociente de varianzas
label='Cociente de varianzas (mm/año)'
outroute='figures\\coc_var.png'
plot(lons_C, lats_C, var_coc, label, 'cividis', outroute)
#..............................................

#..........TEST F DE SNEDECOR..............
# Grados de libertad (Nº elementos - ddof), ddof=1
df1 = dfC - 1  # Grados de libertad para CERRA
df2 = dfW - 1  # Grados de libertad para WRF

#El estadístico F se calcula como la razón de las varianzas, ya calculado.
#Calculamos el valor p para cada celda y vemos en qué regiones las diferencias en la 
#varianza son significativas.
p_value=f.cdf(var_coc, df1, df2)
label=f'Valores p menores que 0.05 (95% de confianza)'
plot(lons_C, lats_C, p_value<0.05, label, 'plasma', 'figures\\p_value.png')
#..........................................
#--------------------------------------------------------------------------------------------------