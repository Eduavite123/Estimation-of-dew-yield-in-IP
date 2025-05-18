# Script para comparar medias en un periodo de una magnitud entre dos series de datos. Las funciones son:
# - Genera un plot (mapas de la PI) que compara lado a lado la media de la magnitud.
# - Calcula y genera un plot de la diferencia relativa de las medias. Realiza también un test U 
#   de Mann-Whitney.
# - Calcula las desviaciones típicas y genera un plot del cociente de las desviaciones. Realiza
#   también un test F de Snedecor.

import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import mannwhitneyu
from scipy.stats import f 
from datetime import datetime

#--------------------------------------------DUAL PLOT--------------------------------------------------------------
#Función que genera el gráfico (mapas de la PI) de comparación lado a lado
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
#La media total del periodo se calcula con cdo (comando timmean)
#--------------------------------------------------------------------------------------------------------------

#----------------------------------------VAR DEW YIELD--------------------------------------------------------
#La varianza total del periodo se calcula con cdo (comando timvar)
#--------------------------------------------------------------------------------------------------------------

#------------------------------------------MAIN--------------------------------------------------
#extraemos latitudes y longitudes para los plots
fn='scripts\\dewyield_results\\CERRA\\dew_yield_CERRA_1991.nc'
cerraf=xr.open_dataset(fn)
lons_C=cerraf.coords['longitude'].values
lats_C=cerraf.coords['latitude'].values
cerraf.close()
#por si no está aplicada la máscara de tierra, la cargamos también
fn='scripts\\data\\land_mask.nc'
cerraf=xr.open_dataset(fn)
land_mask=cerraf['lsm'].values
cerraf.close()
land_mask=np.where(land_mask==0, np.nan, land_mask) #convertimos a NaN los valores de mar

#Selección de archivos que usaremos en este script
#dew yield
file_mdy_CERRA='scripts\\dewyield_results\\mean_dy_CERRA.nc'   #selname 'dew_yield'
file_vdy_CERRA='scripts\\dewyield_results\\var_dy_CERRA.nc'
file_mdy_WRF='scripts\\dewyield_results\\mean_dy_WRF.nc'   #selname 'dew_yield'
file_vdy_WRF='scripts\\dewyield_results\\var_dy_WRF.nc'
#precipitation
file_mpr_CERRA='scripts\\precipitacion\\mean_pr_CERRA.nc'  #selname 'tp'
file_vpr_CERRA='scripts\\precipitacion\\var_pr_CERRA.nc'
file_mpr_WRF='scripts\\precipitacion\\mean_pr_WRF.nc'    #selname 'pr'
file_vpr_WRF='scripts\\precipitacion\\var_pr_WRF.nc'

#archivos y nombres de variables a ejecutar
file1=file_mdy_CERRA #media 1
file2=file_mdy_WRF #media 2
file3=file_vdy_CERRA #varianza 1
file4=file_vdy_WRF #varianza 2
selname1='dew_yield'
selname2='dew_yield'

#nombres de los archivos de figuras y sus títulos y leyendas
label_comp='Producción de rocío (mm/año)' #leyenda comparación lado a lado
title_comp='Producción de rocío entre 1991 y 2020' #título comparación lado a lado
out_name_comp='CERRA_WRF_dy.png' #nombre del archivo de comparación

out_name_dif='dif_relativas_dy.png' #nombre del archivo de diferencias relativas

label_var='Cociente de varianzas (mm/año)' #leyenda cociente de varianzas
out_name_var='coc_var_dy.png' #nombre del archivo de cociente de varianzas

#....................PLOT COMPARACIÓN......................
#cargamos los datos
fn=file1
cerraf=xr.open_dataset(fn)
mean_C=cerraf[selname1].values * land_mask #media de CERRA (mm/año)
cerraf.close()
#######
fn=file2
wrff=xr.open_dataset(fn)
mean_W=wrff[selname2].values * land_mask #media de WRF (mm/año)
wrff.close()

mean_C=mean_C[0,:,:]
mean_W=mean_W[0,:,:]

#Creamos el plot
leyenda=label_comp
title= title_comp
outroute=f'figures\\{out_name_comp}'
plot_dual(lons_C, lats_C, mean_C, mean_W, leyenda, 'CERRA', 'WRF', title, outroute)
#.........................................................

#........DIFERENCIAS RELATIVAS..............
h_difrel=(mean_W-mean_C)/mean_C*100

#Plot de diferencias relativas
label='Diferencia relativa (%)'
outroute=f'figures\\{out_name_dif}'
plot(lons_C, lats_C, np.where(h_difrel<400, h_difrel, np.nan), label, 'rainbow', outroute)
#..............................................

#..........TEST U DE MANN-WHITNEY..............
# Filtrar valores NaN de h_C y h_W para el test de Mann-Whitney
valid_mask = ~np.isnan(mean_C) & ~np.isnan(mean_W)
h_C_valid = mean_C[valid_mask]
h_W_valid = mean_W[valid_mask]

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
#Cargamos los datos de varianzas
fn=file1
cerraf=xr.open_dataset(fn)
var_C=cerraf[selname1].values * land_mask #varianza de CERRA
cerraf.close()
########
fn=file2
wrff=xr.open_dataset(fn)
var_W=wrff[selname2].values *land_mask #varianza de WRF
wrff.close()

var_C=var_C[0,:,:]
var_W=var_W[0,:,:]

#Calculamos el cociente de varianzas
var_coc=var_W/var_C

#Plot de cociente de varianzas
label=label_var
outroute='figures\\{out_name_var}'
plot(lons_C, lats_C, np.where(var_coc<12, var_coc, np.nan), label, 'rainbow', outroute)
#..............................................

#..........TEST F DE SNEDECOR..............
# # Grados de libertad (Nº elementos - ddof), ddof=1
# df1 = dfC - 1  # Grados de libertad para CERRA
# df2 = dfW - 1  # Grados de libertad para WRF

# #El estadístico F se calcula como la razón de las varianzas, ya calculado.
# #Calculamos el valor p para cada celda y vemos en qué regiones las diferencias en la 
# #varianza son significativas.
# p_value=f.cdf(var_coc, df1, df2)
# label=f'Valores p menores que 0.05 (95% de confianza)'
# plot(lons_C, lats_C, p_value<0.05, label, 'plasma', 'figures\\p_value.png')
#..........................................
#--------------------------------------------------------------------------------------------------