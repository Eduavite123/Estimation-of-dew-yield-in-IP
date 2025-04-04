#Script para el estudio de los campos medios y de campo de varianzas de la producción de rocío
#recolectable en el periodo 1991-2020 en la Península Ibérica.

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import mannwhitneyu

#--------------------------------------------PLOT--------------------------------------------------------------
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

#------------------------------------------------MAIN--------------------------------------------------
#Cargamos los datos 
#CERRA
fn=f'scripts\\dewyield_results\\mean_dy_CERRA.nc'
cerraf=xr.open_dataset(fn)
h_C=cerraf['mean_dew_yield'].values
lat=cerraf['latitude'].values
lon=cerraf['longitude'].values
cerraf.close()
#WRF
fn=f'scripts\\dewyield_results\\mean_dy_WRF.nc'
wrff=xr.open_dataset(fn)
h_W=wrff['mean_dew_yield'].values
wrff.close()

#Calculamos la diferencia relativa
h_difrel=(h_W-h_C)/h_C*100
#Hay valores en Francia que son muy ruidosos e innecesarios, eliminamos dichos valores
h_difrel=np.where(h_difrel>250, np.nan, h_difrel)

print('Minimo:', np.nanmin(h_difrel))
print('Maximo:', np.nanmax(h_difrel))

#Plot de diferencias relativas
label='Diferencia relativa (%)'
outroute='figures\\dif_relativas.png'
plot(lon, lat, h_difrel, label, 'plasma', outroute)

# Filtrar valores NaN de h_C y h_W
valid_mask = ~np.isnan(h_C) & ~np.isnan(h_W)
h_C_valid = h_C[valid_mask]
h_W_valid = h_W[valid_mask]

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