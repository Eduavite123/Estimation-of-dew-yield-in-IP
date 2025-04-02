#Script que genera un plot para comparar la producción de rocío por hora en WRF y CERRA entre los años 1991 y 2020.

import xarray as xr 
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs

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

