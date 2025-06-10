import numpy as np
import h5py

f = h5py.File('scripts\\dewyield_results\\dewyield_CERRA_data.mat','r')
data = f.get('daily_dewyield')
mod1 = np.array(data)
mod1=mod1*24*365

mean1=np.nanmean(mod1, axis=2)

land_mask = np.where(mean1>0, 1, np.nan)  # Assuming land_mask is a binary mask where land is 1 and ocean is NaN

print(land_mask.shape)  # Check the shape of the land mask

np.save("scripts\\data\\PI_landmask.npy", land_mask)