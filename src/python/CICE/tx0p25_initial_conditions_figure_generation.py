
import os
import afim
import xarray            as xr
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs



fig = plt.figure(figsize=(20,18))
ax  = fig.add_subplot(1,1,1, projection=ccrs.Robinson())
ax.add_feature(cartopy.feature.COASTLINE.with_scale(res))
plt.scatter(x=G0p25.ulon.values[1::7],
            y=G0p25.ulat.values[1::7],
            color='blue',
            s=1,
            transform=ccrs.PlateCarree())
plt.show()
