#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install rasterio --user')


# In[2]:


import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import glob
import os
import zipfile
import tempfile


# In[3]:


import os
print(os.getcwd())
os.chdir('/scratch/cloudstor/BBC_Sentinel1')


# In[3]:


# this assumes that the downloaded files are unzipped and in this PRODUCT directory, change the path where needed

filenames = glob.glob("s1_data_mawson_roughness/*.zip")
dest_path="s1_data_mawson_roughness_processed"
from rasterio.warp import calculate_default_transform, reproject, Resampling

dst_crs = 'EPSG:3031'
for file in filenames:
    path_to_zip_file, fname = os.path.split(file)
    dst_name=os.path.splitext(fname)[0]
    print(dst_name)
    #unzip the original file
    with tempfile.TemporaryDirectory() as tempdir:
        with zipfile.ZipFile(file, 'r') as zip_ref:
            zip_ref.extractall(tempdir)  
        with rasterio.open(os.path.join(tempdir,"{}.SAFE".format(dst_name))) as src:
            gcps, gcp_crs = src.gcps

            transform,width,height = calculate_default_transform(gcp_crs, dst_crs, src.width, src.height,  gcps=gcps)
            kwargs = src.meta.copy()
            kwargs.update({
                'crs': dst_crs,
                'transform': transform,
                'width':width,
                'height':height,
                'driver':'GTIFF',
                'compress': 'lzw',
                'NBITS':10

            })
            print(kwargs)

            with rasterio.open("{}/{}.tif".format(dest_path,dst_name), 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.nearest)


   


# In[ ]:




