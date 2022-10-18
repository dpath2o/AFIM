#!/usr/bin/env python
# coding: utf-8

# In[15]:


import os
import pprint
from datetime import date
from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt


# In[16]:


pp     = pprint.PrettyPrinter(indent=4)
Dbase  = os.path.join('/','Volumes','SEASONG','PHD')
Ddata  = os.path.join(Dbase,'data','satellite','S1')
#Ddata  = os.path.join(Dbase,'data','satellite','S2')
Dsrc   = os.path.join(Dbase,'src','py','BBC_Sentinel1')
Fftprt = os.path.join(Dsrc,'Mawson.geojson')
ftprt  = geojson_to_wkt(read_geojson(Fftprt))
user   = "dpath2o"
pscd   = "bimgaN-6kixce-kixnep"
query_kwargs = {'platformname': 'Sentinel-1',
                'producttype': 'GRD',
                'date': ('2021-07-01T00:00:00.000Z','2021-08-13T23:59:00.000Z' ),
                'instrumentshortname': 'SAR-C SAR',
                'sensoroperationalmode': 'EW'}
#query_kwargs = {'platformname': 'Sentinel-2',
#                'producttype': 'S2MSI1C',
#                'date': ('NOW-14DAYS', 'NOW')}


# In[17]:


api         = SentinelAPI(user, pscd, 'https://scihub.copernicus.eu/dhus')
products    = api.query(ftprt,**query_kwargs)
products_df = api.to_dataframe(products)
for file in products:
    print('Filename {}'.format(products_df.title[file]))
    product_info = api.get_product_odata(file)
    if product_info['Online']:
        print('file is online')
    else:
        print('files is NOT online')
    api.download(file,Ddata)
    pp.pprint(file)


# In[ ]:




