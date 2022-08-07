#!/usr/bin/env python
# coding: utf-8

# ## LOAD LIBRARIES

# In[1]:


import os
import re
import datetime as dt
import requests
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from wand.image import Image, COMPOSITE_OPERATORS
from wand.drawing import Drawing
from wand.display import display
from turfpy.measurement import boolean_point_in_polygon
from geojson import Point, Polygon, Feature
from bs4 import BeautifulSoup
import urllib


# ## Parameters

# In[7]:


Ubase = 'http://south.seaice.dk/'
SatTy = ('s1a','s1b','s1ab','amsr')
Scene = ('mos','ice','scene')
ImgTy = ('jpg','gif')
Uimgs = ('.s1a.1km.s.mos.jpg',
         '.s1b.1km.s.mos.jpg',
         '.s1ab.1km.s.mos.jpg',
         '.amsr2.s.ice.gif')
Uinfo = '.s.scene.info'
Dbase = os.path.join('/','Volumes','SEASONG','PHD')
Dimgs = os.path.join(Dbase,'data','satellite','sentinel')
t0    = dt.datetime.strptime('2021-03-18', '%Y-%m-%d')
tN    = dt.datetime.strptime('2021-08-14', '%Y-%m-%d')
ts    = [t0 + dt.timedelta(days=x) for x in range(0, (tN-t0).days)]
yrmo  = set(list(map(lambda x: os.path.join(Dimgs,x.strftime('%Y'),x.strftime('%m')) , ts)))
#SEARCH POLYGON
ply   = Polygon([[(73.54555, -70.44110), 
                  (78.28614, -62.76219),
                  (52.24452, -60.65423),
                  (54.16312, -69.13106),
                  (73.54555, -70.44110)]])


# # Search for desired images and download

# In[8]:


#URL (example): http://south.seaice.dk/2016/01/04/20160104074753.s1a.s.scene.info
#
#FILE CONTENTS (example)
#
#IMAGE       "2016/01/04/20160104074753.s1a.s.scene.jpg"
#TITLE       "s1a.s.scene 07:47:53"
#TIMESTAMP	"2016-01-04"
#PIXELRES	40 40 300 300
#PROJECTION	401 0.000000 -90.000000 -70.000000
#FIXPOINT	1.0 1.0 -69.102418 -168.765950
#
#We want the fix-point coordinates, which will use as our point inside a polygon test
#which is to say that if the fix-point location is within a user-defined polygon then 
#the corresponding satellite image will contain a portion of that image inside
#the user-defined polygon. Not the most surgical approach, but it'll at least get the
#job most of the way there
#
p = re.compile(r'-?\d+\.\d+')

def download_image(URL,FILE_LOCAL):
    req = requests.get(URL,stream=True)
    if req.status_code == 200:
        req.raw.decode_content = True
        with open(FILE_LOCAL,'wb') as f:
            shutil.copyfileobj(req.raw, f)
            print('downloaded {} to local {}'.format(URL,FILE_LOCAL))
    else:
        print('NOT downloaded {}'.format(URL))

for t in ts:
    Udir = os.path.join(Ubase,
                         t.strftime('%Y'),
                         t.strftime('%m'),
                         t.strftime('%d'))
    Dls  = urllib.request.urlopen(Udir).read()
    soup = BeautifulSoup(Dls)
    print('\nSearching {} for possible images'.format(Udir))
    for link in soup.findAll('a'):
        Fdown  = ''
        Fin    = link.string
        Fparts = Fin.split('.')
        if len(Fparts)<=4: continue
        Ufile = os.path.join(Udir,Fin)
        if ('amsr2' in Fparts[1]) and (Fparts[2]=='s') and (Fparts[3]=='ice') and (Fparts[4]=='gif'):
            #print('AMSR image found {}'.format(Fin))
            Dimg  = os.path.join(Dimgs,'amsr',
                             t.strftime('%Y'),
                             t.strftime('%m'))
            Path(Dimg).mkdir(parents=True, exist_ok=True)
            Fdown = os.path.join(Dimg,Fin)
            if os.path.isfile(Fdown): 
                print('{} already exists skipping'.format(Fdown))
                continue
            else:
                download_image(Ufile,Fdown)
        elif ('s1' in Fparts[1]) and (Fparts[2]=='1km') and (Fparts[3]=='s') and (Fparts[4]=='mos') and (Fparts[5]=='jpg'):
            #print('Sentinel image found {}'.format(Fin))
            Dimg  = os.path.join(Dimgs,'mosaic',
                             t.strftime('%Y'),
                             t.strftime('%m'))
            Path(Dimg).mkdir(parents=True, exist_ok=True)
            Fdown = os.path.join(Dimg,Fin)
            if os.path.isfile(Fdown): 
                print('{} already exists skipping'.format(Fdown))
                continue
            else:
                download_image(Ufile,Fdown)
        elif Uinfo in Fin:
            #print('Info file found {}'.format(Fin))
            Dimg  = os.path.join(Dimgs,'MawsonCoast',
                                 t.strftime('%Y'),
                                 t.strftime('%m'))
            Path(Dimg).mkdir(parents=True, exist_ok=True)
            Ugob = '{}/{}'.format(Udir,Fin)
            df   = pd.read_csv(Ugob,sep='\t',engine='python',names=['col1','col2'])
            tmp  = [float(i) for i in p.findall(df['col2'][5])]
            pnt  = Feature(geometry=Point((tmp[3],tmp[2])))
            if boolean_point_in_polygon(pnt,ply):
                Fjpg   = '{}.{}.{}.{}.jpg'.format(Fparts[0],Fparts[1],Fparts[2],Fparts[3])
                Fdown1 = os.path.join(Dimg,Fjpg)
                Fdown2 = os.path.join(Dimg,Fin)
                Udown1 = os.path.join(Ubase,
                                      t.strftime('%Y'),
                                      t.strftime('%m'),
                                      t.strftime('%d'),
                                      Fjpg)
                Udown2 = os.path.join(Ubase,
                                      t.strftime('%Y'),
                                      t.strftime('%m'),
                                      t.strftime('%d'),
                                      Fin)
                if (os.path.isfile(Fdown1) and os.path.isfile(Fdown2)): 
                    print('{} and {} already exists, skipping'.format(Fdown1,Fdown2))
                    continue
                download_image(Udown1,Fdown1)
                download_image(Udown2,Fdown2)


# In[ ]:




