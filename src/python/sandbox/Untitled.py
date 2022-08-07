#!/usr/bin/env python
# coding: utf-8

# In[36]:


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


# In[67]:


Ubase = 'http://south.seaice.dk/'
Ujpg  = '.s1a.1km.s.mos.jpg'
Uinfo = '.s.scene.info'
Dbase = os.path.join('/','Volumes','SEASONG','PHD')
Dimgs = os.path.join(Dbase,'data','satellite','S1A','mosaic')
t0    = dt.datetime.strptime('2015-01-01', '%Y-%m-%d')
tN    = dt.datetime.strptime('2021-07-31', '%Y-%m-%d')
ts    = [t0 + dt.timedelta(days=x) for x in range(0, (tN-t0).days)]
yrmo  = list(set(list(map(lambda x: os.path.join(Dimgs,x.strftime('%Y'),x.strftime('%m')) , ts))))
print(yrmo[0])
#SEARCH POLYGON
ply   = Polygon([[(73.54555, -70.44110), 
                  (78.28614, -62.76219),
                  (52.24452, -60.65423),
                  (54.16312, -69.13106),
                  (73.54555, -70.44110)]])


# In[25]:


response = requests.get(target_url)
data = response.text
df = pd.read_csv('/Users/dpath2o/Downloads/20160104060808.s1a.s.scene.info.txt',
                 sep='\t',
                 engine='python',
                 names=['col1','col2'])
print(df['col2'][5])



# In[68]:



p    = re.compile(r'-?\d+\.\d+')
data = urllib.request.urlopen('http://south.seaice.dk/2018/09/09/').read()
soup = BeautifulSoup(data)
for link in soup.findAll('a'):
    if Uinfo in link.string:
        Fin  = link.string
        Ugob = 'http://south.seaice.dk/2018/09/09/{}'.format(Fin)
        df   = pd.read_csv(Ugob,sep='\t',engine='python',names=['col1','col2'])
        tmp  = [float(i) for i in p.findall(df['col2'][5])]
        pnt  = Feature(geometry=Point((tmp[3],tmp[2])))
        print(Fin.split('.'))
        if boolean_point_in_polygon(pnt,ply):
            print('inside')
        #resp = requests.get(Ugob)
        #info = response.text
        


# In[ ]:


#shell command: convert ./2015/07/201507*.jpg -gravity center -background None -compose lighten -layers Flatten ./2015/201507_s1a.1km.s.mos.png
gog = Image(filename ='gog.png')
road = Image(filename ='rd.jpg')

for t in ts:
    Dimg  = os.path.join(Dimgs,
                         t.strftime('%Y'),
                         t.strftime('%m'))

path = "____absolute_dir_path____ (ex. /home/kim/work/)"

dirList=os.listdir(path)
for fname in dirList:
    print fname
    with Image(filename=path+fname) as img:
        print img.size
g = gog.clone()
r = road.clone()
with Drawing() as draw:
    # composite image with color_burn operator
    draw.composite(operator ='color_burn', left = 20, top = 30,
                   width = r.width, height = r.height, image = r)
    draw(g)
    g.save(filename ="colorburn.png")
    display(g)

