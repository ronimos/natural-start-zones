# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import pdal
import numpy as np
import pandas as pd
import json
import tqdm
import matplotlib.pyplot as plt

DATA_PATH = '../data'

in_file = os.path.join(DATA_PATH, 'sisters_crop.laz')
#outfile = os.path.join(DATA_PATH, 'sisters_crop_no_trees.laz')

# itree_filter = {
#     "pipeline":[
#         {
#             "type": "readers.las",
#             "filename":in_file
#         },
#         {
#             "type":"filters.hag_delaunay"
#         },
#         {
#             "type":"filters.sort",
#             "dimension":"HeightAboveGround",
#             "order":"DESC"
#         },
#         {
#             "type":"filters.litree",
#             "min_points":10,
#             "min_height":3.0,
#             "radius":100.0
#         },
#         {
#             "type":"writers.las",
#             "filename":outfile,
#             "minor_version":1.4,
#             "extra_dims":"all"
#         }    
#         ]

# }

# json.dump(itree_filter, open('itree_filter_pipeline.json', 'w'))

#r = pdal.Pipeline(json.dumps(itree_filter))


reader_pipeline = {
    "pipeline":[
        {
            "type": "readers.las",
            "filename":in_file
        },
    ]

}


r = pdal.Pipeline(json.dumps(reader_pipeline))

count = r.execute()
arrays = r.arrays
    
df = pd.DataFrame(arrays[0])[['X','Y','Z']]
print(df.shape)
x_min = 2889242.0
x_max = 2891136.0
y_min = 1670983.0
y_max = 1672982.0 

#df = df.loc[df['Y'].between(1669000, 1672900), :]
df_s = df.query(f'X.between({x_min}, {x_max}) & Y.between({y_min}, {y_max})')
df_s = df_s.loc[df_s.index[::50], :]
print(df_s.shape)

df = pd.DataFrame(arrays[0])[['X','Y','Z']]

x_min = 2889500.0
x_max = 2890610.0
y_min = 1671500.0  
y_max = 1672650.0
df_s34 = df.query(f'X.between({x_min}, {x_max}) & Y.between({y_min}, {y_max})')




df_s = df_s34.loc[df_s34.index[::5], :].dropna()

ax = plt.figure(figsize=(22,22)).add_subplot(projection='3d')

ax.plot_trisurf(df_s['X'], df_s['Y'], df_s['Z'], linewidth=0.2, antialiased=True)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.view_init(elev=90, azim=90)
plt.show()

