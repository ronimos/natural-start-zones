# -*- coding: utf-8 -*-
"""
Created on Thu May 23 09:35:39 2024

@author: Avalanche
"""

import pdal

DATA_PATH = '../data'
files =  DATA_PATH + '/*.laz'
outfile = DATA_PATH + '/Sisters.laz'

print("PDAL pipeline")

pipeline = (
    pdal.Reader.las(filename=files)
    | pdal.Filter.merge()
    | pdal.Writer.las(outfile)
)

print("Executing pipeline to merge files inside directory")
pipeline.execute()


import os
import pdal
import numpy as np
import pandas as pd
import json
import tqdm
import matplotlib.pyplot as plt

DATA_PATH = '../data'

x_min = 2889242.0
x_max = 2891136.0
y_min = 1670983.0
y_max = 1672982.0

in_file = "C:/Users/Avalanche/Documents/projects/Summer2024/Startzones/data/Sisters.laz"
out_file = "C:/Users/Avalanche/Documents/projects/Summer2024/Startzones/data/sisters_crop.laz"
crop_pipeline = (
    pdal.Reader.las(filename=in_file)
    | pdal.Filter.crop(bounds=f"([{x_min},{x_max}],[{y_min},{y_max}])")
    | pdal.Writer.las(out_file)
    )
crop_pipeline.execute()






siters_bounds = dict(x_min = 2889242.0,
                     x_max = 2891136.0,
                     y_min = 1670983.0,
                     y_max = 1672982.0,
                     )
sisters_3_4_bounds = dict(x_min = 2889242.0,
                          x_max = 2890060.0,
                          y_min = 1670983.0,  
                          y_max = 1671660.0,
                          )