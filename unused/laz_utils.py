# -*- coding: utf-8 -*-
"""
Created on Thu May 23 09:35:39 2024

@author: Avalanche
"""

import os
import pdal
import numpy as np
import pandas as pd
import json
import tqdm
import matplotlib.pyplot as plt
import cv2

DATA_PATH = '../data'

def crop_laz(bounds,
             in_file,
             outfile):
    
    in_file = os.path.join(DATA_PATH, in_file)
    outfile = os.path.join(DATA_PATH, outfile)
    bound_str = f"([{bounds['x_min']},{bounds['x_max']}],[{bounds['y_min']},{bounds['y_max']}])"
    crop_pipeline = (
        pdal.Reader.las(filename=in_file)
        | pdal.Filter.crop(bounds=bound_str)
        | pdal.Writer.las(outfile)
        )
    crop_pipeline.execute()


def merge_files(outfile):

    lazfiles =  os.path.join(DATA_PATH,'Star_mtn/*.laz')
    outfile = os.path.join(DATA_PATH,'/Star_mtn/Star_mtn.laz')


    pipeline = (
        pdal.Reader.las(filename=lazfiles)
        | pdal.Filter.merge()
        | pdal.Writer.las(outfile)
    )

    print("Executing pipeline to merge files inside directory")
    pipeline.execute()


def laz_2_tif(in_file,
              outfile):
    
    conver_pipeline = {
        "pipeline":[
            f"{in_file}",
            {
                "type":"writers.gdal",
                f"filename":"{outfile}",
                "resolution":2.0,
                "output_type":"max"
            }
            ]
        }
    r = pdal.Pipeline(json.dumps(conver_pipeline))
    r.execute()

def gen_mask(tif_file):
    pass
    


import os
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib.path import Path

DATA_PATH = r'C:\Users\Avalanche\Documents\projects\Summer2024\Startzones\data'

#global drawing, mask_img, xi, yi, tpoly, mask

class MaskGennerator():
    
    def __init__(self, img_path):
        
        self.drawing = False
        self.img = cv2.imread(img_path)
        self.mask_img = self.img.copy()
        self.xi, self.yi = None, None
        self.tpoly = []
        self.mask_poly = []

    def deaw_mask(self, event, x, y, flags, params):
        """
        A mouse callback function that draws a rectangle on a video frame
        and save its corners as ROI for future analasys
        --------------------------------------------------------------
        """
    
        if event == cv2.EVENT_LBUTTONDOWN:
            self.tpoly.append((x, y))
            self.ix = x
            self.iy = y
            self.drawing = True
            
        if event == cv2.EVENT_MOUSEMOVE:
            if self.drawing==True:
                self.mask_img = self.img.copy()
                _poly = self.tpoly + [(x, y), self.tpoly[0]]
                for p1, p2 in zip(_poly[:-1], _poly[1:]):
                    cv2.line(self.mask_img, p1, p2, (255, 0, 0), 2)
                    
        if event == cv2.EVENT_LBUTTONDBLCLK:
            self.drawing = False
            self.mask_poly = self.tpoly[:]
            self.tpoly = []



    def draw_olygon_mask(self):
        
        self.mask_poly = []
        window_name = 'Draw Mask'
        cv2.namedWindow(window_name)
        cv2.setMouseCallback(window_name, self.deaw_mask)
        self.drawing = False
        self.tpoly = []
        
        while True:
            
            cv2.imshow(window_name, self.mask_img)
            k = cv2.waitKey(30) & 0xFF
            if k == 27:
                cv2.destroyAllWindows()
                break
        return self.mask_poly

    def gen_mask_from_polygon(self,
                              mask_path, 
                              show_mask=True):
        
        heigth, width = self.img.shape[:2]
        x, y = np.meshgrid(range(width), range(heigth))
        x, y = x.flatten(), y.flatten()
        
        points = np.vstack((x,y)).T
        path = Path(self.mask_poly)
        grid = path.contains_points(points)
        mask = grid.reshape((heigth, width))
        if show_mask:
            plt.imshow(mask)
        np.save(mask_path, mask)
        self.mask = mask
mask_g = MaskGennerator(os.path.join(DATA_PATH, 'Sisters_slope.png'))
mask_g.draw_olygon_mask()

mask_path = os.path.join(DATA_PATH, 'masks','sister_2_mask')
mask_g.gen_mask_from_polygon(mask_path)



if __name__ == "__main__":
    
    in_file = os.path.join(DATA_PATH, "sisters_crop.laz")
    outfile = os.path.join(DATA_PATH, "sisters_crop.tif")
    laz_2_tif(in_file, outfile)
    
    

    siters_bounds = dict(x_min = 2889242.0,
                         x_max = 2891136.0,
                         y_min = 1670983.0,
                         y_max = 1672982.0,
                         )
    
    in_file = os.path.join(DATA_PATH, "Sisters.laz")
    outfile = os.path.join(DATA_PATH, "sisters_crop.laz")
    crop_laz(siters_bounds, in_file, outfile)
    
    sisters_3_4_bounds = {'x_min': 2889800.0,
                          'x_max': 2890610.0,
                          'y_min': 1671500.0,
                          'y_max': 1672650.0
                          }

    in_file = os.path.join(DATA_PATH, "Sisters.laz")
    outfile = os.path.join(DATA_PATH, "sisters_234_crop.laz")
    crop_laz(sisters_3_4_bounds, in_file, outfile)
    