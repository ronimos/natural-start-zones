# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:52:22 2024

@author: Avalanche
"""
import os
from glob import glob
import numpy as np
from osgeo import gdal
import richdem as rd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter
from tqdm import tqdm
from scipy.stats import kurtosis#, skew
import warnings
warnings.filterwarnings("ignore")


def get_polygons_from_kml(path):
    
    import xml.etree.ElementTree as ET
    root = ET.parse(path).getroot()
    p_cords = [c.text.strip() for c in list(root.iter()) if 'coordinates' in c.tag]
    polys = []
    for p_c in p_cords:
        p_c = [c.split(',') for c in p_c.split(' ')]
        poly = np.array(p_c)[:, :-1].astype(float)
        polys.append(poly)


def load_tif_as_Array(path='../data/Sisters/Sisters.tif'):
    
    gdal.UseExceptions()
    ds = gdal.Open(path)
    
    band = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    return band


def crop_tif(in_path,
             out_path,
             polygon_path=None,
             minx=None, maxx=None,
             miny=None, maxy=None):
    
    gdal.UseExceptions()
    assert polygon_path is not None or minx is not None, 'Geometry is missing...' 
    
    clip =  gdal.Warp(destNameOrDestDS = out_path, 
                      srcDSOrSrcDSTab  = in_path,
                      cutlineDSName    = polygon_path, 
                      cropToCutline    = True,
                      copyMetadata     = True,
                      dstNodata        = np.nan)
    
    dem = clip.GetRasterBand(1).ReadAsArray()
    clip = None
    return dem


def show_array(a, cmap='winter', title=''):
    
    fig, ax = plt.subplots(figsize=(12,12))
    im = ax.imshow(a, cmap=cmap)
    plt.title(title)
    fig.colorbar(im, ax=ax)
    plt.show()
    
def plot_3d(a, cmap='winter_r', title='', **kwargs):
    
    view_elev = kwargs.get('view elev', 20)
    view_azi = kwargs.get('view azi', 60)
    ax = plt.figure(figsize=(22,22)).add_subplot(projection='3d')
    height, width = a.shape
    x = np.arange(width)
    y = np.arange(height)
    x, y = np.meshgrid(x, y)
    x, y = x.flatten(), y.flatten()
    ax.plot_trisurf(x, y, a.flatten(), linewidth=0.2, antialiased=True, cmap='winter_r')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.view_init(elev=view_elev, azim=view_azi)
    ax.grid(False)
    plt.title(title)
    plt.show()
    
   
def get_slope_attributes(dem_path, *args):
    """
    This function gets a path to a DEM and a list of sloe attributes to 
    generate. It saves a tig files with these attriutes and returns a dict 
    where the keys are the paths to the tif files.  
    
    Parameters
    ----------
    dem_path : TYPE
        path to the dem to generate the attributs from.
    *args : list os str
        A list of the variables to gerate.
        Possibole options are: slope_riserun, slope_percentage, slope_degrees, 
                               slope_radians, aspect, curvature, 
                               planform_curvature, profile_curvature

    Returns
    -------
    attr_paths : list arrays contaning the attributes data
        DESCRIPTION.

    """
    
    dem = rd.LoadGDAL(dem_path)
    dem = rd.FillDepressions(dem, epsilon=False, in_place=False)

    attr_arrays = {}
    for attr in args:
        name = os.path.split(path)[1]
        name = name[:name.rfind(' ')] + f' {attr}'
        try:
            attr_band = rd.TerrainAttribute(dem, attrib=attr)
            attr_arrays[name] = attr_band
        except Exception as e:
            print(e)

        # try:
        #     ds = gdal.DEMProcessing(out_tif, dem_path, attr, computeEdges=True)
        #     ds = None
        #     attr_paths.append(out_tif)
        # except RuntimeError as e:
        #     print(e)
        
    return attr_arrays
    

def plot_attr_vals_probability(attr, title, *paths):
    
    gdal.UseExceptions()
    
    clip_vals = {
                    'slope_degrees':
                    {
                        'min': 30, 
                        'max': 60,
                    },
                    'slope_riserun':
                    {
                        'min': 0.6, 
                        'max': 2,
                    },
                      'aspect':
                    {
                        'min': 0, 
                        'max': 360
                    }
                  }
    
    fig, ax = plt.subplots()
    for path in tqdm(paths, total=len(paths), desc=f'Fetching start zones {attr}'): 
        raster = rd.LoadGDAL(path)
        raster = rd.FillDepressions(raster, epsilon=False, in_place=False)
        raster = raster[~np.isnan(raster)]
#        raster = rd.TerrainAttribute(dem, attrib=attrib)
#        r = raster[~np.isnan(raster)]
#        unique, counts = np.unique(r, return_counts=True)
#        if counts.max()/r.shape[0] > 0.75:
#            raster = r[r != unique[counts.argmax()]]
#        slope = raster[np.where(np.logical_and(raster>=clip_vals[attrib]['min'], raster<=clip_vals[attrib]['max']))]
        if attr == 'aspect':
            raster = (raster + 180) % 360
        ku = round(kurtosis(raster.flatten()), 2)
        label = f'{os.path.split(path)[1][:-4]} (kurtosis: {ku})'
        sns.histplot(slope.flatten(), stat='probability', bins=100, 
                     ax=ax, kde=True, line_kws={'lw': 3}, label=label, 
                     alpha=0.4)
    if attr == 'aspect':
        corected_xticks = (ax.get_xticks() + 180) % 360
        #ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels(corected_xticks.astype(int)) 
    plt.title(title)
    plt.legend()
    





if __name__ == '__main__':
    
    cmaps = {
                'deafult': 'winter',
                'slope_degrees': 'magma_r',
                'slope_riserun': 'magma_r', 
                'slope_percentage': 'magma_r', 
                'slope_degrees': 'magma_r', 
                'slope_radians': 'magma_r',
                'aspect': 'twilight_shifted',
                }

    sisters_dem_path = '../data/Sisters/Sisters.tif'
    dem = crop_tif(sisters_dem_path, '../data/Sisters/Sisters_crp.tif',  'Sisters.kml')
    show_array(dem)
    sisters_dem_path = '../data/Sisters/Sisters_crp.tif'
    sisters_kml = glob('../data/Sisters/Sister ? SZ.kml')
    
    
    for kml in sisters_kml:
        out_path = kml.replace('.kml', '_DEM.tif')
        dem = crop_tif(sisters_dem_path, out_path, kml)
    
    
    for path in glob('../data/Sisters/Sister ? SZ_DEM.tif'):
        
        sister_id = path[path.find(' '): path.rfind(' ')].strip()
        slope_params = get_slope_attributes(path, 'curvature')
        for k, arr in slope_params.items():
            show_array(arr[::-1,::-1], cmap='magma_r', title=f'Sister {sister_id} Slope')
    




"""
slope_riserun, slope_percentage, slope_degrees, slope_radians, aspect, curvature, planform_curvature, profile_curvature
"""





sSZdem = crop_tif('../data/Sisters/Sisters_crp.tif', '../data/Sisters/Sisters_SZ_DEM.tif', '../data/Sisters/Sisters_SZs.kml')
show_array(sSZdem)


ds = gdal.Open('../data/Sisters/Sisters_crp.tif')


slope = gdal.DEMProcessing('../data/Sisters/Sisters_slope_crp.tif', ds, 'slope',computeEdges=True)
slp = slope.GetRasterBand(1).ReadAsArray()
slope=None
slp[slp<0] = np.nan
show_array(slp[::-1,::-1], cmap='magma_r', title='Slope angle')

# aspect = gdal.DEMProcessing('../data/Sisters/Sisters_aspect_crp.tif', ds, 'aspect',computeEdges=True)
# asp = aspect.GetRasterBand(1).ReadAsArray()
# slope=None
# asp[asp<0] = np.nan
# show_array(slp, cmap='twilight_shifted')

hillshade = gdal.DEMProcessing('../data/Sisters/Sisters_hilside_crp.tif', ds, 'hillshade',computeEdges=True)
hs = hillshade.GetRasterBand(1).ReadAsArray()
hillshade = None
show_array(hs, cmap='hsv', title='Hillshade')

trindex = gdal.DEMProcessing('../data/Sisters/Sisters_tri_crp.tif', ds, 'TRI',computeEdges=True)
tri = trindex.GetRasterBand(1).ReadAsArray()
trindex = None
tri[tri < 0] = np.nan
show_array(tri, cmap='jet', title='Terrain Roughness Index')

Roughness = gdal.DEMProcessing('../data/Sisters/Sisters_Roughness_crp.tif', ds, 'Roughness',computeEdges=True)
rgh = Roughness.GetRasterBand(1).ReadAsArray()
Roughness = None
rgh[rgh < 0] = np.nan
show_array(rgh, cmap='jet', title='Roughnes')





sigma=2.0                  # standard deviation for Gaussian kernel
truncate=4.0               # truncate filter at this many sigmas

gdal.UseExceptions()

ds = gdal.Open('../data\\Star_mtn\\Star_mtn_crp.tif')

a = np.array(ds.GetRasterBand(1).ReadAsArray())

a[a==a.min()] = np.nan
a1=a.copy()
a1[np.isnan(a)]=0
VV=gaussian_filter(a1,sigma=sigma,truncate=truncate)

a2=0*a.copy()+1
a2[np.isnan(a)]=0
WW=gaussian_filter(a2,sigma=sigma,truncate=truncate)

Z=VV/WW

plt.imshow(Z, cmap='winter')
plt.colorbar()



