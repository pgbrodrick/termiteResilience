import gdal
import numpy as np
import pandas as pd
import os,sys,subprocess
from scipy import signal
import numpy.matlib

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



bdir = '/Carnegie/DGE/Data/Shared/Labs/Asner/Private/Research/Researcher/Davies/4.Kruger/2.Termites_and_landuse/'

mound_center_csvs = [\
os.path.join('full_landscape_mound_predictions','andover_welv_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','erosionOliver_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L1_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L23_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L45_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L7_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L8_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','nwaswitshaka_cnn_mounds.csv')]

polygon_files = [\
os.path.join(bdir,'polygons','andover_raster.tif'),\
os.path.join(bdir,'polygons','erosion_raster.tif'),\
os.path.join(bdir,'polygons','landuse_1_raster.tif'),\
os.path.join(bdir,'polygons','landuse_2_raster.tif'),\
os.path.join(bdir,'polygons','landuse_456_raster.tif'),\
os.path.join(bdir,'polygons','landuse_7_raster.tif'),\
os.path.join(bdir,'polygons','landuse_8_raster.tif'),\
os.path.join(bdir,'polygons','nwas_raster.tif')]

def crop_image(dat):
  x_has_nodata = np.any(dat != -9999,axis=0).astype(int)
  y_has_nodata = np.any(dat != -9999,axis=1).astype(int)

  trim_x_b = np.argmax(x_has_nodata)
  trim_x_t = np.argmax(x_has_nodata[::-1])
  trim_y_b = np.argmax(y_has_nodata)
  trim_y_t = np.argmax(y_has_nodata[::-1])
  
  dat = dat[trim_y_b:-trim_y_t,trim_x_b:-trim_x_t]
  return(dat)



key_df = pd.read_csv(bdir + '/polygons/polygon_key.csv',sep=',')
print(list(key_df))
output = []
driver = gdal.GetDriverByName('GTiff') 
driver.Register()


for _f in range(0,len(polygon_files)):
  print(('filepair',mound_center_csvs[_f],polygon_files[_f]))

  df = pd.read_csv(mound_center_csvs[_f],sep=',')
  # filter for only mounds < 5 m

  poly_set = gdal.Open(polygon_files[_f],gdal.GA_ReadOnly)
  poly = poly_set.ReadAsArray()
  trans = poly_set.GetGeoTransform()
  
  per_ha_conv = 1.0e4/float(trans[1]*trans[1])
  
  un_poly = np.unique(poly)
  un_poly = un_poly[un_poly != 0]
  un_poly = un_poly[un_poly != 31]

  for m in range(0,len(un_poly)):
    current_poly = un_poly[m]

    lp = poly == un_poly[m]
    
    of_name = 'poly_' + str(un_poly[m])
    outDataset = driver.Create('rk_runs/roi_raster/' + of_name +'.tif',lp.shape[1],lp.shape[0],1,gdal.GDT_Byte)
    outDataset.SetProjection(poly_set.GetProjection())
    outDataset.SetGeoTransform(trans)
    outDataset.GetRasterBand(1).WriteArray(lp,0,0)
    del outDataset

    cmd_str = 'gdal_polygonize.py rk_runs/roi_raster/' + of_name + '.tif rk_runs/roi_shp/' + \
              of_name + '.shp -mask ' + 'rk_runs/roi_raster/' + of_name + '.tif -f \'ESRI Shapefile\''  
     
    subprocess.call(cmd_str,shell=True)

    pd.DataFrame(data=np.transpose(np.vstack([np.array(df['x']),np.array(df['y'])]))).to_csv('rk_runs/roi_points/'+of_name+'.csv',index=False)

    outstr = 'sh sub_rk.csh rk_runs/roi_points/' +  of_name  + '.csv rk_runs/roi_output/' + of_name + '.csv rk_runs/roi_shp/' + of_name + '.shp'
    print(outstr)
    subprocess.call(outstr,shell=True)



















