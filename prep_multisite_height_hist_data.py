import pandas as pd
import numpy as np
from sklearn import preprocessing

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy import stats
import sys, os
from matplotlib.lines import Line2D
import gdal

def get_mound_heights(roi,df,px_x,px_y):

  heights = np.zeros(roi.shape)-1
  heights[px_y,px_x] = np.array(df['height'])
  heights[roi == 0] = -1
  heights[heights == -9999] = -1

  heights = heights[heights != -1]
  return heights



#bdir = '/Carnegie/DGE/Data/Shared/Labs/Asner/Private/Research/Researcher/Davies/4.Kruger/2.Termites_and_landuse/'
bdir = 'conservative_tests'

mound_center_csvs = [\
os.path.join('full_landscape_mound_predictions','andover_welv_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','erosionOliver_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L1_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L23_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L45_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L7_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','L8_cnn_mounds.csv'),\
os.path.join('full_landscape_mound_predictions','nwaswitshaka_cnn_mounds.csv')]

bdir = '/Carnegie/DGE/Data/Shared/Labs/Asner/Private/Research/Researcher/Davies/4.Kruger/2.Termites_and_landuse/'
polygon_files = [\
os.path.join(bdir,'polygons','andover_raster.tif'),\
os.path.join(bdir,'polygons','erosion_raster.tif'),\
os.path.join(bdir,'polygons','landuse_1_raster.tif'),\
os.path.join(bdir,'polygons','landuse_2_raster.tif'),\
os.path.join(bdir,'polygons','landuse_456_raster.tif'),\
os.path.join(bdir,'polygons','landuse_7_raster.tif'),\
os.path.join(bdir,'polygons','landuse_8_raster.tif'),\
os.path.join(bdir,'polygons','nwas_raster.tif')]



key_df = pd.read_csv(bdir + '/polygons/polygon_key.csv',sep=',')
poly_size_numbers = []
poly_size_area = []

un_treat = np.unique(np.array(key_df['treatment']))
un_treat_label = ['Subsistance Agriculture','Communal Grazing', 'Kruger NP','Private Reserve']
height_list = [np.zeros(0) for x in range(len(un_treat))]
for _f in range(0,len(polygon_files)):
 print(('filepair',mound_center_csvs[_f],polygon_files[_f]))

 df = pd.read_csv(mound_center_csvs[_f],sep=',')

 poly_set = gdal.Open(polygon_files[_f],gdal.GA_ReadOnly)
 poly = poly_set.ReadAsArray()

 un_poly = np.unique(poly)
 un_poly = un_poly[un_poly != 0]
 un_poly = un_poly[un_poly != 31]
 for m in range(0,len(un_poly)):
   poly_size_numbers.append(un_poly[m])
   poly_size_area.append(np.sum(poly == un_poly[m]))


 trans = poly_set.GetGeoTransform()
 
 per_ha_conv = 1.0e4/float(trans[1]*trans[1])
 
 px_y = ((np.array(df['y'])-trans[3])/float(trans[5])).astype(int)
 px_x = ((np.array(df['x'])-trans[0])/float(trans[1])).astype(int)

 un_poly = np.unique(poly)
 un_poly = un_poly[un_poly != 0]
 for m in range(0,len(un_poly)):
   current_poly = un_poly[m]
   ret_list = get_mound_heights(poly == un_poly[m],df,px_x,px_y)
   height_list[un_treat.tolist().index(np.array(key_df['treatment'])[np.where(key_df['polygon_number'] == un_poly[m])][0])] = np.append(height_list[un_treat.tolist().index(np.array(key_df['treatment'])[np.where(key_df['polygon_number'] == un_poly[m])][0])],ret_list)

 #  np.savez('munged_dat/height_poly_' + str(un_poly[m]) + '.npz',histogram=ret_list)

np.savez('munged_dat/poly_areas.npz',poly_size_area = poly_size_area,poly_size_numbers=poly_size_numbers)
np.save('munged_dat/height_values.npy',height_list)
