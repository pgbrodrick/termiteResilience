import gdal
import numpy as np
import pandas as pd
import os,sys,subprocess

import matplotlib as mpl
mpl.use('Agg')



bdir = 'path/to/root/raster/data'

ensemble_files = [\
os.path.join(bdir,'ensemble','andover_welv_128_256_ensemble.tif'),\
os.path.join(bdir,'ensemble','erosionOliver_128_256_ensemble.tif'),\
os.path.join(bdir,'ensemble','L1_128_256_ensemble.tif'),\
os.path.join(bdir,'ensemble','L23_128_256_ensemble.tif'),\
os.path.join(bdir,'ensemble','L45_128_256_ensemble.tif'),\
os.path.join(bdir,'ensemble','L7_128_256_ensemble.tif'),\
os.path.join(bdir,'ensemble','L8_128_256_ensemble.tif'),\
os.path.join(bdir,'ensemble','nwaswitshaka_128_256_ensemble.tif')]

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


def calc_density(roi,df,px_x,px_y):

  mounds = np.zeros(roi.shape)
  mounds[px_y,px_x] = 1
  mounds[roi == 0] = 0

  return float(np.sum(mounds)) / float(np.sum(roi == 1))

def calc_mean_mound_height(roi,df,px_x,px_y):

  heights = np.zeros(roi.shape)-1
  heights[px_y,px_x] = np.array(df['height'])
  heights[roi == 0] = -1
  heights[heights == -9999] = -1

  heights = heights[heights != -1]
  return np.mean(heights),np.std(heights)

def calc_n_mounds(roi,px_x,px_y):
  count = np.zeros(roi.shape)
  count[px_y,px_x] = 1
  count[roi == 0] = 0
  return np.sum(count)

def all_calcs(roi,df,px_x,px_y):
  
  mound_density = calc_density(roi,df,px_x,px_y)*per_ha_conv
  mean_mound_height,std_mound_height = calc_mean_mound_height(roi,df,px_x,px_y)
  area = np.sum(roi)
  mound_count = calc_n_mounds(roi,px_x,px_y)

  return [mound_density,mean_mound_height,std_mound_height,area,mound_count]



key_df = pd.read_csv('polygon_info/polygon_key.csv',sep=',')


output_df = pd.DataFrame(data=np.array(key_df['landscape']),columns=['landscape'])
output_df['treatment'] = key_df['treatment']

area                   = np.zeros(len(output_df['treatment']))
n_mounds               = np.zeros(len(output_df['treatment']))
bs_mean_mound_density  = np.zeros(len(output_df['treatment']))
bs_std_mound_density   = np.zeros(len(output_df['treatment']))
mean_mound_height      = np.zeros(len(output_df['treatment']))
std_mound_height       = np.zeros(len(output_df['treatment']))
mound_density          = np.zeros(len(output_df['treatment']))
low_rem_frac           = np.zeros(len(output_df['treatment']))
high_rem_frac          = np.zeros(len(output_df['treatment']))
veg_g1                 = np.zeros(len(output_df['treatment']))
veg_g3                 = np.zeros(len(output_df['treatment']))
veg_b13                = np.zeros(len(output_df['treatment']))


polygon_number = np.array(key_df['polygon_number'])


for _f in range(0,len(ensemble_files)):
  print(('filepair',mound_center_csvs[_f],polygon_files[_f]))

  df = pd.read_csv(mound_center_csvs[_f],sep=',')

  poly_set = gdal.Open(polygon_files[_f],gdal.GA_ReadOnly)
  poly = poly_set.ReadAsArray()
  trans = poly_set.GetGeoTransform()
  
  per_ha_conv = 1.0e4/float(trans[1]*trans[1])
  
  px_y = ((np.array(df['y'])-trans[3])/float(trans[5])).astype(int)
  px_x = ((np.array(df['x'])-trans[0])/float(trans[1])).astype(int)

  un_poly = np.unique(poly)
  un_poly = un_poly[un_poly != 0]


  for m in range(0,len(un_poly)):
    arr_idx = np.where(polygon_number == un_poly[m])[0][0]
    current_poly = un_poly[m]

    roi = poly == un_poly[m]

    area[arr_idx] = np.sum(roi)
    mound_density[arr_idx] = calc_density(roi,df,px_x,px_y)*per_ha_conv
    n_mounds[arr_idx] = calc_n_mounds(roi,px_x,px_y)

    mean, std = calc_mean_mound_height(roi,df,px_x,px_y)

    mean_mound_height[arr_idx] = mean
    std_mound_height[arr_idx] = std


for _p in range(0,len(polygon_number)):

    npzf = np.load('munged_dat/rem_poly_' + str(polygon_number[_p]) + '.0.npz') 
    frequency = npzf['cur_frequency']
    ratio = np.sum(frequency[-10:]) / np.sum(frequency)
    high_rem_frac[_p] = ratio
    ratio = np.sum(frequency[:10]) / np.sum(frequency)
    low_rem_frac[_p] = ratio


df = pd.read_csv('bootstrapped_results/ss_1000.csv')
bs_polygon = np.array(df['polygon'])

for _p in range(0,len(polygon_number)):
    bs_mean_mound_density[_p] = np.mean(np.array(df['mound_density'])[bs_polygon == polygon_number[_p]])
    bs_std_mound_density[_p] = np.std(np.array(df['mound_density'])[bs_polygon == polygon_number[_p]])


for _p in range(0,len(polygon_number)):
    npzf = np.load('munged_dat/tch_rep_poly_' + str(polygon_number[_p]) + '.0.npz') 
    hist = npzf['histogram']
    bins = npzf['bins']
    veg_g1[_p] = np.sum(hist[1:]) / np.sum(hist)
    veg_g3[_p] = np.sum(hist[3:]) / np.sum(hist)
    veg_b13[_p] = np.sum(hist[1:3]) / np.sum(hist)




output_df['area'] = area
output_df['number_of_mounds'] = n_mounds
output_df['bs_mean_mound_density'] = bs_mean_mound_density
output_df['bs_std_mound_density'] = bs_std_mound_density
output_df['mean_mound_height'] = mean_mound_height
output_df['std_mound_height'] = std_mound_height
output_df['mound_density'] = mound_density

output_df['low_rem_ratio'] = low_rem_frac
output_df['high_rem_ratio'] = high_rem_frac
output_df['veg_g1'] = veg_g1
output_df['veg_g3'] = veg_g3
output_df['veg_b13'] = veg_b13


output_df.to_csv('results.csv',sep=',')

