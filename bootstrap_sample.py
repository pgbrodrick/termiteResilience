import argparse
import gdal
import numpy as np
import pandas as pd
import os,sys,subprocess
from scipy import signal
import numpy.matlib

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from tqdm import tqdm


parser = argparse.ArgumentParser(description='Boostrap sample across sites and polygons')
parser.add_argument('-sample_length',default=1000,type=int)
parser.add_argument('-samples_per_treatment',default=2000,type=int)
args = parser.parse_args()


bdir = '/Carnegie/DGE/Data/Shared/Labs/Asner/Private/Research/Researcher/Davies/4.Kruger/2.Termites_and_landuse/'

mound_center_csvs = [\
os.path.join(bdir,'shapefiles','andover_welv.csv'),\
os.path.join(bdir,'shapefiles','erosionOliver.csv'),\
os.path.join(bdir,'shapefiles','L1.csv'),\
os.path.join(bdir,'shapefiles','L23.csv'),\
os.path.join(bdir,'shapefiles','L45.csv'),\
os.path.join(bdir,'shapefiles','L7.csv'),\
os.path.join(bdir,'shapefiles','L8.csv'),\
os.path.join(bdir,'shapefiles','nwaswitshaka.csv')]

polygon_files = [\
os.path.join(bdir,'polygons','andover_raster.tif'),\
os.path.join(bdir,'polygons','erosion_raster.tif'),\
os.path.join(bdir,'polygons','landuse_1_raster.tif'),\
os.path.join(bdir,'polygons','landuse_2_raster.tif'),\
os.path.join(bdir,'polygons','landuse_456_raster.tif'),\
os.path.join(bdir,'polygons','landuse_7_raster.tif'),\
os.path.join(bdir,'polygons','landuse_8_raster.tif'),\
os.path.join(bdir,'polygons','nwas_raster.tif')]

tch_files = [\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','andover_welv_tch'),\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','erosionOliver_tch'),\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','L1_tch'),\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','L23_tch'),\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','L45_tch'),\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','L7_tch'),\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','L8_tch'),\
os.path.join('/lustre/scratch/pbrodrick/termite_tmp','nwaswitshaka_tch')]



def calc_density(roi,px_x,px_y):
  mounds = np.zeros(roi.shape)
  mounds[px_y,px_x] = 1
  mounds[roi == 0] = 0
  return float(np.sum(mounds)) / float(np.sum(roi == 1))

def calc_mean_mound_height(roi,heights,px_x,px_y):
  lh = np.zeros(roi.shape)-1
  lh[px_y,px_x] = heights
  lh[roi == 0] = -1
  lh[lh == -9999] = -1
  lh = lh[lh != -1]
  return np.mean(lh)

def calc_tch_covers(roi,tch):
  bins = np.arange(0,80)
  bins[0] = -10
  hist, bins = np.histogram(tch[roi != 0].flatten(),bins=bins)
  bins = bins[1:]

  g3_fraction = np.sum(hist[bins > 3])/float(np.sum(hist))
  g1_fraction = np.sum(hist[bins > 1])/float(np.sum(hist))
  b13_fraction = np.sum(hist[np.logical_and(bins > 1,bins <=3)])/float(np.sum(hist))

  return g3_fraction,g1_fraction,b13_fraction

def all_calcs(roi,heights,px_x,px_y,tch):
  
  mound_density = calc_density(roi,px_x,px_y)*per_ha_conv
  mean_mound_height = calc_mean_mound_height(roi,heights,px_x,px_y)
  g1,g3,b13 = calc_tch_covers(roi,tch)
  
  poly_size = np.sum(roi)

  return [mound_density,mean_mound_height,poly_size,g1,g3,b13]

key_df = pd.read_csv(bdir + '/polygons/polygon_key.csv',sep=',')
treatment = np.array(key_df['treatment'])
polygon_number = np.array(key_df['polygon_number'])

un_treatments = np.unique(treatment)
un_poly = np.unique(polygon_number)
sizes = np.zeros((len(un_poly),len(un_treatments)))-9999

# predetermine the number of samples to be take from each treatment, from each landscape
if (os.path.isfile('munged_dat/treatment_sizes.npy') == False):
  for _f in range(0,len(polygon_files)):
    poly_set = gdal.Open(polygon_files[_f],gdal.GA_ReadOnly)
    poly = poly_set.ReadAsArray()
    un_poly = np.unique(poly)
    un_poly = un_poly[un_poly != 0]
    for m in range(0,len(un_poly)):
      sizes[int(un_poly[m])-1,un_treatments.tolist().index(treatment[polygon_number == un_poly[m]][0])] = np.sum(poly == un_poly[m])
  np.save('munged_dat/treatment_sizes.npy',sizes)
else:
  sizes = np.load('munged_dat/treatment_sizes.npy')





subsample_size = [args.sample_length,args.sample_length]

sample_numbers = np.zeros(sizes.shape)
sizes[sizes == -9999] = 0
sample_numbers = (sizes / np.matlib.repmat(np.sum(sizes,axis=0),sizes.shape[0],1) * args.samples_per_treatment).astype(int)

np.random.seed(13)

output = []
for _f in range(0,len(polygon_files)):
  print(('filepair',mound_center_csvs[_f],polygon_files[_f]))

  df = pd.read_csv(mound_center_csvs[_f],sep=',')
  # filter for only mounds < 5 m

  poly_set = gdal.Open(polygon_files[_f],gdal.GA_ReadOnly)
  poly = poly_set.ReadAsArray()
  trans = poly_set.GetGeoTransform()

  tch_set = gdal.Open(tch_files[_f],gdal.GA_ReadOnly)
  tch = tch_set.ReadAsArray()
  
  per_ha_conv = 1.0e4/float(trans[1]*trans[1])
  
  px_y = ((np.array(df['y'])-trans[3])/float(trans[5])).astype(int)
  px_x = ((np.array(df['x'])-trans[0])/float(trans[1])).astype(int)
  heights = np.array(df['height'])

  un_poly = np.unique(poly)
  un_poly = un_poly[un_poly != 0]

  for m in range(0,len(un_poly)):
    print('Polygon ' + str(un_poly[m]))
    lp = poly == un_poly[m]
    
    x_has_nodata = np.any(lp != 0,axis=0).astype(int)
    y_has_nodata = np.any(lp != 0,axis=1).astype(int)

    y_list,x_list = np.where(lp == 1)
    subset = np.logical_and.reduce((y_list >= subsample_size[1]/2,y_list < lp.shape[0] - subsample_size[1]/2,x_list >= subsample_size[0]/2,x_list < lp.shape[1] - subsample_size[0]/2))
    x_list = x_list[subset]
    y_list = y_list[subset]

    perm = np.random.permutation(len(x_list))
    y_list = y_list[perm[:sample_numbers[int(un_poly[m])-1,un_treatments.tolist().index(treatment[polygon_number == un_poly[m]][0])]]]
    x_list = x_list[perm[:sample_numbers[int(un_poly[m])-1,un_treatments.tolist().index(treatment[polygon_number == un_poly[m]][0])]]]

    for i in tqdm(range(0,len(x_list)), ncols=80):
      y_b = int(y_list[i]-subsample_size[1]/2)
      y_t = int(y_list[i]+subsample_size[1]/2)
      x_b = int(x_list[i]-subsample_size[0]/2)
      x_t = int(x_list[i]+subsample_size[0]/2)

      subset = np.logical_and.reduce((px_y >= y_b, px_y < y_t,px_x >= x_b, px_x < x_t))

      ret_list = all_calcs(lp[y_b:y_t,x_b:x_t],heights[subset],px_x[subset]-x_b,px_y[subset]-y_b,tch[y_b:y_t,x_b:x_t])

      ret_list.insert(0,int(un_poly[m]))
      ret_list.append(np.array(key_df['treatment'])[np.where(key_df['polygon_number'] == un_poly[m])][0])
      ret_list.append(np.array(key_df['landscape'])[np.where(key_df['polygon_number'] == un_poly[m])][0])
      ret_list.append(i)

      output.append(np.array(ret_list))
  
output = np.array(output)

df = pd.DataFrame(data=output,columns=['polygon','mound_density','mean_mound_height','rep_poly_count','cover_g1','cover_g3','cover_b13','treatment','landscape','rep'])
df.to_csv('bootstrapped_results/ss_' + str(subsample_size[0]) + '.csv',sep=',',index=False)







