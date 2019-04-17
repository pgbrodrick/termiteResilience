


import argparse
import gdal
import numpy as np
import pandas as pd
import os,sys,subprocess
from scipy import signal
from skimage import measure
import numpy.matlib
from tqdm import tqdm


import fiona
import rasterio.features



parser = argparse.ArgumentParser(description='Re-center mound estimates based on local DEM max')
parser.add_argument('prob_file')
parser.add_argument('dem_file')
parser.add_argument('out_file')
parser.add_argument('-boundary_file',default=None)
parser.add_argument('-band',default=2,type=int)
args = parser.parse_args()


print(args.prob_file)
print(args.dem_file)
print(args.out_file)
csv_name = os.path.join(os.path.dirname(args.out_file) , os.path.basename(args.out_file).split('.')[0] + '.csv')
print(csv_name)


prob_set = gdal.Open(args.prob_file,gdal.GA_ReadOnly)
dem_set = gdal.Open(args.dem_file,gdal.GA_ReadOnly)


prob = prob_set.GetRasterBand(args.band).ReadAsArray()
dem = dem_set.GetRasterBand(1).ReadAsArray()
if (args.boundary_file is not None):
    boundary_set = gdal.Open(args.boundary_file,gdal.GA_ReadOnly)
    boundary = boundary_set.ReadAsArray()
    prob[boundary == 0] = 0

thresh = 0.99
mounds = prob > thresh
labeled_mounds = measure.label(mounds,connectivity=2)


un_mounds = np.unique(labeled_mounds)
un_mounds = un_mounds[un_mounds != 0]
print('# mounds: ' + str(len(un_mounds)))

ringmask = \
[\
[0,0,0,0,1,1,1,0,0,0,0],
[0,0,1,1,0,0,0,1,1,0,0],
[0,1,0,0,0,0,0,0,0,1,0],
[0,1,0,0,0,0,0,0,0,1,0],
[1,0,0,0,0,0,0,0,0,0,1],
[1,0,0,0,0,0,0,0,0,0,1],
[1,0,0,0,0,0,0,0,0,0,1],
[0,1,0,0,0,0,0,0,0,1,0],
[0,1,0,0,0,0,0,0,0,1,0],
[0,0,1,1,0,0,0,1,1,0,0],
[0,0,0,0,1,1,1,0,0,0,0]]
ringmask = np.array(ringmask).astype(bool)

ringmask = \
[\
[0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0],
[0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0],
[0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0]]
ringmask = np.array(ringmask).astype(bool)

interior = \
[\
[0,0,0,0,1,1,1,0,0,0,0],
[0,0,1,1,1,1,1,1,1,0,0],
[0,1,1,1,1,1,1,1,1,1,0],
[0,1,1,1,1,1,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1],
[0,1,1,1,1,1,1,1,1,1,0],
[0,1,1,1,1,1,1,1,1,1,0],
[0,0,1,1,1,1,1,1,1,0,0],
[0,0,0,0,1,1,1,0,0,0,0]]
interior = np.array(interior).astype(bool)


interior = \
[\
[0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
[0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
[0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0]]
interior = np.array(interior).astype(bool)


trans = dem_set.GetGeoTransform()
xmin = trans[0]

np.set_printoptions(linewidth=200)
x_loc = []
y_loc = []
height = []
cum_prob = []
max_prob = []
hm_prob = []
area = []


all_x = np.matlib.repmat(np.arange(0,dem_set.RasterXSize,step=1).reshape(1,dem_set.RasterXSize),dem_set.RasterYSize,1).flatten().astype(int)
all_y = np.matlib.repmat(np.arange(0,dem_set.RasterYSize,step=1).reshape(dem_set.RasterYSize,1),1,dem_set.RasterXSize).flatten().astype(int)

labeled_mounds = labeled_mounds.flatten()
mounds = mounds.flatten()
labeled_mounds = labeled_mounds[mounds]
all_x = all_x[mounds]
all_y = all_y[mounds]
sub_dem = dem.flatten()[mounds]
sub_prob = prob.flatten()[mounds]

for n in tqdm(range(len(un_mounds)), ncols=80):

  loc = np.where(labeled_mounds == un_mounds[n])

  peak = [all_y[loc][np.argmax(sub_dem[loc])],all_x[loc][np.argmax(sub_dem[loc])]]
  if (peak[0] > 10 and peak[0] < dem_set.RasterYSize-10 and peak[1] > 10 and peak[1] < dem_set.RasterXSize-10):

    #peak = [loc[0][np.argmax(dem[loc])],loc[1][np.argmax(dem[loc])]]
    p_subset = prob[peak[0]-10:peak[0]+11,peak[1]-10:peak[1]+11]
    p_subset = p_subset[interior]
    l_area = np.sum(p_subset > thresh)
    if (l_area > 7):
        area.append(l_area)
        subset = dem[peak[0]-10:peak[0]+11,peak[1]-10:peak[1]+11]

        hm_prob.append(p_subset[np.argmax(subset[interior])])

        if (not np.any(subset == -9999)):
          mean = np.mean(subset[ringmask])
          height.append(np.max(subset[interior]) - mean)
        else:
          height.append(-9999)
        

        x_loc.append(trans[0] + trans[1]*peak[1])
        y_loc.append(trans[3] + trans[5]*peak[0])
        cum_prob.append(np.sum(sub_prob[loc]))
        max_prob.append(np.max(sub_prob[loc]))
 
    #print((n,len(un_mounds),height[-1]))

dat = np.transpose(np.vstack([x_loc,y_loc,height,cum_prob,max_prob,hm_prob,area]))
df = pd.DataFrame(data=dat,columns=['x','y','height','cumulative_probability','max_probability','probability_at_max_height','area'])

df.to_csv(csv_name,sep=',',index=False)
cmd_str = 'ogr2ogr ' + args.out_file + ' ' + csv_name + ' -a_srs EPSG:32736 -f \'ESRI Shapefile\' -oo X_POSSIBLE_NAMES=x -oo Y_POSSIBLE_NAMES=y'
print(cmd_str)
subprocess.call(cmd_str,shell=True)














