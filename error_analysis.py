
import argparse
import gdal
import numpy as np
import pandas as pd
import os,sys,subprocess
from scipy import signal
from skimage import measure
import numpy.matlib
from tqdm import tqdm




act_df = pd.read_csv(sys.argv[1],sep=',')
pred_df = pd.read_csv(sys.argv[2],sep=',')


pred_xy = np.vstack([pred_df['x'],pred_df['y']]).T
pred_area = np.array(pred_df['area'])
try:
    act_xy = np.vstack([act_df['X'],act_df['Y']]).T
except:
    act_xy = np.vstack([act_df['x'],act_df['y']]).T

valid = np.logical_and(pred_xy[:,0] >= np.min(act_xy[:,0]), pred_xy[:,0] <= np.max(act_xy[:,0]))
valid[pred_xy[:,1] < np.min(act_xy[:,1])] = False
valid[pred_xy[:,1] > np.max(act_xy[:,1])] = False

pred_xy = pred_xy[valid,:]
pred_area = pred_area[valid]



error_dist = 10


for min_area in range(0,40,2):
  tp = 0
  fp = 0
  print(min_area)
  for _i in range(len(pred_xy)):

    #if (np.sum(np.sum(act_xy == pred_xy[_i,:].reshape((1,2)),axis=1) == 2) > 0):
    closest_dist = np.min(np.power(np.power(act_xy[:,0]-pred_xy[_i,0],2) + np.power(act_xy[:,1]-pred_xy[_i,1],2),0.5))
    if (np.sum(closest_dist < error_dist) > 0 and pred_area[_i] > min_area):
        tp += 1

    if (np.sum(closest_dist < error_dist) == 0 and pred_area[_i] > min_area):
        fp += 1

  print(('T:',len(act_xy)))
  print(('TP:',tp))
  print(('FP:',fp))
  print(('TPR:',tp/float(len(act_xy))))
  print(('FPR:',fp/float(len(act_xy))))
  print('\n')




    



