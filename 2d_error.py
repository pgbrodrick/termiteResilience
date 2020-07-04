
import argparse
import gdal
import numpy as np
import pandas as pd
import os,sys,subprocess
from scipy import signal
from skimage import measure
import numpy.matlib
from tqdm import tqdm





bdir = '/path/to/data/base/'

act_df = pd.read_csv(sys.argv[1],sep=',')
pred_df = pd.read_csv(sys.argv[2],sep=',')


pred_xy = np.vstack([pred_df['x'],pred_df['y']]).T
pred_area = np.array(pred_df['area'])
pred_prob = np.array(pred_df['cumulative_probability'])

try:
    act_xy = np.vstack([act_df['X'],act_df['Y']]).T
except:
    act_xy = np.vstack([act_df['x'],act_df['y']]).T

valid = np.logical_and(pred_xy[:,0] >= np.min(act_xy[:,0]), pred_xy[:,0] <= np.max(act_xy[:,0]))
valid[pred_xy[:,1] < np.min(act_xy[:,1])] = False
valid[pred_xy[:,1] > np.max(act_xy[:,1])] = False

pred_xy = pred_xy[valid,:]
pred_area = pred_area[valid]
pred_prob = pred_prob[valid]

closest_dist = np.zeros(len(pred_area))


for _i in range(len(pred_xy)):
    closest_dist[_i] = np.min(np.power(np.power(act_xy[:,0]-pred_xy[_i,0],2) + np.power(act_xy[:,1]-pred_xy[_i,1],2),0.5))



error_dist = 10


prob_range = np.arange(np.min(pred_prob),np.max(pred_prob) + 0.01,0.1)
area_range = np.arange(0,30)
prob_range = [0]
area_range = [0]

output = np.zeros((6,len(prob_range),len(area_range)))
for _min_prob in range(len(prob_range)):
    for _min_area in range(len(area_range)):
        min_prob = prob_range[_min_prob]
        min_area = area_range[_min_area]
        
        tp = np.sum(np.logical_and.reduce((closest_dist < error_dist, pred_area > min_area, pred_prob > min_prob)))
        fp = np.sum(np.logical_and.reduce((closest_dist > error_dist, pred_area > min_area, pred_prob > min_prob)))

        tpr = round(tp/float(len(act_xy)),3)
        fpr = round(fp/float(len(act_xy)),3)

        output[:,_min_prob,_min_area] = np.array([min_prob,min_area,tp,fp,tpr,fpr])

output = np.reshape(output,(output.shape[0],np.prod(output.shape[1:]))).T
        
subset = np.zeros(len(output))
cutoff = 0.9
while (np.sum(subset) < 1):
  subset = output[:,-2] >= cutoff
  final_output = output[subset,:]
  cutoff -= 0.1

  if (cutoff <= 0):
    print('no valid cutoff found')
    quit()

_best = np.argmin(final_output[:,-1])
print('{0}, {1}, {2}'.format(final_output[_best,2],final_output[_best,3], len(act_xy)))









