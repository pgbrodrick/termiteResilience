import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.ensemble import RandomForestRegressor

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
import argparse
import gdal
from scipy import stats
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import fmin_l_bfgs_b,minimize
import multiprocessing
import matplotlib.patches as mpatches
import subprocess
import sys
from scipy import stats



def plot_lhist(hist,edge):
  hist = hist.tolist()
  edge = edge.tolist()
  phist = [0]
  pedge = [edge[0]]
  for _e in range(0,len(edge)-1):
    phist.append(hist[_e])
    phist.append(hist[_e])
  
    pedge.append(edge[_e])
    pedge.append(edge[_e+1])
  
  phist.append(0)
  pedge.append(edge[-1])
  phist = np.array(phist)
  pedge = np.array(pedge)
  return phist,pedge


subsample_size = 1000

df = pd.read_csv('bootstrapped_results/ss_' + str(subsample_size) + '.csv',sep=',')

density = np.array(df['mound_density'])
height = np.array(df['mean_mound_height'])
treatment = np.array(df['treatment'])
reps = np.array(df['rep_poly_count'])
rem = np.array(df['mean_rem'])
tch = np.array(df['mean_tch'])

ss = 0.10
density_bins = np.arange(0.0,2.6,step=ss)
un_treat = np.unique(treatment)
un_treat = un_treat[un_treat != 'nature_reserve']

def calc_stat(value,name_key):

    distros = []
    for n in range(len(un_treat)):
      l_density = value[treatment == un_treat[n]]
       
      h,b = np.histogram(l_density,bins=density_bins)
      b,h = plot_lhist(h,b)
      distros.append(l_density)
    
    ks_stat = np.zeros((len(distros),len(distros)))
    p_value = np.zeros((len(distros),len(distros)))
    for _d in range(len(distros)):
     for __d in range(len(distros)):
       d,p = stats.ks_2samp(distros[_d], distros[__d]) 
       ks_stat[_d,__d] = d
       p_value[_d,__d] = p
    
    df = pd.DataFrame(data=ks_stat,columns=un_treat)
    df.to_csv('stats/' + name_key + '_ks_stat.csv',sep=',',index=False)
    
    df = pd.DataFrame(data=p_value,columns=un_treat)
    df.to_csv('stats/' + name_key + '_ks_pvalue.csv',sep=',',index=False)


calc_stat(density, 'density')
calc_stat(tch, 'tch')
calc_stat(rem, 'rem')


