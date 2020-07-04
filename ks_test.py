import pandas as pd
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
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

un_treat = np.unique(treatment)
un_treat = un_treat[un_treat != 'nature_reserve']

def calc_stat(value,name_key,type=0):

    distros = []
    for n in range(len(un_treat)):
      if (type == 0):
        l_density = value[treatment == un_treat[n]]
      else:
        l_density = value[n]
       
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


ss = 0.10
density_bins = np.arange(0.0,2.6,step=ss)
calc_stat(density, 'density')
calc_stat(tch, 'tch')
ss = 0.05
density_bins = np.arange(0.0,1.1,step=ss)
calc_stat(rem, 'rem')





key_bootstrap_df = pd.read_csv('polygon_info/polygon_key.csv',sep=',')
un_treat = np.unique(np.array(key_bootstrap_df['treatment']))
un_treat_label = ['Subsistance Agriculture','Communal Grazing', 'Kruger NP','Private Reserve']
height_list = np.load('munged_dat/height_values.npy')

height_list = height_list[np.arange(len(height_list)) != un_treat.tolist().index('nature_reserve')]
un_treat = un_treat[un_treat != 'nature_reserve']


ss = 0.10
density_bins = np.arange(0.2,3.0,step=ss)
calc_stat(height_list, 'height',type=1)











