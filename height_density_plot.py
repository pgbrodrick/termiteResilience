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


############################ Figure 3 ################################

fig = plt.figure(figsize=(12,8))
#gs = gridspec.GridSpec(3, 4)
#gs.update(wspace=0.3,hspace=1.0)

subsample_size = 1000
df = pd.read_csv('bootstrapped_results/mc_results_ss_' + str(subsample_size) + '.csv',sep=',')

density = np.array(df['mound_density'])
height = np.array(df['mean_mound_height'])
cover_g1 = np.array(df['cover_g1'])
cover_g3 = np.array(df['cover_g3'])
cover_b13 = np.array(df['cover_b13'])
treatment = np.array(df['treatment'])
reps = np.array(df['rep_poly_count'])

un_treat = np.unique(treatment)
un_treat = un_treat[un_treat != 'nature_reserve']
un_treat_label = ['Subsistance Agriculture','Communal Grazing', 'Kruger NP','Private Reserve']

################# 3a
ax = fig.add_axes([0.07,0.96,0.02,0.02],zorder=1)
plt.text(0,0,'a',fontweight='bold')
plt.axis('off')
ax = fig.add_axes([0.1,0.5,0.4,0.45],zorder=0)

ss = 0.10
density_bins = np.arange(0.0,2.6,step=ss)
for n in range(len(un_treat)):
  
  l_density = density[treatment == un_treat[n]]
  l_reps = reps[treatment == un_treat[n]]
  l_density = l_density[l_reps >= float(int(subsample_size)**2)/2.]

  h,b = np.histogram(l_density,bins=density_bins)
  b,h = plot_lhist(h,b)
  plt.plot(h + float(n) / 15. * ss,b)

plt.xlabel('Termite Mound Density [mounds ha$^{-1}$]')
plt.ylabel('Frequency')

plt.legend(un_treat_label)



################# 3b


key_df = pd.read_csv('polygon_info/polygon_key.csv',sep=',')
un_treat = np.unique(np.array(key_df['treatment']))
print(un_treat)
un_treat_label = ['Subsistance Agriculture','Communal Grazing', 'Kruger NP','Private Reserve']

#ax = plt.subplot(gs[:2, 2:])
ax = fig.add_axes([0.52,0.96,0.02,0.02],zorder=1)
plt.text(0,0,'b',fontweight='bold')
plt.axis('off')
ax = fig.add_axes([0.55,0.5,0.4,0.45],zorder=0)

height_list = np.load('munged_dat/height_values.npy')

ss = 0.10
height_bins = np.arange(0.2,3.0,step=ss)
for n in range(len(un_treat)):
  if (un_treat[n] != 'nature_reserve'):   
    h,b = np.histogram(height_list[n],bins=height_bins)
    b,h = plot_lhist(h/np.max(h),b)
    #b,h = plot_lhist(h,b)
    plt.plot(h + float(n) / 15. * ss,b)
    print(un_treat[n])

plt.xlabel('Termite Mound Height [m]')
plt.ylabel('Relative Frequency')

plt.legend(un_treat_label)



################# 3c


def get_rk_extr(f):
  df = np.array(pd.read_csv(f))
  df = df[:,[0,2]]
  df = df[np.isnan(df[:,1]) == False,:]

  df = df[df[:,0] < 250,:]
  df[:,1] = np.power(df[:,1]/np.pi,0.5)
  return df


if (os.path.isfile('munged_dat/poly_areas.npz') == False):
  poly_size_numbers = []
  poly_size_area = []
  for _f in range(0,len(polygon_files)):
    poly_set = gdal.Open(polygon_files[_f],gdal.GA_ReadOnly)
    poly = poly_set.ReadAsArray()
    un_poly = np.unique(poly)
    un_poly = un_poly[un_poly != 0]
    un_poly = un_poly[un_poly != 31]
    for m in range(0,len(un_poly)):
      poly_size_numbers.append(un_poly[m])
      poly_size_area.append(np.sum(poly == un_poly[m]))
  np.savez('munged_dat/poly_areas.npz',poly_size_area = poly_size_area,poly_size_numbers=poly_size_numbers)
else:
  npzf = np.load('munged_dat/poly_areas.npz')
  poly_size_numbers = npzf['poly_size_numbers'].tolist()
  poly_size_area = npzf['poly_size_area'].tolist()


def bin_vals(x,binsize=5):
  out = []
  for n in range(0,len(x),binsize):
    if (n+binsize >= len(x)):
      break
    out.append(np.nanmean(x[n:n+binsize]))
  if (n != len(x)-1):
    out.append(np.nanmean(x[n:]))
  return(np.array(out))


un_treat = np.unique(treatment)
un_treat = un_treat[un_treat != 'nature_reserve']
un_treat_label = ['Subsistance Agriculture','Communal Grazing', 'Kruger NP','Private Reserve']

colors = ['blue','green','red','orange','grey','purple','cyan','brown','black','yellow']

key_df = pd.read_csv('polygon_info/polygon_key.csv',sep=',')
treatment = np.array(key_df['treatment'])
polygon_number = np.array(key_df['polygon_number'])

ax = fig.add_axes([0.07,0.39,0.02,0.02],zorder=1)
plt.text(0,0,'c',fontweight='bold')
plt.axis('off')

ind = np.zeros(4).astype(int)
for _t in range(0,len(un_treat)):
  loc_poly = polygon_number[treatment == un_treat[_t]]
  for m in range(len(loc_poly)):
   if (loc_poly[m] != 31 and poly_size_area[poly_size_numbers.index(loc_poly[m])]/10000. > 300):
    fname = 'roi_output/poly_'+str(loc_poly[m]) + '.0.csv' 
    rk = get_rk_extr(fname)

    ax = fig.add_axes([0.1 + 0.9/float(len(un_treat)) * _t,0.1,0.9/float(len(un_treat))-0.05,0.26],zorder=0)
    if(ind[_t] == 0):
      plt.title(un_treat_label[_t])
      plt.xlabel('Distance [m]')
    if (_t == 0):
      plt.ylabel('Ripley\'s K (L transformation)')

    x = bin_vals(rk[:,0])
    y = bin_vals(rk[:,1])
    plt.plot([0,250],[0,250],c='red',ls='--',lw=1,zorder=0)
    #plt.plot(x,y,lw=2,c=colors[ind[_t]],zorder=1) 
    plt.plot(x,y,lw=1,c='black',zorder=1) 

    plt.axis('equal')
    plt.xlim([-5,255])
    plt.ylim([-5,255])
    plt.xticks([0,50,100,150,200,250])
    plt.yticks([0,50,100,150,200,250])
    ind[_t] += 1
 

plt.savefig('figs/figure_3.png',dpi=200)
plt.clf()





############################## Figure 4 ###########################################
fig = plt.figure(figsize=(7,6))

color_ref_dict = ['blue','orange','green','red']
color_ref = [mpl.colors.to_rgba(x) for x in color_ref_dict]

treatment = np.array(df['treatment'])
colors = np.ones((len(treatment),3))
for i in range(0,len(un_treat)):
  colors[treatment == un_treat[i],:] = np.array(color_ref[i])[:3]

valid = np.isnan(height)==False
valid[np.logical_and(height > 1.6,density < 0.5)] = False
plt.scatter(height[valid],density[valid],c=colors[valid],s=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(height[valid],density[valid])

plt.xlabel('Mean Mound Height [m]')
plt.ylabel('Mean Mound Density [Mounds ha$^{-1}$]')


legend_elements = []
for _tr in range(len(un_treat)):
  tr = un_treat[_tr]
  subset = np.logical_and.reduce((np.isnan(height) == False,np.isnan(density) == False,treatment == tr))
  subset[np.logical_and(height > 1.6,density < 0.5)] = False
  slope, intercept, r_value, p_value, std_err = stats.linregress(height[subset],density[subset])

  plt.plot([np.min(height[valid]),np.max(height[valid])],[np.min(height[valid])*slope+intercept,np.max(height[valid])*slope+intercept],ls='-',lw=2.8,c='black')
  plt.plot([np.min(height[valid]),np.max(height[valid])],[np.min(height[valid])*slope+intercept,np.max(height[valid])*slope+intercept],ls='-',lw=2,c=color_ref_dict[_tr])
  print((tr,'R-squared = ' + str(r_value**2),'P value = ' + str(p_value)))
  legend_elements.append(Line2D([0], [0], color=color_ref_dict[_tr], label=un_treat_label[_tr] + ', R$^2=$' +  str(round(r_value**2,2)),lw=2,ls='-'))

slope, intercept, r_value, p_value, std_err = stats.linregress(height[valid],density[valid])
plt.legend(handles=legend_elements,fontsize=8)

plt.savefig('figs/figure_4.png',dpi=200,bbox_inches='tight')


