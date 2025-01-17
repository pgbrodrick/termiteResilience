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
from matplotlib.patches import Patch




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



def single_plot(subsample_size):
    bootstrap_df = pd.read_csv('bootstrapped_results/ss_' + str(subsample_size) + '.csv',sep=',')

    
    density = np.array(bootstrap_df['mound_density'])
    height = np.array(bootstrap_df['mean_mound_height'])
    cover_g1 = np.array(bootstrap_df['cover_g1'])
    cover_g3 = np.array(bootstrap_df['cover_g3'])
    cover_b13 = np.array(bootstrap_df['cover_b13'])
    treatment = np.array(bootstrap_df['treatment'])
    reps = np.array(bootstrap_df['rep_poly_count'])
    
    un_treat = np.unique(treatment)
    un_treat = un_treat[un_treat != 'nature_reserve']
    
    bin_step = 0.20
    density_bins = np.arange(0.0,2.6,step=bin_step)
    for n in range(len(un_treat)):
      
      l_density = density[treatment == un_treat[n]]
    
      h,b = np.histogram(l_density,bins=density_bins)
      b,h = plot_lhist(h / np.max(h),b)
      plt.plot(h + float(n) / 15. * bin_step,b,c=treat_colors[n])
    
    plt.xlabel('Termite Mound Density [mounds ha$^{-1}$]')
    plt.ylabel('Relative Frequency')
    
    plt.title('Bootstrap size: \n' + str(subsample_size**2/10000.) +  ' ha')
    #plt.legend(un_treat_label)



fig = plt.figure(figsize=(10,16))

treat_colors = [[255,170,0],[190,232,255],[56,168,0],[239,138,98]]
treat_colors = [(np.array(x)/255.).tolist() for x in treat_colors]
un_treat_label = ['Subsistance\nAgriculture','Communal\nGrazing', 'Kruger NP','Private\nReserve']
########### SI Figure 1
ss_list = [x*100 for x in range(4,20,2)]
edge_buffer = 0.06
n_rows = np.ceil(len(ss_list)/2.)
plot_y_size = (1 - edge_buffer*(1+ n_rows )) / n_rows
plot_x_size = plot_y_size*2

for _s in range(len(ss_list)):

    if (_s % 2 == 0):
        x_pos = edge_buffer
        y_pos = edge_buffer * (n_rows - np.floor(_s / 2)) + plot_y_size * (n_rows - np.floor(_s/2) - 1)
    else:
        x_pos = edge_buffer * 2 + plot_x_size

    ax = fig.add_axes([x_pos, y_pos, plot_x_size, plot_y_size],zorder=1)
    single_plot(ss_list[_s])

    if (_s == 3):
      plt.legend(un_treat_label)

plt.savefig('figs/figure_si_1.png',dpi=200,bbox_inches='tight')



########### SI Figure 2 / 3


def plot_rk(f,size_thresh):
  bootstrap_df = np.array(pd.read_csv(f))
  bootstrap_df = bootstrap_df[:,[0,2]]
  bootstrap_df = bootstrap_df[np.isnan(bootstrap_df[:,1]) == False,:]

  bootstrap_df = bootstrap_df[bootstrap_df[:,0] < 250,:]
  bootstrap_df[:,1] = np.power(bootstrap_df[:,1]/np.pi,0.5)

  rk = bootstrap_df

  plt.xlabel('Distance [m]')
  if (_t == 0):
    plt.ylabel('Ripley\'s K\n(L transformation)')

  x = bin_vals(rk[:,0])
  y = bin_vals(rk[:,1])
  plt.plot([0,250],[0,250],c='red',ls='--',lw=1,zorder=0)
  plt.plot(x,y,lw=1,c='black',zorder=1) 

  #plt.axis('equal')
  plt.xlim([-5,255])
  plt.ylim([-5,255])
  plt.xticks([0,50,100,150,200,250],rotation=90)
  plt.yticks([0,50,100,150,200,250])
  if (size_thresh):
    plt.text(10,200,str(loc_poly[m]))
  else:
    plt.text(10,200,str(loc_poly[m]),color='red')




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




key_bootstrap_df = pd.read_csv('polygon_info/polygon_key.csv',sep=',')
treatment = np.array(key_bootstrap_df['treatment'])
polygon_number = np.array(key_bootstrap_df['polygon_number'])

un_treat = np.unique(treatment)
un_treat = un_treat[un_treat != 'nature_reserve']






plt.figure(figsize=(8,28))
gs1 = gridspec.GridSpec(14,4)
gs1.update(hspace=0.8,wspace=0.5)

ind = np.zeros(4).astype(int)
for _t in range(0,len(un_treat)):
  loc_poly = polygon_number[treatment == un_treat[_t]]
  for m in range(len(loc_poly)):
    #if (loc_poly[m] != 31 and poly_size_area[poly_size_numbers.index(loc_poly[m])]/10000. > 300):
    if (loc_poly[m] != 31 ):

      if (ind[_t] <= 5):
        ax = plt.subplot(gs1[ind[_t],_t])
        fname = 'rk_runs/roi_output/poly_'+str(loc_poly[m]) + '.0.csv' 

        size_thresh = poly_size_area[poly_size_numbers.index(loc_poly[m])]/10000. > 300
        plot_rk(fname,size_thresh)

        if(ind[_t] == 0):
          plt.title(un_treat_label[_t])

        ind[_t] += 1
 

plt.savefig('figs/figure_si_2.png',dpi=200,bbox_inches='tight')
plt.clf()


plt.figure(figsize=(8,28))
gs1 = gridspec.GridSpec(14,4)
gs1.update(hspace=0.8,wspace=0.5)

ind = np.zeros(4).astype(int)
for _t in range(0,len(un_treat)):
  loc_poly = polygon_number[treatment == un_treat[_t]]
  for m in range(len(loc_poly)):
    #if (loc_poly[m] != 31 and poly_size_area[poly_size_numbers.index(loc_poly[m])]/10000. > 300):
    if (loc_poly[m] != 31 ):

      if (ind[_t] > 5):
        ax = plt.subplot(gs1[ind[_t]-5,_t])
        fname = 'rk_runs/roi_output/poly_'+str(loc_poly[m]) + '.0.csv' 

        size_thresh = poly_size_area[poly_size_numbers.index(loc_poly[m])]/10000. > 300
        plot_rk(fname,size_thresh)

        if(ind[_t]-6 == 0):
          plt.title(un_treat_label[_t])

      ind[_t] += 1
 

plt.savefig('figs/figure_si_3.png',dpi=200,bbox_inches='tight')






############# SI Figure 4 ################

# Use the selected 100 ha bootstrap size

bootstrap_df = pd.read_csv('bootstrapped_results/ss_1000.csv',sep=',')
density = np.array(bootstrap_df['mound_density'])
height = np.array(bootstrap_df['mean_mound_height'])
cover_g1 = np.array(bootstrap_df['cover_g1'])
cover_g3 = np.array(bootstrap_df['cover_g3'])
cover_b13 = np.array(bootstrap_df['cover_b13'])
treatment = np.array(bootstrap_df['treatment'])
reps = np.array(bootstrap_df['rep_poly_count'])

ax_s = 0.35
ax_b = 0.075
colors = np.ones((len(treatment),3))
for i in range(0,len(un_treat)):
  colors[treatment == un_treat[i],:] = treat_colors[i]


fig = plt.figure(figsize=(10,10))

def make_legend(r2_values):
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=un_treat_label[x].replace('\n',' ') + ', R$^2$ = ' + str(round(r2_values[x],2)),
                   markerfacecolor=treat_colors[x], markersize=8) for x in range(len(un_treat))]
    return legend_elements


ax = fig.add_axes([ax_b,ax_b,ax_s,ax_s])
plt.scatter(cover_g1,height,c=colors,s=1)
r2_values = []
for _n in range(len(un_treat_label)):
    valid = np.logical_and(np.isnan(height) == False, un_treat[_n] == treatment)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cover_g1[valid],height[valid])
    print((un_treat_label[_n],r_value**2))
    r2_values.append(r_value**2)

ax.legend(handles=make_legend(r2_values), labelspacing=0.2, fontsize=8)
plt.xlabel('Land Cover Fraction with Vegetation > 1 m')
plt.ylabel('Mean Mound Height [m]')
plt.ylim([0,2])

ax = fig.add_axes([ax_b,(ax_b*2+ax_s),ax_s,ax_s])
plt.scatter(cover_g1,density,c=colors,s=1)
r2_values = []
slope, intercept, r_value, p_value, std_err = stats.linregress(cover_g1,density)
print('cover > 1, density')
plt.xlabel('Land Cover Fraction with Vegetation > 1 m')
plt.ylabel('Mound Density [Mounds ha$^{-1}$]')
for _n in range(len(un_treat_label)):
    valid = np.logical_and(np.isnan(density) == False, un_treat[_n] == treatment)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cover_g1[valid],density[valid])
    print((un_treat_label[_n],r_value**2))
    r2_values.append(r_value**2)

plt.ylim([0,2])
ax.legend(handles=make_legend(r2_values), labelspacing=0.2, fontsize=8)

#plt.legend([x.replace('\n',' ') for x in un_treat_label])
#legend = ax.legend([x.replace('\n',' ') for x in un_treat_label], frameon=True)
#lgnd = plt.legend([x.replace('\n',' ') for x in un_treat_label])
#ax.legend(markerscale=5)

ax = fig.add_axes([(ax_b*2+ax_s),ax_b,ax_s,ax_s])
plt.scatter(cover_g3,height,c=colors,s=1)
r2_values = []
slope, intercept, r_value, p_value, std_err = stats.linregress(cover_g3[np.isnan(height) == False],height[np.isnan(height) == False])
print('cover > 3, height')
plt.xlabel('Land Cover Fraction with Vegetation > 3 m')
plt.ylabel('Mean Mound Height [m]')
for _n in range(len(un_treat_label)):
    valid = np.logical_and(np.isnan(height) == False, un_treat[_n] == treatment)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cover_g3[valid],height[valid])
    print((un_treat_label[_n],r_value**2))
    r2_values.append(r_value**2)

plt.ylim([0,2])
ax.legend(handles=make_legend(r2_values), labelspacing=0.2, fontsize=8)

ax = fig.add_axes([(ax_b*2+ax_s),(ax_b*2+ax_s),ax_s,ax_s])
plt.scatter(cover_g3,density,c=colors,s=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(cover_g3,density)
r2_values = []
print('cover > 3, density')
plt.xlabel('Land Cover Fraction with Vegetation > 3 m')
plt.ylabel('Mound Density [Mounds ha$^{-1}$]')
for _n in range(len(un_treat_label)):
    valid = np.logical_and(np.isnan(density) == False, un_treat[_n] == treatment)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cover_g3[valid],density[valid])
    print((un_treat_label[_n],r_value**2))
    r2_values.append(r_value**2)

plt.ylim([0,2])
ax.legend(handles=make_legend(r2_values), labelspacing=0.2, fontsize=8)

ax = fig.add_axes([(ax_b*3+ax_s*2),ax_b,ax_s,ax_s])
plt.scatter(cover_b13,height,c=colors,s=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(cover_b13[np.isnan(height) == False],height[np.isnan(height) == False])
r2_values = []
print('cover (> 1 and <=3), height')
plt.xlabel('Cover (> 1 and <= 3) Fraction')
plt.xlabel('Land Cover Fraction with Vegetation Between 1 and 3 m')
plt.ylabel('Mean Mound Height [m]')
for _n in range(len(un_treat_label)):
    valid = np.logical_and(np.isnan(height) == False, un_treat[_n] == treatment)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cover_b13[valid],height[valid])
    print((un_treat_label[_n],r_value**2))
    r2_values.append(r_value**2)

plt.ylim([0,2])
ax.legend(handles=make_legend(r2_values), labelspacing=0.2, fontsize=8)

ax = fig.add_axes([(ax_b*3+ax_s*2),(ax_b*2+ax_s),ax_s,ax_s])
plt.scatter(cover_b13,density,c=colors,s=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(cover_b13,density)
r2_values = []
print('cover (> 1 and <=3), density')
plt.xlabel('Land Cover Fraction with Vegetation Between 1 and 3 m')
plt.ylabel('Mound Density [Mounds ha$^{-1}$]')
for _n in range(len(un_treat_label)):
    valid = np.logical_and(np.isnan(density) == False, un_treat[_n] == treatment)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cover_b13[valid],density[valid])
    print((un_treat_label[_n],r_value**2))
    r2_values.append(r_value**2)

plt.ylim([0,2])
ax.legend(handles=make_legend(r2_values), labelspacing=0.2, fontsize=8)

plt.savefig('figs/figure_si_4.png',dpi=100,bbox_inches='tight')
plt.clf()






############# Figure SI 5
bootstrap_df = pd.read_csv('bootstrapped_results/ss_1000.csv',sep=',')
rem = np.array(bootstrap_df['mean_rem'])
tch = np.array(bootstrap_df['mean_tch'])
treatment = np.array(bootstrap_df['treatment'])


fig = plt.figure(figsize=(12,5))


################# 5a
ax = fig.add_axes([0.05,0.98,0.02,0.02],zorder=1)
plt.text(0,0,'a',fontweight='bold')
plt.axis('off')
ax = fig.add_axes([0.1,0.1,0.4,0.9],zorder=0)

ss = 0.10
density_bins = np.arange(0.0,2.6,step=ss)
c = 0
for n in range(len(un_treat)):
  
  l_density = tch[treatment == un_treat[n]]

  h,b = np.histogram(l_density,bins=density_bins)
  b,h = plot_lhist(h / np.max(h),b)
  plt.plot(h + float(c) / 15. * ss,b,c=treat_colors[n])
  c+=1

plt.xlabel('Top of Canopy Height [m ha$^{-1}$]')
plt.ylabel('Relative Frequency')

plt.legend([x.replace('\n',' ') for x in un_treat_label], frameon=True, loc='upper right')



################# 5b

ax = fig.add_axes([0.51,0.98,0.02,0.02],zorder=1)
plt.text(0,0,'b',fontweight='bold')
plt.axis('off')
ax = fig.add_axes([0.55,0.1,0.4,0.9],zorder=0)

ss = 0.05
density_bins = np.arange(0.0,1.1,step=ss)
c = 0
for n in range(len(un_treat)):
  
  l_density = rem[treatment == un_treat[n]]

  h,b = np.histogram(l_density,bins=density_bins)
  b,h = plot_lhist(h / np.max(h),b)
  plt.plot(h + float(c) / 15. * ss,b,c=treat_colors[n])
  c+=1

plt.xlabel('Relative Elevation [unitless]')
plt.ylabel('Relative Frequency')

plt.legend([x.replace('\n',' ') for x in un_treat_label], frameon=True, loc='upper right')

plt.savefig('figs/figure_si_5.png',dpi=200,bbox_inches='tight')





