# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:55:36 2019

@author: Rohini
"""

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import numpy as np
fmri = sns.load_dataset("fmri")
import pandas as pd 

ref_hypervolume=0.5614373077

data_borg = np.loadtxt('hypervolume_Borg.txt',delimiter=' ', skiprows=1)
data_new_borg=np.transpose(data_borg)
data_stacked=np.hstack(data_new_borg)
data_stacked=data_stacked/ref_hypervolume
NFE_borg=200
NFE_vector_borg=np.linspace(0,NFE_borg,data_borg.shape[0])
NFE_vector_borg_rep=np.tile(NFE_vector_borg, 50)

ax = sns.lineplot(NFE_vector_borg_rep, data_stacked)
#ax = sns.lineplot(x="timepoint", y="signal", data=fmri)

Borg=pd.DataFrame({'NFE':NFE_vector_borg_rep, 'Hypervolume':data_stacked})



data_MOEAD = np.loadtxt('hypervolume_MOEAD.txt',delimiter=' ', skiprows=1)
data_new_MOEAD=np.transpose(data_MOEAD)
data_stacked_MOEAD=np.hstack(data_new_MOEAD)
data_stacked_MOEAD=data_stacked_MOEAD/ref_hypervolume
NFE=200
NFE_vector_MOEAD=np.linspace(0,NFE,data_MOEAD.shape[0])
NFE_vector_MOEAD_rep=np.tile(NFE_vector_MOEAD, 50)

ax = sns.lineplot(NFE_vector_MOEAD_rep, data_stacked_MOEAD)

MOEAD=pd.DataFrame({'NFE':NFE_vector_MOEAD_rep, 'Hypervolume':data_stacked_MOEAD})


data_NSGAII = np.loadtxt('hypervolume_NSGAII.txt',delimiter=' ', skiprows=1)
data_new_NSGAII=np.transpose(data_NSGAII)
data_stacked_NSGAII=np.hstack(data_new_NSGAII)
data_stacked_NSGAII=data_stacked_NSGAII/ref_hypervolume
NFE=200
NFE_vector_NSGAII=np.linspace(0,NFE,data_NSGAII.shape[0])
NFE_vector_NSGAII_rep=np.tile(NFE_vector_NSGAII, 50)

ax = sns.lineplot(NFE_vector_NSGAII_rep, data_stacked_NSGAII)

NSGAII=pd.DataFrame({'NFE':NFE_vector_NSGAII_rep, 'Hypervolume':data_stacked_NSGAII})



data_NSGAIII = np.loadtxt('hypervolume_NSGAIII.txt',delimiter=' ', skiprows=1)
data_new_NSGAIII=np.transpose(data_NSGAIII)
data_stacked_NSGAIII=np.hstack(data_new_NSGAIII)
data_stacked_NSGAIII=data_stacked_NSGAIII/ref_hypervolume
NFE=200
NFE_vector_NSGAIII=np.linspace(0,NFE,data_NSGAIII.shape[0])
NFE_vector_NSGAIII_rep=np.tile(NFE_vector_NSGAIII, 50)

ax = sns.lineplot(NFE_vector_NSGAIII_rep, data_stacked_NSGAIII)

NSGAIII=pd.DataFrame({'NFE':NFE_vector_NSGAIII_rep, 'Hypervolume':data_stacked_NSGAIII})



data_RVEA = np.loadtxt('hypervolume_RVEA.txt',delimiter=' ', skiprows=1)
data_new_RVEA=np.transpose(data_RVEA)
data_stacked_RVEA=np.hstack(data_new_RVEA)
data_stacked_RVEA=data_stacked_RVEA/ref_hypervolume
NFE=200
NFE_vector_RVEA=np.linspace(0,NFE,data_RVEA.shape[0])
NFE_vector_RVEA_rep=np.tile(NFE_vector_RVEA, 50)

RVEA=pd.DataFrame({'NFE':NFE_vector_RVEA_rep, 'Hypervolume':data_stacked_RVEA})

ax = sns.lineplot(NFE_vector_RVEA_rep, data_stacked_RVEA)

set_RVEA=pd.DataFrame({'NFE':NFE_vector_RVEA_rep, 'Hypervolume':data_stacked_RVEA})


concatenated = pd.concat([Borg.assign(Algorithm='Borg'), NSGAII.assign(Algorithm='NSGA-II'),NSGAIII.assign(Algorithm='NSGA-III'),RVEA.assign(Algorithm='RVEA'),MOEAD.assign(Algorithm='MOEA/D')])
flatui = ["Blue", "Orange", "Green", "Purple", "Red"]

sns.set_palette(flatui)
    
ax = sns.lineplot(x='NFE',y='Hypervolume',data=concatenated,hue='Algorithm')    
sns.set(font_scale=1.5)
#ax.set_title('Hypervolume Runtime Dynamics',fontsize=15)
ax.set_ylabel('Hypervolume',size=14)
plt.setp(ax.get_legend().get_title(), fontsize='18')
plt.setp(ax.get_legend().get_texts(), fontsize='14') # for legend text
ax.set_xlabel('Number of Function Evaluations ($10^3$)', size=14)

ax.set_ylim(0,1)

fig.savefig("Hypervolume_Confidence_Interval.png")

