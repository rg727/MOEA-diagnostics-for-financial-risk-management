import numpy as np
#import string
#import brewer2mpl
#import seaborn as sns
import matplotlib.pyplot as plt
#import matplotlib.cbook as cbook
#sns.set(style='ticks', palette='Set2_r')

# cbar=brewer2mpl.get_map('PRGn', 'diverging', 11).mpl_colormap
#labels=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q']
labels=['Borg','MOEAD','NSGAII','NSGAIII','RVEA']
fig,ax=plt.subplots(1)

y= np.linspace(0,100,101)
x=np.arange(0,6)
#z=np.loadtxt('./hv_files50MC/attainment_vectors.txt')
z=np.loadtxt('./Attainment/attainment_vectors.txt')
z=z[::-1]

#print np.size(z,0), np.size(z,1), np.size(x,0),x
#print z[:,0]
best=np.loadtxt('Attainment/best_epsilon.txt') 
colormap='RdBu'
#plt.pcolormesh(fig, ax, np.loadtxt('attainment_vectors.txt'), 
               # xticklabels=string.uppercase[:100])
plt.pcolormesh(x,y,z, cmap=colormap) #, vmin=z_min, vmax=z_max)

color='#636363'
plt.title('Epsilon Indicator', fontsize=20)
# set the limits of the plot to the limits of the data

#plt.colorbar()
cbar=plt.colorbar()
cbar.set_label("Probability of attainment",size=15)
font_size = 16 # Adjust as appropriate.
cbar.ax.tick_params(labelsize=font_size)
#cbar.set_label(size=18)

#x2=np.linspace(0.5,16.5,17)
#plt.scatter(x2,best*100.5, s=80, c="w")
#x2=np.linspace(0.5,1.5,2)
x2=np.linspace(0.5,4.5,5)
plt.scatter(x2,best*100.5, s=80, c="w")
plt.axis([x.min(), x.max(), y.min(), y.max()])
ax.set_xticks(x + 0.5, minor=False)
# fontsize = 50
# fontweight = 'bold'
#fontproperties = {'family':'sans-serif','sans-serif':['Helvetica'],'weight' : fontweight, 'size' : fontsize}

ax.set_xticklabels(labels[:])
ax.tick_params(labelsize=12, color=color)
ax.set_ylabel('% of best generational distance',size=15)



#x=[.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5]

#ax.scatter(x, np.loadtxt('./Attainment/best.txt'),s=50, c='k')
plt.savefig('./Attainment/pcolormesh_prettyplotlib_labels.png')
#plt.savefig('pcolormesh_prettyplotlib_labels.pdf')