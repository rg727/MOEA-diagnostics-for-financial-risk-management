from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


 # read in reference sets and negate all objectives
reference_set = 1*np.loadtxt('portfolio.ref')

reference_set[:,1] = -reference_set[:,1]
reference_set[:,0] = -reference_set[:,0]
fig = plt.figure()
#plt.figure(figsize=(20,10))
ax = fig.add_subplot(111, projection='3d')



z = reference_set[:,0]
x = reference_set[:,1]
y = reference_set[:,2]
color = reference_set[:,3]



f=ax.scatter(x, y, z, c=color, marker='o',cmap='magma')
ax.set_facecolor('white')
#ax.set_zlabel('Average Annualized \n Adjusted Revenue ($M)',fontsize=14)
#ax.set_xlabel('Average Minimum \n Adjusted Revenue ($M)',fontsize=14)
#ax.set_ylabel('Average Maximum\n Hedging Complexity',fontsize=14)
ax.tick_params(axis = 'both', which = 'major', labelsize = 18)

ax.xaxis.labelpad = 20
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 40

ax.tick_params(axis='both', which='major', pad=15)


#cbar = plt.colorbar(f)
#cbar.set_label('Average Maximum Contingency Fund Balance ($M)',fontsize=14)
#cbar.ax.tick_params(labelsize=12)

ax.scatter(max(x),min(y),max(z) , c='black', s=500, linewidth=0, marker='*',)

cbar = plt.colorbar(f)
cbar.set_label('Average Maximum Contingency Fund Balance ($M)',fontsize=14)
cbar.ax.tick_params(labelsize=16)

fig.savefig("Best_Approximate_Pareto_Set.svg")

##########################################################



#######################################################################
#Best reference set for each algorithm
reference_set = 1*np.loadtxt('RVEA_portfolio.set')

reference_set[:,1] = -reference_set[:,1]
reference_set[:,0] = -reference_set[:,0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=18, azim=44)




z = reference_set[:,0]
x = reference_set[:,1]
y = reference_set[:,2]
color = reference_set[:,3]



f=ax.scatter(x, y, z, c=color, marker='o',cmap='magma')
ax.set_facecolor('white')
#ax.set_zlabel('Average Annualized \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_xlabel('Average Minimum \n Adjusted Revenue ($M)',fontsize=12)
#x.set_ylabel('Average Maximum\n Hedging Complexity',fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 10)


ax.xaxis.labelpad = 20
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 20


#cbar = plt.colorbar(f)
#cbar.set_label('Average Maximum Contingency Fund Balance ($M)',fontsize=14)
#cbar.ax.tick_params(labelsize=14)

ax.set_title("RVEA",fontsize=16)

fig.savefig("RVEA_REFERENCE_SET.pdf")




plt.show()

#########################################################################

#######################################################################
#Best reference set for each algorithm
reference_set = 1*np.loadtxt('MOEAD_portfolio.set')

reference_set[:,1] = -reference_set[:,1]
reference_set[:,0] = -reference_set[:,0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=18, azim=44)


z = reference_set[:,0]
x = reference_set[:,1]
y = reference_set[:,2]
color = reference_set[:,3]



f=ax.scatter(x, y, z, c=color, marker='o',cmap='magma')
ax.set_facecolor('white')
#ax.set_zlabel('Average Annualized \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_xlabel('Average Minimum \n Adjusted Revenue ($M)',fontsize=12)
#x.set_ylabel('Average Maximum\n Hedging Complexity',fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 10)


ax.xaxis.labelpad = 20
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 20


#cbar = plt.colorbar(f)
#cbar.set_label('Average Maximum Contingency Fund Balance ($M)',fontsize=14)
#cbar.ax.tick_params(labelsize=14)

ax.set_title("MOEA/D",fontsize=16)

fig.savefig("MOEAD_REFERENCE_SET.pdf")




plt.show()

#################################################################
#Best reference set for each algorithm
reference_set = 1*np.loadtxt('NSGAII_portfolio.set')

reference_set[:,1] = -reference_set[:,1]
reference_set[:,0] = -reference_set[:,0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=18, azim=44)


z = reference_set[:,0]
x = reference_set[:,1]
y = reference_set[:,2]
color = reference_set[:,3]



f=ax.scatter(x, y, z, c=color, marker='o', cmap='magma')
ax.set_facecolor('white')
#ax.set_zlabel('Average Annualized \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_xlabel('Average Minimum \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_ylabel('Average Maximum\n Hedging Complexity',fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 10)

ax.xaxis.labelpad = 20
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 20


#cbar = plt.colorbar(f)
#cbar.set_label('Average Maximum Contingency Fund Balance ($M)',fontsize=14)
#cbar.ax.tick_params(labelsize=14)


ax.set_title("NSGA-II",fontsize=16)

fig.savefig("NSGAII_REFERENCE_SET.pdf")




plt.show()

#############################################################################
#Best reference set for each algorithm
reference_set = 1*np.loadtxt('Borg_portfolio.set')

reference_set[:,1] = -reference_set[:,1]
reference_set[:,0] = -reference_set[:,0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=18, azim=44)


z = reference_set[:,0]
x = reference_set[:,1]
y = reference_set[:,2]
color = reference_set[:,3]



f=ax.scatter(x, y, z, c=color, marker='o', cmap='magma')
ax.set_facecolor('white')
#ax.set_zlabel('Average Annualized \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_xlabel('Average Minimum \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_ylabel('Average Maximum\n Hedging Complexity',fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 10)

ax.xaxis.labelpad = 20
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 20


#cbar = plt.colorbar(f)
#cbar.set_label('Average Maximum Contingency Fund Balance ($M)',fontsize=14)
#cbar.ax.tick_params(labelsize=14)


ax.set_title("Borg",fontsize=16)

fig.savefig("Borg_REFERENCE_SET.pdf")



plt.show()
##############################################################################
#Best reference set for each algorithm
reference_set = 1*np.loadtxt('NSGAIII_portfolio.set')

reference_set[:,1] = -reference_set[:,1]
reference_set[:,0] = -reference_set[:,0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=18, azim=44)


z = reference_set[:,0]
x = reference_set[:,1]
y = reference_set[:,2]
color = reference_set[:,3]



f=ax.scatter(x, y, z, c=color, marker='o', cmap='magma')
ax.set_facecolor('white')
#ax.set_zlabel('Average Annualized \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_xlabel('Average Minimum \n Adjusted Revenue ($M)',fontsize=12)
#ax.set_ylabel('Average Maximum\n Hedging Complexity',fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 10)

ax.xaxis.labelpad = 20
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 20
ax.set_zlim(125, 127.5)


#cbar = plt.colorbar(f)
#cbar.set_label('Average Maximum Contingency Fund Balance ($M)',fontsize=14)
#cbar.ax.tick_params(labelsize=14)


ax.set_title("NSGA-III",fontsize=16)

fig.savefig("NSGAIII_REFERENCE_SET.pdf")



plt.show()