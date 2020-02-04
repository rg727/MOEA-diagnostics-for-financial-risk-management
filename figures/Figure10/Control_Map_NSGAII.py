import matplotlib.mlab as mlab # matlab compatibility functions
from matplotlib.backends import backend_agg as agg # raster backend
import numpy        # numeric library
import pandas       # data analysis library

table = pandas.read_table("FINAL_TABLE_NSGAII.txt", 
                          sep=' ', header=0, 
                          names=["x", "y", "z"])
#ticks = numpy.arange(0,10,0.1)

#x_ticks=numpy.arange(10,251)

ref_hypervolume=0.6118455239
x_ticks=numpy.linspace(10,250,50)

y_ticks=numpy.arange(5000,200000,50)

z=table.z.values/ref_hypervolume
grid = mlab.griddata(table.x.values, table.y.values, z, 
                     x_ticks, y_ticks,interp='linear')

# plot
fig = matplotlib.figure.Figure() # create the figure
agg.FigureCanvasAgg(fig)         # attach the rasterizer
ax = fig.add_subplot(1, 1, 1)    # make axes to plot on

#ax.set_xlabel("Initial Population Size")
#ax.set_ylabel("Number of Functional Evaluations")


cmap = matplotlib.cm.get_cmap("RdBu") # get the "hot" color map
contourset = ax.contourf(x_ticks, y_ticks, grid, 50, cmap=cmap,vmin=0.4, vmax=1)
#contourset.set_clim(vmin=0, vmax=1)
#bounds = [0, 50,100]
#norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

#cbar = fig.colorbar(contourset,norm=norm,boundaries=bounds,orientation='horizontal')
#v = numpy.linspace(0, 1, 15, endpoint=True)
#cbar = fig.colorbar(contourset)

#cbar.set_ticks([0,50,100])
sm.set_array([40, 100])
cbar = fig.colorbar(sm, ticks=[40, 50, 60, 70, 80, 90,100])
cbar.set_label('% of reference set hypervolume')
#fig.axes[-1].set_ylabel("% of reference set hypervolume") # last axes instance is the colorbar

ax.set_xlim(15,240)
#ax.set_yticks(list(numpy.arange(5000,205000, 20000)))
ax.set_ylim(5000,200000)
ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
ax.set_title("NSGAII",fontsize=18)
#contourset.clim(0,1)

fig.set_size_inches(7, 5)
#matplotlib.rcParams.update({'font.size': 8})

fig.savefig("NSGAII_Control_Map.png")

# vim:ts=4:sw=4:expandtab:ai:colorcolumn=68:number:fdm=indent