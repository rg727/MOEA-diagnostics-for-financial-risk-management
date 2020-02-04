import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import pandas
import seaborn.apionly as sns

C1 = np.loadtxt('portfolio.ref')



def makeFigure4():
    '''Makes Parallel Axis plots for Hetch Hetchy DPS Formulation)'''
    
    sns.set_style("dark")
    
    # load thinned reference sets from each problem formulation
    C1 = np.loadtxt('portfolio.ref')
    C2=pandas.read_csv('NSGAII_portfolio.set',delimiter=' ',header=None)
    C2.columns = ['Annual', 'Min', 'Hedge', 'Complexity']
    C2=C2[C2.Min<-115]
    np.savetxt('NSGAII_portfolio_req_1.txt', C2, delimiter=',')
    C2=np.loadtxt('NSGAII_portfolio_req_1.txt',delimiter=',')
    #C3 = np.loadtxt('./../EV/EV_thinned.csv',delimiter=',',skiprows=1)
    #C4 = np.loadtxt('./../EVSDH/EVSDH_thinned.csv',delimiter=',',skiprows=1)
    
    # set plotting characteristics
    formulations = [C1]
    formulations_1=[C2]
    labels = [['Annualized\nAdjusted Revenue\n($M)','Minimum\nAdjusted Revenue\n($M)','Maximum\nHedge Complexity','Maximum\nFund Balance\n($M)']]
        #['WP1 Hydro\n(Gwh/day)','WP1 Deficit$\mathregular{^2}$\n(m$\mathregular{^3}\!$/s)$\mathregular{^2}$','WP1 Flood\n(m above 11.25 m)','WP1 Recovery\n(days)'],\
        #['EV Hydro\n(Gwh/day)','EV Deficit$\mathregular{^2}$\n(m$\mathregular{^3}\!$/s)$\mathregular{^2}$','WP1 Flood\n(m above 11.25 m)','EV Recovery\n(days)'],\
        #['EV Hydro\n(Gwh/day)','EV Deficit$\mathregular{^2}$\n(m$\mathregular{^3}\!$/s)$\mathregular{^2}$','WP1 Flood\n(m above 11.25 m)','EV Recovery\n(days)','SD Hydro\n(Gwh/day)']]
    #cmaps = ['Reds_r','Blues_r','Greens_r','Purples_r']
    cmaps = ['Reds_r']
    titles = ['NSGA-II']
    precision = [[0,0,2,1]]
    
    # make 2 x 2 subplot with parallel axes for each problem formulation
    plot(formulations, formulations_1, labels, precision, cmaps, titles, 'NSGAII_first_requirement.pdf')
    
    return None
    
def plot(formulations, formulations_1,labels, precision, cmaps, titles, filename):
    fig = plt.figure()
    shadeIndex = [0,0,0,0,0,0,0,0,0,0]
    for i in range(1):
        ax = fig.add_subplot(1,1,1)
        table = pandas.DataFrame(formulations[i],columns=[labels[i]])
        mins = np.min(formulations[i],0)
        maxs = np.max(formulations[i],0)
        # round number of significant digits shown on objective labels
        for j in range(len(labels[i])):
            if precision[i][j] != 0:
                labels[i][j] = str(np.round(mins[j],precision[i][j])) + '\n' + labels[i][j]
            else:
                labels[i][j] = str(int(mins[j]))+ '\n' + labels[i][j]
            # don't show negative sign on maximization objectives
            if mins[j] < 0:
                labels[i][j] = labels[i][j][1:]
        parallel_coordinate(ax, table, mins, maxs, 'gray', shadeIndex[i], titles[i], labels[i], precision[i])
       
    for i in range(1):
        ax = fig.add_subplot(1,1,1)
        table_1 = pandas.DataFrame(formulations_1[i],columns=[labels[i]])
        mins_1 = np.min(formulations_1[i],0)
        maxs_1 = np.max(formulations_1[i],0)
        # round number of significant digits shown on objective labels
        for j in range(len(labels[i])):
            if precision[i][j] != 0:
                labels[i][j] = str(np.round(mins[j],precision[i][j])) + '\n' + labels[i][j]
            else:
                labels[i][j] = str(int(mins[j]))+ '\n' + labels[i][j]
            # don't show negative sign on maximization objectives
            if mins[j] < 0:
                labels[i][j] = labels[i][j][1:]
        parallel_coordinate_2(ax, table_1, mins, maxs, 'Oranges', shadeIndex[i], titles[i], labels[i], precision[i])
         
        
    fig.set_size_inches([13,10])
    fig.tight_layout()
    fig.savefig(filename)
    
    
    return None
    
def parallel_coordinate(ax, table, mins, maxs, cmap, shadeIndex, \
    title, xlabels, precision):
        
    toplabels = []
    # round number of significant digits shown on objective labels
    for i in range(len(xlabels)):
        if precision[i] != 0:
            toplabels.append(str(np.round(maxs[i],precision[i])))
        else:
            toplabels.append(str(int(maxs[i])))
        if maxs[i] < 0:
            # don't show negative sign on maximization objectives
            toplabels[i] = toplabels[i][1:]
        
    cmap = cmap
    scaled = table.copy()
    index = 0
    for column in table.columns:
        scaled[column] = (table[column] - mins[index]) / (maxs[index] - mins[index])
        index = index + 1  
    
    for solution in scaled.iterrows():
        ys = solution[1]
        xs = range(len(ys))
        ax.plot(xs, ys, cmap, linewidth=1)
    scaled.to_pickle("./dummy.pkl")    
        
    
    ax.set_title(title, y=1.1)
    #ax.set_xticks(np.arange(0,np.shape(table)[1],1))
    ax.set_xlim([0,np.shape(table)[1]-1])
    ax.set_ylim([0,1])
    ax.set_xticklabels(xlabels)
    ax.tick_params(axis='y',which='both',labelleft='off',left='off',right='off')
    ax.tick_params(axis='x',which='both',top='off',bottom='off')
    
    # make subplot frames invisible
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #ax.plot(xs, scaled.iloc[567,:], c='r', linewidth=6)
    #print(scaled.iloc[567,:])

    
    # draw in axes
    for i in np.arange(0,np.shape(table)[1],1):
        ax.plot([i,i],[0,1],c='k')
    
    # create twin y axis to put x tick labels on top
    ax2 = ax.twiny()
    ax2.set_xticks(np.arange(0,np.shape(table)[1],1))
    ax2.set_xlim([0,np.shape(table)[1]-1])
    ax2.set_ylim([0,1])
    ax2.set_xticklabels(toplabels)
    ax2.tick_params(axis='y',which='both',labelleft='off',left='off',right='off')
    ax2.tick_params(axis='x',which='both',top='off',bottom='off')
    
    # make subplot frames invisible
    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    return ax,

def parallel_coordinate_2(ax, table, mins, maxs, cmap, shadeIndex, \
    title, xlabels, precision):
        
    toplabels = []
    # round number of significant digits shown on objective labels
    for i in range(len(xlabels)):
        if precision[i] != 0:
            toplabels.append(str(np.round(maxs[i],precision[i])))
        else:
            toplabels.append(str(int(maxs[i])))
        if maxs[i] < 0:
            # don't show negative sign on maximization objectives
            toplabels[i] = toplabels[i][1:]
        
    cmap = matplotlib.cm.get_cmap(cmap)
    scaled = table.copy()
    index = 0
    for column in table.columns:
        scaled[column] = (table[column] - mins[index]) / (maxs[index] - mins[index])
        index = index + 1
    
    for solution in scaled.iterrows():
        ys = solution[1]
        xs = range(len(ys))
        ax.plot(xs, ys, c=cmap(ys[shadeIndex]), linewidth=2)
    
    ax.set_title(title,y=1.1)
    ax.set_xticks(np.arange(0,np.shape(table)[1],1))
    ax.set_xlim([0,np.shape(table)[1]-1])
    ax.set_ylim([0,1])
    ax.set_ylabel("← Preference")
    #ax.set_xticklabels(xlabels,fontsize=14)
    ax.tick_params(axis='y',which='both',labelleft='off',left='off',right='off')
    ax.tick_params(axis='x',which='both',top='off',bottom='off')
    
    # make subplot frames invisible
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
#    # draw in axes
#    for i in np.arange(0,np.shape(table)[1],1):
#        ax.plot([i,i],[0,1],c='k')
    
#    # create twin y axis to put x tick labels on top
#    ax2 = ax.twiny()
#    ax2.set_xticks(np.arange(0,np.shape(table)[1],1))
#    ax2.set_xlim([0,np.shape(table)[1]-1])
#    ax2.set_ylim([0,1])
#    ax2.set_xticklabels(toplabels)
#    ax2.tick_params(axis='y',which='both',labelleft='off',left='off',right='off')
#    ax2.tick_params(axis='x',which='both',top='off',bottom='off')
    
#    # make subplot frames invisible
#    ax2.spines["top"].set_visible(False)
#    ax2.spines["bottom"].set_visible(False)
#    ax2.spines["left"].set_visible(False)
#    ax2.spines["right"].set_visible(False)
#    
    return ax

makeFigure4()

