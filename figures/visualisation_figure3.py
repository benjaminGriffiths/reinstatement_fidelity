# ----------- Data Visualisation ----------- #
# Written by Benjamin J. Griffiths
# Created on Friday 26th October, 2018
# ------------------------------------------ #

# %% -- import modules ------------------------------------------------------ #
import matplotlib as mpl
import numpy as np
import pandas
import ptitprince as pt
import seaborn as sns
from matplotlib import pyplot
from matplotlib import transforms

# %% -- define functions ---------------------------------------------------- #
# plot raincloud
def custom_rainplot(data,colour,axes,fontname,labels,ylim,offset,pvalue):
    
    # get transform data
    trans   = axes.transData
    offset  = transforms.ScaledTranslation(offset,0,f.dpi_scale_trans)
    
    # plot violin
    axes=pt.half_violinplot(data = data,bw = "scott",inner = None, scale = "count",
                          width = 0.5, linewidth = 1, cut = 1, palette = colour,
                          ax = axes, edgecolor = [0,0,0])
    
    # plot single points
    axes=sns.swarmplot(data = data, edgecolor =[0,0,0], size = 1.5, 
                       transform = trans + offset, palette = colour,
                       ax = axes)
    
    # plot mean and confidence intervals
    axes=sns.boxplot(data = data, palette = colour, width = 0.1, ax = axes, linewidth = 1, fliersize = 1, whis = 3)
    
    # plot significance
    sig_offset = ylim[1]-(ylim[1]*0.05)
    for i in np.arange(0,np.size(pvalue)):
        if pvalue[i] < 0.001:
            pyplot.scatter(np.array([-0.05,0,0.05])+i,[sig_offset,sig_offset,sig_offset],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.01:
            pyplot.scatter(np.array([-0.025,0.025])+i,[sig_offset,sig_offset],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.05:
            pyplot.scatter(np.array([0])+i,[sig_offset],s=3,c='black',marker='*',edgecolors=None)
    
    # add horizontal line
    axes.axhline(y=0, xmin=-1, xmax=3, color=[0,0,0], linewidth = 1)
    
    # aesthetics
    axes.set_ylabel(labels['ylabel'],fontname=fontname,fontsize=7,labelpad=5,fontweight='light')   # add Y axis label
    axes.set_ylim(ylim)                  # set Y axis range to 0 -> 1
    axes.set_xlim(-0.65,-0.5+len(data.columns))                  # set Y axis range to 0 -> 1
    axes.tick_params(axis='x',          # change X tick parameters
                   which='both',          # affect both major and minor ticks
                   bottom=False,          # turn off bottom ticks
                   labelbottom=True,  # keep bottom labels
                   pad=2.5,
                   width=1)   
    axes.tick_params(axis='y',          # change X tick parameters
                       pad=3,
                       width=1,
                       length=2.5)
    axes.set_yticks(labels['yticks'])
    axes.set_xticklabels(labels['xticklabel'],fontname=fontname,fontweight='light',fontsize=6)
    axes.set_yticklabels(labels['yticklabel'],fontname=fontname,fontweight='light',fontsize=6)

    # change axes
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['left'].set_linewidth(1)
              
def custom_timeseriesplot(data,variables,axes,colour,labels,xlim,ylim,xtick,vertical,horizontal):
    
    sns.lineplot(x=variables['x'],
             y=variables['y'],
             data=data,
             ci='sem',
             hue=variables['condition'],
             palette=colour,
             ax = axes,
             linewidth=1)
        
    # add horizontal line
    if vertical:
        ax.axvline(x=0, linewidth = 1, color = [0,0,0], linestyle='--')
    if horizontal:
        ax.axhline(y=0, linewidth = 1, color = [0,0,0], linestyle='-')
    
    # aesthetics
    ax.set_ylabel(labels['ylabel'],fontname='Calibri',fontsize=6,labelpad=-5,fontweight='light')   # add Y axis label
    ax.set_xlabel(labels['xlabel'],fontname='Calibri',fontsize=6,labelpad=3,fontweight='light')   # add Y axis label
    ax.set_ylim(ylim)                  # set Y axis range to 0 -> 1
    ax.set_xlim(xlim)                  # set Y axis range to 0 -> 1  
    ax.set_yticks([ylim[0],0,ylim[1]])
    ax.set_xticks(xtick)
    ax.set_yticklabels([ylim[0],0,ylim[1]],fontname='Calibri',fontweight='light',fontsize=5)
    ax.set_xticklabels(xtick,fontname='Calibri',fontweight='light',fontsize=5)
    ax.tick_params(axis='both',          # change X tick parameters
                       pad=3,
                       length=2.5)
    # change axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # remove legend
    ax.get_legend().remove()
    
# %% -- define key parameters ----------------------------------------------- #
# get current directory
wdir = 'E:/bjg335/projects/reinstatement_fidelity/'

# set context
sns.set_context("paper")

# set plot defaults
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.color'] = [0,0,0]
mpl.rcParams['ytick.color'] = [0,0,0]
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['axes.linewidth'] = 1

# %% -- prepare data -------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "data/fig3_data/group_task-all_eeg-cluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_raincloud = data_raincloud.drop(labels = ['forgotten','per_noBold','ret_noBold','ret_noConf'], axis = 1)

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(5/2.54) # 4inches 
f.set_figwidth(4.2/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Greens",n_colors=7)
colour = [colour[1],colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'Power-Similarity Correlation (z)',
          'xticklabel':['Perception','Retrieval'],
          'yticks':[-0.4,-0.2,0,0.2,0.4],
          'yticklabel':['-0.4','-0.2','0','0.2','0.4']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'calibri',labels,[-0.4,0.4],0.125,[1,1])
   
# save image
pyplot.savefig(wdir + "/figures/fig3a.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% --- correlation plot
# define raincloud data filename
data_fname = wdir + "data/fig3_data/subj-01_task-ers_xy.csv"

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.3/2.54) # 4inches 
f.set_figwidth(2.3/2.54) # 12inches
f.set_dpi(300)

# load raincloud data
data_corr = pandas.read_csv(data_fname,
                       delimiter=",",
                       header=None,
                       names=['X','Y'])

# define colour scheme
colour = sns.color_palette("Greens",n_colors=7)
colour = [colour[4]]

# draw scatter
sns.scatterplot(x='X',y='Y',data=data_corr,linewidth=0,alpha = 0.7,ax=ax,**{'s':3,'edgecolors':colour,'c':colour})
sns.regplot(x='X',y='Y',data=data_corr,ax=ax,scatter=False,color=colour[0])

# aesthetics
ax.set_ylabel('Similarity Index (z)',fontname='Calibri',fontsize=6,labelpad=5,fontweight='light')   # add Y axis label
ax.set_xlabel('Alpha/Beta Power (z)',fontname='Calibri',fontsize=6,labelpad=3,fontweight='light')   # add Y axis label

ax.set_ylim(-2,2.6)                  # set Y axis range to 0 -> 1
ax.set_xlim(-2.5,2.1)                  # set Y axis range to 0 -> 1  
ax.set_yticks([])
ax.set_xticks([])

# change axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# save image
pyplot.savefig(wdir + "/figures/fig3b.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% Corr 2
# define raincloud data filename
data_fname = wdir + "data/fig3_data/subj-02_task-ers_xy.csv"

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.3/2.54) # 4inches 
f.set_figwidth(2.3/2.54) # 12inches
f.set_dpi(300)

# load raincloud data
data_corr = pandas.read_csv(data_fname,
                       delimiter=",",
                       header=None,
                       names=['X','Y'])

# define colour scheme
colour = sns.color_palette("Greens",n_colors=7)
colour = [colour[4]]

# draw scatter
sns.scatterplot(x='X',y='Y',data=data_corr,linewidth=0,alpha = 0.7,ax=ax,**{'s':3,'edgecolors':colour,'c':colour})
sns.regplot(x='X',y='Y',data=data_corr,ax=ax,scatter=False,color=colour[0])

# aesthetics
ax.set_ylabel('Similarity Index (z)',fontname='Calibri',fontsize=6,labelpad=5,fontweight='light')   # add Y axis label
ax.set_xlabel('Alpha/Beta Power (z)',fontname='Calibri',fontsize=6,labelpad=3,fontweight='light')   # add Y axis label

ax.set_ylim(-2.3,2.6)                  # set Y axis range to 0 -> 1
ax.set_xlim(-2.6,3.1)                  # set Y axis range to 0 -> 1  
ax.set_yticks([])
ax.set_xticks([])

# change axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# save image
pyplot.savefig(wdir + "/figures/fig3c.jpg",bbox_inches='tight',transparent=True,dpi='figure')
