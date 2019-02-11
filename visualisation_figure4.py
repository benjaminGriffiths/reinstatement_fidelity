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
    axes=sns.boxplot(data = data, palette = colour, width = 0.1, ax = axes, linewidth = 1, fliersize = 1)
    
    # plot significance
    for i in np.arange(0,np.size(pvalue)):
        if pvalue[i] < 0.001:
            pyplot.scatter(np.array([-0.05,0,0.05])+i,[0.19,0.19,0.19],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.01:
            pyplot.scatter(np.array([-0.025,0.025])+i,[0.19,0.19],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.05:
            pyplot.scatter(np.array([0])+i,[0.19],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
    
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

# %% -- define key parameters ----------------------------------------------- #
# get current directory
wdir = 'E:/bjg335/projects/reinstatement_fidelity/'

# set context
sns.set_context("paper")

# set plot defaults
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.color'] = [0,0,0]
mpl.rcParams['ytick.color'] = [0,0,0]
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.linewidth'] = 3

# %% -- prepare data -------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "data/fig1_data/percept_betas.csv"

# ----- RSA Plot ----- # 
# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_raincloud = data_raincloud.drop(labels = ['FrontCent','LeftTemp'], axis = 1)

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(6.3/2.54) # 4inches 
f.set_figwidth(3.5/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Reds",n_colors=7)
colour = [colour[2],colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'Encoding-Encoding Similarity Index (z)',
          'xticklabel':[''],
          'yticks':[-0.05,0,0.05,0.1,0.15],
          'yticklabel':['-0.05','0','0.05','0.1','0.15']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'calibri',labels,[-0.05,0.15],0.125,[0.0001])
   
# save image
pyplot.savefig(wdir + "/figures/fig4a.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- Power Plot ----- # 
# define raincloud data filename
data_fname = wdir + "data/fig2_data/group_task-rf_eeg-cluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_raincloud = data_raincloud.drop(labels = 'retrieval', axis = 1)

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(6.3/2.54) # 4inches 
f.set_figwidth(3.5/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'Post-Stim. > Pre-Stim. Alpha/Beta Power (z)',
          'xticklabel':[''],
          'yticks':[-0.5,-0.25,0,0.25],
          'yticklabel':['-0.5','-0.25','0','0.25']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'calibri',labels,[-0.5,0.25],0.125,[0.0001])
   
# save image
pyplot.savefig(wdir + "/figures/fig4a.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- Correlation Plot ----- # 
# define raincloud data filename
data_fname = wdir + "data/fig3_data/group_task-all_eeg-cluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_raincloud = data_raincloud.drop(labels = ['retrieval','forgotten','partial'], axis = 1)

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(6.3/2.54) # 4inches 
f.set_figwidth(3.5/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Greens",n_colors=7)
colour = [colour[2],(0.7,0.7,0.7)]

# define labels
labels = {'title':'',
          'ylabel':'Power-Similarity Correlation (z)',
          'xticklabel':[''],
          'yticks':[-0.4,-0.2,0,0.2,0.4],
          'yticklabel':['-0.4','-0.2','0','0.2','0.4']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'calibri',labels,[-0.4,0.4],0.125,[0.033])
   
# save image
pyplot.savefig(wdir + "/figures/fig3a.jpg",bbox_inches='tight',transparent=True,dpi='figure')
