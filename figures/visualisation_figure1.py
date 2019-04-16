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
    axes=sns.boxplot(data = data, palette = colour, width = 0.1, ax = axes, linewidth = 1, fliersize = 1, whis = 5)
    
    # plot significance
    sig_offset = ylim[1]-(ylim[1]*0.05)
    for i in np.arange(0,np.size(pvalue)):
        if pvalue[i] < 0.001:
            pyplot.scatter(np.array([-0.05,0,0.05])+i,[sig_offset,sig_offset,sig_offset],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.01:
            pyplot.scatter(np.array([-0.025,0.025])+i,[sig_offset,sig_offset],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.05:
            pyplot.scatter(np.array([0])+i,[sig_offset],s=3,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
    
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
                   pad=-2,
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
data_fname = wdir + "data/fig1_data/ers_betas.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# read in percept surface images
img = [pyplot.imread(wdir + "/data/fig1_surface/left_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/left_medial.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/right_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/right_medial.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/both_ventral.tiff")]

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(5/2.54) # 4inches 
f.set_figwidth(4.2/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Reds",n_colors=7)
colour = [colour[2],colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'Similarity Index (z)',
          'xticklabel':['Left Fusiform','Right Fusiform'],
          'yticks':[-0.1,0,0.1,0.2],
          'yticklabel':['-0.1','0','0.1','0.2']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'calibri',labels,[-0.1,0.2],0.125,[1,1])
   
# save image
pyplot.savefig(wdir + "/figures/fig1a.tif",bbox_inches='tight',transparent=True,dpi='figure')
