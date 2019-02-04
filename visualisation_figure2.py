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
def custom_rainplot(data,colour,axes,fontname,labels,ylim,offset):
    
    # get transform data
    trans   = axes.transData
    offset  = transforms.ScaledTranslation(offset,0,f.dpi_scale_trans)
    
    # plot violin
    axes=pt.half_violinplot(data = data,bw = "scott",inner = None, scale = "count",
                          width = 0.5, linewidth = 3, cut = 1, palette = colour,
                          ax = axes, edgecolor = [0,0,0])
    
    # plot single points
    axes=sns.swarmplot(data = data, edgecolor =[0,0,0], size = 5, 
                       transform = trans + offset, palette = colour,
                       ax = axes)
    
    # plot mean and confidence intervals
    axes=sns.boxplot(data = data, palette = colour, width = 0.1, ax = axes, linewidth = 3)
    
    # add horizontal line
    axes.axhline(y=0, xmin=-1, xmax=3, color=[0,0,0], linewidth = 3)
    
    # aesthetics
    axes.set_title(labels['title'],pad=25,fontname=fontname,fontsize=28,fontweight='light')
    axes.set_ylabel(labels['ylabel'],fontname=fontname,fontsize=24,labelpad=15,fontweight='light')   # add Y axis label
    axes.set_ylim(ylim)                  # set Y axis range to 0 -> 1
    axes.set_xlim(-0.65,-0.5+len(data.columns))                  # set Y axis range to 0 -> 1
    axes.tick_params(axis='x',          # change X tick parameters
                   which='both',          # affect both major and minor ticks
                   bottom=False,          # turn off bottom ticks
                   labelbottom=True,  # keep bottom labels
                   pad=10,
                   width=3)   
    axes.tick_params(axis='y',          # change X tick parameters
                       pad=10)
    axes.set_yticks(labels['yticks'])
    axes.set_xticklabels(labels['xticklabel'],fontname=fontname,fontweight='light')
    axes.set_yticklabels(labels['yticklabel'],fontname=fontname,fontweight='light')

    # change axes
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['left'].set_linewidth(3)
    
               
# plot brain surface plots
def custom_surfaceplot(img,axes,colour,title,xrange):
    
    # plot surface images
    axes[0,0].imshow(img[0],interpolation="bicubic")
    axes[0,0].axis('off')
    
    axes[1,0].imshow(img[1],interpolation="bicubic")
    axes[1,0].axis('off')
    
    axes[0,1].imshow(img[2],interpolation="bicubic")
    axes[0,1].axis('off')
    
    axes[1,1].imshow(img[3],interpolation="bicubic")
    axes[1,1].axis('off')
    
    axes[0,2].imshow(img[4],interpolation="bicubic")
    axes[0,2].axis('off')
    
    # plot colourbar
    cX = np.tile(np.linspace(0,1,128),[12,1])
    axes[1,2].imshow(X=cX,
                 cmap=colour,
                 aspect='equal')
    axes[1,2].plot(np.array([0,0]),np.array([0,12]),color='black',linewidth=3)
    axes[1,2].plot(np.array([128,128]),np.array([0,12]),color='black',linewidth=3)
    axes[1,2].plot(np.array([0,128]),np.array([0,0]),color='black',linewidth=3)
    axes[1,2].plot(np.array([0,128]),np.array([12,12]),color='black',linewidth=3)
    axes[1,2].set_yticks([])
    axes[1,2].set_xticks([0,128])
    axes[1,2].set_xlim(-20,148)
    axes[1,2].set_ylim(-48,12)
    axes[1,2].set_xticklabels(xrange,fontname='Calibri',fontweight='light')
    axes[1,2].set_xlabel(title,labelpad=-18,fontname='Calibri',fontweight='light')
    axes[1,2].tick_params(axis='x',
                        pad=-48,
                        bottom=False,
                        top=False,
                        labelbottom=True)
    axes[1,2].spines['top'].set_visible(False)
    axes[1,2].spines['right'].set_visible(False)
    axes[1,2].spines['bottom'].set_visible(False)
    axes[1,2].spines['left'].set_visible(False)


def custom_timeseriesplot(data,variables,axes,colour,labels):
    
    # plot time series
    sns.lineplot(data=data, x=variables['x'], y=variables['y'], hue=variables['hue'],
                 ci = 70, palette=colour, ax = axes, linewidth=3)
    
    axes.legend(labels['legend'],frameon=False)
    
    # add horizontal line
    axes.axhline(y=0, xmin=-1, xmax=3, color=[0,0,0], linewidth = 3)
    axes.axvline(x=0, ymin=-0.4, ymax=1, color=[0,0,0], linewidth = 3)
    
    # aesthetics
    axes.set_title(labels['title'],pad=25,fontname='Calibri',fontsize=28,fontweight='light')
    axes.set_ylabel(labels['ylabel'],fontname='Calibri',fontsize=24,labelpad=10,fontweight='light')   # add Y axis label
    axes.set_xlabel(labels['xlabel'],fontname='Calibri',fontsize=24,labelpad=10,fontweight='light')   # add Y axis label
    axes.set_ylim(-0.3,0.5)                  # set Y axis range to 0 -> 1
    axes.set_xlim(-0.5,2)                  # set Y axis range to 0 -> 1  
    axes.set_yticks([-0.2,0,0.2,0.4])
    axes.set_yticklabels([-0.2,0,0.2,0.4],fontname='Calibri',fontweight='light')
    axes.set_xticklabels([-0.5,0,0.5,1,1.5,2],fontname='Calibri',fontweight='light')
    axes.tick_params(axis='both',          # change X tick parameters
                       pad=25)
    # change axes
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)


# %% -- define key parameters ----------------------------------------------- #
# get current directory
wdir = 'E:/bjg335/projects/reinstatement_fidelity/'

# set context
sns.set_context("talk")

# set plot defaults
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.color'] = [0,0,0]
mpl.rcParams['ytick.color'] = [0,0,0]
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.linewidth'] = 3

# %% -- prepare data -------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "data/fig2_data/group_task-rf_eeg-cluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# read in percept surface images
img_percept = [pyplot.imread(wdir + "/data/fig2_surface/encoding_left_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/encoding_left_medial.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/encoding_right_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/encoding_right_medial.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/encoding_both_ventral.tiff")]

# read in rse surface images
img_retrieval = [pyplot.imread(wdir + "/data/fig2_surface/retrieval_left_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/retrieval_left_medial.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/retrieval_right_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/retrieval_right_medial.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/retrieval_both_ventral.tiff")]

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(6) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'Change in Alpha/Beta Power (z)',
          'xticklabel':['Perception','Retrieval'],
          'yticks':[-0.6,-0.4,-0.2,0,0.2,0.4],
          'yticklabel':['-0.6','-0.4','-0.2','0','0.2','0.4']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'Calibri',labels,[-0.6,0.35],0.5)
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.jpg",bbox_inches='tight',transparent=True,dpi='figure')
  
# ----- Percept Surface Plot ----- # 
# create figure
f,ax = pyplot.subplots(2,3)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(11) # 12inches
f.set_dpi(300)
f.suptitle('',fontname='Calibri',fontsize=28,fontweight='light')

# plot surface
custom_surfaceplot(img_percept,ax,'Blues_r','Pre > Post\nPower (t)',['-7','-2'])

# save image
pyplot.savefig(wdir + "/figures/fig2b.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- RSE Surface Plot ----- # 
# create figure
f,ax = pyplot.subplots(2,3)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(11) # 12inches
f.set_dpi(300)
f.suptitle('',fontname='Calibri',fontsize=28,fontweight='light')

# plot surface
custom_surfaceplot(img_retrieval,ax,'Blues_r','Hits > Misses\nPower (t)',['-4','-2'])
