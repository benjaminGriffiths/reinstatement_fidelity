# ----------- Data Visualisation ----------- #
# Written by Benjamin J. Griffiths
# Created on Friday 26th October, 2018
# ------------------------------------------ #

# %% -- import modules ------------------------------------------------------ #
import matplotlib as mpl
import numpy as np
import os
import pandas
import ptitprince as pt
import scipy.io
import scipy.stats
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
    
    
# plot rdms
def custom_grandRDMplot(data,colour,axes,fontname,labels,drawX,drawY):
    
    # plot RDM
    im = axes.imshow(X=data,
                 cmap=colour,
                 aspect='equal')
    axes.set_title(labels['title'],pad=12,fontname=fontname,fontsize=20,fontweight='light')
    axes.tick_params(axis='both', pad=0, width = 0)
    
    if drawX: 
        axes.set_xticks([0,1,2,3])
        axes.set_xticklabels(labels['xticks'],fontname=fontname,fontsize=12,rotation=45,fontweight='light')
        axes.set_xlabel(labels['xlabel'],fontname=fontname,fontsize=14,labelpad=10,fontweight='light')
    else:
        axes.set_xticks([])   
        
    if drawY: 
        axes.set_yticks([0,1,2,3])
        axes.set_yticklabels(labels['yticks'],fontname=fontname,fontsize=12,rotation=45,fontweight='light')
        axes.set_ylabel(labels['ylabel'],fontname=fontname,fontsize=14,labelpad=10,fontweight='light')
    else:
        axes.set_yticks([])   
    
    return im
   
               
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
wdir = os.path.dirname(os.path.realpath(__file__))

# move up one directory
#wdir = wdir[:-13]

# set context
sns.set_context("talk")

# set plot defaults
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.color'] = [0,0,0]
mpl.rcParams['ytick.color'] = [0,0,0]
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.linewidth'] = 3

# %% -- Figure 1. fMRI-RSA -------------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "/data/searchlight_betas.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# define RDM filename
data_fname = wdir + "/data/searchlight_rdm.mat"

# load RDM data
data_rdm  = scipy.io.loadmat(data_fname)

# get average RDM
grand_rdm = np.mean(data_rdm['sRDMs'],axis=0)

# read all surface images
img = [pyplot.imread(wdir + "/data/fig1_surface/left_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/left_medial.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/right_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/right_medial.tiff"),
       pyplot.imread(wdir + "/data/fig1_surface/both_ventral.tiff")]

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(11) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Reds",n_colors=6)
colour = colour[1:4]

# define labels
labels = {'title':'Similarity per Cluster',
          'ylabel':'Similarity Index\n(Same > Different; z)',
          'xticklabel':['Left Fusiform','Right Fusiform','Right Parieto-Occipital\nFissure'],
          'yticks':[-0.05,0,0.05,0.10,0.15,0.20],
          'yticklabel':['-0.05','0','0.05','0.10','0.15','0.20']}

# plot
custom_rainplot(data_raincloud,colour,ax,'Calibri',labels,[-0.05,0.20],0.55)

# save image
pyplot.savefig(wdir + "/figures/fig1a.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(8) # 12inches
f.set_dpi(300)

# create grid
grid = pyplot.GridSpec(2,11,figure=f)
ax    = [0,0,0,0,0]
ax[0] = pyplot.subplot(grid[0,0:4])
ax[1] = pyplot.subplot(grid[1,0:4])
ax[2] = pyplot.subplot(grid[0,5:9])
ax[3] = pyplot.subplot(grid[1,5:9])
ax[4] = pyplot.subplot(grid[:,10:])

# define RDM labels
rdm_label = ['Left Fusiform','Right Fusiform','Right POF','Combined']

# cycle through each RDM
for i in np.arange(0,np.size(grand_rdm,0)):
       
    # restrict RDM to visual encoding-retrieval
    roi_rdm = grand_rdm[i,0:4,8:12]
    
    # rank rdm between 0 and 1
    roi_rdm = np.reshape(scipy.stats.rankdata(roi_rdm)-1,[4,4])/((4*4)-1)
        
    # decide whether to draw X labels
    if i == 1 or i == 3: 
        drawX = True
    else:
        drawX = False
        
    # decide whether to draw Y labels
    if i <= 1: 
        drawY = True
    else:
        drawY = False
        
    # define labels
    labels = {'title'  : rdm_label[i],
              'xticks' : ['Video A','Video B','Video C','Video D'],
              'xlabel' : 'Retrieval',
              'yticks' : ['Video A','Video B','Video C','Video D'],
              'ylabel' : 'Encoding'}
    
    # plot rdm
    im = custom_grandRDMplot(roi_rdm,'OrRd_r',ax[i],'Calibri',labels,drawX,drawY)

# add colourbar
cbar = f.colorbar(im,cax=ax[4],
           ticks=[0,1],
           fraction=0.1,
           pad=0.075)
cbar.set_label('Ranked Distance\n(Cross-Validated Mahalanobis)',fontname='Calibri',fontsize=24,fontweight='light')
cbar.ax.tick_params(width=0)

# save image
pyplot.savefig(wdir + "/figures/fig1b.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# create second figure for surface plots
f,ax = pyplot.subplots(2,3)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(11) # 12inches
f.set_dpi(300)
f.suptitle('Searchlight-based RSA',fontname='Calibri',fontsize=28,fontweight='light')

custom_surfaceplot(img,ax,'OrRd_r','Similarity Index\n(Same > Different; t)',['2','7'])

# save image
pyplot.savefig(wdir + "/figures/fig1c.jpg",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 2. Source EEG RSE -------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "/data/sourceEEG.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# read in surface images
img = [pyplot.imread(wdir + "/data/fig2_surface/left_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/left_medial.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/right_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/right_medial.tiff"),
       pyplot.imread(wdir + "/data/fig2_surface/both_ventral.tiff")]

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(5) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],(0.8,0.8,0.8)]

# define labels
labels = {'title':'Retrieval Success Effect',
          'ylabel':'Alpha/Beta Power (Hits > Misses; z)',
          'xticklabel':[''],
          'yticks':[-0.5,-0.25,0,0.25],
          'yticklabel':['-0.5','-0.25','0','0.25']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'Calibri',labels,[-0.5,0.25],0.65)
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.jpg",bbox_inches='tight',transparent=True,dpi='figure')
  
# ----- Surface Plot ----- # 
# create figure
f,ax = pyplot.subplots(2,3)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(11) # 12inches
f.set_dpi(300)
f.suptitle('Alpha/Beta (8-30Hz) Retrieval Success Effect',fontname='Calibri',fontsize=28,fontweight='light')

# plot surface
custom_surfaceplot(img,ax,'Blues_r','Hits > Misses\nPower (z)',['-3.5','-0.5'])

# save image
pyplot.savefig(wdir + "/figures/fig2c.jpg",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 3. Combined Analysis ----------------------------------------- #
# define raincloud data filename
rc_fname = wdir + "/data/RSA-EEG_correlation.csv"
  
# load time series data
data_raincloud = pandas.read_csv(rc_fname,
                       delimiter=",")

# read all surface images
img = [pyplot.imread(wdir + "/data/fig3_surface/left_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig3_surface/left_medial.tiff"),
       pyplot.imread(wdir + "/data/fig3_surface/right_lateral.tiff"),
       pyplot.imread(wdir + "/data/fig3_surface/right_medial.tiff"),
       pyplot.imread(wdir + "/data/fig3_surface/both_ventral.tiff")]

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(8) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Greens",n_colors=7)
colour = [colour[2],(0.8,0.8,0.8)]

# define labels
labels = {'title':'ρ(Alpha/Beta Power, Similarity Index)',
          'ylabel':"Pearson's Correlation (r)",
          'xticklabel':['Hits','Misses'],
          'yticks':[-0.6,-0.4,-0.2,0,0.2,0.4,0.6],
          'yticklabel':['-0.6','-0.4','-0.2','0','0.2','0.4','0.6']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'Calibri',labels,[-0.6,0.6],0.5)

# save image
pyplot.savefig(wdir + "/figures/fig3a.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- Surface Plot ----- # 
# prepare figure
f,ax = pyplot.subplots(2,3)
f.set_figheight(20/3) # 4inches 
f.set_figwidth(11) # 12inches
f.set_dpi(300)
f.suptitle('ρ(Alpha/Beta Power, BOLD)',fontname='Calibri',fontsize=28,fontweight='light')

# plot surface
custom_surfaceplot(img,ax,'Greens_r','T-Statistic',[-5,0])

# save image
pyplot.savefig(wdir + "/figures/fig3c.jpg",bbox_inches='tight',transparent=True,dpi='figure')