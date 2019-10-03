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
    sig_offset = ylim[1]-(ylim[1]*0.05)
    for i in np.arange(0,np.size(pvalue)):
        if pvalue[i] < 0.001:
            pyplot.scatter(np.array([-0.05,0,0.05])+i,[sig_offset,sig_offset,sig_offset],s=1,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.01:
            pyplot.scatter(np.array([-0.025,0.025])+i,[sig_offset,sig_offset],s=1,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.05:
            pyplot.scatter(np.array([0])+i,[sig_offset],s=2,c='black',marker='*',edgecolors=None)
    
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
    
               
def custom_timeseriesplot(data,variables,axes,colour,labels,xlim,ylim,xtick,legend,vertical,horizontal=False):
    
    sns.lineplot(x=variables['x'],
             y=variables['y'],
             data=data,
             hue=variables['hue'],
             hue_order=[1,0],
             ci='sem',
             palette=colour,
             ax = axes,
             linewidth=1)

    if legend:
        ax.legend(labels['legend'],frameon=False,fontsize=5,bbox_to_anchor=(0., 1.02, 1., .102), 
                  loc=3, ncol=2, mode="expand", borderaxespad=0.)
    else:
        ax.get_legend().remove()
        
        
    # add horizontal line
    if vertical:
        ax.axvline(x=0, linewidth = 1, color = [0,0,0], linestyle='--')
    if horizontal:
        ax.axhline(y=0, linewidth = 1, color = [0,0,0], linestyle='-')
    
    # aesthetics
    ax.set_ylabel(labels['ylabel'],fontname='Calibri',fontsize=6,labelpad=0,fontweight='light')   # add Y axis label
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
data_fname = wdir + "data/fig2_data/group_task-rf_eeg-wavecluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_vis = data_raincloud.drop(labels = ['aud_enc_pow','aud_ret_pow','aud_rse_pow'], axis = 1)
data_aud = data_raincloud.drop(labels = ['vis_enc_pow','vis_ret_pow','vis_rse_pow'], axis = 1)

# %% ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.7/2.54) # 4inches 
f.set_figwidth(8.3/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Power (z)',
          'xticklabel':['Percept. ERD','Ret. ERD','Ret. Success'],
          'yticks':[-0.6,-0.3,0,0.3],
          'yticklabel':['-0.6','-0.3','0','0.3']}

# plot raincloud
custom_rainplot(data_vis,colour,ax,'Calibri',labels,[-0.6,0.3],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.7/2.54) # 4inches 
f.set_figwidth(8.3/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Power (z)',
          'xticklabel':['Percept. ERD','Ret. ERD','Ret. Success'],
          'yticks':[-0.6,-0.3,0,0.3],
          'yticklabel':['-0.6','-0.3','0','0.3']}

# plot raincloud
custom_rainplot(data_aud,colour,ax,'Calibri',labels,[-0.6,0.3],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')
  

# %% -- IRASA DATA -------------- #

# %% -- prepare data -------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "data/fig2_data/group_task-rf_eeg-irasacluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_pow = data_raincloud.drop(labels = ['aud_enc_pow','aud_ret_pow','aud_rse_pow',
                                         'vis_enc_slp','vis_ret_slp','vis_rse_slp',
                                         'aud_enc_slp','aud_ret_slp','aud_rse_slp',
                                         'vis_enc_int','vis_ret_int','vis_rse_int',
                                         'aud_enc_int','aud_ret_int','aud_rse_int'], axis = 1)
    
# restrict to retrieval data
data_slp = data_raincloud.drop(labels = ['aud_enc_pow','aud_ret_pow','aud_rse_pow',
                                         'vis_enc_pow','vis_ret_pow','vis_rse_pow',
                                         'aud_enc_slp','aud_ret_slp','aud_rse_slp',
                                         'vis_enc_int','vis_ret_int','vis_rse_int',
                                         'aud_enc_int','aud_ret_int','aud_rse_int'], axis = 1)
    
# restrict to retrieval data
data_int = data_raincloud.drop(labels = ['aud_enc_pow','aud_ret_pow','aud_rse_pow',
                                         'vis_enc_pow','vis_ret_pow','vis_rse_pow',
                                         'aud_enc_slp','aud_ret_slp','aud_rse_slp',
                                         'vis_enc_slp','vis_ret_slp','vis_rse_slp',
                                         'aud_enc_int','aud_ret_int','aud_rse_int'], axis = 1)

# --- plot power
# rescale power
pow_sign = np.sign(data_pow)
data_pow = np.log10(abs(data_pow / (10**7))+1) * pow_sign
    
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.7/2.54) # 4inches 
f.set_figwidth(8.3/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Log-Norm. Power [a.u.]',
          'xticklabel':['Percept. ERD','Ret. ERD','Ret. Success'],
          'yticks':[-1,-0.5,0,0.5],
          'yticklabel':['-1','-0.5','0','0.5']}

# plot raincloud
custom_rainplot(data_pow,colour,ax,'Calibri',labels,[-1,0.5],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')


# --- plot slope
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.7/2.54) # 4inches 
f.set_figwidth(8.3/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Aperiodic Slope',
          'xticklabel':['Percept. ERD','Ret. ERD','Ret. Success'],
          'yticks':[-40,-20,0,20],
          'yticklabel':['-40','-20','0','20']}

# plot raincloud
custom_rainplot(data_slp,colour,ax,'Calibri',labels,[-40,20],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')


# --- plot int
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.7/2.54) # 4inches 
f.set_figwidth(8.3/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Aperiodic Intercept (t)',
          'xticklabel':['Percept. ERD','Ret. ERD','Ret. Success'],
          'yticks':[-100,0,100,200],
          'yticklabel':['-100','0','100','200']}

# plot raincloud
custom_rainplot(data_int,colour,ax,'Calibri',labels,[-100,200],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')


