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
data_fname = wdir + "data/fig2_data/group_task-rf_eeg-cluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_raincloud = data_raincloud.drop(labels = 'encoding', axis = 1)


# -- prep frequency
# load frequency data
datatmp = pandas.read_csv(wdir + "data/fig2_data/group_task-memory_eeg-freqseries.csv",
                                 delimiter=',',
                                 header=None)

# create new structure for frequency data
data_frequency = pandas.DataFrame(data=np.reshape(datatmp.values,[datatmp.size]),columns=['signal'])

# create subject number array
data_frequency = data_frequency.assign(subj=pandas.Series(np.repeat(np.arange(0,21),[150])).values)
    
# create condition array
data_frequency = data_frequency.assign(condition=pandas.Series(np.tile(np.append(np.zeros(75),np.ones(75)),[21])).values)

# create frequency array
data_frequency = data_frequency.assign(frequency=pandas.Series(np.tile(np.linspace(3,40,75),[42])).values)


# -- prep frequency diff
# get frames for hits and misses
data_freqA = data_frequency[data_frequency['condition']==1];
data_freqB = data_frequency[data_frequency['condition']==0];
data_freqA = data_freqA.reset_index()
data_freqB = data_freqB.reset_index()

# create new frame
data_freqdiff = data_freqA
data_freqdiff['signal'] = data_freqB['signal']-data_freqA['signal']


# -- prep time
# load frequency data
datatmp = pandas.read_csv(wdir + "data/fig2_data/group_task-memory_eeg-timeseries.csv",
                                 delimiter=',',
                                 header=None)

# create new structure for frequency data
data_time = pandas.DataFrame(data=np.reshape(datatmp.values,[datatmp.size]),columns=['signal'])

# create subject number array
data_time = data_time.assign(subj=pandas.Series(np.repeat(np.arange(0,21),[122])).values)
    
# create condition array
data_time = data_time.assign(condition=pandas.Series(np.tile(np.append(np.zeros(61),np.ones(61)),[21])).values)

# create frequency array
data_time = data_time.assign(frequency=pandas.Series(np.tile(np.linspace(-1,2,61),[42])).values)


# -- prep frequency diff
# get frames for hits and misses
data_timeA = data_time[data_time['condition']==1];
data_timeB = data_time[data_time['condition']==0];
data_timeA = data_timeA.reset_index()
data_timeB = data_timeB.reset_index()

# create new frame
data_timediff = data_timeA
data_timediff['signal'] = data_timeB['signal']-data_timeA['signal']


# %% ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(6.2/2.54) # 4inches 
f.set_figwidth(3.1/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Power (Rem. > Forgot.; z)',
          'xticklabel':[''],
          'yticks':[-0.5,-0.25,0,0.25],
          'yticklabel':['-0.5','-0.25','0','0.25']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'Calibri',labels,[-0.5,0.25],0.15,[0.018])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.jpg",bbox_inches='tight',transparent=True,dpi='figure')
  
# %% ----- Frequency Individual Series ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(4/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [(0.7,0.7,0.7),colour[5]]

# define labels and variables
labels = {'legend':['Forgot.','Rem.'],
          'ylabel':'Power (z)',
          'xlabel':'Frequency (Hz.)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_frequency,variables,ax,colour,labels,[3,40],[-0.25,0.25],[5,10,15,20,25,30,35,40],True,False)

# save image
pyplot.savefig(wdir + "/figures/fig2b.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- Frequency Difference Series ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(4/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[5],colour[5]]

# define labels and variables
labels = {'legend':[''],
          'ylabel':'Power\n(Rem. > Forg.; z)',
          'xlabel':'Frequency (Hz.)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_freqdiff,variables,ax,colour,labels,[3,40],[-0.25,0.1],[5,10,15,20,25,30,35,40],False,True,True)

# save image
pyplot.savefig(wdir + "/figures/fig2c.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- Time Series ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(4/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [(0.7,0.7,0.7),colour[5]]

# define labels and variables
labels = {'legend':['Forgot.','Rem.'],
          'ylabel':'Power (z)',
          'xlabel':'Time (s)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_time,variables,ax,colour,labels,[-0.5,2],[-0.25,0.25],[-0.5,0,0.5,1,1.5,2],True,True)

# save image
pyplot.savefig(wdir + "/figures/fig2d.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- Time Series Individual ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(4/2.54) # 12inches
f.set_dpi(300)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[5],colour[5]]

# define labels and variables
labels = {'legend':[''],
          'ylabel':'Power\n(Rem. > Forg.; z)',
          'xlabel':'Time (s)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_timediff,variables,ax,colour,labels,[-0.5,2],[-0.3,0.2],[-0.5,0,0.5,1,1.5,2],False,True,True)

# save image
pyplot.savefig(wdir + "/figures/fig2e.jpg",bbox_inches='tight',transparent=True,dpi='figure')