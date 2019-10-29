# ----------- Data Visualisation ----------- #
# Written by Benjamin J. Griffiths
# Created on Friday 26th October, 2018
# ------------------------------------------ #

# %% -- import modules ------------------------------------------------------ #
import matplotlib as mpl
import numpy as np
import scipy
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
    
    # drop outliers from violin
    data_violin = data.copy()
    for i in np.arange(0,3):
        qs  = data_violin.iloc[:,i].quantile([.25,.75])
        qs = qs.tolist()
        iqr = scipy.stats.iqr(data_violin.iloc[:,i])
        idx = (data_violin.iloc[:,i] > (qs[1] + 1.5*iqr)) | (data_violin.iloc[:,i] < (qs[0] - 1.5*iqr))
        idx = [x for x, val in enumerate(idx) if val] 
        if idx:
            data_violin.iloc[idx,i] = float('nan')
    
    # plot violin
    axes=pt.half_violinplot(data = data_violin,bw = "scott",inner = None, scale = "count",
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
    
               
def custom_timeseriesplot(data,variables,axes,colour,labels,xlim,ylim,xtick,vertical,horizontal=False):
    
    sns.lineplot(x=variables['x'],
             y=variables['y'],
             data=data,
             hue=variables['hue'],
             ci='sem',
             palette=colour,
             ax = axes,
             linewidth=1)
        
    # add horizontal line
    if vertical:
        ax.axvline(x=0, linewidth = 1, color = [0,0,0], linestyle='--')
    if horizontal:
        ax.axhline(y=0, linewidth = 1, color = [0,0,0], linestyle='-')
    
    # aesthetics
    ax.set_ylabel(labels['ylabel'],fontname='Calibri',fontsize=6,labelpad=3,fontweight='light')   # add Y axis label
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
    ax.get_legend().remove()
    
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

# %% ----- WAVELET RAINCLOUD ----- # 
# define raincloud data filename
data_fname = wdir + "data/fig2_data/group_task-rf_eeg-wavecluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# restrict to retrieval data
data_raincloud = data_raincloud.drop(labels = ['vis_ret_pow','aud_ret_pow','aud_rse_pow'], axis = 1)

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.4/2.54) # 4inches 
f.set_figwidth(5.4/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Power (z)',
          'xticklabel':['Vis. ERD','Aud. ERD','Vis. RSE'],
          'yticks':[-0.8,-0.4,0,0.4],
          'yticklabel':['-0.8','-0.4','0','0.4']}

# -- PLOT
custom_rainplot(data_raincloud,colour,ax,'Calibri',labels,[-0.8,0.4],0.125,[1])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- plot irasa spectra
# define raincloud data filename
data = wdir + "data/fig2_data/group_task-irasa_eeg_osci.csv"
data = pandas.read_csv(data,delimiter=",",header=None)
data.columns = ['pow','subj','task','epoch']
data['freq'] = np.tile(np.linspace(5,25,26),252)

# rescale power
data['pow'] = data['pow'] / (10**7)

# get task sets
dat    = [[],[],[]]
dat[0] = data.loc[data['task']==1]
dat[1] = data.loc[data['task']==4]
dat[2] = data.loc[data['task']==3]

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[5],(0.5,0.5,0.5)]

# cycle through tasks
for t in np.arange(0,3):
  
    # define labels and varibales
    variables = {'x':'freq','y':'pow','hue':'epoch'}
    if t == 0:       
        labels = {'xlabel':'Freq. (Hz)','ylabel':'Periodic Power (a.u.)'}
    else:
        labels = {'xlabel':'Freq. (Hz)','ylabel':''}
    
    # get maxabs range
    datgrp = dat[t].groupby(['subj','epoch']).mean()
    maxabs = np.round(max(np.abs(datgrp['pow'])) * 0.15) * 10
    
    # create figure combined
    f,ax = pyplot.subplots(1,1)
    f.set_figheight(1.8/2.54) # 4inches 
    f.set_figwidth(2.4/2.54) # 12inches
    f.set_dpi(1000)        
    custom_timeseriesplot(dat[t],variables,ax,colour,labels,[5,25],[-30,30],[5,10,15,20,25],False,True)
    
    # save image
    pyplot.savefig(wdir + "/figures/supp7_fi" + str(t) + ".tif",bbox_inches='tight',transparent=True,dpi='figure')
    
    # get difference
    datdiff = dat[t][dat[t]['epoch']==1].copy()
    datdiff['pow'] = dat[t]['pow'][dat[t]['epoch']==1].as_matrix() - dat[t]['pow'][dat[t]['epoch']==2].as_matrix()
    
    # get maxabs range
    datgrp = datdiff.groupby(['subj']).mean()
    maxabs = np.round(max(np.abs(datgrp['pow'])) * 2)
    
    # create figure combined
    f,ax = pyplot.subplots(1,1)
    f.set_figheight(1.8/2.54) # 4inches 
    f.set_figwidth(2.4/2.54) # 12inches
    f.set_dpi(1000)        
    custom_timeseriesplot(datdiff,variables,ax,[colour[0]],labels,[5,25],[-12,12],[5,10,15,20,25],False,True)
    
    # save image
    pyplot.savefig(wdir + "/figures/supp7_fj" + str(t) + ".tif",bbox_inches='tight',transparent=True,dpi='figure')
    

# %% ----- plot irasa spectra
# define raincloud data filename
data = wdir + "data/fig2_data/group_task-irasa_eeg_frac.csv"
data = pandas.read_csv(data,delimiter=",",header=None)
data.columns = ['pow','subj','task','epoch']
data['freq'] = np.tile(np.linspace(5,25,26),252)

# rescale power
data['pow'] = data['pow'] / (10**7)

# get task sets
dat    = [[],[],[]]
dat[0] = data.loc[data['task']==1]
dat[1] = data.loc[data['task']==4]
dat[2] = data.loc[data['task']==3]

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[5],(0.5,0.5,0.5)]

# define labels and varibales
labels = {'xlabel':'Time (s)','ylabel':'Alpha/Beta Power (t)'}
variables = {'x':'freq','y':'pow','hue':'epoch'}

# cycle through tasks
for t in np.arange(0,3):

    # define labels and varibales
    variables = {'x':'freq','y':'pow','hue':'epoch'}
    if t == 0:       
        labels = {'xlabel':'Freq. (Hz)','ylabel':'Aperiodic Power (a.u.)'}
    else:
        labels = {'xlabel':'Freq. (Hz)','ylabel':''}
    
    # get maxabs range
    datgrp = dat[t].groupby(['subj','epoch']).mean()
    maxabs = np.round(max(np.abs(datgrp['pow'])) * 0.15) * 10
    
    # create figure combined
    f,ax = pyplot.subplots(1,1)
    f.set_figheight(1.8/2.54) # 4inches 
    f.set_figwidth(2.4/2.54) # 12inches
    f.set_dpi(1000)        
    custom_timeseriesplot(dat[t],variables,ax,colour,labels,[5,25],[0,120],[5,10,15,20,25],False,True)
    
    # save image
    pyplot.savefig(wdir + "/figures/supp7_fg" + str(t) + ".tif",bbox_inches='tight',transparent=True,dpi='figure')
    
    # get difference
    datdiff = dat[t][dat[t]['epoch']==1].copy()
    datdiff['pow'] = dat[t]['pow'][dat[t]['epoch']==1].as_matrix() - dat[t]['pow'][dat[t]['epoch']==2].as_matrix()
    
    # get maxabs range
    datgrp = datdiff.groupby(['subj']).mean()
    maxabs = np.round(max(np.abs(datgrp['pow'])) * 2)
    
    # create figure combined
    f,ax = pyplot.subplots(1,1)
    f.set_figheight(1.8/2.54) # 4inches 
    f.set_figwidth(2.4/2.54) # 12inches
    f.set_dpi(1000)        
    custom_timeseriesplot(datdiff,variables,ax,[colour[0]],labels,[5,25],[-8,4],[5,10,15,20,25],False,True)
   
    # save image
    pyplot.savefig(wdir + "/figures/supp7_fh" + str(t) + ".tif",bbox_inches='tight',transparent=True,dpi='figure')
    
# %% ----- TIME-SERIES ----- # 
# define filenames
filename = [wdir+"data/fig2_data/group_task-time_eeg-visenc.csv",
            wdir+"data/fig2_data/group_task-time_eeg-audenc.csv",
            wdir+"data/fig2_data/group_task-time_eeg-visrse.csv"]

# intialise
df = pandas.DataFrame(columns=['pow','time','subj','task'])

# cycle through each file
for f in np.arange(0,2):
    
    # load time-series data
    dat = np.genfromtxt(filename[f],delimiter=",")
    
    # define colour scheme
    colour = sns.color_palette("Blues",n_colors=7)
    colour = [colour[5]]
    
    # create panda object
    pdat = pandas.DataFrame(dat)
    pdat.columns = ['pow','time','subj']
    pdat['task'] = np.tile(f+1,np.size(pdat,0))
    
    # add to group
    df = df.append(pdat)
   
# get condition difference
df['condition'] = np.ones(np.size(df,0))

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.7/2.54) # 4inches 
f.set_figwidth(2/2.54) # 12inches
f.set_dpi(1000)

# define variables
vari = {'x':'time','y':'pow','hue':'condition'}

# define labels
labels = {'title':'','ylabel':'Alpha/Beta Power (z)','xlabel':'Freq. (Hz)',
          'yticks':[-0.25,0,0.25],'yticklabel':['-0.2','0','0.2']}

# plot
custom_timeseriesplot(df,vari,ax,colour,labels,xlim=[-1,3],ylim=[-0.25,0.25],xtick=[-1,0,1,2,3],vertical=True,horizontal=True)
 
# save image
pyplot.savefig(wdir + "/figures/fig2b.tif",bbox_inches='tight',transparent=True,dpi='figure')

# ---generate per condition plot ---
# define filenames
filename = [wdir+"data/fig2_data/group_task-time_eeg-visenc.csv",
            wdir+"data/fig2_data/group_task-time_eeg-audenc.csv",
            wdir+"data/fig2_data/group_task-time_eeg-visrse.csv"]

# cycle through each file
for fi in filename:
    print(fi)

    # load time-series data
    dat = np.genfromtxt(fi,delimiter=",")
    
    # define colour scheme
    colour = sns.color_palette("Blues",n_colors=7)
    
    # create panda object
    pdat = pandas.DataFrame(dat)
    if np.size(pdat,1) == 3:
        pdat.columns = ['pow','time','subj']
        pdat['condition'] = np.ones(np.size(pdat,0))
        colour = [colour[5]]
    else:
        pdat.columns = ['pow','time','subj','condition']
        colour = [colour[5],(0.8,0.8,0.8)]
           
    # create figure
    f,ax = pyplot.subplots(1,1)
    f.set_figheight(1.7/2.54) # 4inches 
    f.set_figwidth(2/2.54) # 12inches
    f.set_dpi(1000)
    
    # define variables
    vari = {'x':'time','y':'pow','hue':'condition'}
    
    # define labels
    labels = {'title':'','ylabel':'Alpha/Beta Power (z)','xlabel':'Time (s)',
              'yticks':[-0.3,0,0.3],'yticklabel':['-0.3','0','0.3']}
    
    # plot
    custom_timeseriesplot(pdat,vari,ax,colour,labels,xlim=[-1,3],ylim=[-0.3,0.3],xtick=[-1,0,1,2,3],vertical=True,horizontal=True)
    
    # save image
    pyplot.savefig(wdir + "/figures/supp_3t_" + fi[77:83] + ".tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% ----- FREQ-SERIES ----- # 
# define filenames
filename = [wdir+"data/fig2_data/group_task-freq_eeg-visenc.csv",
            wdir+"data/fig2_data/group_task-freq_eeg-audenc.csv",
            wdir+"data/fig2_data/group_task-freq_eeg-visrse.csv"]

# intialise
df = pandas.DataFrame(columns=['pow','freq','subj','condition','task'])

# cycle through each file
for f in np.arange(0,2):
    
    # load time-series data
    dat = np.genfromtxt(filename[f],delimiter=",")
    
    # define colour scheme
    colour = sns.color_palette("Blues",n_colors=7)
    colour = [colour[5]]
    
    # create panda object
    pdat = pandas.DataFrame(dat)
    pdat.columns = ['pow','freq','subj','condition']
    pdat['task'] = np.tile(f+1,np.size(pdat,0))
    
    # add to group
    df = df.append(pdat)
   
# get condition difference
pow  = df['pow'].loc[df['condition']==1].as_matrix() - df['pow'].loc[df['condition']==2].as_matrix()
df   = df.loc[df['condition']==1]
df['pow'] = pow

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.7/2.54) # 4inches 
f.set_figwidth(2/2.54) # 12inches
f.set_dpi(1000)

# define variables
vari = {'x':'freq','y':'pow','hue':'condition'}

# define labels
labels = {'title':'','ylabel':'Alpha/Beta Power (z)','xlabel':'Freq. (Hz)',
          'yticks':[-0.4,0,0.2],'yticklabel':['-0.25','0','0.25']}

# plot
custom_timeseriesplot(df,vari,ax,colour,labels,xlim=[3,30],ylim=[-0.4,0.2],xtick=[3,10,20,30],vertical=False,horizontal=True)
 
# save image
pyplot.savefig(wdir + "/figures/fig2c.tif",bbox_inches='tight',transparent=True,dpi='figure')

# ---generate per condition plot ---
# define filenames
filename = [wdir+"data/fig2_data/group_task-freq_eeg-visenc.csv",
            wdir+"data/fig2_data/group_task-freq_eeg-audenc.csv",
            wdir+"data/fig2_data/group_task-freq_eeg-visrse.csv"]

# cycle through each file
for fi in filename:
    print(fi)

    # load time-series data
    dat = np.genfromtxt(fi,delimiter=",")
    
    # define colour scheme
    colour = sns.color_palette("Blues",n_colors=7)
    colour = [colour[5]]
    
    # create panda object
    pdat = pandas.DataFrame(dat)
    pdat.columns = ['pow','freq','subj','condition']
    
    # get condition difference
    pow  = pdat['pow'].loc[pdat['condition']==1].as_matrix() - pdat['pow'].loc[pdat['condition']==2].as_matrix()
    pdat = pdat.loc[pdat['condition']==1]
    pdat['pow'] = pow
    
    # create figure
    f,ax = pyplot.subplots(1,1)
    f.set_figheight(1.7/2.54) # 4inches 
    f.set_figwidth(2/2.54) # 12inches
    f.set_dpi(1000)
    
    # define variables
    vari = {'x':'freq','y':'pow','hue':'condition'}
    
    # define ylimit
    pavg = pdat.groupby('freq').mean()
    yl = np.floor(min(pavg['pow'])*13)/10
    
    # define labels
    labels = {'title':'','ylabel':'Alpha/Beta Power (z)','xlabel':'Freq. (Hz)',
              'yticks':[yl,0,-yl],'yticklabel':[str(yl),'0',str(-yl)]}
    
    # plot
    custom_timeseriesplot(pdat,vari,ax,colour,labels,xlim=[3,30],ylim=[yl,-yl],xtick=[3,10,20,30],vertical=False,horizontal=True)
    
    # save image
    pyplot.savefig(wdir + "/figures/supp_3f_" + fi[77:83] + ".tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% -- IRASA DATA -------------- #

# %% -- prepare data -------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "data/fig2_data/group_task-rf_eeg-irasacluster.csv"

# load raincloud data
data_raincloud = pandas.read_csv(data_fname,
                       delimiter=",")

# reorder columns
data_pow = data_raincloud[['vis_enc_pow','aud_enc_pow','vis_rse_pow']]
data_slp = data_raincloud[['vis_enc_slp','aud_enc_slp','vis_rse_slp']]
data_int = data_raincloud[['vis_enc_int','aud_enc_int','vis_rse_int']]

# --- plot power
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/2.54) # 4inches 
f.set_figwidth(4.8/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Periodic Power (z x 10^7)',
          'xticklabel':['Vis. ERD','Aud. ERD','Ret. Success'],
          'yticks':[-4,-2,0,2],
          'yticklabel':['-4','2','0','2']}

# plot raincloud
data_pow = data_pow/10**7
custom_rainplot(data_pow,colour,ax,'Calibri',labels,[-1,0.5],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/supp_7a.tif",bbox_inches='tight',transparent=True,dpi='figure')


# --- plot slope
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/2.54) # 4inches 
f.set_figwidth(4.8/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Aperiodic Slope (t)',
          'xticklabel':['Vis. ERD','Aud. ERD','Ret. Success'],
          'yticks':np.arange(-30,40,15),
          'yticklabel':['-30','','0','','30']}

# plot raincloud
custom_rainplot(data_slp,colour,ax,'Calibri',labels,[-30,30],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/supp_7b.tif",bbox_inches='tight',transparent=True,dpi='figure')


# --- plot int
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/2.54) # 4inches 
f.set_figwidth(4.8/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[2],colour[3],colour[4]]

# define labels
labels = {'title':'',
          'ylabel':'Aperiodic Intercept (t)',
          'xticklabel':['Vis. ERD','Aud. ERD','Ret. Success'],
          'yticks':np.arange(-100,250,100),
          'yticklabel':['-100','0','100','200']}

# plot raincloud
custom_rainplot(data_int,colour,ax,'Calibri',labels,[-100,200],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/supp_7c.tif",bbox_inches='tight',transparent=True,dpi='figure')


