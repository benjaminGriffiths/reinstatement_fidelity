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
import scipy
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
    for i in np.arange(0,np.size(data_violin,1)):
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
    axes=sns.boxplot(data = data, palette = colour, width = 0.1, ax = axes, linewidth = 1, fliersize = 1, whis = 1.5)
    
    # plot significance
    sig_offset = ylim[1]-(ylim[1]*0.05)
    for i in np.arange(0,np.size(pvalue)):
        if pvalue[i] < 0.001:
            pyplot.scatter(np.array([-0.05,0,0.05])+i,[sig_offset,sig_offset,sig_offset],s=3,c=[0.2,0.2,0.2],marker='*',edgecolors=None)
        elif pvalue[i] < 0.01:
            pyplot.scatter(np.array([-0.025,0.025])+i,[sig_offset,sig_offset],s=3,c=[0.2,0.2,0.2],marker='*',edgecolors=None)
        elif pvalue[i] < 0.05:
            pyplot.scatter(np.array([0])+i,[sig_offset],s=3,c=[0.2,0.2,0.2],marker='*',edgecolors=None)
        elif pvalue[i] < 0.1:
            pyplot.scatter(np.array([0])+i,[sig_offset],s=3,c=[0.2,0.2,0.2],marker='x',edgecolors=None)
    
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
              
def custom_timeseriesplot(data,variables,axes,colour,labels,xlim,ylim,xtick,ytick,vertical,horizontal):
    
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
    ax.set_ylabel(labels['ylabel'],fontname='Calibri',fontsize=6,labelpad=3,fontweight='light')   # add Y axis label
    ax.set_xlabel(labels['xlabel'],fontname='Calibri',fontsize=6,labelpad=3,fontweight='light')   # add Y axis label
    ax.set_ylim(ylim)                  # set Y axis range to 0 -> 1
    ax.set_xlim(xlim)                  # set Y axis range to 0 -> 1  
    ax.set_yticks(ytick)
    ax.set_xticks(xtick)
    ax.set_yticklabels(ytick,fontname='Calibri',fontweight='light',fontsize=5)
    ax.set_xticklabels(xtick,fontname='Calibri',fontweight='light',fontsize=5)
    axes.tick_params(axis='y',          # change X tick parameters
                       pad=2,
                       width=1,
                       length=2.5)
    axes.tick_params(axis='x',          # change X tick parameters
                   which='both',          # affect both major and minor ticks
                   bottom=True,          # turn off bottom ticks
                   labelbottom=True,  # keep bottom labels
                   pad=2,
                   width=1,
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
data_A = wdir + "data/fig3_data/group_task-rf_eegfmri-parametricwavecluster.csv"
data_B = wdir + "data/fig3_data/group_task-rf_eegfmri-parametricaudiowavecluster.csv"

# load raincloud data
data_A = pandas.read_csv(data_A,delimiter=",")
data_B = pandas.read_csv(data_B,delimiter=",")

# rename columns
data_A.columns = ['ve_bold','vr_bold','ve_conf','vr_conf','ve_postpow','vr_postpow','ve_prepow','vr_prepow']
data_B.columns = ['ae_bold','ae_conf','ae_postpow','ae_prepow']

# create dataframes
data_post = pandas.concat([data_A[['ve_postpow','vr_postpow']].copy(),data_B[['ae_postpow']].copy()],axis=1)
data_other = pandas.concat([data_A[['ve_prepow','vr_prepow','ve_bold','vr_bold','ve_conf','vr_conf']].copy(),data_B[['ae_prepow','ae_bold','ae_conf']].copy()],axis=1)

# define colour scheme
colour = sns.color_palette("Greens",n_colors=8)
colour = [colour[5],colour[3],colour[1]]

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.8/2.54) # 4inches 
f.set_figwidth(5.5/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Pow. (t)',
          'xticklabel':['Vis. Percept.','Vis. Retrieval','Aud. Percept'],
          'yticks':[-3,-1.5,0,1.5,3],
          'yticklabel':['-3','','0','','3']}

# plot raincloud
custom_rainplot(data_post,colour,ax,'calibri',labels,[-3,3],0.125,[0.055,0,0.006])
   
# save image
pyplot.savefig(wdir + "/figures/fig3a.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/2.54) # 4inches 
f.set_figwidth(9/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Regressor Weight (t)',
          'xticklabel':['Vis/Per','Vis/Ret','Aud/Per','Vis/Per','Vis/Ret','Aud/Per','Vis/Per','Vis/Ret','Aud/Per'],
          'yticks':[-5,-2.5,0,2.5,5],
          'yticklabel':['-5','','0','','5']}

# plot raincloud
custom_rainplot(data_other,colour,ax,'calibri',labels,[-5,5],0.1,[1,1,1,1,1,1,1,1,0.04])
   
# save image
pyplot.savefig(wdir + "/figures/fig3b.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- PLOT FREQSERIES
# define raincloud data filename
data = wdir + "data/fig3_data/group_task-rf_eegfmri-freqgeneralisation.csv"
data = pandas.read_csv(data,delimiter=",",header=None)
data.columns = ['pow','subj','freq','cond']
data_ve = data.loc[data['cond']==1]
data_vr = data.loc[data['cond']==2]
data_ae = data.loc[data['cond']==3]
data['cond'] = np.ones(np.shape(data['cond']))

# define colour scheme
colour = sns.color_palette("Greens",n_colors=8)
colour = [colour[5]]

# define labels and varibales
labels = {'xlabel':'Frequency (Hz.)','ylabel':'Post-Stim. Power (t)'}
variables = {'x':'freq','y':'pow','condition':'cond'}

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.8/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)
custom_timeseriesplot(data,variables,ax,colour,labels,[3,40],[-.3,.2],[10,20,30,40],[-.3,0,.2],True,True)
pyplot.savefig(wdir + "/figures/fig3c.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.0/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)
custom_timeseriesplot(data_ve,variables,ax,colour,labels,[3,40],[-.4,.4],[10,20,30,40],[-.4,0,.4],True,True)
pyplot.savefig(wdir + "/figures/supp5a.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
labels = {'xlabel':'Frequency (Hz)','ylabel':''}
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.0/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)
custom_timeseriesplot(data_ae,variables,ax,colour,labels,[3,40],[-.4,.4],[10,20,30,40],[],True,True)
pyplot.savefig(wdir + "/figures/supp5b.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.0/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)
custom_timeseriesplot(data_vr,variables,ax,colour,labels,[3,40],[-.4,.4],[10,20,30,40],[],True,True)
pyplot.savefig(wdir + "/figures/supp5c.jpg",bbox_inches='tight',transparent=True,dpi='figure')


# %% ----- PLOT TimeSERIES
# define raincloud data filename
data = wdir + "data/fig3_data/group_task-rf_eegfmri-timegeneralisation.csv"
data = pandas.read_csv(data,delimiter=",",header=None)
data.columns = ['pow','subj','time','cond']
data_ve = data.loc[data['cond']==1]
data_vr = data.loc[data['cond']==2]
data_ae = data.loc[data['cond']==3]
data['cond'] = np.ones(np.shape(data['cond']))

# define colour scheme
colour = sns.color_palette("Greens",n_colors=8)
colour = [colour[5]]

# define labels and varibales
labels = {'xlabel':'Time (s)','ylabel':'Alpha/Beta Power (t)'}
variables = {'x':'time','y':'pow','condition':'cond'}

# create figure combined
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.8/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)        
custom_timeseriesplot(data,variables,ax,colour,labels,[-1,1.5],[-.4,.3],[-1,-0.5,0,.5,1,1.5],[-.4,0,.3],True,True)
pyplot.savefig(wdir + "/figures/fig3d.jpg",bbox_inches='tight',transparent=True,dpi='figure')

f,ax = pyplot.subplots(1,1)
f.set_figheight(2.0/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)        
custom_timeseriesplot(data_ve,variables,ax,colour,labels,[-1,1.5],[-.5,.5],[-1,0,1.5],[-.5,0,.5],True,True)
pyplot.savefig(wdir + "/figures/supp5d.jpg",bbox_inches='tight',transparent=True,dpi='figure')

labels = {'xlabel':'Time (s)','ylabel':''}
# create figure seperate
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.0/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)        
custom_timeseriesplot(data_ae,variables,ax,colour,labels,[-1,1.5],[-.5,.5],[-1,0,1.5],[],True,True)
pyplot.savefig(wdir + "/figures/supp5e.jpg",bbox_inches='tight',transparent=True,dpi='figure')

f,ax = pyplot.subplots(1,1)
f.set_figheight(2.0/2.54) # 4inches 
f.set_figwidth(2.4/2.54) # 12inches
f.set_dpi(1000)        
custom_timeseriesplot(data_vr,variables,ax,colour,labels,[-1,1.5],[-.5,.5],[-1,0,1.5],[],True,True)
pyplot.savefig(wdir + "/figures/supp5f.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- PLOT IRASA
# define raincloud data filename
data = wdir + "data/fig3_data/group_task-rf_eegfmri-irasacluster.csv"

# load raincloud data
data = pandas.read_csv(data,delimiter=",")

# create dataframes
data_post = data[['visenc_postpow','visret_postpow','audenc_postpow']].copy()
data_slp  = data[['visenc_postslp','visret_postslp','audenc_postslp']].copy()
data_int  = data[['visenc_postint','visret_postint','audenc_postint']].copy()

# define colour scheme
colour = sns.color_palette("Greens",n_colors=8)
colour = [colour[5],colour[3],colour[1]]

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.8/2.54) # 4inches 
f.set_figwidth(5.5/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Pow. (t)',
          'xticklabel':['Vis. Percept.','Vis. Retrieval','Aud. Percept'],
          'yticks':[-4,-2,0,2,4],
          'yticklabel':['-4','','0','','4']}

# plot raincloud
custom_rainplot(data_post,colour,ax,'calibri',labels,[-4,4],0.125,[])
   
# save image
pyplot.savefig(wdir + "/figures/supp7_d.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.8/2.54) # 4inches 
f.set_figwidth(5.5/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Slope (t)',
          'xticklabel':['Vis. Percept.','Vis. Retrieval','Aud. Percept'],
          'yticks':[-4,-2,0,2,4],
          'yticklabel':['-4','','0','','4']}

# plot raincloud
custom_rainplot(data_slp,colour,ax,'calibri',labels,[-4,4],0.125,[])
   
# save image
pyplot.savefig(wdir + "/figures/supp7_e.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.8/2.54) # 4inches 
f.set_figwidth(5.5/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Intercept (t)',
          'xticklabel':['Vis. Percept.','Vis. Retrieval','Aud. Percept'],
          'yticks':[-4,-2,0,2,4],
          'yticklabel':['-4','','0','','4']}

# plot raincloud
custom_rainplot(data_int,colour,ax,'calibri',labels,[-4,4],0.125,[])
   
# save image
pyplot.savefig(wdir + "/figures/supp7_f.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% ---- PLOT POWER DISTRIBUTION
# define raincloud data filename
data = wdir + "data/fig2_data/group_task-bimodal_eeg-power.csv"
data = pandas.read_csv(data,delimiter=",",header=None)
data.columns = ['count','bin','condition','subj']

# define colour scheme
colour = sns.color_palette("Greens",n_colors=8)
colour = [colour[6],colour[4],colour[2]]

# cycle through each subject
for s in np.arange(1,22):
    
    # create figure
    f,ax = pyplot.subplots(1,1)
    f.set_figheight(2/2.54) # 4inches 
    f.set_figwidth(2/2.54) # 12inches
    f.set_dpi(1000)
    
    dats = data[data.subj==s]
    sns.lineplot(data=dats,x='bin',y='count',hue='condition',axes=ax,linewidth=1,palette=colour,)
    ax.axvline(x=0, linewidth = 1, color = [0,0,0], linestyle='--')
    ax.set_xlabel('Power (z)',fontname='Calibri',fontweight='light',fontsize=5)
    ax.set_xticks([-4,-2,0,2,4])
    ax.set_xticklabels(['-4','','0','','4'],fontname='Calibri',fontweight='light',fontsize=5)
    ax.set_xlim([-4,4])
    ax.set_yticks([0,5,10,15,20,25,30])
    ax.set_yticklabels(['0','','10','','20','','30'],fontname='Calibri',fontweight='light',fontsize=5)
    ax.set_ylim([-0,30])
    ax.set_ylabel('Count',fontname='Calibri',fontweight='light',fontsize=5)
    ax.tick_params(axis='both',pad=3,length=2.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_legend().remove()
        
    # save image
    pyplot.savefig(wdir + "/figures/fig3e_s" + str(s) + ".jpg",bbox_inches='tight',transparent=True,dpi='figure')

# %% ---- PLOT BIMODAL REGRESSION
# define raincloud data filename
data_A = wdir + "data/fig3_data/group_task-rf_eegfmri-medianwavecluster.csv"
data_B = wdir + "data/fig3_data/group_task-rf_eegfmri-medianaudiowavecluster.csv"

# load raincloud data
data_A = pandas.read_csv(data_A,delimiter=",")
data_B = pandas.read_csv(data_B,delimiter=",")

# rename columns
data_A.columns = ['ve_bold','vr_bold','ve_conf','vr_conf','ve_postpow','vr_postpow','ve_prepow','vr_prepow']
data_B.columns = ['ae_bold','ae_conf','ae_postpow','ae_prepow']

# create dataframes
data_post = pandas.concat([data_A[['ve_postpow','vr_postpow']].copy(),data_B[['ae_postpow']].copy()],axis=1)
data_pre = pandas.concat([data_A[['ve_prepow','vr_prepow']].copy(),data_B[['ae_prepow']].copy()],axis=1)

# define colour scheme
colour = sns.color_palette("Greens",n_colors=8)
colour = [colour[5],colour[3],colour[1]]

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/2.54) # 4inches 
f.set_figwidth(4.8/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Post-Stim. Power (t)',
          'xticklabel':['Vis. Percept.','Vis. Retrieval','Aud. Percept'],
          'yticks':[-4,-2,0,2,4],
          'yticklabel':['-4','','0','','4']}

# plot raincloud
custom_rainplot(data_post,colour,ax,'calibri',labels,[-4,4],0.125,[0.2,0.04,0.17])
   
# save image
pyplot.savefig(wdir + "/figures/fig3f.jpg",bbox_inches='tight',transparent=True,dpi='figure')

# ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/2.54) # 4inches 
f.set_figwidth(4.8/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Pre-Stim. Power (t)',
          'xticklabel':['Vis. Percept.','Vis. Retrieval','Aud. Percept'],
          'yticks':[-4,-2,0,2,4],
          'yticklabel':['-4','','0','','4']}

# plot raincloud
custom_rainplot(data_pre,colour,ax,'calibri',labels,[-4,4],0.125,[0.8,0.7,0.04])
   
# save image
pyplot.savefig(wdir + "/figures/fig3g.jpg",bbox_inches='tight',transparent=True,dpi='figure')
