'''
Author: Thais Silva Abelha
Email: thais.silva.abelha@gmail.com

This work was done during a undergraduated research "Comparing different methods of position
reconstruction considering 1D readout of GEM detectors"
'''

import matplotlib.pyplot as plt
from matplotlib import container
import numpy as np
import pandas as pd

#my libraries
from strips_integrator import Strips_Integrator
from monte_carlo import Monte_Carlo
#%%

'''
linear = weighted average considering the charge as a weight
quadratic = weighted average considering the charge*charge as a weighet
logarithmic = weighted average considering the ln(charge) as a weight
Q = charge
s=standard deviation or sigma σ of the gaussian
u=mean of the gaussian or electron cloud center

prefix "s" stands for standard deviation, which is represented by σ

RMS = Root Mean Square 
'''
def Plot_Cluster(first_evt_cluster,u,s,Labels,seed,threshold):
    fig, ax = plt.subplots()
    plt.axhline(seed,label="Seed", linestyle=(0,(5,10)),color="green",linewidth=1.0)
    plt.axhline(threshold,label="Threshold",color="red",linewidth=1.0)
    ax.plot(first_evt_cluster['Position_strip'],first_evt_cluster['Charge_strip_noise'],"o", markersize=4,label=(Labels["Label_Charges"]))
    ax.plot(first_evt_cluster['Position_Cluster'],first_evt_cluster['Charge_Cluster'],"*",markersize=7,label=Labels["Label_Cluster"])
    fig.suptitle(Labels["Title"])
    ax.set_xlabel(Labels["xlabel"])
    ax.set_ylabel(Labels["ylabel"])
    ax.legend(loc='upper right')
    ax.grid()
    ax.set_ylim(-0.05,0.35)
    textstr = '\n'.join((
    r'$\mu=%.2f$' % (u, )+"mm",
    r'$\sigma=%.2f$' % (s, )+"mm"
    ))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.08, 0.95, textstr, transform=ax.transAxes, fontsize=11,verticalalignment='top', bbox=props)
    plt.show()
    
def Dataframe(strip_width,pitch,slist,snoise,N,Position_strips,electron_cloud_centers,seed,threshold):
    '''For the case when    pitch=0.39mm, the electron cloud center is going to be reconstructed approximately -0.8mm<x<0.8mm
    Also the electron cloud center is going to be reconstructed with steps 10 times smaller than the strip pitch'''
    
    #You can comment "Plot_Cluster()" It is just a Graph of control then it is going to run faster
    #Plot_Cluster() is used just to see if verything is fine with the algorithm
    
    Labels={"Title":"Cluster - Event 1/1000",
            "xlabel":"Position of the strip center (mm)",
            "ylabel":"Charge collected by each strip (mm)",
            "Label_Charges":"Charge collected by strips",
            "Label_Cluster":"Cluster"}
    
    Labels_pt={"Title":"Cluster - Evento 1/1000",
            "xlabel":"Posição do centro da fita (mm)",
            "ylabel":"Carga coletada por cada fita",
            "Label_Charges":"Carga coletada pelas fitas",
            "Label_Cluster":"Cluster"}
    Data=[]
    for s in slist:
    
        for u in electron_cloud_centers:
            Charge_strips=Strips_Integrator(strip_width,pitch,snoise,Position_strips).Charge_Strip(s,u)
            Charge,sCharge,E_linear,sE_linear,E_quadratic,sE_quadratic,E_logarithmic,sE_logarithmic,first_evt_cluster=Monte_Carlo(snoise,N,u,Charge_strips,Position_strips,seed,threshold)
            Data.append([s,u,Charge,sCharge,E_linear,sE_linear,E_quadratic,sE_quadratic,E_logarithmic,sE_logarithmic])
            #Plot_Cluster(first_evt_cluster,u,s,Labels)
            L=[s,u,first_evt_cluster['Position_strip'],first_evt_cluster['Charge_strip_noise']]

    df = pd.DataFrame(Data, columns=["s","x","Q","sQ","E_linear","sE_linear","E_quadratic","sE_quadratic","E_logarithmic","sE_logarithmic"])
    return df
    
def Plot_Charge_or_Error(Name,slist, df, markers, markersizes,ymin,ymax,Title,Ylabel):
    for s in slist:
        t=r'$\sigma $ '+' '+str(s)+"mm"
        plt.errorbar(df.loc[df['s'].isin([s]),'x'], df.loc[df['s'].isin([s]),Name], marker=markers[slist.index(s)],markersize=markersizes[slist.index(s)], linestyle='' ,label=t,xerr=0, yerr=df.loc[df['s'].isin([s]),'sE_linear'])
    
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    
    ax.legend(handles, labels,loc='upper right')
    plt.title(Title)
    plt.ylabel(Ylabel) 
    plt.ylim(ymin,ymax)
    plt.xlabel(r'$\mu$ (mm)')
    plt.grid(True)
    plt.show()

    
def RMS(v):
    return np.sqrt((v*v).sum()/len(v))

def RMS_uncertainty_propagation(v,sv):
    return np.sqrt(sum((v*sv)**2))/(len(v)*RMS(v))

def RMS_Dataframe(slist,df):
    Data=[]
    for s in slist:
        Data.append([
                     s,
                     RMS(df.loc[df['s'].isin([s]),'E_linear']), 
                     RMS_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_linear'],df.loc[df['s'].isin([s]),'sE_linear']), 
                     RMS(df.loc[df['s'].isin([s]),'E_quadratic']), 
                     RMS_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_quadratic'],df.loc[df['s'].isin([s]),'sE_quadratic']), 
                     RMS(df.loc[df['s'].isin([s]),'E_logarithmic']), 
                     RMS_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_logarithmic'],df.loc[df['s'].isin([s]),'sE_logarithmic'])
                     ])
        
    RMS_df = pd.DataFrame(Data,columns=['s','E_linear','sE_linear','E_quadratic','sE_quadratic','E_logarithmic','sE_logarithmic'])
    return RMS_df

def Plot_RMS(RMS_df,Labels):
    plt.title(Labels["Title"]) 
    plt.errorbar(RMS_df['s'],RMS_df['E_linear'],marker=markers[0], markersize=markersizes[0]*1.5,linestyle='' ,label="Linear weight",xerr=0, yerr=RMS_df['sE_linear'])
    plt.errorbar(RMS_df['s'],RMS_df['E_quadratic'],marker=markers[1],  markersize=markersizes[1]*1.5, linestyle='' ,label="Squared weight",xerr=0, yerr=RMS_df['sE_quadratic'])
    plt.errorbar(RMS_df['s'],RMS_df['E_logarithmic'],marker=markers[2],  markersize=markersizes[2]*1.5,linestyle='' ,label="Logarithmic weight",xerr=0, yerr=RMS_df['sE_logarithmic'])
    plt.ylabel(Labels["ylabel"]) 
    
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    ax.legend(handles, labels,loc='upper right')
    plt.xlabel(r'$\sigma$ (mm)')  
    plt.grid(True)
    plt.show()

#%%
'''For this study we have modeled a 1 cm readout, from -5 mm to 5 mm  where the coordinate 0 is in the middle of the
strips pattern. The clouds were generated in the range −0.8 mm to 0.8 mm, so the cloud didn’t
suffer edge effect. The clouds of electron colected by conductive strips were modeled as normalized Gaussian functions centered
in μ, which we call u, with standard deviation σ, which we call s.
The strips have 0.2 mm width and ~0.39 mm pitch, which means the distance between two consecutive strips.'''

#%%############################################################################
#Chose the parameters for the 1D readout geometry
###############################################################################
strip_width=0.2#(mm). The typical strip width is 0.2mm
pitch=100/256#(mm). The typical strip pitch is 0.39mm (100mm/256strips)

###############################################################################
#Defining the strip certer positions
#The region R is where we are going to calculate the charge collected by the strips is:
#   first_strip_pos< R< last_strip_pos
#The strip_centers list define the position of the strip center 
###############################################################################
first_strip_pos=-5#(mm)
last_strip_pos=5#(mm)
strip_centers=np.arange(first_strip_pos,last_strip_pos,pitch)

###############################################################################
#Defining positions of the electron cloud centers
#The region R2 is where we are going to calculate the charge collected by the strips is:
#   first_cloud_pos< R2 <last_cloud_pos
#The electron_cloud_centers list define the position where the electron cloud centers are going to be
###############################################################################
first_cloud_pos=-2*pitch#(mm)
last_cloud_pos=2*pitch#(mm)
electron_cloud_centers=np.arange(first_cloud_pos,last_cloud_pos,pitch/10)

###############################################################################
#Parameters for the clusterization
#the seed  value must be less than the maximum charge that a strip can collect considering the electron cloud choosen
#A charge value collected by a strip that is lower than the threshold is not considered part of an electron cloud by the algorithm
###############################################################################
seed=0.08
threshold=0
###############################################################################7

###############################################################################
#list of standard deviations for the electron clouds modeled as gaussians
slist=[0.2,0.25,0.3,0.35,0.4]
#sigma for the white noise which is drawn from a gaussian distribution
snoise=0.01
N=1000   

df=Dataframe(strip_width,pitch,slist,snoise,N,strip_centers,electron_cloud_centers,seed,threshold)

#%%############################################################################
#Graphs that depend on the electron cloud position relative to the readout strips
###############################################################################
#List with the indexes of the cloud widths (sigmas or standar deviations) we are going to plot in the Graph of the Errors
markers=['s','^','P','*','x']
markersizes=[4,4,4,14/3,4]
###############################################################################
ymax=0.6
ymin=0.45
Title={"Q":"Collected charge"}
Title_pt={"Q":'Carga coletada pelas fitas'}
Ylabel={"Q":'Fraction of collected charge'}
Ylabel_pt={"Q":'Fração de carga coletada'}
Plot_Charge_or_Error('Q',slist, df, markers, markersizes,ymin,ymax,Title["Q"],Ylabel['Q'])

###############################################################################
ymin=-0.045
ymax=0.045
#ymin=-0.055
#ymax=0.3
Title={"E_linear":"Position reconstruction using a linear weight",
       "E_quadratic":"Position reconstruction using a squared weight",
       "E_logarithmic": "Position reconstruction using a logarithmic weight"}
Title_pt={ "E_linear":"Reconstrução da posição usando peso linear",
          "E_quadratic":"Reconstrução da posição usando peso quadrático",
          "E_logarithmic": "Reconstrução da posição usando peso logaritmico"}
Ylabel={"E":'Error (mm)'}
Ylabel_pt={"E":'Erro (mm)'}
Plot_Charge_or_Error("E_linear",slist, df, markers, markersizes,ymin,ymax,Title["E_linear"],Ylabel['E'])
Plot_Charge_or_Error("E_quadratic",slist, df, markers, markersizes,ymin,ymax,Title["E_quadratic"],Ylabel['E'])
Plot_Charge_or_Error("E_logarithmic",slist, df, markers, markersizes,ymin,ymax,Title["E_logarithmic"],Ylabel['E'])

#%%############################################################################
#RMS Graph that depend on the electron cloud width (standard deviation of the gaussian)
###############################################################################
RMS_df=RMS_Dataframe(slist,df)

Labels={"Title":"Comparing methods for position reconstruction",
        "ylabel":'RMS of the Error (mm)'}

Labels_pt={"Title":"Comparing methods for position reconstruction",
           "ylabel":'Erro do RMS (mm)'}

Plot_RMS(RMS_df,Labels)
##############################################################################