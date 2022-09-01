'''
Author: Thais Silva Abelha
Email: thais.silva.abelha@gmail.com

This work was done during a undergraduated research "Comparing different methods of position
reconstruction considering 1D readout of GEM detectors"
'''

import matplotlib.pyplot as plt
from matplotlib import container
import numpy as np
from cloud_collection_simulation import Cloud_Collection_Simulation
import pandas as pd

'''
WAC = weighted average considering the charge as a weight
WACC = weighted average considering the charge*charge as a weighet
WAlnC = weighted average considering the ln(charge) as a weight
Q = charge
s=standard deviation or sigma σ of the gaussian
u=mean of the gaussian or electron cloud center

prefix "s" stands for standard deviation or sigma σ

RMS = Root Mean Square 
'''
def Plot_Cluster(first_evt_cluster,u,s,Labels):
    fig, ax = plt.subplots()
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
    
def Dataframe(strip_width,pitch,slist,signoiselist,Position_strips,electron_cloud_centers):
    '''For the case when    pitch=0.39mm, the electron cloud center is going to be reconstructed approximately -0.8mm<x<0.8mm
    Also the electron cloud center is going to be reconstructed with steps 10 times smaller than the strip pitch'''
    
    #You can comment "Plot_Cluster()" It is just a Graph of control then it is going to run faster
    #Plot_Cluster() is used just to see if verything is fine with the algorithm
    
    Labels={"Title":"Cluster - Event 1/1000",
            "xlabel":"Position of the strip center (mm)",
            "ylabel":"Charge collected by each strip (mm)",
            "Label_Charges":"Charge collected by strips",
            "Label_Cluster":"Cluster"}
    
    Labels={"Title":"Cluster - Evento 1/1000",
            "xlabel":"Posição do centro da fita (mm)",
            "ylabel":"Carga coletada por cada fita",
            "Label_Charges":"Carga coletada pelas fitas",
            "Label_Cluster":"Cluster"}
    Data=[]
    for s in slist:
        sig_noise=snoiselist[slist.index(s)]
        for u in electron_cloud_centers:
            Charge,sCharge,E_WAC,sE_WAC,E_WACC,sE_WACC,E_WAlnC,sE_WAlnC,first_evt_cluster=Cloud_Collection_Simulation(strip_width,pitch, s, u,sig_noise,Position_strips).Monte_Carlo()
            Data.append([s,u,Charge,sCharge,E_WAC,sE_WAC,E_WACC,sE_WACC,E_WAlnC,sE_WAlnC])
            Plot_Cluster(first_evt_cluster,u,s,Labels)
    df = pd.DataFrame(Data, columns=["s","x","Q","sQ","E_WAC","sE_WAC","E_WACC","sE_WACC","E_WAlnC","sE_WAlnC"])
    return df
    
def Plot_Charge_Error(Name,slist, df, markers, markersizes,ymin,ymax,Title,Ylabel):
    for s in slist:
        t=r'$\sigma $ '+' '+str(s)+"mm"
        plt.errorbar(df.loc[df['s'].isin([s]),'x'], df.loc[df['s'].isin([s]),Name], marker=markers[slist.index(s)],markersize=markersizes[slist.index(s)], linestyle='' ,label=t,xerr=0, yerr=df.loc[df['s'].isin([s]),'sE_WAC'])
    
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
                     RMS(df.loc[df['s'].isin([s]),'E_WAC']), 
                     RMS_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_WAC'],df.loc[df['s'].isin([s]),'sE_WAC']), 
                     RMS(df.loc[df['s'].isin([s]),'E_WACC']), 
                     RMS_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_WACC'],df.loc[df['s'].isin([s]),'sE_WACC']), 
                     RMS(df.loc[df['s'].isin([s]),'E_WAlnC']), 
                     RMS_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_WAlnC'],df.loc[df['s'].isin([s]),'sE_WAlnC'])
                     ])
        
    RMS_df = pd.DataFrame(Data,columns=['s','E_WAC','sE_WAC','E_WACC','sE_WACC','E_WAlnC','sE_WAlnC'])
    return RMS_df

def Plot_RMS(RMS_df,Labels):
    plt.title(Labels["Title"]) 
    plt.errorbar(RMS_df['s'],RMS_df['E_WAC'],marker=markers[0], markersize=markersizes[0]*1.5,linestyle='' ,label="Linear weight",xerr=0, yerr=RMS_df['sE_WAC'])
    plt.errorbar(RMS_df['s'],RMS_df['E_WACC'],marker=markers[1],  markersize=markersizes[1]*1.5, linestyle='' ,label="Squared weight",xerr=0, yerr=RMS_df['sE_WACC'])
    plt.errorbar(RMS_df['s'],RMS_df['E_WAlnC'],marker=markers[2],  markersize=markersizes[2]*1.5,linestyle='' ,label="Logarithmic weight",xerr=0, yerr=RMS_df['sE_WAlnC'])
    plt.ylabel(Labels["ylabel"]) 
    
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    ax.legend(handles, labels,loc='upper right')
    plt.xlabel(r'$\sigma$ (mm)')  
    plt.grid(True)
    plt.show()

#%%
'''For this study we have modeled a 10mm readout, from -5mm to 5mm  where the coordinate 0 is the middle of the
strips pattern, in this case. The clouds were generated in the range −0.8mm to 0.8mm, so the cloud didn’t
suffer edge effect. The charge clouds were modeled as normalized Gaussian functions centered
in μ with standard deviation σ'''
###############################################################################
#Chose the parameters for the 1D readout geometry
###############################################################################
strip_width=0.2 #The typical strip width is 0.2mm
#The typical strip pitch is 0.39mm (100mm/256strips)
pitch=0.39
#list of standard deviations for the electron clouds modeled as gaussians
slist=[0.2,0.25,0.3,0.35,0.4]
#sigma for the white noise which is drawn from a gaussian distribution
snoiselist=[0.01,0.01,0.01,0.01,0.01]
###############################################################################
#Calculating the Collected Charge and the Errors of Position Reconstruction
#Region we are going to calculate the charge collected by the strips and strip width
###############################################################################
first_strip_pos=-5#(mm)
last_strip_pos=5#(mm)
first_cloud_pos=-2*pitch#(mm)
last_cloud_pos=2*pitch#(mm)

###############################################################################
#Defining the strip certer positions
#And positions of the electron cloud centers
###############################################################################
strip_centers=np.arange(first_strip_pos,last_strip_pos,pitch)
electron_cloud_centers=np.arange(first_cloud_pos,last_cloud_pos,pitch/10)
    
'''Now we are just organazing the data in a DataFrame. We are calculating the data for the colected charge and 
the Errors whith linear, quadratic and logarithmic weights, considering different sigmas for the gaussian that 
stands for the electron cloud'''

df=Dataframe(strip_width,pitch,slist,snoiselist,strip_centers,electron_cloud_centers)
#%%############################################################################
#Graphs that depend on electron cloud position relative to the readout strips
###############################################################################
#List with the indexes of the cloud widths (sigmas or standar deviations) we are going to plot in the Graph of the Errors
s_indexes_toplot=[0,2,4] 
markers=['s','^','P','*','x']
markersizes=[4,4,4,14/3,4]
ymax=0.6
ymin=0.45

Title={"Q":"Collected charge",
       "E_WAC":"Position reconstruction using a linear weight",
       "E_WACC":"Position reconstruction using a squared weight",
       "E_WAlnC": "Position reconstruction using a logarithmic weight"}

Title_pt={"Q":'Carga coletada pelas fitas',
          "E_WAC":"Reconstrução da posição usando peso linear",
          "E_WACC":"Reconstrução da posição usando peso quadrático",
          "E_WAlnC": "Reconstrução da posição usando peso logaritmico"}

Ylabel={"Q":'Fraction of collected charge',"E":'Error (mm)'}
Ylabel_pt={"Q":'Fração de carga coletada',"E":'Erro (mm)'}

Plot_Charge_Error('Q',slist, df, markers, markersizes,ymin,ymax,Title["Q"],Ylabel['Q'])

ymin=-0.045
ymax=0.045
#ymin=-0.055
#ymax=0.3

Plot_Charge_Error("E_WAC",slist, df, markers, markersizes,ymin,ymax,Title["E_WAC"],Ylabel['E'])
Plot_Charge_Error("E_WACC",slist, df, markers, markersizes,ymin,ymax,Title["E_WACC"],Ylabel['E'])
Plot_Charge_Error("E_WAlnC",slist, df, markers, markersizes,ymin,ymax,Title["E_WAlnC"],Ylabel['E'])

#%%############################################################################
#RMS Graph that depend on the electron cloud width (standard deviation of the gaussian)
###############################################################################
RMS_df=RMS_Dataframe(slist,df)

Labels={"Title":"Comparing methods for position reconstruction",
        "ylabel":'RMS of the Error (mm)'}

Labels_pt={"Title":"Comparing methods for position reconstruction",
           "ylabel":'Erro do RMS (mm)'}

Plot_RMS(RMS_df,Labels)