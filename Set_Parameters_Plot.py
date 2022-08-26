#%%
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
s=standard deviation or sigma σ of the gaussian
u=mean of the gaussian or electron cloud center

prefix "s" stands for standard deviation or sigma σ

RMS = Root Mean Square 
'''

def Dataframe(strip_width,pitch,slist,signoiselist,Position_strips):
    '''For the case when    pitch=0.39mm, the electron cloud center is going to be reconstructed approximately -0.8mm<x<0.8mm
    Also the electron cloud center is going to be reconstructed with steps 10 times smaller than the strip pitch'''
    #u0=-2*pitch  #First position of the electron Cloud relative to the readout strips (u = mean value of the gaussian)
    u0=-2*pitch
    du=pitch/10  #Step for the electron cloud positions
    
    Data=[]
    for s in slist:
        sig_noise=snoiselist[slist.index(s)]
        u=u0
        while u<2*pitch:
            Charge,sCharge,E_WAC,sE_WAC,E_WACC,sE_WACC,E_WAlnC,sE_WAlnC=Cloud_Collection_Simulation(strip_width,pitch, s, u,sig_noise,Position_strips).Monte_Carlo()
            Data.append([s,u,Charge,sCharge,E_WAC,sE_WAC,E_WACC,sE_WACC,E_WAlnC,sE_WAlnC])
            u=u+du
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
    
def rms(v):
    return np.sqrt((v*v).sum()/v.shape[0])

def rms_uncertainty_propagation(v,sv):
    return np.sqrt(sum((v*sv)**2))/(len(v)*rms(v))

def rms_Dataframe(slist,df):
    Data=[]
    for s in slist:
        Data.append([s,
                     rms(df.loc[df['s'].isin([s]),'E_WAC']), 
                     rms_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_WAC'],df.loc[df['s'].isin([s]),'sE_WAC']), 
                     rms(df.loc[df['s'].isin([s]),'E_WACC']), 
                     rms_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_WACC'],df.loc[df['s'].isin([s]),'sE_WACC']), 
                     rms(df.loc[df['s'].isin([s]),'E_WAlnC']), 
                     rms_uncertainty_propagation(df.loc[df['s'].isin([s]),'E_WAlnC'],df.loc[df['s'].isin([s]),'sE_WAlnC'])
                     ])
        
    rms_df = pd.DataFrame(Data,columns=['s','E_WAC','sE_WAC','E_WACC','sE_WACC','E_WAlnC','sE_WAlnC'])
    return rms_df

'''
for this study we have modeled a 10mm readout, from -5mm to 5mm  where the coordinate 0 is the middle of the
strips pattern, in this case. The clouds were generated in the range −0.8mm to 0.8mm, so the cloud didn’t
suffer edge effect. The charge clouds were modeled as normalized Gaussian functions centered
in μ with standard deviation σ

'''
###############################################################################
#Input
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
###############################################################################
#Region we are going to calculate the charge collected by the strips and strip width
first_strip_pos=-5
last_strip_pos=5

#Defining the strip positions
Position_strips=[]
p=first_strip_pos
i=0
while p<(last_strip_pos):
    Position_strips.append(p)
    p=p+pitch
    i=i+1
    
'''Now we are just organazing the data in a Matrix to save run time. We are calculating the data for the colected charge and 
the Errors whith linear, quadratic and logarithmic weights, considering different sigmas for the gaussian that 
stands for the electron cloud'''
###############################################################################
#Matrix with the data
###############################################################################
#x,Charges_Matrix,sCharges_Matrix,E_WAC_Matrix,sE_WAC_Matrix,E_WACC_Matrix,sE_WACC_Matrix,E_WAlnC_Matrix,sE_WAlnC_Matrix=Charge_ErrorPosition_Matrix(strip_width,pitch,slist,snoiselist,Position_strips)
###############################################################################
df=Dataframe(strip_width,pitch,slist,snoiselist,Position_strips)
#%%
#y axis limits
ymax=0.6
ymin=0.45
s_indexes_toplot=[0,2,4] #List with the indexes of the cloud widths (sigmas) we are going to plot

markers=['s','^','P','*','x']
markersizes=[4,4,4,14/3,4]

Title={"Q":"Collected charge",
       "E_WAC":"Position reconstruction using a linear weight",
          "E_WACC":"Position reconstruction using a squared weight",
          "E_WAlnC": "Position reconstruction using a logarithmic weight"}
Titulo={"Q":'Carga coletada pelas fitas',
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

rms_df=rms_Dataframe(slist,df)

#RMS Graph
#En
plt.title("Comparing methods for position reconstruction") 
plt.errorbar(rms_df['s'],rms_df['E_WAC'],marker=markers[0], markersize=markersizes[0]*1.5,linestyle='' ,label="Linear weight",xerr=0, yerr=rms_df['sE_WAC'])
plt.errorbar(rms_df['s'],rms_df['E_WACC'],marker=markers[1],  markersize=markersizes[1]*1.5, linestyle='' ,label="Squared weight",xerr=0, yerr=rms_df['sE_WACC'])
plt.errorbar(rms_df['s'],rms_df['E_WAlnC'],marker=markers[2],  markersize=markersizes[2]*1.5,linestyle='' ,label="Logarithmic weight",xerr=0, yerr=rms_df['sE_WAlnC'])
plt.ylabel('Error RMS (mm)') #En

'''
#Pt
plt.title("Comparando diferente métodos para reconstrução da posição") #Pt
plt.errorbar(slist,RMS_WAC, marker=markers[0], markersize=markersizes[0]*1.5, linestyle='' , label="Peso linear", xerr=0, yerr=sRMS_WAC)
plt.errorbar(slist,RMS_WACC, marker=markers[1], markersize=markersizes[1]*1.5, linestyle='' , label="Peso quadrático", xerr=0, yerr=sRMS_WACC)
plt.errorbar(slist,RMS_WAlnC, marker=markers[2], markersize=markersizes[2]*1.5, linestyle='' , label="Peso logarítmico", xerr=0, yerr=sRMS_WAlnC)
plt.ylabel('RMS do Erro (mm)') 
'''

ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels,loc='upper right')
plt.xlabel(r'$\sigma$ (mm)')  
plt.grid(True)
plt.show()