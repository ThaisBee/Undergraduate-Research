'''
Author: Thais Silva Abelha
Email: thais.silva.abelha@gmail.com

This work was done during a undergraduated research "Comparing different methods of position
reconstruction considering 1D readout of GEM detectors"
'''

import numpy as np
import math
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import matplotlib.cm as cm

from models.significativos import dois_significativos
from models.monte_carlo import Monte_Carlo
from models.strips_integrator import Strips_Integrator

#%%
def x_projection(df):
    xprojection=df.groupby('x')['E'].sum().reset_index()
    Bin=xprojection.iat[1,0]-xprojection.iat[0,0]
    xprojection['E']=xprojection['E']/(xprojection['E']*Bin).sum()    
    return xprojection['x'],xprojection['E']

def interpolation(x_coordinate,charge_density):
    f=interp1d(x_coordinate, charge_density,kind='cubic')
    n=50
    x=np.linspace(x_coordinate[0],x_coordinate[len(x_coordinate)-1],num=n*len(x_coordinate),endpoint=True)
    E=f(x)
    return x,E

def Dataframe(strip_width,pitch,s,snoise,N,Position_strips,x_coordinate,charge_fraction,electron_cloud_centers,seed,threshold):
    '''For the case when pitch=0.39mm, the electron cloud center is going to be reconstructed approximately -0.8mm<x<0.8mm
    Also the electron cloud center is going to be reconstructed with steps 10 times smaller than the strip pitch'''

    Data=[]
    for u in electron_cloud_centers:
        x_new_coordinate=x_coordinate+u
        Charge_strips=Strips_Integrator(strip_width,pitch,snoise,Position_strips).ChargeBins_Strip(x_new_coordinate,charge_fraction)
        Charge,sCharge,E_linear,sE_linear,E_quadratic,sE_quadratic,E_logarithmic,sE_logarithmic,first_evt_cluster=Monte_Carlo(snoise,N,u,Charge_strips,Position_strips,seed,threshold)
        Data.append([s,u,Charge,sCharge,E_linear,sE_linear,E_quadratic,sE_quadratic,E_logarithmic,sE_logarithmic])
        #Plot_Cluster(first_evt_cluster,u,s,seed,threshold)

    df = pd.DataFrame(Data, columns=["s","x","Q","sQ","E_linear","sE_linear","E_quadratic","sE_quadratic","E_logarithmic","sE_logarithmic"])
    return df

def Func_Normal(x, x0, sigma):
    return (1/(sigma*np.sqrt(2*math.pi))) * np.exp(-(x - x0)**2 / (2.0 * sigma**2))

def Plot_Normal(x_coordinate,charge_fraction,Position_strips,popt,pcov):
    inc=np.sqrt(np.diagonal(pcov))
       
    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0,height_ratios=[2,1])
    axs = gs.subplots(sharex=True)
    title="Distribuição das cargas da nuvem"
    fig.suptitle(title)
    
    axs[0].plot(x_coordinate, charge_fraction, 'b+:', label='Nuvem de elétrons, Garfield ++',color='blue')
    axs[0].plot(x_coordinate, Func_Normal(x_coordinate, *popt), 'r-',label='fit:'+r' $\mu$='+dois_significativos(popt[0],inc[0]).forma2()+"mm, "+r'$\sigma$='+dois_significativos(popt[1],inc[1]).forma2()+"mm")
    axs[0].set_ylabel('Carga (U.A.)')
    axs[0].set_ylim(0,3)
    axs[0].legend(loc='upper right')
    axs[0].xaxis.set_major_locator(MultipleLocator(0.4))
    axs[0].xaxis.set_minor_locator(MultipleLocator(0.2))
    
    Residuals=charge_fraction-Func_Normal(x_coordinate,*popt)
    axs[1].plot(x_coordinate,Residuals,linestyle="None",marker="*",markersize=1,color='b')
    axs[1].set_xlabel('Coordenadas em x (mm)')
    axs[1].set_ylabel('Resíduos (U.A.)')
    axs[1].set_xlim(-1.2,1.2)
    axs[1].grid(axis="y") 
    
    rectangles=[k for k in Position_strips if k >-1.2 and k<1.2]
    trans0 = transforms.blended_transform_factory(axs[0].transData, axs[0].transAxes)
    trans1 = transforms.blended_transform_factory(axs[1].transData, axs[1].transAxes)
    
    for i in rectangles:
        rect0 = mpatches.Rectangle((i, 0), width=0.2, height=1, transform=trans0,color='red', alpha=0.5)
        rect1 = mpatches.Rectangle((i, 0), width=0.2, height=1, transform=trans1,color='red', alpha=0.5)

        axs[0].add_patch(rect0)
        axs[1].add_patch(rect1)

    plt.show()
    return 

def Plot_Hist2d(x,y,z,nbins):
    fig, ax = plt.subplots()
    counts, xedges, yedges, im = ax.hist2d(x, y,bins=nbins,weights=z,cmap=cm.jet)
    fig.colorbar(im, ax=ax,label="Elétrons")
    ax.set_ylim(-1,0.85)
    ax.set_xlim(-1,0.85)
    ax.set_xlabel('Coordeadas em x (mm)')
    ax.set_ylabel('Coordenadas em y (mm)')
    plt.title("Nuvem de elétrons depois da multiplicação \n Fóton de 8 KeV "+r"$Ar/CO_2$ (70/30)")
    plt.plot(0,0,"*",color="black")
    plt.show()
    
def Plot_Cluster(first_evt_cluster,u,s,seed,threshold):
    fig, ax = plt.subplots()
    plt.axhline(seed,label="Seed", linestyle=(0,(5,10)),color="green",linewidth=1.0)
    plt.axhline(threshold,label="Threshold",color="red",linewidth=1.0)
    #print(first_evt_cluster['Position_strip'])
    #print(first_evt_cluster['Charge_strip_noise'])
    ax.plot(first_evt_cluster['Position_strip'],first_evt_cluster['Charge_strip_noise'],"o", markersize=4,label="Charge collected by strips")
    ax.plot(first_evt_cluster['Position_Cluster'],first_evt_cluster['Charge_Cluster'],"*",markersize=7,label="Cluster")
    fig.suptitle("Cluster - Event 1/1000")
    ax.set_xlabel("Position of the strip center (mm)")
    ax.set_ylabel("Charge collected by each strip (charge fraction)")
   
    '''ax.plot(first_evt_cluster['Position_strip'],first_evt_cluster['Charge_strip_noise'],"o", markersize=4,label="Carga coletada pelas strips")
    ax.plot(first_evt_cluster['Position_Cluster'],first_evt_cluster['Charge_Cluster'],"*",markersize=7,label="Cluster (nuvem de elétrons)")
    fig.suptitle("Clusterização")
    ax.set_xlabel("Posição do centro da fita (mm)")
    ax.set_ylabel("Carga coletada por cada fita (mm)")'''
    
    ax.grid()
    ax.set_ylim(-0.05,0.35)

    textstr = '\n'.join((
    r'$\mu=%.2f$' % (u, )+"mm",
    r'$\sigma=%.2f$' % (s, )+"mm"
    ))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.08, 0.95, textstr, transform=ax.transAxes, fontsize=11,verticalalignment='top', bbox=props)
    ax.legend(loc='upper right')

    plt.show()
                   
def Plot_Errors(df,label):
    plt.errorbar(df['x'],df['E_linear'],marker='o',markersize=3, linestyle='',label=label['E_linear'],xerr=0, yerr=df['sE_linear'])
    plt.errorbar(df['x'],df['E_quadratic'],marker='o',markersize=3, linestyle='',label=label['E_quadratic'],xerr=0, yerr=df['sE_quadratic'])
    plt.errorbar(df['x'],df['E_logarithmic'],marker='o',markersize=3, linestyle='',label=label['E_logarithmic'],xerr=0, yerr=df['sE_logarithmic'])
    plt.title(label['Title'])
    plt.legend(loc='upper right')
    plt.ylabel(label['ylabel'])
    plt.xlabel(r'$\mu$ (mm)')
    plt.grid()
    plt.show()

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
###############################################################################

#%%############################################################################
#Reading file with the electron cloud data which is a Garfiled++ output
###############################################################################
file="data/evt1.txt"
electron_cloud_df=pd.read_csv(file,sep='\s+',header=None)
electron_cloud_df.rename({0: 'x', 1: 'y',2:'E'}, axis=1, inplace=True)
electron_cloud_df[['x','y']]=electron_cloud_df[['x','y']]/1000


#%%############################################################################
#Input data: Electron Cloud -> Electron cloud projection in x direction
# -> Interpolation of the data
###############################################################################
x_coordinate,charge_fraction=x_projection(electron_cloud_df)
x_coordinate_interpolate,charge_fraction_interpolate=interpolation(x_coordinate,charge_fraction)

'''Since we are going to fit a gaussian function, we can try those expressions for getting an approximation
for the mean and sigma values of the fit'''

mean = sum(x_coordinate * charge_fraction) / sum(charge_fraction)
sigma = np.sqrt(sum(charge_fraction * (x_coordinate - mean)**2) / sum(charge_fraction))
popt,pcov = curve_fit(Func_Normal, x_coordinate, charge_fraction, p0=[mean, sigma])

###############################################################################
#we fit on the electron cloud generated by Garfield++ software, just too see how gaussian it is
#And the s value is the standard deviation of the gaussian function
#We modeled noise which comes from the eletronics as a gaussian distribution 
#and its standard deviation is represented by snoise
##############################################################################
s=popt[1]# popt[1] is the fit result for sigma 
snoise=0.01
N=1000
###############################################################################
#Plotting the original 2D data and its projection in x direction
###############################################################################
Plot_Hist2d(electron_cloud_df['x'], electron_cloud_df['y'], electron_cloud_df['E'], len(x_coordinate))
Plot_Normal(x_coordinate, charge_fraction,strip_centers,popt,pcov)

#%%############################################################################
#Calculating the Collected Charge considering the noise, using Monte Carlo's method
#A clustarization algorithm, made especially for this problem, separetes the electron cloud signal from the noise
#and the Errors of Position Reconstruction are calculetade for 1000 events
###############################################################################
df=Dataframe(strip_width,pitch,s,snoise,N,strip_centers,x_coordinate_interpolate,charge_fraction_interpolate,electron_cloud_centers,seed,threshold)

###############################################################################
#%% Plot of the errors of the position of the electron cloud center calculated by different methods
#to understand how this errors are calculated go to errors.py
##############################################################################
'''Errors of position reconstructions considering different positions of the electron cloud center
relative to the readout strips'''

labels={'E_linear':'Linear weight',
       'E_quadratic':'Squared weight',
       'E_logarithmic':'Logarithmic weight',
       'Title': 'Electron Cloud generated by Garfield++',
       'ylabel': 'Errors'
       }
   
labels_pt={'E_linear':'Peso linear',
       'E_quadratic':'Peso quadrático',
       'E_logarithmic':'Peso logarítmico',
       'Title':'Nuvem de elétrons geradas pelo Garfield++',
       'ylabel':'Erros (mm)'
       }

Plot_Errors(df,labels_pt)

