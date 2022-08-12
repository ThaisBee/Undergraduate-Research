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

'''
WAC = weighted average considering the charge as a weight
WACC = weighted average considering the charge*charge as a weighet
WAlnC = weighted average considering the ln(charge) as a weight

prefix "s" means standard deviation or sigma σ

RMS = Root Mean Square 
'''

def Cargas(first_strip_pos,last_strip_pos,strip_width,pitch,s,sig_noise):
    
    u0=-2*pitch  #First position of the electron Cloud relative to the readout strips (u = mean value of the gaussian)
    du=pitch/10  #Step for the electron cloud positions
    n=int(abs(u0*2)/du)+1
    x=np.zeros(n) #list of all electron cloud center positions

    #Charge Reconstructions and its uncertanty
    Charge_list=np.zeros(n)
    sCharge_list=np.zeros(n)
    
    #Errors od position reconstructions and their uncertanties
    E_WAC_list=np.zeros(n)
    sE_WAC_list=np.zeros(n)
    E_WACC_list=np.zeros(n)
    sE_WACC_list=np.zeros(n)
    E_WAlnC_list=np.zeros(n)
    sE_WAlnC_list=np.zeros(n)
    
    u=u0
    for i in range(n):
        Charge_list[i],sCharge_list[i],E_WAC_list[i],sE_WAC_list[i],E_WACC_list[i],sE_WACC_list[i],E_WAlnC_list[i],sE_WAlnC_list[i]=Cloud_Collection_Simulation(first_strip_pos,last_strip_pos,strip_width,pitch, s, u,sig_noise).Monte_Carlo()
        x[i]=u
        u=u+du
    
    return x,Charge_list,sCharge_list,E_WAC_list,sE_WAC_list,E_WACC_list,sE_WACC_list,E_WAlnC_list,sE_WAlnC_list

def RMS_func(E_WAC_Matrix,sE_WAC_Matrix,E_WACC_Matrix,sE_WACC_Matrix,E_WAlnC_Matrix,sE_WAlnC_Matrix,slist):
    RMS_WAC=np.zeros(len(slist))
    RMS_WACC=np.zeros(len(slist))
    RMS_WAlnC=np.zeros(len(slist))
    sRMS_WAC=np.zeros(len(slist))
    sRMS_WACC=np.zeros(len(slist))
    sRMS_WAlnC=np.zeros(len(slist))
    
    for i in range(len(slist)):
        RMS_WAC[i]=np.sqrt(sum(E_WAC_Matrix[i]*E_WAC_Matrix[i])/len(E_WAC_Matrix[i]))
        RMS_WACC[i]=np.sqrt(sum(E_WACC_Matrix[i]*E_WACC_Matrix[i])/len(E_WACC_Matrix[i]))
        RMS_WAlnC[i]=np.sqrt(sum(E_WAlnC_Matrix[i]*E_WAlnC_Matrix[i])/len(E_WAlnC_Matrix[i]))
        
        sRMS_WAC[i]=sum((E_WAC_Matrix[i]*sE_WAC_Matrix[i])**2)/(RMS_WAC[i]*len(E_WAC_Matrix[i]))
        sRMS_WACC[i]=sum((E_WACC_Matrix[i]*sE_WACC_Matrix[i])**2)/(RMS_WACC[i]*len(E_WACC_Matrix[i]))
        sRMS_WAlnC[i]=sum((E_WAlnC_Matrix[i]*sE_WAlnC_Matrix[i])**2)/(RMS_WAlnC[i]*len(E_WAlnC_Matrix[i]))
        
    return RMS_WAC, sRMS_WAC, RMS_WACC, sRMS_WACC, RMS_WAlnC, sRMS_WAlnC

def Mean_func(E_WAC_Matrix,sE_WAC_Matrix,E_WACC_Matrix,sE_WACC_Matrix,E_WAlnC_Matrix,sE_WAlnC_Matrix,slist):
    Mean_WAC=np.zeros(len(slist))
    Mean_WACC=np.zeros(len(slist))
    Mean_WAlnC=np.zeros(len(slist))
    
    sMean_WAC=np.zeros(len(slist))
    sMean_WACC=np.zeros(len(slist))
    sMean_WAlnC=np.zeros(len(slist))
    
    for i in range(len(slist)):
        Mean_WAC[i]=(np.mean(E_WAC_Matrix[i]))
        Mean_WACC[i]=(np.mean(E_WACC_Matrix[i]))
        Mean_WAlnC[i]=(np.mean(E_WAlnC_Matrix[i]))
        
        sMean_WAC[i]=np.std(E_WAC_Matrix[i],ddof=1)/np.sqrt(len(E_WAC_Matrix[i]))
        sMean_WACC[i]=np.std(E_WACC_Matrix[i],ddof=1)/np.sqrt(len(E_WACC_Matrix[i]))
        sMean_WAlnC[i]=np.std(E_WAlnC_Matrix[i],ddof=1)/np.sqrt(len(E_WAlnC_Matrix[i]))

    return Mean_WAC, sMean_WAC, Mean_WACC, sMean_WACC, Mean_WAlnC, sMean_WAlnC

'''
for this study we have modeled a 10 mm readout, from -5 mm to 5 mm  where the coordinate 0 is the middle of the
strips pattern, in this case. The clouds were generated in the range −0.8 mm to 0.8 mm, so the cloud didn’t
suffer edge effect. The charge clouds were modeled as normalized Gaussian functions centered
in μ with standard deviation σ

'''
#######################################################################################
#Input
######################################################################################
#region we are going to calculate the charge collected by the strips and strip width
first_strip_pos=-5
last_strip_pos=5
#######################################################################################
strip_width=0.2 #The typical strip width is 0.2mm
#The typical strip pitch is 0.39mm (100mm/256strips)
pitch=0.39
#list of standard deviations for the electron clouds modeled as gaussians
slist=[0.2,0.25,0.3,0.35,0.4]
#sigma for the white noise which is drawn from a gaussian distribution
snoiselist=[0.01,0.01,0.01,0.01,0.01]
###########################################################################
#Calculating the Collected Charge and the Errors of Position Reconstruction
############################################################################

Charges_Matrix=[]
sCharges_Matrix=[]
E_WAC_Matrix=[]
sE_WAC_Matrix=[]
E_WACC_Matrix=[]
sE_WACC_Matrix=[]
E_WAlnC_Matrix=[]
sE_WAlnC_Matrix=[]

for i in range(len(slist)):
    s=slist[i]
    sig_noise=snoiselist[i]
    x,Q,sQ,E_WAC,sE_WAC,E_WACC,sE_WACC,E_WAlnC,sE_WAlnC=Cargas(first_strip_pos, last_strip_pos, strip_width, pitch, s,sig_noise)    
    
    Charges_Matrix.append(Q)
    sCharges_Matrix.append(sQ)
    E_WAC_Matrix.append(E_WAC)
    sE_WAC_Matrix.append(sE_WAC)
    E_WACC_Matrix.append(E_WACC)
    sE_WACC_Matrix.append(sE_WACC) 
    E_WAlnC_Matrix.append(E_WAlnC)
    sE_WAlnC_Matrix.append(sE_WAlnC)


#y axis limits
ymax=0.6
ymin=0.45
L=[0,2,4] #List with the indexes of the cloud widths (sigmas) we are going to plot

markers=['s','^','P','*','x']
markersizes=[4,6*2/3,6*2/3,7*2/3,6*2/3]

for i in L:
    t=r'$\sigma $ '+' '+str(slist[i])+" mm"
    plt.errorbar(x, Charges_Matrix[i],marker=markers[i],markersize=markersizes[i], linestyle='' ,label=t,xerr=0, yerr=sCharges_Matrix[i])
    
plt.title('Collected charge')

#Legends without error bar
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels,loc='upper right')

plt.ylabel('Fraction of collected charge')
plt.xlabel(r'$\mu$ (mm)')
plt.grid(True)
plt.show()

ymin=-0.045
ymax=0.045
for i in L:
    t=r'$\sigma $ '+' '+str(slist[i])+" mm"
    plt.errorbar(x, E_WAC_Matrix[i],marker=markers[i],markersize=markersizes[i], linestyle='' ,label=t,xerr=0, yerr=sE_WAC_Matrix[i])
t="Position reconstruction using a linear weight"#'Weighted average of the charges position'
plt.title(t)

ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels,loc='upper right')
plt.ylabel('Error (mm)')
#plt.ylabel('Residuos (mm)')
plt.ylim((ymin,ymax))
plt.xlabel(r'$\mu$ (mm)')
plt.grid(True)
plt.show()

for i in L:
    t=r'$\sigma $ '+' '+str(slist[i])+" mm"
    plt.errorbar(x, E_WACC_Matrix[i],marker=markers[i],markersize=markersizes[i], linestyle='' ,label=t,xerr=0, yerr=sE_WACC_Matrix[i])
    
t="Position reconstruction using a squared weight"#'Weighted average of the '+r'$charges^{2}$'+' position '
plt.title(t)
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels,loc='upper right')
plt.ylabel('Error (mm)')
plt.xlabel(r'$\mu$ (mm)')
plt.ylim((ymin,ymax))
plt.grid(True)
plt.show()

for i in L:
    t=r'$\sigma $ '+' '+str(slist[i])+" mm"
    plt.errorbar(x, E_WAlnC_Matrix[i],marker=markers[i],markersize=markersizes[i], linestyle='' ,label=t,xerr=0, yerr=sE_WAlnC_Matrix[i])
plt.title("Position reconstruction using a logarithmic weight")

ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels,loc='upper right')
plt.ylabel('Error (mm)')
plt.xlabel(r'$\mu$ (mm)')
plt.ylim((ymin,ymax))
plt.grid(True)
plt.show()

#RMS Graph
RMS_WAC, sRMS_WAC, RMS_WACC, sRMS_WACC, RMS_WAlnC, sRMS_WAlnC=RMS_func(E_WAC_Matrix,sE_WAC_Matrix, E_WACC_Matrix,sE_WACC_Matrix, E_WAlnC_Matrix, sE_WAlnC_Matrix, slist)
plt.title("Comparing methods for position reconstruction")
plt.errorbar(slist,RMS_WAC,marker=markers[0], markersize=markersizes[0]*1.5,linestyle='' ,label="Linear weight",xerr=0, yerr=sRMS_WAC)
plt.errorbar(slist,RMS_WACC,marker=markers[1],  markersize=markersizes[1]*1.5, linestyle='' ,label="Squared weight",xerr=0, yerr=sRMS_WACC)
plt.errorbar(slist,RMS_WAlnC,marker=markers[2],  markersize=markersizes[2]*1.5,linestyle='' ,label="Logarithmic weight",xerr=0, yerr=sRMS_WAlnC)
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels,loc='upper right')
plt.ylabel('Error RMS (mm)')
plt.xlabel(r'$\sigma$ (mm)')  
plt.grid(True)
plt.show()