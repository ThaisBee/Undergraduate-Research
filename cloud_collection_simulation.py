import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
from cluster import Cluster
from errors import Errors
import pandas as pd

''''
This simulates all electron clouds and the noise.
It also gives the error of the electron cloud center compared to the electron cloud 
calculated using the three different methods.

'''

class Cloud_Collection_Simulation:

    def __init__(self,strip_width,pitch,sigma,u,sig_noise,Position_strips,N=1000,seed=0.04):
        self.strip_width=strip_width
        self.pitch=pitch
        self.s=sigma
        self.u=u
        self.sig=sig_noise
        self.Position_strips=Position_strips
        self.seed=seed
        self.N=N
        
    def Gauss(self,x,s,u):
        return (1/(s*np.sqrt(2*math.pi)))*math.exp( -((x-u)**2) / (2*s**2)  )
    
    def IntegraStrip(self,p):
        I = quad(self.Gauss, p-self.strip_width/2, p+self.strip_width/2, args=(self.s,self.u))
        return I[0]

    def Charge_StripPos(self):
        Charges_strips=[self.IntegraStrip(p) for p in self.Position_strips]
        return Charges_strips,self.Position_strips
    
    def Monte_Carlo(self):
        Data=[]                           
        Charge_strip,Position_strip=self.Charge_StripPos()
        
        for i in range(self.N):
            Noise=np.random.normal(0,self.sig,len(Charge_strip))
            
            Charge_strip_noise=Charge_strip+Noise
            Charge_Cluster,Position_Cluster=Cluster(self.seed,Charge_strip_noise,Position_strip).Find_Cluster()
            Data.append([sum(Charge_Cluster),Errors(self.u, Charge_Cluster,Position_Cluster).E_WAC(),Errors(self.u, Charge_Cluster,Position_Cluster).E_WACC(),Errors(self.u, Charge_Cluster,Position_Cluster).E_WAlnC()])

            if i==0:
                first_evt_cluster={"Position_strip":Position_strip,
                                   "Charge_strip_noise":Charge_strip_noise,
                                   "Position_Cluster":Position_Cluster,
                                   "Charge_Cluster":Charge_Cluster }
                
        df=pd.DataFrame(Data, columns=["Q","E_WAC","E_WACC","E_WAlnC"])
        
        '''
        Considering N events we are going to calaculate Error means and Stardard deviations of the Errors menas
        '''
        
        Charge_mean=np.mean(df['Q'])
        E_WAC=np.mean(df['E_WAC'])
        E_WACC=np.mean(df['E_WACC'])
        E_WAlnC=np.mean(df['E_WAlnC'])
        
        sCharge=np.std(df['Q'],ddof=1)/np.sqrt(self.N)
        sE_WAC=np.std(df['E_WAC'],ddof=1)/np.sqrt(self.N)
        sE_WACC=np.std(df['E_WACC'],ddof=1)/np.sqrt(self.N)   
        sE_WAlnC=np.std(df['E_WAlnC'],ddof=1)/np.sqrt(self.N)
        
        return Charge_mean,sCharge,E_WAC,sE_WAC,E_WACC,sE_WACC,E_WAlnC,sE_WAlnC,first_evt_cluster