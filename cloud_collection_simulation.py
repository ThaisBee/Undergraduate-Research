import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad

from cluster import Cluster
from errors import Errors

''''
This simulates all electron clouds and the noise.
It also gives the error of the electron cloud center compared to the electron cloud 
calculated using the three different methods.

'''

class Cloud_Collection_Simulation:

    def __init__(self,first_strip_pos,last_strip_pos,strip_width,pitch,sigma,u,sig_noise,N=1000,seed=0.04):
        self.first_strip_pos=first_strip_pos
        self.last_strip_pos=last_strip_pos
        self.strip_width=strip_width
        self.pitch=pitch
        self.s=sigma
        self.u=u
        self.sig=sig_noise
        self.seed=seed
        self.N=N
        
    def Gauss(self,x,s,u):
        return (1/(s*np.sqrt(2*math.pi)))*math.exp( -((x-u)**2) / (2*s**2)  )
    
    def IntegraStrip(self,p):
        I = quad(self.Gauss, p-self.strip_width/2, p+self.strip_width/2, args=(self.s,self.u))
        return I[0]
    
    def Charge_StripPos(self):
        p=self.first_strip_pos
        Charges_strips=[]
        Positions_strips=[]
        i=0
        
        while p<(self.last_strip_pos):
            Charges_strips.append(self.IntegraStrip(p))
            Positions_strips.append(p)
            p=p+self.pitch
            i=i+1
        return Charges_strips,Positions_strips
    
    def Monte_Carlo(self):
        self.N=1000
        #charges
        List_TotalCharge=np.zeros(self.N)
        List_TotalCharge_Cluster=np.zeros(self.N)
        #recontrução de posição
        E_WAC_List=np.zeros(self.N)
        E_WACC_List=np.zeros(self.N)
        E_WAlnC_List=np.zeros(self.N)
        
        ns=[]
        pN=self.N//1
        
        v=pN-1
        while v<(self.N):
            ns.append(v)
            v=v+pN
                           
        Charge_strip,Position_strip=self.Charge_StripPos()
        
        for i in range(self.N):
            Noise=np.random.normal(0,self.sig,len(Charge_strip))
            
            Charge_strip_noise=Charge_strip+Noise
            Sum_Charge=sum(Charge_strip_noise)
            List_TotalCharge[i]=Sum_Charge
            cont=0

            Charge_Cluster,Position_Cluster=Cluster(self.seed,Charge_strip_noise,Position_strip).Find_Cluster()
            
            '''
            Here we plot the clusterization of the first event, for each electron width.
            It is done to have a graph of control, to see if everything works as it is
            supposed to work.
            '''
            
            if i==0:
                plt.plot(Position_strip,Charge_strip_noise,"o", markersize=4,label="Charge collected by strips")
                plt.plot(Position_Cluster,Charge_Cluster,"*",markersize=7,label="Cluster")
                plt.legend(loc='upper right')
                plt.title('Charge cloud centered on '+r'$\mu$='+"%2.2f"%(self.u)+' mm position, '+r'$\sigma$ '+str(self.s)+" mm")#+ r' $\mu$' +"%2.2f"%(self.u)+" mm e ")
                plt.xlabel("Position of the strip center (mm)")
                plt.ylabel("Charge collected by each strip (mm)")
                plt.grid()
                plt.ylim(-0.05,0.35)
                plt.show()
                
            somaCharge_Cluster=sum(Charge_Cluster)
            List_TotalCharge_Cluster[i]=somaCharge_Cluster
                    
            for k in Noise:
                cont=cont+1

            Charge_Cluster=np.asarray(Charge_Cluster)
            Position_Cluster=np.asarray(Position_Cluster)
            
            E_WAC_List[i]=Errors(self.u, Charge_Cluster,Position_Cluster).E_WAC()
            E_WACC_List[i]=Errors(self.u, Charge_Cluster,Position_Cluster).E_WACC()
            E_WAlnC_List[i]=Errors(self.u, Charge_Cluster,Position_Cluster).E_WAlnC()

        Charge_mean=np.mean(List_TotalCharge_Cluster)
        sCharge=np.std(List_TotalCharge_Cluster,ddof=1)/np.sqrt(self.N)
        
        '''
        Errors considering the three different methods
        '''
        
        E_MPmean_C=np.mean(E_WAC_List)
        E_MPmean_CC=np.mean(E_WACC_List)
        E_MPmean_lnC=np.mean(E_WAlnC_List)
        
        E_MPs_C=np.std(E_WAC_List,ddof=1)/np.sqrt(self.N)
        E_MPs_CC=np.std(E_WACC_List,ddof=1)/np.sqrt(self.N)   
        E_MPs_lnC=np.std(E_WAlnC_List,ddof=1)/np.sqrt(self.N)
        
        return Charge_mean,sCharge,E_MPmean_C,E_MPs_C,E_MPmean_CC,E_MPs_CC,E_MPmean_lnC,E_MPs_lnC
