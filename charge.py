#Python library
import numpy as np
import math
import matplotlib.pyplot as plt
#My library
from scipy.integrate import quad
from cluster import Cluster
from errors import Errors

class Charge:

    def __init__(self,first_strip_pos,last_strip_pos,strip_width,pitch,sigma,u,sig_noise,N=1000,seed=0.04):
        self.first_strip_pos=first_strip_pos
        self.last_strip_pos=last_strip_pos
        self.strip_width=strip_width
        self.pitch=pitch
        self.s=sigma
        self.u=u
        self.sig=sig_noise#ruido
        self.seed=seed
        self.N=N
        
    def Gauss(self,x,s,u):
        return (1/(s*np.sqrt(2*math.pi)))*math.exp( -((x-u)**2) / (2*s**2)  )
    
    def IntegraStrip(self,p):#strip com ruido
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
        E_MPC_List=np.zeros(self.N)
        E_MPCC_List=np.zeros(self.N)
        E_MPlnC_List=np.zeros(self.N)
        
        ns=[]
        pN=self.N//1
        
        v=pN-1
        while v<(self.N):
            ns.append(v)
            v=v+pN
                           
        A,P=self.Charge_StripPos()
        
        for i in range(self.N):
            Noise=np.random.normal(0,self.sig,len(A))
            
            Q=A+Noise
            somaQ=sum(Q)
            List_TotalCharge[i]=somaQ
            cont=0

            Q_Cluster,P_Cluster=Cluster(self.seed,Q,P).Find_Cluster()
            
            if i==0:
                plt.plot(P,Q,"o", markersize=4,label="Charge collected by strips")
                plt.plot(P_Cluster,Q_Cluster,"*",markersize=7,label="Cluster")
                plt.legend(loc='upper right')
                plt.title('Charge cloud centered on '+r'$\mu$='+"%2.2f"%(self.u)+' mm position, '+r'$\sigma$ '+str(self.s)+" mm")#+ r' $\mu$' +"%2.2f"%(self.u)+" mm e ")
                plt.xlabel("Position of the strip center (mm)")
                plt.ylabel("Charge collected by each strip (mm)")
                plt.grid()
                plt.ylim(-0.05,0.35)
                plt.show()
                
            somaQ_Cluster=sum(Q_Cluster)
            List_TotalCharge_Cluster[i]=somaQ_Cluster
                    
            for k in Noise:
                cont=cont+1

            Q_Cluster=np.asarray(Q_Cluster)
            P_Cluster=np.asarray(P_Cluster)
            
            E_MPC_List[i]=Errors(self.u, Q_Cluster,P_Cluster).E_MPC()
            E_MPCC_List[i]=Errors(self.u, Q_Cluster,P_Cluster).E_MPCC()
            E_MPlnC_List[i]=Errors(self.u, Q_Cluster,P_Cluster).E_MPlnC()

        Qmean=np.mean(List_TotalCharge_Cluster)
        Qs=np.std(List_TotalCharge_Cluster,ddof=1)/np.sqrt(self.N)
        
        E_MPmean_C=np.mean(E_MPC_List)
        E_MPmean_CC=np.mean(E_MPCC_List)
        E_MPmean_lnC=np.mean(E_MPlnC_List)
        
        E_MPs_C=np.std(E_MPC_List,ddof=1)/np.sqrt(self.N)
        E_MPs_CC=np.std(E_MPCC_List,ddof=1)/np.sqrt(self.N)   
        E_MPs_lnC=np.std(E_MPlnC_List,ddof=1)/np.sqrt(self.N)
        
        return Qmean,Qs,E_MPmean_C,E_MPs_C,E_MPmean_CC,E_MPs_CC,E_MPmean_lnC,E_MPs_lnC
