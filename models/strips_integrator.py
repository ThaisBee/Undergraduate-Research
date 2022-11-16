from scipy.integrate import quad
import numpy as np
import math

class Strips_Integrator:

    def __init__(self,strip_width,pitch,sig_noise,Position_strips,N=1000):
        self.strip_width=strip_width
        self.pitch=pitch
        self.sig=sig_noise
        self.Position_strips=Position_strips
        self.N=N
        
    def Normal(self,x,s,u):
        return (1/(s*np.sqrt(2*math.pi)))*math.exp( -((x-u)**2) / (2*s**2)  )
    
    def IntegraStrip(self,p,s,u):
        I = quad(self.Normal, p-self.strip_width/2, p+self.strip_width/2, args=(s,u))
        return I[0]

    def Charge_Strip(self,s,u):
        Charges_strips=[self.IntegraStrip(p,s,u) for p in self.Position_strips]
        return Charges_strips
    
    def IntegraBins_Strip(self,p,x_coordinate,charge_fraction):
        strip_begin=p-self.strip_width/2
        strip_end=p+self.strip_width/2
        last_index_x=len(x_coordinate)-1
        '''
        The charge signal was generated in a region that does not correspond to all the strips in the readout, 
        only part of them was considered.
         So strips that do not include the sign will be zero.
        '''
        if (strip_end)<x_coordinate[0]:
            return 0
        
        if (strip_begin)>x_coordinate[last_index_x]:
            return 0
        ''' 
        From here we find the p-centered strip with charge collected from the generated electron cloud 
        '''
        i=0
        I=0
        while i<(last_index_x) and x_coordinate[i] < (strip_end): 
            
            if x_coordinate[i]>(strip_begin) and x_coordinate[i]<(strip_end):
                I=I+charge_fraction[i]*(x_coordinate[i+1]-x_coordinate[i])  
            i=i+1
        return I
        
    def ChargeBins_Strip(self,x_coordinate,charge_fraction):
        Charges_strips=[self.IntegraBins_Strip(p,x_coordinate,charge_fraction) for p in self.Position_strips]
        return Charges_strips
 

