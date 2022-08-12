import numpy as np
'''
A clustering algorithm is used to identify the strips corresponding to charges that came from the
electron cloud. The algorithm first looks for a strip whose signal is above a defined threshold,
0.08 AU in this work, to find a strip that definitely contains a signal. Then, starting from this
seed strip, the algorithm scans in both direction to construct a cluster. The scan stops when it
finds a strip that collected a negative charge, which means just noise. So we are going to find
the positions I1 and I2. Which I1 represents the first strip of the cluster and I2 the last one
'''
class Cluster:
    def __init__(self,seed,Charges,Positions,zero=0):
        self.seed=seed
        self.Q=Charges
        self.P=Positions
        self.zero=zero
        
    def find_seed(self):
        i=0
        Qi=self.Q[i]
         
        while Qi<self.seed:
            i=i+1
            Qi=self.Q[i]
        seed_index=i
        
        return seed_index

    def scans_left(self,I1):
        Qleft=self.Q[I1]
        while Qleft>self.zero and I1!=0:
                I1=I1-1
                Qleft=self.Q[I1]
               
        if Qleft>self.zero:
            return I1
        else:
            I1=I1+1 
            return I1
    
    def scans_right(self,I2): 
        Qright=self.Q[I2] 
        last_index=len(self.Q)-1
        
        while Qright>self.zero and I2<last_index:
            I2=I2+1
            Qright=self.Q[I2]
        
        if Qright>self.zero:
            return I2
        
        else: 
            I2=I2-1
            return I2
    
    def Find_Cluster(self):
        seed_index=self.find_seed()
        I1=seed_index
        I2=seed_index
        last_index=len(self.Q)-1

        if seed_index>0 and seed_index<last_index:
            I1=self.scans_left(seed_index-1)
            I2=self.scans_right(seed_index+1)
            
        if seed_index==0:
            I1=0
            I2=self.scans_right(seed_index+1)
  
        if seed_index==last_index:
            I2=last_index 
            I1=self.scans_left(I2)
            
        if I1!=I2:
            I2=I2+1
            P_Cluster=self.P[I1:I2]
            Q_Cluster=self.Q[I1:I2]
            
        else:
            P_Cluster=np.zeros(1)
            Q_Cluster=np.zeros(1)
            P_Cluster[0]=self.P[I1]
            Q_Cluster[0]=self.Q[I1]
            
        return Q_Cluster,P_Cluster
   