import numpy as np 

'''
WAC = weighted average considering the charge as a weight
WACC = weighted average considering the charge*charge as a weighet
WAlnC = weighted average considering the ln(charge) as a weight

These weights were chosen in order to evaluate the impact of lower and higher charges
 in the charge and position reconstructions. The logarithmic weight gives more 
 importance to lower charge values, the squared weight to the higher values, and the
 linear weight is between them. 
'''
class Errors:
    def __init__(self,u,Charges,Positions):
        self.u=u #Gaussian Center
        self.StripCharges=Charges
        self.StripPositions=Positions
        
    def E_WAC(self):
        WA=sum(self.StripCharges*self.StripPositions)/sum(self.StripCharges)
        return WA-self.u
    
    def E_WACC(self):
        WA=sum(self.StripCharges*self.StripCharges*self.StripPositions)/sum(self.StripCharges*self.StripCharges)
        return WA-self.u
    
    def E_WAlnC(self):
        charges=np.zeros(len(self.StripCharges))
        
        for i in range(len(self.StripCharges)):
            charges[i]=np.log(self.StripCharges[i]+1)
    
        WA=sum(charges*self.StripPositions)/sum(charges)
        return WA-self.u