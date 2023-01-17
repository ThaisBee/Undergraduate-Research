import numpy as np

"""
WAC = weighted average considering the charge as a weight
WACC = weighted average considering the charge*charge as a weighet
WAlnC = weighted average considering the ln(charge) as a weight

These weights were chosen in order to evaluate the impact of lower and higher charges 
widths in the charge and position reconstructions. The logarithmic weight gives more 
importance to lower charge values, the squared weight to the higher values, and the
linear weight is between them. 
"""


class Position_Reconstruction:
    def __init__(self, Charges: list, Positions: list):
        self.StripCharges = Charges
        self.StripPositions = Positions

    def linear(self) -> float:
        WA = sum(self.StripCharges * self.StripPositions) / sum(self.StripCharges)
        return WA

    def quadratic(self) -> float:
        WA = sum((self.StripCharges**2) * self.StripPositions) / sum(
            self.StripCharges**2
        )
        return WA

    def logarithmic(self) -> float:
        charges = [np.log(k + 1) for k in self.StripCharges]
        WA = sum(charges * self.StripPositions) / sum(charges)
        return WA
