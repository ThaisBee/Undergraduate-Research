from scipy.integrate import quad
import numpy as np
import math


class Strips_Integrator:
    def __init__(self, strip_width: float, pitch: float, sig_noise: float, Position_strips: list, N: int=1000):
        self.strip_width = strip_width
        self.pitch = pitch
        self.sig = sig_noise
        self.Position_strips = Position_strips
        self.N = N

    def Normal(self, x: float, s: float, u: float) -> float:
        return (1 / (s * np.sqrt(2 * math.pi))) * math.exp(
            -((x - u) ** 2) / (2 * s**2)
        )

    def IntegraStrip(self, p: float, s: float, u: float) -> float:
        I = quad(
            self.Normal, p - self.strip_width / 2, p + self.strip_width / 2, args=(s, u)
        )
        return I[0]

    def Charge_Strip(self, s: float, u: float) -> float:
        Charges_strips = [self.IntegraStrip(p, s, u) for p in self.Position_strips]
        return Charges_strips

    def IntegraBins_Strip(self, p: float, x_coordinates: list, charge_fractions: list) -> float:
        strip_begin = p - self.strip_width / 2
        strip_end = p + self.strip_width / 2
        last_index_x = len(x_coordinates) - 1
        """
        The charge signal was generated in a region that does not correspond to all the strips in the readout, 
        only part of them was considered.
         So strips that do not include the sign will be zero.
        """
        if (strip_end) < x_coordinates[0]:
            return 0

        if (strip_begin) > x_coordinates[last_index_x]:
            return 0
        """ 
        From here we find the p-centered strip with charge collected from the generated electron cloud 
        """
        i = 0
        I = 0
        while i < (last_index_x) and x_coordinates[i] < (strip_end):

            if x_coordinates[i] > (strip_begin) and x_coordinates[i] < (strip_end):
                I = I + charge_fractions[i] * (x_coordinates[i + 1] - x_coordinates[i])
            i = i + 1
        return I

    def ChargeBins_Strip(self, x_coordinates: list, charge_fractions: list) -> list:
        Charges_strips = [
            self.IntegraBins_Strip(p, x_coordinates, charge_fractions)
            for p in self.Position_strips
        ]
        return Charges_strips
