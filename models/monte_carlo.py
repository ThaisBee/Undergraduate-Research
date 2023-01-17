import numpy as np
import pandas as pd

# my libraries
from models.cluster import Cluster
from models.position_reconstruction import Position_Reconstruction


class Monte_Carlo:
    def __init__(
        self,
        sig_noise: float,
        N: int,
        u: int,
        Charge_strip: list,
        Position_strip: list,
        seed: float,
        threshold: float,
    ):
        self.sig_noise = sig_noise
        self.N = N
        self.u = u
        self.Charge_strip = Charge_strip
        self.Position_strip = Position_strip
        self.seed = seed
        self.threshold = threshold
        self.df = None
        self.first_evt = {}

    def compute(self) -> None:

        Data = []
        for i in range(self.N):
            Noise = np.random.normal(0, self.sig_noise, len(self.Charge_strip))
            Charge_strip_noise = self.Charge_strip + Noise
            Charge_Cluster, Position_Cluster = Cluster(
                self.seed, self.threshold, Charge_strip_noise, self.Position_strip
            ).Find_Cluster()

            Data.append(
                [
                    sum(Charge_Cluster),
                    Position_Reconstruction(Charge_Cluster, Position_Cluster).linear()
                    - self.u,
                    Position_Reconstruction(
                        Charge_Cluster, Position_Cluster
                    ).quadratic()
                    - self.u,
                    Position_Reconstruction(
                        Charge_Cluster, Position_Cluster
                    ).logarithmic()
                    - self.u,
                ]
            )

            if i == 0:
                self.first_evt = {
                    "Position_strip": self.Position_strip,
                    "Charge_strip_noise": Charge_strip_noise,
                    "Position_Cluster": Position_Cluster,
                    "Charge_Cluster": Charge_Cluster,
                }
        self.df = pd.DataFrame(
            Data,
            columns=[
                "Q",
                "linear_weight_error",
                "quadratic_weight_error",
                "logarithmic_weight_error",
            ],
        )


    @property
    def charge(self) -> tuple:
        charge_mean = np.mean(self.df["Q"])
        std_Charge = np.std(self.df["Q"], ddof=1) / np.sqrt(self.N)
        return charge_mean, std_Charge

    @property
    def linear_weight_error(self) -> tuple:
        linear_weight_error = np.mean(self.df["linear_weight_error"])
        std_linear_weight_error = np.std(
            self.df["linear_weight_error"], ddof=1
        ) / np.sqrt(self.N)
        return linear_weight_error, std_linear_weight_error

    @property
    def quadratic_weight_error(self) -> tuple:
        quadratic_weight_error = np.mean(self.df["quadratic_weight_error"])
        std_quadratic_weight_error = np.std(
            self.df["quadratic_weight_error"], ddof=1
        ) / np.sqrt(self.N)
        return quadratic_weight_error, std_quadratic_weight_error

    @property
    def logarithmic_weight_error(self) -> tuple:
        logarithmic_weight_error = np.mean(self.df["logarithmic_weight_error"])
        std_logarithmic_weight_error = np.std(
            self.df["logarithmic_weight_error"], ddof=1
        ) / np.sqrt(self.N)
        return logarithmic_weight_error, std_logarithmic_weight_error

    @property
    def first_evt_cluster(self) -> dict:
        return self.first_evt
