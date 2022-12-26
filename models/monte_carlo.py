import numpy as np
import pandas as pd

# my libraries
from models.cluster import Cluster
from models.position_reconstruction import Position_Reconstruction


def Monte_Carlo(sig_noise, N, u, Charge_strip, Position_strip, seed, threshold):
    Data = []
    for i in range(N):
        Noise = np.random.normal(0, sig_noise, len(Charge_strip))
        Charge_strip_noise = Charge_strip + Noise
        Charge_Cluster, Position_Cluster = Cluster(
            seed, threshold, Charge_strip_noise, Position_strip
        ).Find_Cluster()

        Data.append(
            [
                sum(Charge_Cluster),
                Position_Reconstruction(Charge_Cluster, Position_Cluster).linear() - u,
                Position_Reconstruction(Charge_Cluster, Position_Cluster).quadratic()
                - u,
                Position_Reconstruction(Charge_Cluster, Position_Cluster).logarithmic()
                - u,
            ]
        )

        if i == 0:
            first_evt_cluster = {
                "Position_strip": Position_strip,
                "Charge_strip_noise": Charge_strip_noise,
                "Position_Cluster": Position_Cluster,
                "Charge_Cluster": Charge_Cluster,
            }
    df = pd.DataFrame(Data, columns=["Q", "E_linear", "E_quadratic", "E_logarithmic"])

    """
    Considering N events we are going to calaculate Error means and Stardard deviations of the Position_Reconstruction menas
    """

    Charge_mean = np.mean(df["Q"])
    E_linear = np.mean(df["E_linear"])
    E_quadratic = np.mean(df["E_quadratic"])
    E_logarithmic = np.mean(df["E_logarithmic"])

    sCharge = np.std(df["Q"], ddof=1) / np.sqrt(N)
    sE_linear = np.std(df["E_linear"], ddof=1) / np.sqrt(N)
    sE_quadratic = np.std(df["E_quadratic"], ddof=1) / np.sqrt(N)
    sE_logarithmic = np.std(df["E_logarithmic"], ddof=1) / np.sqrt(N)

    return (
        Charge_mean,
        sCharge,
        E_linear,
        sE_linear,
        E_quadratic,
        sE_quadratic,
        E_logarithmic,
        sE_logarithmic,
        first_evt_cluster,
    )
