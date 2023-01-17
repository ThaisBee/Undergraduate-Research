from models.strips_integrator import Strips_Integrator
from models.monte_carlo import Monte_Carlo
import pandas as pd

def create_a_dataframe(
    STRIP_WIDTH,
    PITCH,
    std_deviation_of_electron_clouds,
    STD_DEVIATION_OF_THE_NOISE ,
    NUMBER_OF_ELECTRON_CLOUDS,
    Position_strips,
    electron_cloud_centers,
    SEED,
    THRESHOLD,
):
    """For the case when    PITCH=0.39mm, the electron cloud center is going to be reconstructed approximately -0.8mm<x<0.8mm
    Also the electron cloud center is going to be reconstructed with steps 10 times smaller than the strip PITCH"""

    # You can comment "Plot_Cluster()" It is just a Graph of control then it is going to run faster
    # Plot_Cluster() is used just to see if verything is fine with the algorithm

    Data = []
    for s in std_deviation_of_electron_clouds:

        for u in electron_cloud_centers:
            Charge_strips = Strips_Integrator(
                STRIP_WIDTH, PITCH, STD_DEVIATION_OF_THE_NOISE , Position_strips
            ).Charge_Strip(s, u)
            monte_carlo_obj = Monte_Carlo(
                STD_DEVIATION_OF_THE_NOISE , NUMBER_OF_ELECTRON_CLOUDS, u, Charge_strips, Position_strips, SEED, THRESHOLD
            )
            monte_carlo_obj.compute()
            Charge, sCharge = monte_carlo_obj.charge
            E_linear, sE_linear = monte_carlo_obj.linear_weight_error
            E_quadratic, sE_quadratic = monte_carlo_obj.quadratic_weight_error
            E_logarithmic, sE_logarithmic = monte_carlo_obj.logarithmic_weight_error
            first_evt_cluster = monte_carlo_obj.first_evt_cluster

            Data.append(
                [
                    s,
                    u,
                    Charge,
                    sCharge,
                    E_linear,
                    sE_linear,
                    E_quadratic,
                    sE_quadratic,
                    E_logarithmic,
                    sE_logarithmic,
                ]
            )
            # Plot_Cluster(first_evt_cluster,u,s,Labels)
            L = [
                s,
                u,
                first_evt_cluster["Position_strip"],
                first_evt_cluster["Charge_strip_noise"],
            ]

    df = pd.DataFrame(
        Data,
        columns=[
            "s",
            "x",
            "Q",
            "sQ",
            "E_linear",
            "sE_linear",
            "E_quadratic",
            "sE_quadratic",
            "E_logarithmic",
            "sE_logarithmic",
        ],
    )
    return df