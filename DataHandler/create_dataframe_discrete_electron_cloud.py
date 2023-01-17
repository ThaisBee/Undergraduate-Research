from models.strips_integrator import Strips_Integrator
from models.monte_carlo import Monte_Carlo
import pandas as pd

def create_a_dataframe(
    STRIP_WIDTH,
    PITCH,
    STD_DEVIATION_OF_ELECTRON_CLOUD_FIT,
    STD_DEVIATION_OF_THE_NOISE,
    NUMBER_OF_ELECTRON_CLOUDS,
    Position_strips,
    x_coordinate,
    charge_fraction,
    electron_cloud_centers,
    SEED,
    THRESHOLD,
):
    """For the case when PITCH=0.39mm, the electron cloud center is going to be reconstructed approximately -0.8mm<x<0.8mm
    Also the electron cloud center is going to be reconstructed with steps 10 times smaller than the strip PITCH"""

    Data = []
    for u in electron_cloud_centers:
        x_new_coordinate = x_coordinate + u
        Charge_strips = Strips_Integrator(
            STRIP_WIDTH, PITCH, STD_DEVIATION_OF_THE_NOISE, Position_strips
        ).ChargeBins_Strip(x_new_coordinate, charge_fraction)
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
                STD_DEVIATION_OF_ELECTRON_CLOUD_FIT,
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
        # plot_cluster(first_evt_cluster,u,s,SEED,THRESHOLD)

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
