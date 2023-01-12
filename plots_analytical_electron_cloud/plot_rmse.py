import matplotlib.pyplot as plt
from matplotlib import container
import numpy as np
import pandas as pd


def rmse(v):
    return np.sqrt((v * v).sum() / len(v))


def rmse_dataframe(std_deviation_of_electron_clouds, df):
    Data = []
    for s in std_deviation_of_electron_clouds:
        Data.append(
            [
                s,
                rmse(df.loc[df["s"].isin([s]), "E_linear"]),
                rmse(df.loc[df["s"].isin([s]), "E_quadratic"]),
                rmse(df.loc[df["s"].isin([s]), "E_logarithmic"]),
            ]
        )

    rmse_df = pd.DataFrame(
        Data,
        columns=[
            "s",
            "E_linear",
            "E_quadratic",
            "E_logarithmic",
        ],
    )
    return rmse_df


def plot_rmse(rmse_df, Labels, markers, markersizes):
    plt.title(Labels["Title"])
    plt.errorbar(
        rmse_df["s"],
        rmse_df["E_linear"],
        marker=markers[0],
        markersize=markersizes[0] * 1.5,
        linestyle="",
        label="Linear weight",
        xerr=0,
    )
    plt.errorbar(
        rmse_df["s"],
        rmse_df["E_quadratic"],
        marker=markers[1],
        markersize=markersizes[1] * 1.5,
        linestyle="",
        label="Squared weight",
        xerr=0,
    )
    plt.errorbar(
        rmse_df["s"],
        rmse_df["E_logarithmic"],
        marker=markers[2],
        markersize=markersizes[2] * 1.5,
        linestyle="",
        label="Logarithmic weight",
        xerr=0,
    )
    plt.ylabel(Labels["ylabel"])

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    handles = [
        h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles
    ]
    ax.legend(handles, labels, loc="upper right")
    plt.xlabel(r"$\sigma$ (mm)")
    plt.grid(True)
    plt.show()
