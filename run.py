"""
Author: Thais Silva Abelha
Email: thais.silva.abelha@gmail.com

This work was done during a undergraduated research "Comparing different methods of position
reconstruction considering 1D readout of GEM detectors"
"""

import numpy as np
import pandas as pd

from constants import (
    STRIP_WIDTH,
    PITCH,
    FIRST_STRIP_POSITION,
    LAST_STRIP_POSITION,
    FIRST_CLOUD_POSITION,
    LAST_CLOUD_POSITION,
    SEED,
    THRESHOLD,
    STD_DEVIATION_OF_THE_NOISE,
    NUMBER_OF_ELECTRON_CLOUDS,
)
from plots_analytical_electron_cloud.plot_for_charges_and_errors import (
    plot_for_carges_and_errors,
)
from plots_analytical_electron_cloud.plot_rmse import rmse, rmse_dataframe, plot_rmse
from DataHandler.create_dataframe_analytical_electron_cloud import create_a_dataframe

#%%

"""
linear = weighted average considering the charge as a weight
quadratic = weighted average considering the charge*charge as a weighet
logarithmic = weighted average considering the ln(charge) as a weight
Q = charge
s=standard deviation or sigma σ of the gaussian
u=mean of the gaussian or electron cloud center

prefix "s" stands for standard deviation, which is represented by σ

rmse = Root Mean Square 
"""

#%%
"""For this study we have modeled a 1 cm readout, from -5 mm to 5 mm  where the coordinate 0 is in the middle of the
strips pattern. The clouds were generated in the range −0.8 mm to 0.8 mm, so the cloud didn’t
suffer edge effect. The clouds of electron colected by conductive strips were modeled as normalized Gaussian functions centered
in μ, which we call u, with standard deviation σ, which we call s.
The strips have 0.2 mm width and ~0.39 mm PITCH, which means the distance between two consecutive strips."""


strip_centers = np.arange(FIRST_STRIP_POSITION, LAST_STRIP_POSITION, PITCH)
electron_cloud_centers = np.arange(
    FIRST_CLOUD_POSITION, LAST_CLOUD_POSITION, PITCH / 10
)

###############################################################################
# list of standard deviations for the electron clouds modeled as gaussians
std_deviation_of_electron_clouds = [0.2, 0.25, 0.3, 0.35, 0.4]

df = create_a_dataframe(
    STRIP_WIDTH,
    PITCH,
    std_deviation_of_electron_clouds,
    STD_DEVIATION_OF_THE_NOISE,
    NUMBER_OF_ELECTRON_CLOUDS,
    strip_centers,
    electron_cloud_centers,
    SEED,
    THRESHOLD,
)

#%%############################################################################
# Graphs that depend on the electron cloud position relative to the readout strips
###############################################################################
markers = ["s", "^", "P", "*", "x"]
markersizes = [4, 4, 4, 14 / 3, 4]
ymax = 0.6
ymin = 0.45
Title = {"Q": "Collected charge"}
Title_pt = {"Q": "Carga coletada pelas fitas"}
Ylabel = {"Q": "Fraction of collected charge"}
Ylabel_pt = {"Q": "Fração de carga coletada"}
plot_for_carges_and_errors(
    "Q",
    std_deviation_of_electron_clouds,
    df,
    markers,
    markersizes,
    ymin,
    ymax,
    Title["Q"],
    Ylabel["Q"],
)

###############################################################################
ymin = -0.045
ymax = 0.045
# ymin=-0.055
# ymax=0.3
Title = {
    "E_linear": "Position reconstruction using a linear weight",
    "E_quadratic": "Position reconstruction using a squared weight",
    "E_logarithmic": "Position reconstruction using a logarithmic weight",
}
Title_pt = {
    "E_linear": "Reconstrução da posição usando peso linear",
    "E_quadratic": "Reconstrução da posição usando peso quadrático",
    "E_logarithmic": "Reconstrução da posição usando peso logaritmico",
}
Ylabel = {"E": "Error (mm)"}
Ylabel_pt = {"E": "Erro (mm)"}
plot_for_carges_and_errors(
    "E_linear",
    std_deviation_of_electron_clouds,
    df,
    markers,
    markersizes,
    ymin,
    ymax,
    Title["E_linear"],
    Ylabel["E"],
)
plot_for_carges_and_errors(
    "E_quadratic",
    std_deviation_of_electron_clouds,
    df,
    markers,
    markersizes,
    ymin,
    ymax,
    Title["E_quadratic"],
    Ylabel["E"],
)
plot_for_carges_and_errors(
    "E_logarithmic",
    std_deviation_of_electron_clouds,
    df,
    markers,
    markersizes,
    ymin,
    ymax,
    Title["E_logarithmic"],
    Ylabel["E"],
)

#%%############################################################################
# rmse Graph that depend on the electron cloud width (standard deviation of the gaussian)
###############################################################################
rmse_df = rmse_dataframe(std_deviation_of_electron_clouds, df)

Labels = {
    "Title": "Comparing methods for position reconstruction",
    "ylabel": "rmse of the Error (mm)",
}

Labels_pt = {
    "Title": "Comparing methods for position reconstruction",
    "ylabel": "Erro do rmse (mm)",
}
markers = ["s", "^", "P", "*", "x"]
markersizes = [4, 4, 4, 14 / 3, 4]

plot_rmse(rmse_df, Labels, markers, markersizes)
