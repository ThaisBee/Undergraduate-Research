import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
from matplotlib.ticker import MultipleLocator

from data_handler.significativos import dois_significativos

import numpy as np


def Func_Normal(x, x0, sigma):
    return (1 / (sigma * np.sqrt(2 * math.pi))) * np.exp(
        -((x - x0) ** 2) / (2.0 * sigma**2)
    )


def plot_normal(x_coordinate, charge_fraction, Position_strips, popt, pcov):
    inc = np.sqrt(np.diagonal(pcov))

    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0, height_ratios=[2, 1])
    axs = gs.subplots(sharex=True)
    title = "Distribuição das cargas da nuvem"
    fig.suptitle(title)

    axs[0].plot(
        x_coordinate,
        charge_fraction,
        "b+:",
        label="Nuvem de elétrons, Garfield ++",
        color="blue",
    )
    axs[0].plot(
        x_coordinate,
        Func_Normal(x_coordinate, *popt),
        "r-",
        label="fit:"
        + r" $\mu$="
        + dois_significativos(popt[0], inc[0]).forma2()
        + "mm, "
        + r"$\sigma$="
        + dois_significativos(popt[1], inc[1]).forma2()
        + "mm",
    )
    axs[0].set_ylabel("Carga (U.A.)")
    axs[0].set_ylim(0, 3)
    axs[0].legend(loc="upper right")
    axs[0].xaxis.set_major_locator(MultipleLocator(0.4))
    axs[0].xaxis.set_minor_locator(MultipleLocator(0.2))

    Residuals = charge_fraction - Func_Normal(x_coordinate, *popt)
    axs[1].plot(
        x_coordinate, Residuals, linestyle="None", marker="*", markersize=1, color="b"
    )
    axs[1].set_xlabel("Coordenadas em x (mm)")
    axs[1].set_ylabel("Resíduos (U.A.)")
    axs[1].set_xlim(-1.2, 1.2)
    axs[1].grid(axis="y")

    rectangles = [k for k in Position_strips if k > -1.2 and k < 1.2]
    trans0 = transforms.blended_transform_factory(axs[0].transData, axs[0].transAxes)
    trans1 = transforms.blended_transform_factory(axs[1].transData, axs[1].transAxes)

    for i in rectangles:
        rect0 = mpatches.Rectangle(
            (i, 0), width=0.2, height=1, transform=trans0, color="red", alpha=0.5
        )
        rect1 = mpatches.Rectangle(
            (i, 0), width=0.2, height=1, transform=trans1, color="red", alpha=0.5
        )

        axs[0].add_patch(rect0)
        axs[1].add_patch(rect1)

    plt.show()
    return
