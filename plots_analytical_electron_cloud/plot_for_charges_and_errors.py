import matplotlib.pyplot as plt
from matplotlib import container


def plot_for_carges_and_errors(
    Name,
    std_deviation_of_electron_clouds,
    df,
    markers,
    markersizes,
    ymin,
    ymax,
    Title,
    Ylabel,
):
    for s in std_deviation_of_electron_clouds:
        t = r"$\sigma $ " + " " + str(s) + "mm"
        plt.errorbar(
            df.loc[df["s"].isin([s]), "x"],
            df.loc[df["s"].isin([s]), Name],
            marker=markers[std_deviation_of_electron_clouds.index(s)],
            markersize=markersizes[std_deviation_of_electron_clouds.index(s)],
            linestyle="",
            label=t,
            xerr=0,
            yerr=df.loc[df["s"].isin([s]), "sE_linear"],
        )

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    handles = [
        h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles
    ]

    ax.legend(handles, labels, loc="upper right")
    plt.title(Title)
    plt.ylabel(Ylabel)
    plt.ylim(ymin, ymax)
    plt.xlabel(r"$\mu$ (mm)")
    plt.grid(True)
    plt.show()
