import matplotlib as plt


def plot_hist2d(x, y, z, nbins):
    fig, ax = plt.subplots()
    counts, xedges, yedges, im = ax.hist2d(x, y, bins=nbins, weights=z, cmap=cm.jet)
    fig.colorbar(im, ax=ax, label="Elétrons")
    ax.set_ylim(-1, 0.85)
    ax.set_xlim(-1, 0.85)
    ax.set_xlabel("Coordeadas em x (mm)")
    ax.set_ylabel("Coordenadas em y (mm)")
    plt.title(
        "Nuvem de elétrons depois da multiplicação \n Fóton de 8 KeV "
        + r"$Ar/CO_2$ (70/30)"
    )
    plt.plot(0, 0, "*", color="black")
    plt.show()
