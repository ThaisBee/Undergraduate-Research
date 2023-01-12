import matplotlib.pyplot as plt


def plot_cluster(
    first_evt_cluster, u, STD_DEVIATION_OF_ELECTRON_CLOUD_FIT, SEED, THRESHOLD
):
    fig, ax = plt.subplots()
    plt.axhline(
        SEED, label="Seed", linestyle=(0, (5, 10)), color="green", linewidth=1.0
    )
    plt.axhline(THRESHOLD, label="Threshold", color="red", linewidth=1.0)
    # print(first_evt_cluster['Position_strip'])
    # print(first_evt_cluster['Charge_strip_noise'])
    ax.plot(
        first_evt_cluster["Position_strip"],
        first_evt_cluster["Charge_strip_noise"],
        "o",
        markersize=4,
        label="Charge collected by strips",
    )
    ax.plot(
        first_evt_cluster["Position_Cluster"],
        first_evt_cluster["Charge_Cluster"],
        "*",
        markersize=7,
        label="Cluster",
    )
    fig.suptitle("Cluster - Event 1/1000")
    ax.set_xlabel("Position of the strip center (mm)")
    ax.set_ylabel("Charge collected by each strip (charge fraction)")

    ax.plot(
        first_evt_cluster["Position_strip"],
        first_evt_cluster["Charge_strip_noise"],
        "o",
        markersize=4,
        label="Carga coletada pelas strips",
    )
    ax.plot(
        first_evt_cluster["Position_Cluster"],
        first_evt_cluster["Charge_Cluster"],
        "*",
        markersize=7,
        label="Cluster (nuvem de elétrons)",
    )
    fig.suptitle("Clusterização")
    ax.set_xlabel("Posição do centro da fita (mm)")
    ax.set_ylabel("Carga coletada por cada fita (mm)")

    ax.grid()
    ax.set_ylim(-0.05, 0.35)

    textstr = "\n".join((r"$\mu=%.2f$" % (u,) + "mm", r"$\sigma=%.2f$" % (s,) + "mm"))
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(
        0.08,
        0.95,
        textstr,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=props,
    )
    ax.legend(loc="upper right")

    plt.show()
