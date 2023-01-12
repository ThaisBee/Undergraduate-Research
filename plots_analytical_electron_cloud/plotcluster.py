import matplotlib.pyplot as plt


def Plot_Cluster(first_evt_cluster, u, s, Labels, SEED, THRESHOLD):
    fig, ax = plt.subplots()
    plt.axhline(
        SEED, label="Seed", linestyle=(0, (5, 10)), color="green", linewidth=1.0
    )
    plt.axhline(THRESHOLD, label="Threshold", color="red", linewidth=1.0)
    ax.plot(
        first_evt_cluster["Position_strip"],
        first_evt_cluster["Charge_strip_noise"],
        "o",
        markersize=4,
        label=(Labels["Label_Charges"]),
    )
    ax.plot(
        first_evt_cluster["Position_Cluster"],
        first_evt_cluster["Charge_Cluster"],
        "*",
        markersize=7,
        label=Labels["Label_Cluster"],
    )
    fig.suptitle(Labels["Title"])
    ax.set_xlabel(Labels["xlabel"])
    ax.set_ylabel(Labels["ylabel"])
    ax.legend(loc="upper right")
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
    plt.show()
