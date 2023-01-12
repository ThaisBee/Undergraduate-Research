import matplotlib.pyplot as plt


def plot_errors(df, label):
    plt.errorbar(
        df["x"],
        df["E_linear"],
        marker="o",
        markersize=3,
        linestyle="",
        label=label["E_linear"],
        xerr=0,
        yerr=df["sE_linear"],
    )
    plt.errorbar(
        df["x"],
        df["E_quadratic"],
        marker="o",
        markersize=3,
        linestyle="",
        label=label["E_quadratic"],
        xerr=0,
        yerr=df["sE_quadratic"],
    )
    plt.errorbar(
        df["x"],
        df["E_logarithmic"],
        marker="o",
        markersize=3,
        linestyle="",
        label=label["E_logarithmic"],
        xerr=0,
        yerr=df["sE_logarithmic"],
    )
    plt.title(label["Title"])
    plt.legend(loc="upper right")
    plt.ylabel(label["ylabel"])
    plt.xlabel(r"$\mu$ (mm)")
    plt.grid()
    plt.show()
