# UNIVERSIDADE FEDERAL DA PARAÍBA
# CENTRO DE TECNOLOGIA
# DEPARTAMENTO DE ENGENHARIA MECÂNICA
#
# DISCENTE THIAGO NEY EVARISTO RODRIGUES
# ORIENTADOR DR. JACQUES CÉSAR DOS SANTOS
#
# TRABALHO DE CONCLUSÃO DE CURSO
#
# ANÁLISE NUMÉRICA DE CONDUÇÃO TRANSIENTE COM TERMO FONTE VARIÁVEL EM
# VARETAS COMBUSTÍVEIS DE REATORES NUCLEARES PELO MÉTODO DAS LINHAS

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def solve_heat_pde(t, theta, r, Bi, Gast, c1):
    """Calculates the time derivative d(theta)/dt for the heat equation PDE."""
    nr = len(theta)
    dr = (r[-1] - r[0]) / (nr - 1)
    dr2 = dr**2

    # Pre-allocate derivative vectors
    dtheta_dr = np.zeros(nr)
    d2theta_dr2 = np.zeros(nr)

    # Calculate spatial derivatives for interior points
    i = np.arange(1, nr)
    dtheta_dr[i] = (theta[i] - theta[i - 1]) / dr

    i = np.arange(1, nr - 1)
    d2theta_dr2[i] = (theta[i + 1] - 2 * theta[i] + theta[i - 1]) / dr2

    # Apply boundary conditions
    dtheta_dr[0] = 0.0  # Center (symmetry)
    d2theta_dr2[0] = 2 * (theta[1] - theta[0]) / dr2  # Center (symmetry)
    d2theta_dr2[-1] = (
        2 * theta[-2] - 2 * (1 + Bi * dr) * theta[-1]
    ) / dr2  # Edge (convection)

    # Define the new source term, which varies with t
    G = Gast * (1 + c1 * t)

    # Calculate the time derivative for all points
    dtheta_dt = (1 / r) * dtheta_dr + d2theta_dr2 + G

    # Correct the equation for the central point (r=0)
    dtheta_dt[0] = 2 * d2theta_dr2[0] + G

    return dtheta_dt


if __name__ == "__main__":
    # Input Parameters
    nr = 9601  # Points in spatial grid
    nt = 100  # Points in temporal grid
    r0 = 0.001  # Beginning of the r axis (small value to avoid division by zero)
    rl = 1.0  # End of the r axis (Length L)
    t0 = 0.0  # Start time
    tl = 1.5  # End time
    Bi = 15.0  # Biot Number
    Gast = 32.4  # Heat source term constant
    c1 = 1  # New constant for the source term

    # Reference Data Curves
    curve01 = np.array(
        [
            [0.00159, 9.19405],
            [0.12440, 9.03756],
            [0.19777, 8.84194],
            [0.29665, 8.45070],
            [0.39394, 7.90297],
            [0.50718, 7.08138],
            [0.61563, 6.06416],
            [0.71451, 4.96870],
            [0.79904, 3.91236],
            [0.905901, 2.38654],
            [0.99681, 0.978091],
        ]
    )
    curve02 = np.array(
        [
            [0.00159, 11.97180],
            [0.09729, 11.81530],
            [0.20096, 11.58060],
            [0.30463, 11.07200],
            [0.39713, 10.40690],
            [0.51675, 9.23318],
            [0.60128, 8.25509],
            [0.709729, 6.76839],
            [0.805423, 5.24257],
            [0.910686, 3.16901],
            [0.995215, 1.52582],
        ]
    )
    curve03 = np.array(
        [
            [0.00000, 16.39280],
            [0.09091, 16.23630],
            [0.19936, 15.88420],
            [0.31579, 15.02350],
            [0.40829, 14.04540],
            [0.50080, 12.87170],
            [0.59330, 11.42410],
            [0.69697, 9.46792],
            [0.807018, 7.04225],
            [0.909091, 4.34272],
            [0.995215, 1.99531],
        ]
    )
    curve04 = np.array(
        [
            [0.00159, 20.97030],
            [0.10367, 20.69640],
            [0.21053, 20.18780],
            [0.29346, 19.44440],
            [0.40510, 17.99690],
            [0.50558, 16.31460],
            [0.59490, 14.51490],
            [0.701754, 11.93270],
            [0.805423, 8.99844],
            [0.896332, 6.06416],
            [0.995215, 2.58216],
        ]
    )

    # Numerical Calculation
    r = np.linspace(r0, rl, nr)
    t_eval = np.linspace(t0, tl, nt)
    theta0 = Gast * (1 - r**2) / 4 + Gast / (2 * Bi)

    # Solve the system of ODEs
    sol = solve_ivp(
        lambda t, theta: solve_heat_pde(t, theta, r, Bi, Gast, c1),
        t_span=[t0, tl],
        y0=theta0,
        method="BDF",  # Method for stiff problems
        t_eval=t_eval,
    )
    theta = sol.y.T  # Transpose solution to shape (time, space)

    # Plotting

    font_axes = {"weight": "normal", "size": 12}
    font_legend = {"size": 12}

    # Create meshgrid for 3D plots
    T_mesh, R_mesh = np.meshgrid(t_eval, r)

    # Figure 1: 3D Surface Plot
    fig1 = plt.figure(figsize=(10, 7), dpi=600)
    ax1 = fig1.add_subplot(111, projection="3d")
    ax1.plot_surface(T_mesh, R_mesh, theta.T, cmap="viridis", edgecolor="none")
    ax1.set_xlabel("Tempo [s]", fontdict=font_axes)
    ax1.set_ylabel("Raio", fontdict=font_axes)
    ax1.set_zlabel("Temperatura", fontdict=font_axes)
    plt.tick_params(axis="both", which="major", labelsize=12)

    # Figure 2: 2D Heatmap
    fig2, ax2 = plt.subplots(figsize=(10, 7), dpi=600)
    c = ax2.pcolormesh(T_mesh, R_mesh, theta.T, cmap="viridis", shading="auto")
    ax2.set_xlabel("Tempo [s]", fontdict=font_axes)
    ax2.set_ylabel("Raio", fontdict=font_axes)
    plt.tick_params(axis="both", which="major", labelsize=12)
    cbar = fig2.colorbar(c, ax=ax2)
    cbar.set_label("Temperatura", fontdict=font_axes)
    cbar.ax.tick_params(labelsize=12)

    # Figure 3: Comparison with reference curves
    vt = np.array([0, 0.5, 1, 1.5])
    v_indices = [np.argmin(np.abs(t_eval - time_point)) for time_point in vt]

    colors = ["#EDB120", "#0072BD", "#D95319", "#7E2F8E"]
    curves = [curve01, curve02, curve03, curve04]

    plt.figure(figsize=(12, 8), dpi=600)
    for i in range(len(vt)):
        plt.plot(
            r,
            theta[v_indices[i], :],
            color=colors[i],
            label=f"t = {vt[i]}",
        )
        plt.plot(
            curves[i][:, 0],
            curves[i][:, 1],
            "x",
            color=colors[i],
            markersize=8,
            label=f"t = {vt[i]} (Soares2017)",
        )

    plt.legend(prop=font_legend)
    plt.xlabel("Raio", fontdict=font_axes)
    plt.ylabel("Temperatura", fontdict=font_axes)
    plt.grid(True)
    plt.tick_params(axis="both", which="major", labelsize=12)

    plt.show()

    # # Comparative Analysis
    # vt = np.array([0.5])
    # v_indices = [np.argmin(np.abs(t_eval - time_point)) for time_point in vt]
    # values = theta[v_indices[0], :]
    # np.save("b40.npy", values)
