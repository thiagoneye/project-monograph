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


def solve_heat_pde(t, theta, r, Bi, Gast, c3):
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

    # Define the new source term, which varies with r and t
    G = Gast * (r**2) * np.exp(c3 * t)

    # Calculate the time derivative for all points
    dtheta_dt = (1 / r) * dtheta_dr + d2theta_dr2 + G

    # Correct the equation for the central point (r=0)
    dtheta_dt[0] = 2 * d2theta_dr2[0] + G[0]

    return dtheta_dt


if __name__ == "__main__":
    # Input Parameters
    nr = 9601  # Points in spatial grid
    nt = 100  # Points in temporal grid
    r0 = 0.001  # Beginning of the r axis (small value to avoid division by zero)
    rl = 1.0  # End of the r axis (Length L)
    t0 = 0.0  # Start time
    tl = 0.5  # End time
    Bi = 15.0  # Biot Number
    Gast = 32.4  # Heat source term constant
    c3 = 4.0  # New constant for the source term

    # Reference Data Curves
    curve01 = np.array(
        [
            [0.00609, 9.22656],
            [0.09589, 9.11830],
            [0.21461, 8.84952],
            [0.30289, 8.47462],
            [0.41400, 7.81486],
            [0.49620, 7.20896],
            [0.59970, 6.28267],
            [0.69863, 5.24980],
            [0.791476, 4.11037],
            [0.896499, 2.72184],
            [0.987823, 1.29799],
        ]
    )
    curve02 = np.array(
        [
            [0.00000, 7.75111],
            [0.11263, 7.64244],
            [0.18722, 7.51667],
            [0.29985, 7.17689],
            [0.40792, 6.69497],
            [0.51294, 6.07088],
            [0.61339, 5.32243],
            [0.71233, 4.46734],
            [0.80670, 3.50566],
            [0.891933, 2.43748],
            [0.977169, 1.31596],
        ]
    )
    curve03 = np.array(
        [
            [0.00000, 5.51111],
            [0.11263, 5.50911],
            [0.18722, 5.50778],
            [0.30137, 5.48798],
            [0.41400, 5.41486],
            [0.49011, 5.32462],
            [0.60122, 5.03820],
            [0.71081, 4.55625],
            [0.81126, 3.80780],
            [0.908021, 2.73959],
            [0.980213, 1.56480],
        ]
    )
    curve04 = np.array(
        [
            [0.00152, 9.31553],
            [0.11263, 9.43800],
            [0.18722, 9.56112],
            [0.29833, 9.89692],
            [0.40944, 10.23270],
            [0.48402, 10.42700],
            [0.59665, 10.42490],
            [0.70472, 9.92525],
            [0.79148, 8.89260],
            [0.91411, 5.95726],
            [0.946073, 4.94336],
        ]
    )

    # Numerical Calculation
    r = np.linspace(r0, rl, nr)
    t_eval = np.linspace(t0, tl, nt)
    theta0 = Gast * (1 - r**2) / 4 + Gast / (2 * Bi)

    # Solve the system of ODEs
    sol = solve_ivp(
        lambda t, theta: solve_heat_pde(t, theta, r, Bi, Gast, c3),
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
    vt = np.array([0, 0.05, 0.25, 0.5])
    v_indices = [np.argmin(np.abs(t_eval - time_point)) for time_point in vt]

    colors = ["#EDB120", "#0072BD", "#D95319", "#7E2F8E"]
    curves = [curve01, curve02, curve03, curve04]
    legends = ["t = 0", "t = 0.05", "t = 0.25", "t = 0.5"]

    plt.figure(figsize=(12, 8), dpi=600)
    for i in range(len(vt)):
        plt.plot(
            r,
            theta[v_indices[i], :],
            color=colors[i],
            label=f"{legends[i]}",
        )
        plt.plot(
            curves[i][:, 0],
            curves[i][:, 1],
            "x",
            color=colors[i],
            markersize=8,
            label=f"{legends[i]} (Soares2017)",
        )

    plt.legend(prop=font_legend)
    plt.xlabel("Raio", fontdict=font_axes)
    plt.ylabel("Temperatura", fontdict=font_axes)
    plt.grid(True)
    plt.tick_params(axis="both", which="major", labelsize=12)

    plt.show()
