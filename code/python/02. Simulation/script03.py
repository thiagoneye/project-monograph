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
from mpl_toolkits.mplot3d import Axes3D  # Required for 3D plotting


def solve_heat_pde(t, theta, r, Bi, Gast):
    nr = len(theta)
    dr = (r[-1] - r[0]) / (nr - 1)
    dr2 = dr**2

    # Pre-allocate derivative vectors
    dtheta_dr = np.zeros(nr)
    d2theta_dr2 = np.zeros(nr)

    # Vectorized Calculation of Spatial Derivatives
    i = np.arange(1, nr)
    dtheta_dr[i] = (theta[i] - theta[i - 1]) / dr
    i = np.arange(1, nr - 1)
    d2theta_dr2[i] = (theta[i + 1] - 2 * theta[i] + theta[i - 1]) / dr2

    # Boundary Conditions
    dtheta_dr[0] = 0.0
    d2theta_dr2[0] = 2 * (theta[1] - theta[0]) / dr2
    d2theta_dr2[-1] = (2 * theta[-2] - 2 * (1 + Bi * dr) * theta[-1]) / dr2

    # Partial Differential Equation
    # Source term varying linearly with time
    Glin = Gast * t

    # Calculate the time derivative for all points
    dtheta_dt = (1 / r) * dtheta_dr + d2theta_dr2 + Glin

    # Correct the equation for the central point (r=r0) using the limit form
    dtheta_dt[0] = 2 * d2theta_dr2[0] + Glin

    return dtheta_dt


if __name__ == "__main__":
    # Inputs
    nr = 9601  # Points in spatial grid
    nt = 100  # Points in temporal grid
    r0 = 0.001  # Beginning of the r axis
    rl = 1.0  # End of the r axis (Length L)
    t0 = 0.0  # Start time
    tl = 1.5  # End time
    Bi = 15.0  # Biot number
    Gast = 32.4  # Heat source term

    # Reference Curves
    curve01 = np.array(
        [
            [0.00174, 9.22535],
            [0.11130, 9.08451],
            [0.19478, 8.92019],
            [0.29739, 8.49765],
            [0.40174, 7.91080],
            [0.50957, 7.11268],
            [0.58957, 6.38498],
            [0.709565, 5.14085],
            [0.808696, 3.92019],
            [0.909565, 2.48826],
            [0.998261, 1.17371],
        ]
    )

    curve02 = np.array(
        [
            [0.00348, 3.61502],
            [0.08174, 3.59155],
            [0.20348, 3.52113],
            [0.28348, 3.40376],
            [0.40348, 3.14554],
            [0.48348, 2.93427],
            [0.60348, 2.53521],
            [0.68174, 2.23005],
            [0.80000, 1.69014],
            [0.91652, 1.00939],
            [0.99130, 0.53991],
        ]
    )

    curve03 = np.array(
        [
            [0.00174, 7.30047],
            [0.08174, 7.25352],
            [0.20174, 7.06573],
            [0.28174, 6.83099],
            [0.40000, 6.36150],
            [0.51652, 5.70423],
            [0.59130, 5.18779],
            [0.70087, 4.29577],
            [0.80696, 3.26291],
            [0.90783, 2.08920],
            [0.972174, 1.24413],
        ]
    )

    curve04 = np.array(
        [
            [0.00348, 11.80750],
            [0.12174, 11.64320],
            [0.20174, 11.43190],
            [0.31826, 10.82160],
            [0.39304, 10.28170],
            [0.49913, 9.34272],
            [0.60348, 8.19249],
            [0.69739, 6.90141],
            [0.78783, 5.49296],
            [0.89913, 3.47418],
            [0.97913, 1.90141],
        ]
    )

    # Calculations
    r = np.linspace(r0, rl, nr)  # r axis
    t = np.linspace(t0, tl, nt)  # Time
    theta0 = Gast * (1 - r**2) / 4 + Gast / (2 * Bi)  # Initial temperature

    # Solve the system of ODEs
    sol = solve_ivp(
        lambda t, theta: solve_heat_pde(t, theta, r, Bi, Gast),
        t_span=[t0, tl],
        y0=theta0,
        method="BDF",  # Good for stiff problems, similar to ode15s
        t_eval=t,
    )
    # Transpose solution to have shape (nt, nr) for easier plotting
    theta = sol.y.T

    # Plots

    font_axes = {"weight": "normal", "size": 12}
    font_legend = {"size": 12}

    # Create meshgrid for 3D plots
    T_mesh, R_mesh = np.meshgrid(t, r)

    # Figure 1: 3D Surface Plot
    fig1 = plt.figure(figsize=(10, 7), dpi=600)
    ax1 = fig1.add_subplot(111, projection="3d")
    ax1.plot_surface(T_mesh, R_mesh, theta.T, cmap="viridis", edgecolor="none")
    ax1.set_xlabel("Tempo [s]", fontdict=font_axes)
    ax1.set_ylabel("Raio", fontdict=font_axes)
    ax1.set_zlabel("Temperatura", fontdict=font_axes)
    plt.tick_params(axis="both", which="major", labelsize=12)

    # Figure 2: 2D Heatmap (Top-down view of the 3D surface)
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
    # Find the indices of the time points closest to the values in vt
    v_indices = [np.argmin(np.abs(t - time_point)) for time_point in vt]

    colors = ["#EDB120", "#0072BD", "#D95319", "#7E2F8E"]
    curves = [curve01, curve02, curve03, curve04]

    plt.figure(figsize=(12, 8), dpi=600)
    for i in range(len(vt)):
        # Plot simulation result at the specific time slice
        plt.plot(r, theta[v_indices[i], :], color=colors[i], label=f"t = {vt[i]}")
        # Plot reference data
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

    # Display all figures
    plt.show()
