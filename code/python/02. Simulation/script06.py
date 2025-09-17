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


def solve_heat_pde(t, theta, r, Bi, Gast, c2, c3):
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
    G = Gast * (1 + c2 * r**2) * np.exp(c3 * t)

    # Calculate the time derivative for all points
    # A small epsilon is added to r to prevent division by zero at r=0.
    # The value at r=0 will be correctly overwritten in the next step.
    dtheta_dt = (1 / (r + 1e-12)) * dtheta_dr + d2theta_dr2 + G

    # Correct the equation for the central point (r=0)
    dtheta_dt[0] = 2 * d2theta_dr2[0] + G[0]

    return dtheta_dt


if __name__ == "__main__":
    # Input Parameters
    nr = 400 # Points in spatial grid #TODO  9601
    nt = 100  # Points in temporal grid
    r0 = 0.001  # Beginning of the r axis
    rl = 1.0  # End of the r axis (Length L)
    t0 = 0.0  # Start time
    tl = 1.5  # End time
    Bi = 40.0  # Biot Number #TODO
    Gast = 32.4  # Heat source term constant
    c2 = 1.0  # New constant for the source term
    c3 = 1.0  # New constant for the source term

    # Reference Data Curves
    curve01 = np.array(
        [
            [0.00609, 9.10682],
            [0.09589, 9.00000],
            [0.18722, 8.83976],
            [0.30441, 8.35905],
            [0.39422, 7.87834],
            [0.50837, 7.02374],
            [0.59209, 6.27596],
            [0.70167, 5.12760],
            [0.808219, 3.81899],
            [0.910198, 2.37685],
            [0.984779, 1.25519],
        ]
    )
    curve02 = np.array(
        [
            [0.00000, 9.37389],
            [0.09741, 9.26706],
            [0.21918, 9.05341],
            [0.30898, 8.73294],
            [0.39878, 8.30564],
            [0.51750, 7.58457],
            [0.60274, 6.89021],
            [0.71385, 5.79525],
            [0.802131, 4.70030],
            [0.890411, 3.25816],
            [0.998478, 1.38872],
        ]
    )
    curve03 = np.array(
        [
            [0.00000, 11.91100],
            [0.09741, 11.80420],
            [0.18265, 11.67060],
            [0.33181, 11.10980],
            [0.39878, 10.73590],
            [0.52816, 9.72107],
            [0.62557, 8.70623],
            [0.69406, 7.85163],
            [0.79300, 6.30267],
            [0.90563, 3.95252],
            [0.99696, 1.81602],
        ]
    )
    curve04 = np.array(
        [
            [0.00000, 15.70330],
            [0.10198, 15.56970],
            [0.20396, 15.30270],
            [0.28767, 14.87540],
            [0.39726, 14.10090],
            [0.51294, 12.95250],
            [0.59817, 11.80420],
            [0.70929, 9.98813],
            [0.799087, 8.09199],
            [0.893455, 5.50148],
            [0.996956, 2.35015],
        ]
    )

    # Numerical Calculation
    r = np.linspace(r0, rl, nr)
    t_eval = np.linspace(t0, tl, nt)
    theta0 = Gast * (1 - r**2) / 4 + Gast / (2 * Bi)

    # Solve the system of ODEs
    sol = solve_ivp(
        lambda t, theta: solve_heat_pde(t, theta, r, Bi, Gast, c2, c3),
        t_span=[t0, tl],
        y0=theta0,
        method="BDF",  # Method for stiff problems
        t_eval=t_eval,
    )
    theta = sol.y.T  # Transpose solution to shape (time, space)

    # Plotting

    font_axes = {"size": 12}
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
    vt = np.array(
        [0, 0.05, 0.24, 0.5]
    )  # Note: Different time points from the legend in MATLAB
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

    # plt.show()

    # Comparative Analysis
    vt = np.array([0.05])
    v_indices = [np.argmin(np.abs(t_eval - time_point)) for time_point in vt]
    values = theta[v_indices[0], :]
    np.save("b40_1.npy", values)

