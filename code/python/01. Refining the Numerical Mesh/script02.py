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
import time
import pandas as pd


def solve_heat_pde(t, theta, r, Bi, Gast):
    nr = len(theta)
    r0 = r[0]
    rl = r[-1]
    dr = (rl - r0) / (nr - 1)
    dr2 = dr**2

    # Pre-allocate derivative vectors
    dtheta_dr = np.zeros(nr)
    d2theta_dr2 = np.zeros(nr)

    # Vectorized Calculation of Spatial Derivatives

    i = np.arange(1, nr - 1)
    dtheta_dr[i] = (theta[i] - theta[i - 1]) / dr
    d2theta_dr2[i] = (theta[i + 1] - 2 * theta[i] + theta[i - 1]) / dr2

    # Boundary Conditions

    d2theta_dr2[0] = 2 * (theta[1] - theta[0]) / dr2
    dtheta_dr[-1] = (theta[-1] - theta[-2]) / dr  # Backward difference
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
    r0 = 0.001  # Beginning of the r axis
    rl = 1.0  # End of the r axis (Length L)
    t0 = 0.0  # Start time
    tl = 1.5  # End time
    Bi = 15.0  # Biot number
    Gast = 32.4  # Heat source term

    nt = 101  # Points in temporal mesh (fixed)
    number_of_meshes = 8  # Number of spatial meshes for convergence study
    mesh = np.zeros(
        number_of_meshes, dtype=int
    )  # Pre-allocation for spatial mesh points

    # Calculations

    # Initial Parameters for convergence study
    mesh[0] = 76  # Initial number of points in the spatial mesh
    # Convergence matrix: (time_points, coarsest_spatial_points, mesh_refinement_level)
    convergence_matrix = np.zeros((nt, mesh[0], number_of_meshes))

    # Generate progressively finer spatial meshes
    for idx in range(1, number_of_meshes):
        mesh[idx] = (mesh[idx - 1] - 1) * 2 + 1

    # Time vector (calculated only once)
    t_eval = np.linspace(t0, tl, nt)

    # Simulation for each spatial mesh
    start_time = time.time()

    for idx in range(number_of_meshes):
        nr = mesh[idx]  # Current number of points in spatial mesh

        # Spatial discretization for the current mesh
        r = np.linspace(r0, rl, nr)

        # Initial condition for the current mesh
        theta0 = Gast * (1 - r**2) / 4 + Gast / (2 * Bi)

        # Solve the system of ODEs
        sol = solve_ivp(
            lambda t, theta: solve_heat_pde(t, theta, r, Bi, Gast),
            t_span=[t0, tl],
            y0=theta0,
            method="BDF",  # Method suitable for stiff problems
            t_eval=t_eval,
        )
        # Transpose result to have shape (nt, nr)
        theta = sol.y.T

        # Sample the results at locations corresponding to the coarsest mesh
        # This allows for a direct comparison between meshes
        index_vector = np.round(np.linspace(0, nr - 1, mesh[0])).astype(int)
        convergence_matrix[:, :, idx] = theta[:, index_vector]

    end_time = time.time()

    # Spatial Discretization Error Calculation

    # Calculate the difference between solutions of consecutive meshes
    error_matrix = np.diff(convergence_matrix, axis=2)

    # Find the maximum absolute error across the entire domain (r, t)
    maximum_error = np.max(np.abs(error_matrix), axis=(0, 1))

    # Create and Display DataFrame

    # The meshes corresponding to the calculated errors are the finer ones
    refined_meshes = mesh[1:]

    # Calculate the spatial mesh size (step size, dr) for each refined mesh
    mesh_sizes = (rl - r0) / (refined_meshes - 1)

    # Create a dictionary with the results
    results_data = {
        "Points in Mesh": refined_meshes,
        "Mesh Size": mesh_sizes,
        "Maximum Error": maximum_error,
    }

    # Create a pandas DataFrame
    results_df = pd.DataFrame(results_data)

    # Set display options to show all rows and format numbers
    pd.set_option("display.max_rows", None)
    pd.options.display.float_format = "{:,.8f}".format

    # Print the DataFrame
    print(results_df)
