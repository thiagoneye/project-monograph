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
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time


def solve_heat_pde(t, theta, r, dr, Bi, Gast):
    nr = len(theta)
    dr2 = dr**2

    # Pre-allocate derivative vectors
    dtheta_dr = np.zeros(nr)
    d2theta_dr2 = np.zeros(nr)

    # Vectorized Calculation of Spatial Derivatives

    # Interior points (2nd order central differences)
    i = np.arange(1, nr - 1)
    dtheta_dr[i] = (theta[i + 1] - theta[i - 1]) / (2 * dr)
    d2theta_dr2[i] = (theta[i + 1] - 2 * theta[i] + theta[i - 1]) / dr2

    # Boundary Conditions
    
    d2theta_dr2[0] = 2 * (theta[1] - theta[0]) / dr2
    # dtheta_dr[0] is zero and is not needed for the final equation at the center.
    dtheta_dr[-1] = (theta[-1] - theta[-2]) / dr
    d2theta_dr2[-1] = (2 * theta[-2] - 2 * (1 + Bi * dr) * theta[-1]) / dr2

    # Partial Differential Equation

    # Source term varying linearly with time
    Glin = Gast * t

    # Calculate the time derivative for all points (except the center, initially)
    # Note: 1/r is problematic at r=0, so we handle that case separately.
    dtheta_dt = (1 / r) * dtheta_dr + d2theta_dr2 + Glin

    # Correct the equation for the central point (r=r0) using the limit form
    dtheta_dt[0] = 2 * d2theta_dr2[0] + Glin

    return dtheta_dt


if __name__ == "__main__":
    # Inputs
    r0 = 0.001  # Start of the r-axis (avoids division by zero)
    rl = 1.0  # End of the r-axis (Length L)
    t0 = 0.0  # Start time
    tl = 1.5  # End time
    Bi = 15.0  # Biot number
    Gast = 32.4  # Heat source term (base constant)

    nr = 86  # Points in spatial mesh
    number_of_meshes = 13  # Number of temporal meshes for convergence study
    mesh = np.zeros(
        number_of_meshes, dtype=int
    )  # Pre-allocation for temporal mesh points

    # Calculations

    # Initial parameters for the convergence study
    mesh[0] = 11  # Initial number of points in the temporal mesh
    # Convergence matrix: (time_points, spatial_points, mesh_refinement_level)
    convergence_matrix = np.zeros((mesh[0], nr, number_of_meshes))

    # Generate progressively finer temporal meshes
    for idx in range(1, number_of_meshes):
        mesh[idx] = (mesh[idx - 1] - 1) * 2 + 1

    # Spatial discretization (done only once)
    r = np.linspace(r0, rl, nr)
    dr = (rl - r0) / (nr - 1)

    # Initial condition: Steady-state temperature profile
    theta0 = Gast * (1 - r**2) / 4 + Gast / (2 * Bi)

    # Simulation for each temporal mesh
    start_time = time.time()

    for idx in range(number_of_meshes):
        nt = mesh[idx]  # Current number of points in temporal mesh
        t_eval = np.linspace(t0, tl, nt)  # Time vector for solver output

        # Solve the system of ODEs
        # We use a lambda function to pass extra parameters to solve_heat_pde
        sol = solve_ivp(
            lambda t, theta: solve_heat_pde(t, theta, r, dr, Bi, Gast),
            t_span=[t0, tl],
            y0=theta0,
            method="BDF",  # A good method for stiff problems, similar to ode15s
            t_eval=t_eval,
        )
        # The result from solve_ivp is in shape (nr, nt), so we transpose it
        theta = sol.y.T

        # Sample the results to compare with the coarsest mesh
        index_vector = np.round(np.linspace(0, nt - 1, mesh[0])).astype(int)
        convergence_matrix[:, :, idx] = theta[index_vector, :]

    end_time = time.time()

    # Temporal Discretization Error Calculation

    # Calculate the difference between solutions of consecutive meshes
    # The axis=2 specifies that the difference is taken along the mesh refinement dimension
    error_matrix = np.diff(convergence_matrix, axis=2)

    # Find the maximum absolute error across the entire domain (r, t) for each refinement level
    maximum_error = np.max(np.abs(error_matrix), axis=(0, 1))

    # Create and Display DataFrame

    # The meshes corresponding to the calculated errors are the finer ones
    refined_meshes = mesh[1:]

    # Calculate the mesh size (time step, dt) for each refined mesh
    mesh_sizes = (tl - t0) / (refined_meshes - 1)

    # Create a dictionary with the results
    results_data = {
        "Pontos na Malha": refined_meshes,
        "Tamanho da Malha": mesh_sizes,
        "Erro Máximo": maximum_error,
    }

    # Create a pandas DataFrame
    results_df = pd.DataFrame(results_data)

    # Set display options to show all rows and format numbers
    pd.set_option("display.max_rows", None)
    pd.options.display.float_format = "{:,.8}".format

    # Print the DataFrame
    print(results_df)
