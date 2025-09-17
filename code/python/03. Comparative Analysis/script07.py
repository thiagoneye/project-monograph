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
import matplotlib.pyplot as plt

# Input Parameters
# Define the radial coordinate for the plot.
# nr is set to 4 to match the length of the data arrays.
Bi = 15
r0 = 0  # Beginning of the r axis
rl = 1  # End of the r axis (Length L)
nr = 400  # Points in spatial grid

# Create the radial coordinate vector
r = np.linspace(r0, rl, nr)

# Data Arrays
# Pre-calculated temperature profile data.
curve01 = np.load("c000.npy")
curve02 = np.load("c025.npy")
curve03 = np.load("c050.npy")
curve04 = np.load("c075.npy")
curve05 = np.load("c100.npy")

# Plotting
plt.figure(figsize=(12, 8), dpi=600)

# Plot each curve
plt.plot(r, curve01, label="$c_1 = 0$")
plt.plot(r, curve02, label="$c_1 = 0.25$")
plt.plot(r, curve03, label="$c_1 = 0.5$")
plt.plot(r, curve04, label="$c_1 = 0.75$")
plt.plot(r, curve05, label="$c_1 = 1$")

# Add plot details
font_axes = {"size": 12}
font_legend = {"size": 12}

plt.legend(prop=font_legend)
plt.xlabel("Raio", fontdict=font_axes)
plt.ylabel("Temperatura", fontdict=font_axes)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=12)

# Display the plot
plt.show()
