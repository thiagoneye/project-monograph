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

c1 = 1

r0 = 0  # Início do eixo r
rl = 1  # Fim do eixo r (Comprimento L)
nr = 400  # Pontos na malha espacial

r = np.linspace(r0, rl, nr)

curve01 = np.load("b15.npy")
curve02 = np.load("b40.npy")

# Para as matrizes, usamos np.array()
curve03 = np.array(
    [
        [0.00681, 11.91590],
        [0.09295, 11.79140],
        [0.20584, 11.47800],
        [0.29079, 11.05880],
        [0.39900, 10.32420],
        [0.50020, 9.37895],
        [0.59906, 8.26524],
        [0.690922, 7.02507],
        [0.800208, 5.30106],
        [0.900169, 3.42948],
        [0.994309, 1.55777],
    ]
)

curve04 = np.array(
    [
        [0.00093, 11.15790],
        [0.08939, 11.03350],
        [0.21160, 10.69920],
        [0.30818, 10.17490],
        [0.40591, 9.46118],
        [0.49664, 8.62098],
        [0.58271, 7.63332],
        [0.70946, 5.93073],
        [0.82106, 4.08044],
        [0.89429, 2.67146],
        [0.99772, 0.54732],
    ]
)

plt.figure(figsize=(12, 8), dpi=600)
plt.plot(r, curve01, color="#EDB120", label="Bi = 15")
plt.plot(
    curve03[:, 0], curve03[:, 1], "x", color="#EDB120", label="Bi = 15 (Soares2017)"
)
plt.plot(r, curve02, color="#0072BD", label="Bi = 40")
plt.plot(
    curve04[:, 0], curve04[:, 1], "x", color="#0072BD", label="Bi = 40 (Soares2017)"
)

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
