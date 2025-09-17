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

c2 = 1
c3 = 1

r0 = 0  # Início do eixo r
rl = 1  # Fim do eixo r (Comprimento L)
nr = 400  # Pontos na malha espacial

r = np.linspace(r0, rl, nr)

curve01 = np.load("b15_1.npy")
curve02 = np.load("b40_1.npy")

# Para as matrizes, usamos np.array()
curve03 = np.array(
    [
        [0.00000, 9.39024],
        [0.09359, 9.30894],
        [0.20820, 9.10569],
        [0.29127, 8.83469],
        [0.40063, 8.33333],
        [0.50578, 7.68293],
        [0.60463, 6.89702],
        [0.697161, 6.00271],
        [0.802313, 4.729],
        [0.91062, 2.98103],
        [0.995794, 1.47696],
    ]
)

curve04 = np.array(
    [
        [0.00105, 8.71274],
        [0.10200, 8.60434],
        [0.20400, 8.44173],
        [0.29232, 8.15718],
        [0.40799, 7.60163],
        [0.504732, 6.99187],
        [0.599369, 6.23306],
        [0.70347, 5.20325],
        [0.808623, 3.90244],
        [0.901157, 2.37127],
        [0.997897, 0.609756],
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
