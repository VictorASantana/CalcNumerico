# Arquivo referente ao exercicio 4
# %%
from xml.etree.ElementTree import tostring
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# Exercicio 4.1
# Definicao e plotagem das redes G1 e G2
# =========================================================

# Rede G1
#
#   0  1   2   3  4   5
# 5 a             
# 4    a   b   c  
# 3       cab     d 
# 2    c   b   a  d  
# 1 c      b      ad 
# 0 c      b      d   d
#

n = 6

def plotagemRede1():
    # linha 1
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([5, 4, 3, 2, 1])
    plt.scatter(x, y, marker='o')
    plt.plot(x, y)

    # linha 2
    x = np.array([3, 3, 3, 3, 3])
    y = np.array([4, 3, 2, 1, 0])
    plt.scatter(x, y, marker='o')
    plt.plot(x, y)

    # linha 3
    x = np.array([1, 1, 2, 3, 4])
    y = np.array([0, 1, 2, 3, 4])
    plt.scatter(x, y, marker='o')
    plt.plot(x, y)

    # linha 4
    x = np.array([5, 5, 5, 5, 6])
    y = np.array([3, 2, 1, 0, 0])
    plt.scatter(x, y, marker='o')
    plt.plot(x, y)

    # especificacoes do grafico
    plt.title(r"$Rede G_{1}$")
    plt.grid(visible=True, axis='both', alpha=0.2)

    for i in range(0, n-1):
        plt.text(x[i] + 0.05, y[i] + 0.07, 'v')
    
    plt.show()

def calculaMatrizAdjacencia():
    a = 1

plotagemRede1()

# %%
