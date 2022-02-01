#Arquivo referente ao Exercicio 1

import numpy as np

def calcula_parada(x_k, x_converge, epsilon):
    if np.linalg.norm(x_k - x_converge) < epsilon:
        return True
    else:
        return False


