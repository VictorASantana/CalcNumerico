# Arquivo referente ao Exercicio 2

# Comandos uteis
# .ndim : dimensao do array/matriz
# .shape : indica o formato da matriz n x m
# .[r, c] : saida é o elemento contido no posicao (r, c) *conta a partir de zero
# .[:, c] : todos os elementos da coluna c, possivel de atribuir um valor a coluna c
# .[r, a:b:c] : acessa elementos da coluna r, indo do elemento a até o b, pulando de c em c
# np.zeros((n,m)) : inicializa uma matriz de zeros de dimensao n x m
# np.ones((n, m)) : inicializa uma matriz de uns de dimensao n x m
# np.full((n, m), a): inicializa uma matriz de dimensao n x m de valor a
# np.full_like (A, b): inicializa uma matriz de mesmo tamanho de A com todos os valores b
# np.random.randint()
# np.norm: retorna norma/modulo do vetor

import numpy as np


#A função calcula parada é usada para verificação do critério de parada, que é atingido quando ||x_k - x*|| < epsilon
def calcula_parada(x_k, x_converge, epsilon):
    if np.linalg.norm(x_k - x_converge) < epsilon:
        return True
    else:
        return False

#Calcula o sistema (3) - método SOR
def calcula_SOR(matriz_A, x_k):
    #Uma vez completa essa função, meio caminho estará andado
    x_k_mais_1 = 0
    return x_k_mais_1


#Calcula a aproximação para \lambda_n^{-1} na k-ésima operação
def calcula_autovalor(x_k, x_k_mais_1):
    mu_k = np.dot(x_k.T, x_k_mais_1)/np.dot(x_k.T, x_k)
    return mu_k