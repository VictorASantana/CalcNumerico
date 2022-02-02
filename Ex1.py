# Arquivo referente ao Exercicio 1

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

# from matplotlib import *
import numpy as np

itmax = 70
epsilon = 10**(-15)

n = 3

A = np.array([[-2, -4, 2], [-2, 1, 2], [4, 2, 5]])

x0 = np.random.rand(n, 1)

B = np.random.rand(2, 2)

# achar autovetor do maior autovalor


def encontra_autovetor_dominante(matriz_A):
    autovalores, autovetores = np.linalg.eig(matriz_A)

    index = encontra_index_autovalor_dominante(autovalores)

    return autovetores[index]


def encontra_index_autovalor_dominante(autovalores):
    index_autovalor_dominante = 0
    autovalor_dominante = 0

    for i in range(autovalores.size):
        modulo_autovalor = np.abs(autovalores[i])
        if (modulo_autovalor > autovalor_dominante):
            autovalor_dominante = modulo_autovalor
            index_autovalor_dominante = i

    return index_autovalor_dominante


def calcula_Xk(A, x0):
    Xk = x0

    produto = np.matmul(A, Xk)
    Xk = produto / np.linalg.norm(produto)

    return Xk


def calcula_Uk(Xk, A):
    Xk_tranposta = Xk.T
    produto_A_Xk = np.matmul(A, Xk)

    Uk = np.matmul(Xk_tranposta, produto_A_Xk) / np.matmul(Xk_tranposta, Xk)

    return Uk


def metodo_das_potencias(matriz_A, vetor_x0):
    autovalor_Uk = None
    autovetor_Xk = vetor_x0
    i = 0

    autovetor_dominante = encontra_autovetor_dominante(matriz_A)

    erro_autovetor = np.linalg.norm(vetor_x0 - autovetor_dominante)

    while(erro_autovetor > epsilon and i < itmax):
        autovetor_Xk = calcula_Xk(matriz_A, autovetor_Xk)
        autovalor_Uk = calcula_Uk(autovetor_Xk, matriz_A)

        erro_autovetor = np.linalg.norm(autovetor_Xk - autovetor_dominante)
        i += 1

    return autovalor_Uk


Uk_final = metodo_das_potencias(A, x0)

print("O autovalor Uk obtido é:", Uk_final[0][0])
