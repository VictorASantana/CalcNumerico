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
from matplotlib.pyplot import xkcd
import numpy as np

itmax = 30
epsilon = 10**(-15)

n = 8

A = np.array([[-2, -4, 2], [-2, 1, 2], [4, 2, 5]])

x0 = np.random.rand(n, 1)

B = np.random.rand(2, 2)

autovalores, autovetores = np.linalg.eigh(A)

# achar autovetor do maior autovalor


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


def metodo_das_potencias(A, x0):
    Uk = None
    Xk = x0
    i = 0

    erro_autovetor = np.linalg.norm(x0 - autovetores[0])

    while(erro_autovetor > epsilon and i < itmax):
        Xk = calcula_Xk(A, Xk)
        Uk = calcula_Uk(Xk, A)

        erro_autovetor = np.linalg.norm(Xk - autovetores[0])
        i += 1
