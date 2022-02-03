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
# %%
from enum import auto
import matplotlib.pyplot as plt
import numpy as np

itmax = 30
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


def encontra_autovalor_dominante(matriz_A):
    autovalores = (np.linalg.eig(matriz_A))[0]
    autovalor_dominante = 0

    for i in range(autovalores.size):
        modulo_autovalor = np.abs(autovalores[i])
        if (modulo_autovalor > autovalor_dominante):
            autovalor_dominante = modulo_autovalor

    return autovalor_dominante


def encontra_lambda_2(matriz_A):
    autovalores = (np.linalg.eig(matriz_A))[0]

    index_autovalor_dominante = encontra_index_autovalor_dominante(autovalores)

    autovalores_sem_dominante = np.delete(
        autovalores, index_autovalor_dominante)

    lambda_2 = 0

    for i in range(autovalores_sem_dominante.size):
        modulo_autovalor = np.abs(autovalores_sem_dominante[i])
        if (modulo_autovalor > lambda_2):
            lambda_2 = modulo_autovalor

    return lambda_2


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


y_erro_autovalor = []
y_erro_autovetor = []
x_iteracoes = []


def calcula_erros_assintoticos(matriz_A):
    lambda_1 = encontra_autovalor_dominante(matriz_A)
    lambda_2 = encontra_lambda_2(matriz_A)

    erros_assintoticos_por_iteracao = []

    for i in range(len(x_iteracoes)):
        erro_assintotico = (lambda_2/lambda_1)**i
        erros_assintoticos_por_iteracao.append(erro_assintotico)

    return erros_assintoticos_por_iteracao


def plotagem_grafico_erros(matriz_A):
    fig = plt.figure(figsize=(8, 6))

    y_erros_assintoticos_simples = calcula_erros_assintoticos(matriz_A)

    y_erros_assintoticos_quadraticos = []
    for erro in y_erros_assintoticos_simples:
        erro_quadratico = erro**2
        y_erros_assintoticos_quadraticos.append(erro_quadratico)

    plt.plot(x_iteracoes, y_erro_autovalor, color='black')
    plt.plot(x_iteracoes, y_erro_autovetor, color='green')
    plt.plot(x_iteracoes,
             y_erros_assintoticos_simples, color='blue')
    plt.plot(x_iteracoes,
             y_erros_assintoticos_quadraticos, color='red')
    plt.yscale("log")
    plt.xlabel("Iterações")
    plt.ylabel("Erro L2")

    plt.show()


def calcula_erro_autovetor(vetor_xk, autovetor_dominante):
    erro_autovetor = np.linalg.norm(vetor_xk - autovetor_dominante)

    return erro_autovetor


def calcula_erro_autovalor(autovalor_uk, autovalor_dominante):
    erro_autovalor = np.abs(autovalor_uk - autovalor_dominante)

    return erro_autovalor


def metodo_das_potencias(matriz_A, vetor_x0):
    autovalor_Uk = None
    autovetor_Xk = vetor_x0
    i = 0

    autovetor_dominante = encontra_autovetor_dominante(matriz_A)
    autovalor_dominante = encontra_autovalor_dominante(matriz_A)

    erro_autovetor = calcula_erro_autovetor(x0, autovetor_dominante)

    while(erro_autovetor > epsilon and i < itmax):
        autovetor_Xk = calcula_Xk(matriz_A, autovetor_Xk)
        autovalor_Uk = calcula_Uk(autovetor_Xk, matriz_A)

        erro_autovetor = calcula_erro_autovetor(
            autovetor_Xk, autovetor_dominante)
        erro_autovalor = calcula_erro_autovalor(
            autovalor_Uk[0][0], autovalor_dominante)

        y_erro_autovetor.append(erro_autovetor)
        y_erro_autovalor.append(erro_autovalor)
        x_iteracoes.append(i)

        i += 1

    return autovalor_Uk[0][0]


Uk_final = metodo_das_potencias(A, x0)

plotagem_grafico_erros(A)


# %%
