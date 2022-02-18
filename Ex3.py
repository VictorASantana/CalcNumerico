import numpy as np

#Exercicio 1
#Calcular o valor exato dos autovalores

#Normalizar os autovetores => autovetor = autovetor/||autovetor||





#Calcular autovalores pelo método QR
def calcula_autovaloresQR(matriz_A):
    return 0

#Calcula autovetores pelo método QR
def calcula_autovetoresQR(matriz_A):
    return 0

#Devolve a i-ésima coluna de uma matriz A
def coluna(matriz_A, i, n):
    vetor_coluna = np.zeros((n, 1))
    j = 0
    while j < n:
        vetor_coluna[j] = matriz_A[j][i]
        j = j + 1
    return vetor_coluna

#Retorna o vetor_k da i-ésima iteração
def calcula_vetor(vetor_coluna, i, n):
    vetor_k = np.zeros((n, 1))
    j = 0
    while j < n:
        if j == i:
            vetor_k[j] = vetor_coluna[j] + np.linalg.norm(vetor_coluna)
        else:
            vetor_k[j] = vetor_coluna[j]
        j = j + 1
    return vetor_k

#Main
matriz_A = np.array([[6, -2, -1],[-2, 6, -1], [-1, -1, 5]])
vetor = coluna(matriz_A, 0, 3)
print(vetor)
print(np.linalg.norm(vetor))
vetor_i = calcula_vetor(vetor, 0, 3)
print(vetor_i)
autovalores, autovetores = np.linalg.eig(matriz_A)
autovetor_normalizado = np.zeros((3,3))