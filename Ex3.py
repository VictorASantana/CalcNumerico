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

#Calcula a aplicação de vetor_k nas colunas com índice maior que i
def calcula_iteracao(vetor_coluna, vetor_k, n):
    vetor_coluna_novo = np.zeros((n, 1))
    j = 0
    print(vetor_coluna)
    print(vetor_k)
    while j < n:
        vetor_coluna_novo[j] = vetor_coluna[j] - 2 * (prod_escalar(vetor_coluna, vetor_k, n))/(prod_escalar(vetor_k, vetor_k, n)) * vetor_k[j]
        j = j + 1
    return vetor_coluna_novo

#Calcula o produto escalar dos vetores A e B
def prod_escalar(A, B, n):
    prod = 0
    i = 0
    while i < n:
        prod = prod + A[i]*B[i]
        i = i + 1
    return prod

#Main
matriz_A = np.array([[6, -2, -1],[-2, 6, -1], [-1, -1, 5]])
vetor = coluna(matriz_A, 0, 3)
print(vetor)
vetor_i = calcula_vetor(vetor, 0, 3)
print(vetor_i)
vetor = coluna(matriz_A, 1, 3)
vetor_novo = calcula_iteracao(vetor, vetor_i, 3)
print(vetor_novo)
autovalores, autovetores = np.linalg.eig(matriz_A)
autovetor_normalizado = np.zeros((3,3))