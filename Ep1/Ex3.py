# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-=
# Exercício Computacional desenvolvido por:                                 =
#                                                                           =
# Nome: Pedro Henrique Rodrigues de Viveiros                                =
# NUSP: 11804035                                                            =
#                                                                           =
# Nome: Victor de Almeida Santana                                           =
# NUSP: 11806718                                                            =
#                                                                           =
# Entrega em 20/02/2022                                                     =
# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-=

import numpy as np

#Exercicio 1
#Calcular o valor exato dos autovalores

#Normalizar os autovetores => autovetor = autovetor/||autovetor||

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
        if j < i:
            vetor_coluna[j] = 0
        j = j + 1
    j = 0
    while j < n:
        if j == i:
            if vetor_coluna[j] >= 0:
                vetor_k[j] = vetor_coluna[j] + np.linalg.norm(vetor_coluna)
            else:
                vetor_k[j] = vetor_coluna[j] - np.linalg.norm(vetor_coluna)
        elif j > i:
            vetor_k[j] = vetor_coluna[j]
        j = j + 1
    return vetor_k

#Calcula a aplicação de vetor_k nas colunas com índice maior que i
def calcula_iteracao(vetor_coluna, vetor_k, n):
    vetor_coluna_novo = np.zeros((n, 1))
    j = 0
    escalar_coluna = prod_escalar(vetor_coluna, vetor_k, n)
    escalar_k = prod_escalar(vetor_k, vetor_k, n)
    #print(vetor_coluna)
    while j < n:
        elemento_coluna = vetor_coluna[j]
        operador = (-2 * (escalar_coluna/escalar_k) * vetor_k[j])
        vetor_coluna_novo[j] = elemento_coluna + operador
        j = j + 1
    return vetor_coluna_novo

#Calcula a matriz R da i-ésima iteração
def calcula_R(matriz_A, n, index):
    matriz_R = np.zeros((n,n))
    i = 0
    j = 0
    #Atribui os coeficientes da matriz_A para a matriz_R
    while i < n:
        while j < n:
            matriz_R[i][j] = matriz_A[i][j]
            j = j + 1
        j = 0
        i = i + 1
    i = index
    j = 0
    k = 0
    #Calcula o v_i da i-ésima iteração
    while i < n-1:
        v_coluna = coluna(matriz_A, i, n)
        #print("Vetor coluna: ")
        #print(v_coluna)
        v_i = calcula_vetor(v_coluna, i, n)
        #print("Vetor i: ")
        #print(v_i)
        #A i-ésima coluna da matriz é subtraida em v_i
        while j < n:
            #print("Elementos: ")
            #print(matriz_R[j][i])
            #print(v_i[j])
            matriz_R[j][i] = matriz_R[j][i] - v_i[j]
            j = j + 1
        k = i + 1
        j = 0
        #Calcula-se a Hx = x - 2*\frac{x \cdot v}{v \cdot v}v
        while k < n:
            j = 0
            v_coluna = coluna(matriz_R, k, n)
            v_coluna = calcula_iteracao(v_coluna, v_i, n)
            while j < n:
                matriz_R[j][k] = v_coluna[j]
                j = j + 1
            k = k + 1
            j = 0
        i = i + 1
        j = 0
    return matriz_R

#Calcula a matriz_Q da i-ésima iteração
def calcula_Q(matriz_A, matriz_R):
    return np.dot(matriz_A, np.linalg.inv(matriz_R))

#Calcula os autovalores pelo metodo QR
def calcula_autovaloresQR(matriz_A, n, autovalores, epsilon):
    i = 0
    j = 0
    k = 0
    global num_iteracoes
    R = np.zeros((n,n))
    Q = np.zeros((n, n))
    diag = np.zeros((n,1))
    while True:
        while i < n:
            while j < n:
                R[i][j] = matriz_A[i][j]
                j = j + 1
            j = 0
            i = i + 1
        i = 0
        j = 0
        while i < n-1:
            R = calcula_R(R, n, i)
            Q = calcula_Q(matriz_A, R)
            vetor_Q.append(Q)
            num_iteracoes = num_iteracoes + 1
            i = i + 1
        #print("Matriz q: ")
        #print(Q)
        matriz_A = np.dot(R, Q)
        diag = diagonal(matriz_A, n)
        if parada(diag, n, autovalores, epsilon):
            break
        elif k > 70:
            break
        i = 0
        k = k + 1
    return diag

#Calcula a matriz V_K
def calcula_Vk(vetor_Q, n, num_int):
    i = 1
    matriz_V = vetor_Q[0]
    while i < num_int :
        matriz_V = np.dot(matriz_V, vetor_Q[i])
        i = i + 1
    return matriz_V

#Devolve um vetor com os elementos da diagonal da matriz
def diagonal(matriz_A, n):
    diagonal = np.zeros((n, 1))
    i = 0
    j = 0
    while i < n:
        while j < n:
            if i == j:
                diagonal[i] = matriz_A[i][j]
            j = j + 1
        i = i + 1
        j = 0
    return diagonal

#Criterio de parada
def parada(autovalores_calculados, n, autovalores, epsilon):
    fim = False
    i = 0
    while i < n:
        if np.abs(autovalores[i] - autovalores_calculados[i]) < epsilon:
            fim = True
        i = i + 1
    return fim

#Calcula o produto escalar dos vetores A e B
def prod_escalar(A, B, n):
    prod = 0
    i = 0
    while i < n:
        prod = prod + A[i]*B[i]
        i = i + 1
    return prod

#Normaliza os autovetores
def normaliza(matriz_autovetores, n):
    autovet_normalizado = np.zeros((n,n))
    vetor_coluna = np.zeros((n,1))
    i = 0
    j = 0
    while i < n:
        vetor_coluna = coluna(matriz_autovetores, i, n)
        while j < n:
            autovet_normalizado[i][j] = matriz_autovetores[i][j]/np.linalg.norm(vetor_coluna)
            j = j + 1
        j = 0
        i = i + 1
    return autovet_normalizado

#Main
matriz_A = np.array([[5, -2], [-2, 8]])
n = 2
autovalores, autovetores = np.linalg.eig(matriz_A)
#Teste
epsilon = 0.00000000001
vetor_Q = []
num_iteracoes = 0
avalores_calculado = calcula_autovaloresQR(matriz_A, n, autovalores, epsilon)
#Calculando autovetores
#matriz_V = autovetores_sistema(matriz_A, n, autovalores)
matriz_V = calcula_Vk(vetor_Q, n, num_iteracoes)
#matriz_V = normaliza(matriz_V, n)

#Main
matriz_A = np.array([[6, -2, -1],[-2, 6, -1], [-1, -1, 5]])
n = 3
autovalores, autovetores = np.linalg.eig(matriz_A)
#Exercicio 1
vetor_Q = []
num_iteracoes = 0
avalores_calculado = calcula_autovaloresQR(matriz_A, n, autovalores, epsilon)
print("Autovalores: ")
print(autovalores)
print("Autovalores calculados: ")
print(avalores_calculado)
#Calculando autovetores
matriz_V = calcula_Vk(vetor_Q, n, num_iteracoes)
#matriz_V = normaliza(matriz_V, n)
print(autovetores)
print(matriz_V)

#Ex 2
matriz_B = np.array([[1, 1], [-3, 1]])
print("Iniciando o exercício 3.2: ")
print("O exercicio 3.2 compara o cálculo dos autovalores de uma matriz B a partir do polinômio característico e do método QR: ")
print("Matriz B: ")
print(matriz_B)
autovalores_B, autovetores_B = np.linalg.eig(matriz_B)
n = 2
autoval_B = calcula_autovaloresQR(matriz_B, n, autovalores_B, epsilon)
print("Autovalores da matriz B: ")
print(autovalores_B)
print("Autovalores calculados B pelo método QR: ")
print(autoval_B)

#Ex 3
matriz_C = np.array([[3, -3], [0.333333, 5]])
print("Iniciando o exercício 3.3: ")
print("O exercicio 3.3 compara o cálculo dos autovalores de uma matriz C a partir do polinômio característico e do método QR: ")
print("Matriz C: ")
print(matriz_C)
autovalores_C, autovetores_C = np.linalg.eig(matriz_C)
n = 2
autoval_C = calcula_autovaloresQR(matriz_C, n, autovalores_C, epsilon)
print("Autovalores da matriz C: ")
print(autovalores_C)
print("Autovalores calculados de C pelo método QR: ")
print(autoval_C)

#Ex 4
#matriz exercicio 1.a
n = 10
matriz_B = np.random.rand(n, n)
vetor_x0 = np.random.rand(n, 1)
matriz_A = np.matmul(matriz_B, matriz_B.T)
print("Iniciando o exercício 3.4: ")
print("O exercicio 3.4 compara o cálculo dos autovalores de uma matriz C a partir do polinômio característico e do método QR: ")
print("Matriz A: ")
print(matriz_A)
autovalores_A, autovetores_A = np.linalg.eig(matriz_A)
autoval_A = calcula_autovaloresQR(matriz_A, n, autovalores_A, epsilon)
print("Autovalores da matriz A: ")
print(autovalores_A)
print("Autovalores calculados de A pelo método QR: ")
print(autoval_A)

#matriz exercicio 1.b.i
n = 7
matriz_B = np.random.rand(n, n)

# # λ1 = 95 e λ2 = 92
matriz_D = np.array(
[[2, 0, 0, 0, 0, 0, 0],
[0, 13, 0, 0, 0, 0, 0],
[0, 0, 92, 0, 0, 0, 0],
[0, 0, 0, 32, 0, 0, 0],
[0, 0, 0, 0, 76, 0, 0],
[0, 0, 0, 0, 0, 95, 0],
[0, 0, 0, 0, 0, 0, 22]]
)

vetor_x0 = np.random.rand(n, 1)

matriz_interm = np.matmul(matriz_B, matriz_D)

inv_matriz_B = np.linalg.inv(matriz_B)
matriz_A = np.matmul(matriz_interm, inv_matriz_B)
print("Matriz A: ")
print(matriz_A)
autovalores_A, autovetores_A = np.linalg.eig(matriz_A)
autoval_A = calcula_autovaloresQR(matriz_A, n, autovalores_A, epsilon)
print("Autovalores da matriz A: ")
print(autovalores_A)
print("Autovalores calculados de A pelo método QR: ")
print(autoval_A)

# matriz exercicio 1.b.ii
n = 7
matriz_B = np.random.rand(n, n)

# # λ1 = 92 e λ2 = 13
matriz_D = np.array(
[[2, 0, 0, 0, 0, 0, 0],
[0, 13, 0, 0, 0, 0, 0],
[0, 0, 92, 0, 0, 0, 0],
[0, 0, 0, 10, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 8, 0],
[0, 0, 0, 0, 0, 0, 11]]
)

vetor_x0 = np.random.rand(n, 1)

matriz_interm = np.matmul(matriz_B, matriz_D)

inv_matriz_B = np.linalg.inv(matriz_B)
matriz_A = np.matmul(matriz_interm, inv_matriz_B)
print("Matriz A: ")
print(matriz_A)
autovalores_A, autovetores_A = np.linalg.eig(matriz_A)
autoval_A = calcula_autovaloresQR(matriz_A, n, autovalores_A, epsilon)
print("Autovalores da matriz A: ")
print(autovalores_A)
print("Autovalores calculados de A pelo método QR: ")
print(autoval_A)
