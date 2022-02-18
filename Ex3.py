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
            i = i + 1
        matriz_A = np.dot(R, Q)
        diag = diagonal(matriz_A, n)
        if parada(diag, n, autovalores, epsilon):
            break
        elif k > 50:
            break
        i = 0
        k = k + 1
    return diag

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


#Main
matriz_A = np.array([[6, -2, -1],[-2, 6, -1], [-1, -1, 5]])
#matriz_A = np.array([[5, -2], [-2, 8]])
n = 3
autovalores, autovetores = np.linalg.eig(matriz_A)
autovetor_normalizado = np.zeros((3,3))
print(autovetores)
#Exercicio 1
epsilon = 0.000000001
avalores_calculado = calcula_autovaloresQR(matriz_A, n, autovalores, epsilon)
print("Autovalores: ")
print(autovalores)
print("Autovalores calculados: ")
print(avalores_calculado)