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
def calcula_SOR(matriz_A, n, vetor_b):
    #vetor_b = np.ones((n,1))
    vetor_x = np.zeros((n,1))
    soma_atual = 0
    soma_anterior = 0
    elemento_diagonal = 0
    omega = 1 + np.random.rand()
    i = 0
    k = 0
    #A iteração abaixo deve ser realizada num número limitado de vezes
    while True:
        for vetor_atual in matriz_A:
            for j in range(vetor_atual.size):
                if i > j:
                    soma_atual = soma_atual + vetor_atual[j]*vetor_x[j]
                elif i < j:
                    soma_anterior = soma_anterior + vetor_atual[j]*vetor_x[j]
                else:
                    elemento_diagonal = vetor_atual[j]
            print("soma atual: ")
            print(soma_atual)
            print("\n soma anterior: ")
            print(soma_anterior)
            vetor_x[i] = (1 - omega)*vetor_x[i] + (omega/elemento_diagonal)*(vetor_b[i] - soma_atual - soma_anterior)
            print("\n b[i]: ")
            print(vetor_b[i])
            print(vetor_x[i])
            soma_atual = 0
            soma_anterior = 0
            i = i + 1
        #vetor_b = vetor_x/np.linalg.norm(vetor_x)
        k = k + 1
        i = 0
        if k > 5:
            break
    print(vetor_x)




#Calcula a aproximação para \lambda_n^{-1} na k-ésima operação
def calcula_autovalor(x_k, x_k_mais_1):
    mu_k = np.dot(x_k.T, x_k_mais_1)/np.dot(x_k.T, x_k)
    return mu_k

#Se M < 1, então o método GS converge para a solução
def criteiro_de_Sassenfeld(matriz_A, n):
    vetor_beta = np.zeros((n,1))
    beta_atual = 0
    beta_maximo = 0
    i = 0
    for vetor_atual in matriz_A:
        for j in range(vetor_atual.size):
            if i > j:
                if vetor_atual[j] > 0:
                    beta_atual = beta_atual + vetor_atual[j]*vetor_beta[j]
                else:
                    beta_atual = beta_atual - vetor_atual[j] * vetor_beta[j]
            elif i < j:
                if vetor_atual[j] > 0:
                    beta_atual = beta_atual + vetor_atual[j]
                else:
                    beta_atual = beta_atual - vetor_atual[j]
            else:
                elemento_diagonal = vetor_atual[j]
        vetor_beta[i] = beta_atual/elemento_diagonal
        beta_atual = 0
        if beta_maximo < vetor_beta[i]:
            beta_maximo = vetor_beta[i]
        i = i + 1
    print("O valor de beta_maximo é: ")
    print(beta_maximo)

#Função main
n = 3
matriz_A = np.array([[4, 1, 1],[-2, 5, 1],[3, 1, 6]])
vetor_b = np.array([5, 0, -6.5])
calcula_SOR(matriz_A, n, vetor_b)





