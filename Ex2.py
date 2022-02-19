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
#%%
import numpy as np
import matplotlib.pyplot as plt

#A função calcula parada é usada para verificação do critério de parada, que é atingido quando ||x_k - x*|| < epsilon
def calcula_parada(x_k, x_converge, epsilon):
    if np.linalg.norm(x_k - x_converge) < epsilon:
        return True
    else:
        return False

def calcula_autovetor_dominante(matriz_A):
    autovalores, autovetores_array = np.linalg.eig(matriz_A)
    autovetor_dominante = np.zeros(shape=(n, 1))
    index = encontra_index_autovalor_dominante(autovalores)
    for i in range(0, n):
        autovetor_dominante[i][0] = autovetores_array[i][index]
    print(autovetor_dominante)
    return autovetor_dominante

def calcula_autovalor_dominante(matriz_A):
    autovalores = (np.linalg.eig(matriz_A))[0]
    autovalor_dominante = 0
    for i in range(autovalores.size):
        modulo_autovalor = np.abs(autovalores[i])
        if (modulo_autovalor > autovalor_dominante):
            autovalor_dominante = modulo_autovalor
    return autovalor_dominante

def encontra_index_autovalor_dominante(autovalores):
    index_autovalor_dominante = 0
    autovalor_dominante = 0
    for i in range(autovalores.size):
        modulo_autovalor = np.abs(autovalores[i])
        if (modulo_autovalor > autovalor_dominante):
            autovalor_dominante = modulo_autovalor
            index_autovalor_dominante = i
    return index_autovalor_dominante

#Calcula o sistema (3) - método SOR
def calcula_SOR(matriz_A, n, epsilon, num_interacoes):
    vetor_b = np.ones((n,1))
    vetor_x = np.empty((n,1))
    soma_atual = 0
    soma_anterior = 0
    elemento_diagonal = 0
    omega = 1
    erro_valor = 0
    erro_vetor = 0
    autovalor_dominante = calcula_autovalor_dominante(np.linalg.inv(matriz_A))
    autovetor_dominante = calcula_autovetor_dominante(np.linalg.inv(matriz_A))
    print("Dominante: ")
    print(autovetor_dominante)
    i = 0
    global k
    #A iteração abaixo deve ser realizada num número limitado em num_interações vezes
    while True:
        for vetor_atual in matriz_A:
            for j in range(vetor_atual.size):
                if i > j:
                    soma_atual = soma_atual + vetor_atual[j]*vetor_x[j]
                elif i < j:
                    soma_anterior = soma_anterior + vetor_atual[j]*vetor_x[j]
                else:
                    elemento_diagonal = vetor_atual[j]
            vetor_x[i] = (1 - omega)*vetor_b[i] + omega*(vetor_b[i] - soma_atual - soma_anterior)/elemento_diagonal
            soma_atual = 0
            soma_anterior = 0
            i = i + 1
        mu_k = calcula_autovalor(vetor_b, vetor_x)
        erro_valor = erro_autovalor(mu_k, autovalor_dominante)
        erros_autovalor.append(erro_valor)
        erro_vetor = erro_autovetor((vetor_x / np.linalg.norm(vetor_x)), autovetor_dominante)
        erros_autovetor.append(erro_vetor)
        if calcula_parada(vetor_x, vetor_b, epsilon):
            print("Calculado: ")
            print(vetor_x / np.linalg.norm(vetor_x))
            break
        if k > num_interacoes:
            print("Calculado: ")
            print(vetor_x/np.linalg.norm(vetor_x))
            break
        for l in range(vetor_x.size):
            vetor_b[l] = vetor_x[l] / np.linalg.norm(vetor_x)
        k = k + 1
        i = 0
    return mu_k

#Calcula a aproximação para \lambda_n^{-1} na k-ésima operação
def calcula_autovalor(x_k, x_k_mais_1):
    x_k_transposta = x_k.T
    mu_k = np.matmul(x_k_transposta, x_k_mais_1)/np.matmul(x_k_transposta, x_k)
    return mu_k[0][0]

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
    if beta_maximo < 1:
        print("O método SOR converge pelo critério de Sassenfeld")
        return True
    else:
        print("O critério de Sassenfeld não permite afirmar se o método SOR converge ou não.")
        return False

#Calcula o erro do autovalor a cada iteração
def erro_autovalor(autovalor_uk, autovalor_dominante):
    erro_autovalor = np.abs(autovalor_uk - autovalor_dominante)

    return erro_autovalor

#Calcula o erro do autovetor a cada iteração
def erro_autovetor(vetor_xk, autovetor_dominante):

    modulo_autovetor_dominante = np.zeros(shape=(n, 1))
    modulo_xk = np.zeros(shape=(n, 1))

    for i in range(0, n):
        modulo_autovetor_dominante[i][0] = np.abs(autovetor_dominante[i][0])
        modulo_xk[i][0] = np.abs(vetor_xk[i][0])
    print("Valor dominante: ")
    print(modulo_autovetor_dominante)
    print("Valor comparado sem modulo: ")
    print(vetor_xk)
    print("Valor comparado com modulo: ")
    print(modulo_xk)
    sub = modulo_xk - modulo_autovetor_dominante

    erro_autovetor = np.linalg.norm(sub)

    return erro_autovetor

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

def calcula_erros_assintoticos(matriz_A, n_iteracoes):
    lambda_1 = calcula_autovalor_dominante(matriz_A)
    lambda_2 = encontra_lambda_2(matriz_A)
    erros_assintoticos_por_iteracao = []
    for i in range(0, n_iteracoes):
        erro_assintotico = (lambda_2/lambda_1)**i
        erros_assintoticos_por_iteracao.append(erro_assintotico)
    return erros_assintoticos_por_iteracao

def plotagem_grafico_erros(matriz_A):
    fig = plt.figure(figsize=(8, 6))

    y_erros_assintoticos_simples = calcula_erros_assintoticos(
        matriz_A, (k+1))

    y_erros_assintoticos_quadraticos = []
    for erro in y_erros_assintoticos_simples:
        erro_quadratico = erro**2
        y_erros_assintoticos_quadraticos.append(erro_quadratico)

    array_iteracoes = np.array(range(0, (k+1)))

    plt.plot(array_iteracoes, erros_autovalor,
             color='black', label="Erro autovalor")
    plt.plot(array_iteracoes, erros_autovetor,
             color='green', label="Erro autovetor")
    plt.plot(array_iteracoes,
             y_erros_assintoticos_simples, color='blue', label=r"$\|λ_{2}/λ_{1}\|^{k}$")
    plt.plot(array_iteracoes,
             y_erros_assintoticos_quadraticos, color='red', label=r"$\|λ_{2}/λ_{1}\|^{2k}$")
    plt.yscale("log")
    plt.xlabel("Iterações")
    plt.ylabel("Erro L2")
    plt.legend(loc="lower left")

    # print(y_erro_autovalor)

    plt.show()


#Função main
epsilon = 10**(-15)
n = 7
k = 0
n_iteracoes = 0
B = np.random.rand(n, n)
identidade = np.identity(n)
matriz_A = B + B.T + n*identidade
erros_autovalor = []
erros_autovetor = []
print(matriz_A)
print(np.linalg.eig(matriz_A))
print(np.linalg.eig(np.linalg.inv(matriz_A)))
print("\n Calculado: \n")
if criteiro_de_Sassenfeld(matriz_A, n):
    n_iteracoes = 50
    mu = calcula_SOR(matriz_A, n, epsilon, n_iteracoes)
else:
    n_iteracoes = 70
    mu = calcula_SOR(matriz_A, n, epsilon, n_iteracoes)
print("\nautovalor dominante: ")
print(mu)
plotagem_grafico_erros(np.linalg.inv(matriz_A))


# matriz exercicio 1.a
    # n = 10
    # matriz_B = np.random.rand(n, n)
    # vetor_x0 = np.random.rand(n, 1)
    # matriz_A = np.matmul(matriz_B, matriz_B.T)

# matriz exercicio 1.b.i
    # n = 7
    # matriz_B = np.random.rand(n, n)

    # # λ1 = 95 e λ2 = 92
    # matriz_D = np.array(
    #     [[2, 0, 0, 0, 0, 0, 0],
    #     [0, 13, 0, 0, 0, 0, 0],
    #     [0, 0, 92, 0, 0, 0, 0],
    #     [0, 0, 0, 32, 0, 0, 0],
    #     [0, 0, 0, 0, 76, 0, 0],
    #     [0, 0, 0, 0, 0, 95, 0],
    #     [0, 0, 0, 0, 0, 0, 22]]
    # )
    
    # vetor_x0 = np.random.rand(n, 1)

    # matriz_interm = np.matmul(matriz_B, matriz_D)

    # inv_matriz_B = np.linalg.inv(matriz_B)
    # matriz_A = np.matmul(matriz_interm, inv_matriz_B)

# matriz exercicio 1.b.ii
    # n = 7
    # matriz_B = np.random.rand(n, n)

    # # λ1 = 92 e λ2 = 13
    # matriz_D = np.array(
    #     [[2, 0, 0, 0, 0, 0, 0],
    #     [0, 13, 0, 0, 0, 0, 0],
    #     [0, 0, 92, 0, 0, 0, 0],
    #     [0, 0, 0, 10, 0, 0, 0],
    #     [0, 0, 0, 0, 1, 0, 0],
    #     [0, 0, 0, 0, 0, 8, 0],
    #     [0, 0, 0, 0, 0, 0, 11]]
    # )
    
    # vetor_x0 = np.random.rand(n, 1)

    # matriz_interm = np.matmul(matriz_B, matriz_D)

    # inv_matriz_B = np.linalg.inv(matriz_B)
    # matriz_A = np.matmul(matriz_interm, inv_matriz_B)

#%%
