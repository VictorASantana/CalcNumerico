# Arquivo referente ao Exercicio 1
# %%
import matplotlib.pyplot as plt
import numpy as np

# ===========================================================================
# Dados para criterio de parada do Metodo das Potencias
# ===========================================================================

itmax = 70
epsilon = np.sqrt(10) * 10**(-15)

# ===========================================================================
# Calculo dos valores e vetores de parametro para o calculo dos erros
# ===========================================================================

# calculo do autovetor associado ao autovalor dominante da matriz A
def encontra_autovetor_dominante(matriz_A):
    autovalores, autovetores_array = np.linalg.eig(matriz_A)

    autovetor_dominante = np.zeros(shape=(len(matriz_A[0]), 1))

    index = encontra_index_autovalor_dominante(autovalores)

    n = len(matriz_A[0])
    for i in range(0, n):
        autovetor_dominante[i][0] = autovetores_array[i][index]

    return autovetor_dominante

# retorna index do autovalor dominante dentro do array retornado
# pela funcao numpy.linalg.eig
def encontra_index_autovalor_dominante(autovalores):
    index_autovalor_dominante = 0
    autovalor_dominante = 0

    for i in range(autovalores.size):
        modulo_autovalor = np.abs(autovalores[i])
        if (modulo_autovalor > autovalor_dominante):
            autovalor_dominante = modulo_autovalor
            index_autovalor_dominante = i

    return index_autovalor_dominante

# retorna autovalor dominante da matriz A
def encontra_autovalor_dominante(matriz_A):
    autovalores = (np.linalg.eig(matriz_A))[0]
    autovalor_dominante = 0

    for i in range(autovalores.size):
        modulo_autovalor = np.abs(autovalores[i])
        if (modulo_autovalor > autovalor_dominante):
            autovalor_dominante = modulo_autovalor

    return autovalor_dominante

# encontra o segundo autovalor de maior módulo da matriz A
# usado para calculo dos erros assintoticos
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

# ===========================================================================
# Erro assintotico a ser usado no grafico
# ===========================================================================

def calcula_erros_assintoticos(matriz_A, n_iteracoes):
    lambda_1 = encontra_autovalor_dominante(matriz_A)
    lambda_2 = encontra_lambda_2(matriz_A)

    erros_assintoticos_por_iteracao = []

    for i in range(0, n_iteracoes):
        erro_assintotico = (lambda_2/lambda_1)**i
        erros_assintoticos_por_iteracao.append(erro_assintotico)

    return erros_assintoticos_por_iteracao

# ===========================================================================
# Calculo do autovetor e autovalor dominantes 
# a cada iteracao do Metodo das Potencias
# ===========================================================================

# calculo do vetor Xk que tende ao autovetor associado
# ao autovalor dominante da matriz A
def calcula_Xk(matrizA, Xk):

    produto = np.matmul(matrizA, Xk)
    Xk = produto / np.linalg.norm(produto)

    return Xk


# calculo do valor Uk que tende ao autovalor dominante lambda 1
def calcula_Uk(Xk, matrizA):
    Xk_tranposta = Xk.T
    produto_A_Xk = np.matmul(matrizA, Xk)

    Uk = np.matmul(Xk_tranposta, produto_A_Xk) / np.matmul(Xk_tranposta, Xk)

    return Uk[0][0]

# ===========================================================================
# Calculo dos erros a cada iteracao do Metodo das Potencias
# ===========================================================================

# calculo de erro do autovetor calculado em relacao ao teorico
# obtido pela funcao numpy.linalg.eig
def calcula_erro_autovetor(vetor_xk, autovetor_dominante):

    sub = vetor_xk - autovetor_dominante

    erro_autovetor = np.linalg.norm(sub)

    return erro_autovetor

# calculo de erro do autovalor calculado em relacao ao teorico
# obtido pela funcao numpy.linalg.eig
def calcula_erro_autovalor(autovalor_uk, autovalor_dominante):
    erro_autovalor = np.abs(autovalor_uk - autovalor_dominante)

    return erro_autovalor

# ===========================================================================
# Funcoes e constantes usadas para plotagem do grafico via MatPlotLib
# ===========================================================================

# arrays para a plotagem de graficos
y_erro_autovalor = []
y_erro_autovetor = []

# funcao que realiza a plotagem dos graficos requeridos
def plotagem_grafico_erros(matriz_A):
    fig = plt.figure(figsize=(8, 6))

    n_iteracoes = len(y_erro_autovalor)
    y_erros_assintoticos_simples = calcula_erros_assintoticos(
        matriz_A, n_iteracoes)

    y_erros_assintoticos_quadraticos = []
    for erro in y_erros_assintoticos_simples:
        erro_quadratico = erro**2
        y_erros_assintoticos_quadraticos.append(erro_quadratico)

    array_iteracoes = np.array(range(0, n_iteracoes))

    plt.plot(array_iteracoes, y_erro_autovalor,
             color='black', label="Erro autovalor")
    plt.plot(array_iteracoes, y_erro_autovetor,
             color='green', label="Erro autovetor")
    plt.plot(array_iteracoes,
             y_erros_assintoticos_simples, color='blue', label=r"$\|λ_{2}/λ_{1}\|^{k}$")
    plt.plot(array_iteracoes,
             y_erros_assintoticos_quadraticos, color='red', label=r"$\|λ_{2}/λ_{1}\|^{2k}$")
    plt.yscale("log")
    plt.xlabel("Iterações")
    plt.ylabel("Erro L2")
    plt.legend(loc="lower left")

    # para salvar a figura no diretorio atual
    #plt.savefig("grafico.png")

    plt.show()

# ===========================================================================
# Implementacao do Metodo das Potencias geral a partir das demais funcoes
# ===========================================================================

def metodo_das_potencias(matriz_A, vetor_x0):
    autovalor_Uk = 0
    autovetor_Xk = vetor_x0

    autovetor_dominante = encontra_autovetor_dominante(matriz_A)
    autovalor_dominante = encontra_autovalor_dominante(matriz_A)

    erro_autovetor = calcula_erro_autovetor(vetor_x0, autovetor_dominante)

    i = 0
    while(erro_autovetor > epsilon and i < itmax):
        autovetor_Xk = calcula_Xk(matriz_A, autovetor_Xk)
        autovalor_Uk = calcula_Uk(autovetor_Xk, matriz_A)

        erro_autovetor = calcula_erro_autovetor(
            autovetor_Xk, autovetor_dominante)
        erro_autovalor = calcula_erro_autovalor(
            autovalor_Uk, autovalor_dominante)

        y_erro_autovetor.append(erro_autovetor)
        y_erro_autovalor.append(erro_autovalor)

        i += 1

    return autovalor_Uk, autovetor_Xk

# ===========================================================================
# Funcoes para exibicao dos resultados
# ===========================================================================

def resultado_item_1():
    n = 10
    matriz_B = np.random.rand(n, n)
    vetor_x0 = np.random.rand(n, 1)
    matriz_A = np.matmul(matriz_B, matriz_B.T)

    Uk_final = metodo_das_potencias(matriz_A, vetor_x0)[0]
    plotagem_grafico_erros(matriz_A)
    print("O valor do autovalor aproximado é de: ", Uk_final)

def resultado_item_2_1():
    n = 7
    matriz_B = np.random.rand(n, n)

    # λ1 = 95 e λ2 = 92
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

    Uk_final = metodo_das_potencias(matriz_A, vetor_x0)[0]
    plotagem_grafico_erros(matriz_A)
    print("O valor do autovalor aproximado é de: ", Uk_final)

def resultado_item_2_2():
    n = 7
    matriz_B = np.random.rand(n, n)

    # λ1 = 92 e λ2 = 13
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

    Uk_final = metodo_das_potencias(matriz_A, vetor_x0)[0]
    plotagem_grafico_erros(matriz_A)
    print("O valor do autovalor aproximado é de: ", Uk_final)

# ===========================================================================
# Funcao de menu para chamada das funcoes via terminal
# ===========================================================================

def menu():
    operacao = 1

    print("Escolha o modo de operacao do algoritmo")
    print("1) A = B * Bt")
    print("""    --> Matriz B de formato 10x10
    --> Matriz B composta de coeficientes aleatorioes entre 0 e 1
    --> Bt representa a transposta de B\n""")

    print("2) A = B * D * B^(-1)")
    print("""    --> Matriz D diagonal com esta composta de coeficiente estritamente positivos distintos
    --> Matriz B de formato 7x7
    --> Matriz B composta de coeficientes aleatorioes entre 0 e 1
    --> B^(-1) representa a inversa de B
    --> λ1 e λ2 sao relativamente proximos\n""")

    print("3) A = B * D * B^(-1)")
    print("""    --> Matriz D diagonal com esta composta de coeficiente estritamente positivos distintos
    --> Matriz B de formato 7x7
    --> Matriz B composta de coeficientes aleatorioes entre 0 e 1
    --> B^(-1) representa a inversa de B
    --> λ1 e λ2 sao relativamente distantes\n""")

    operacao = int(input("Modo de operacao: "))

    if(operacao == 1):
        resultado_item_1()
    elif(operacao == 2):
        resultado_item_2_1()
    elif(operacao == 3):
        resultado_item_2_2()
    else:
        print("Modo de operacao escolhido invalido.")

# Descomentar somente para testagem do exercicio 1
#menu()

# %%
