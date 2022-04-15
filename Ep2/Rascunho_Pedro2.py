# TCC sobre o EP: https://www.ime.usp.br/~map/tcc/2014/Cassiano%20Reinert%20Novais%20dos%20Santos.pdf

#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
import numpy as np
import matplotlib.pyplot as plt

# calculo de uij com veteorizacao
def calculaUijVetorizado(sigma, N, M, K, L, T):
    # u_j+1 = A*u_j + u_j
    # A = Δtau/Δx^2 * sigma^2 / 2 * [[-2, 1 ...],
    #                               [1, -2, 1, ...]
    #                               ...
    #                               [... 1, -2, 1]
    #                               [...     1, -2]]

    deltaTau = T/M
    deltaX = 2*L/N

    # calculo da matriz A
    const = (deltaTau / (deltaX) ** 2) * (sigma ** 2 / 2)
    A = np.zeros(shape=(N+1, N+1))

    for i in range(0, N+1):
        A[i][i] = -2 * const
        if (i < N):
            A[i + 1][i] = 1 * const
            A[i][i + 1] = 1 * const

    # definicao da matriz de u
    u = np.zeros(shape=(N+1, M+1))

    # calculo da primeira iteracao u0
    u_inicial = np.zeros(shape=(N+1, 1))
    for i in range(0, N+1):
        xi = i * deltaX - L
        u_inicial[i][0] = K * np.maximum(np.exp(xi) - 1, 0)

    u[:, [0]] = u_inicial

    # calcula demais iteracoes ate N-1
    u_atual = u_inicial
    for j in range(1, M+1):
        u_prox = np.matmul(A, u_atual) + u_atual

        tauJ = j * deltaTau
        u_prox[0][0] = 0
        u_prox[N][0] = K * np.exp(L + (sigma ** 2) * (tauJ / 2))

        # salva resultado na matriz
        u[:, [j]] = u_prox

        # atual recebe proximo para iteracao seguinte
        u_atual = u_prox

    return u

def calculaVij(u, T, M, N, r):
    V = np.zeros(shape=(N+1, N+1))
    
    for j in range(M+1):
        tauJ = j * T / M
        V[:, [j]] = u[:, [j]] * np.exp(-1 * r * tauJ)

    return V

def escolheMelhorVij(M, N, L, S, K, r, sigma, T, t):
    i_xProximo = 0
    j_tauProximo = 0

    tauAnalitico = T - t
    x_Analitico = np.log(S/K) + (r - sigma**2 / 2) * tauAnalitico

    difX = np.abs(-L - x_Analitico)
    for i in range(1, N+1):
        xi = i * (2*L/N) - L
        dif = np.abs(xi - x_Analitico)
        if dif < difX:
            difX = dif
            i_xProximo = i

    difTau = np.abs(tauAnalitico)
    for j in range(1, M+1):
        tauJ = j * T / M
        dif = np.abs(tauJ - tauAnalitico)
        if dif < difTau:
            difTau = dif
            j_tauProximo = j

    return i_xProximo, j_tauProximo

def plotagemGrafico(vetorVij, vetorS):
    plt.plot(vetorS, vetorVij, color='blue')
    plt.xlabel("USD/BRL no vencimento")
    plt.ylabel("Lucro/Prejuizo (comprador)")

    plt.savefig("grafico_comprador_cenarioFicticio.png")

    plt.show()

#Sobre o cálculo de S_t
#Calcula V_i_j pela equação do calor
#A iteração X permite obter S_t
#Constrói-se um intervalo que abarque S_t, tal como S_0, que deve possuir espaçamentos equidistantes
#A partir dos valores do intervalo de S, calculam-se valores para S_t, de modo que seja possível construir um gráfico que os relacione

#Cálculo do Prêmio (By Pedro Bacic):
#Calcula-se V(S0, t = 0) e multiplica-se pela Quantidade de Ativos = Prêmio

def calculaPremio(sigma, S0, K, N, M, L, T, r, t, quantOpcoes):
    #Declaração de variáveis

    u = calculaUijVetorizado(sigma, N, M, K, L, T)
    V = calculaVij(u, T, M, N, r)
    
    # decide qual é o melhor Vij a se retornar
    i_xProximo, j_tauProximo = escolheMelhorVij(M, N, L, S0, K, r, sigma, T, t)

    return V[i_xProximo][j_tauProximo] * quantOpcoes

def imprimeMenu():
    print("""
    ********** EP2 de MAP3122 **********
    *** Opcoes no Mercado Financeiro ***
        
    Cenarios disponiveis
    1. Cenario ficticio
    2. Cenario de cambio
    3. Cenario real        
    """)
    cenario = int(input("Escolha o cenario desejado: "))

    if cenario == 1:
        print("""
        1) Precificacao da opcao de compra
        2) Analise do lucro/prejuizo da opcao em 6 meses
        3) Analise com volatilidade a 2 por cento 
        4) Analise com volatilidade e taxa de juros a 10 por cento
        """)
        
        questao = int(input("Escolha o metodo: "))

        if questao == 1:
            premio = calculaPremio(0.01, 1, 1, 10000, 50, 10, 1, 0.01, 0, 1000)
            print("O premio calculado é de R$" + str(premio))
        elif questao == 2:
            # analise do lucro
            print("so pra nao aparecer erro rs")
        elif questao == 3:
            # analise do lucro com parametros diferentes
            print("so pra nao aparecer erro rs")
        elif questao == 4:
            # analise do lucro com parametros diferentes 2
            print("so pra nao aparecer erro rs")
        else:
            print("Metodo invalido!")

    elif cenario == 2:
        volatilidade = 0.1692
        taxaJuros = 0.1075
        valorAtual = 5.6376
        quantidadeOpcoes = 100000
        T = 3/12
        precoExecucao = 5.7

        premio = calculaPremio(volatilidade, valorAtual, precoExecucao, 10000, 50, 10, T, taxaJuros, 0, quantidadeOpcoes)

        print(premio)
    elif cenario == 3:
        # mudar dados para realidade 
        volatilidade = 0.1692
        taxaJuros = 0.1075
        valorAtual = 5.6376
        quantidadeOpcoes = 100000
        T = 3/12
        precoExecucao = 5.7


    else:   
        print("Cenario invalido!")

imprimeMenu()
