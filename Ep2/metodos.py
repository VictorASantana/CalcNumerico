# TCC sobre o EP: https://www.ime.usp.br/~map/tcc/2014/Cassiano%20Reinert%20Novais%20dos%20Santos.pdf

#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
# %%
import numpy as np
import matplotlib.pyplot as plt

# ================================================================================
# Calculo de Uij 
# ================================================================================

def calculaUijVetorizado(sigma, N, M, K, L, T):
    u = np.zeros(shape=(1, N+1))
    u_prox = np.zeros(shape=(1, N+1))
    deltaX = 2 * L / N
    deltaTau = T / M

    for i in range(N+1):
        xi = i * deltaX - L
        u[0][i] = K * np.maximum(np.exp(xi)-1, 0)

    const = (deltaTau / (deltaX) ** 2) * (sigma ** 2 / 2)
    u_atual = u
    
    for j in range(1, M+1):
        uij_prox = np.roll(u_atual, -1)
        uij_prox[0][N] = 0
        uij_ant = np.roll(u_atual, 1)
        uij_ant[0][N] = 0
        u_prox = u_atual + const * (uij_ant - 2*u_atual + uij_prox)
        u = np.concatenate((u, u_prox), axis=0)
        u_atual = u_prox

    return u.T

# calculo de uij iterativamente
def uIterativo(sigma, n, m, K, L, T):
    #Declaração de variáveis
    DeltaT = T/m
    DeltaX = 2*L/n
    u_i_j = np.zeros(shape=(n+1, m+1))

    for i in range (0, n+1):
        xi = i * DeltaX - L
        for j in range(0, m+1):
            tauJ = j * DeltaT
            if i == 0:
                u_i_j[i][j] = 0
            elif i == n:
                u_i_j[i][j] =  K * np.exp(L + (sigma ** 2) * (tauJ / 2))
            elif j == 0:
                u_i_j[i][j] = K * np.maximum(np.exp(xi) - 1, 0)
                u_i_j[i][j + 1] = u_i_j[i][j] + ((DeltaT/DeltaX) * (sigma**2/2))*(u_i_j[i-1][j] - 2*u_i_j[i][j] + u_i_j[i+1][j])
            else:
                if j < m:
                    u_i_j[i][j+1] = u_i_j[i][j] + ((DeltaT/DeltaX) * (sigma**2/2))*(u_i_j[i-1][j] - 2*u_i_j[i][j] + u_i_j[i+1][j])

    return u_i_j

# ================================================================================
# Obtencao de Vij 
# ================================================================================

# calcula matriz Vij
def calculaVij(u, T, M, N, r):
    V = np.zeros(shape=(N+1, N+1))
    
    for j in range(M+1):
        tauJ = j * T / M
        V[:, [j]] = u[:, [j]] * np.exp(-1 * r * tauJ)

    return V

# calcula melhor valor para Vij de acordo com os
# valores de xi e tauJ mais proximos das variaveis analiticas
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

# ================================================================================
# Plotagem do gráfico para Análise de Lucro e Prejuízo
# ================================================================================

# realiza a plotagem do grafico e salva no arquivo analiseLucroPrejuizo.png
def plotagemGraficoLucroPrejuizo(vetorLucroPrejuizo, vetorS):
    x = vetorS[0]
    y = vetorLucroPrejuizo[0]

    plt.figure()

    plt.plot(x, y)

    plt.xlabel("USD/BRL no vencimento")
    plt.ylabel("Lucro/Prejuizo (comprador)")

    plt.savefig("analiseLucroPrejuizo.png")

    plt.show()

# gera intervalos de S para plotagem do grafico
def geraIntervaloS(S0):
    deltaS = 0.05 * S0
    Smin = 0.5* S0
    nIteracoes = 25

    vetorS = np.zeros(shape=(1, nIteracoes))

    for i in range (nIteracoes):
        vetorS[0][i] = Smin + i * deltaS

    return vetorS

# gera vetor de Vij que servira de base para o eixo y 
# do grafico de analise de lucro e prejuizo
def geraVetorVij(vetorS, sigma, N, M, L, K, T, r, t, vetorizado):
    tamanhoVetor = vetorS.shape[1]
    vetorVij = np.zeros(shape=(1, tamanhoVetor))

    for i in range(tamanhoVetor):
        vetorVij[0][i] = calculaOpcao(sigma, vetorS[0][i], K, N, M, L, T, r, t, vetorizado)
        print("Carregando...")

    return vetorVij

# calcula vetor que sera o eixo y do grafico 
# de analise de lucro e prejuizo
def geraVetorLucroPrejuizo(vetorVij, quantidadeOpcoes, premio):
    tamanhoVetor = vetorVij.shape[1]
    vetorLucroPrejuizo = np.zeros(shape=(1, tamanhoVetor))

    for i in range(tamanhoVetor):
        vetorLucroPrejuizo[0][i] = quantidadeOpcoes * vetorVij[0][i] - premio

    return vetorLucroPrejuizo

# calculo de V(S,t) com interpolacao
def V_S_t(x_i, x_proximo, V_i_j, x_ideal, V_proximo_j):
    return ((x_proximo - x_ideal)*V_i_j - (x_i - x_ideal)*V_proximo_j)/(x_proximo - x_i)

# realiza chamadas das funcoes necessarias para analise do Lucro e Prejuizo
def analiseLucroPrejuizo(sigma, N, M, L, K, T, r, t, S0, quantidadeOpcoes, vetorizado):
    vetorS = geraIntervaloS(S0)
    vetorVij = geraVetorVij(vetorS, sigma, N, M, L, K, T, r, t, vetorizado)
    premio = calculaPremio(sigma, S0, K, N, M, L, T, r, t, quantidadeOpcoes, vetorizado)
    vetorLucroPrejuizo = geraVetorLucroPrejuizo(vetorVij, quantidadeOpcoes, premio)

    plotagemGraficoLucroPrejuizo(vetorLucroPrejuizo, vetorS)

# calcula o premio com base no Vij associado
def calculaPremio(sigma, S0, K, N, M, L, T, r, t, quantOpcoes, vetorizado):
    melhorVij = calculaOpcao(sigma, S0, K, N, M, L, T, r, t, vetorizado)

    return melhorVij * quantOpcoes

# base para precificacao de opcao
def calculaOpcao(sigma, S, K, N, M, L, T, r, t, vetorizado):
    # Declaração de variáveis
    if vetorizado == 1:
        u = calculaUijVetorizado(sigma, N, M, K, L, T)
    else:
        u = uIterativo(sigma, N, M, K, L, T)

    V = calculaVij(u, T, M, N, r)

    # decide qual é o melhor Vij a se retornar
    i_xProximo, j_tauProximo = escolheMelhorVij(M, N, L, S, K, r, sigma, T, t)

    return V[i_xProximo][j_tauProximo]

# considera a interpolação para um resultado mais preciso
def calculaOpcaoInterpolacao(sigma, S, K, N, M, L, T, r, t, vetorizado):
    # Declaração de variáveis
    DeltaX = 2*L/N

    if vetorizado == 1:
        u = calculaUijVetorizado(sigma, N, M, K, L, T)
    else:
        u = uIterativo(sigma, N, M, K, L, T)

    V = calculaVij(u, T, M, N, r)

    tauAnalitico = T - t
    x_Analitico = np.log(S / K) + (r - sigma ** 2 / 2) * tauAnalitico
    
    # decide qual é o melhor Vij a se retornar
    i_xProximo, j_tauProximo = escolheMelhorVij(M, N, L, S, K, r, sigma, T, t)

    xi = i_xProximo * DeltaX - L
    xi_proximo = (i_xProximo + 1) * DeltaX - L

    if(xi > x_Analitico):
        xi_proximo = xi
        xi = (i_xProximo - 1) * DeltaX - L

    V_interpolacao = V_S_t(xi, xi_proximo, V[i_xProximo][j_tauProximo], x_Analitico, V[i_xProximo][j_tauProximo + 1])

    return V_interpolacao

# %%
