#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
import numpy as np

def itera(n, m, K, x, sigma, DeltaT, DeltaX, L, r, S0, T, t, tau_t):
    V_i_j = np.zeros(shape=(n+1, m+1))
    S_i_j = np.zeros(shape=(n+1, m+1))
    x_i = np.zeros(shape=(n+1, 1))
    intervalo = 0
    maxvalor = 0
    minvalor = 0
    nintervalos = 10
    matriz_J = np.zeros(shape=(n+1, m+1))
    # para i = 0, u = 0
    matriz_J[0][0] = 0
    matriz_J[1][0] = 0
    tau_j = np.zeros(shape=(m+1, 1))
    max = maximo(np.exp(x) - 1, 0)
    tau_escolhido = 0
    x_escolhido = 0
    diferenca_tau = 0
    diferenca_x = 0
    i_escolhido = 0
    j_escolhido = 0

    # para j = 0, u = K*max(e^{x_i} - 1,0)
    for j in range(1, m+1):
        matriz_J[0][j] = K * max

    # para i = n, u = Ke^{L + \sigma^2\Tau_j/2}

    for i in range(1, n-1):
        for j in range(1, m):
            tau_j[j] = calcula_tau(DeltaT, j)
            matriz_J[n][j] = K * np.exp(L + sigma ** 2 * tau_j[j]/2)
            matriz_J[i][j+1] = matriz_J[i][j] + (DeltaT/DeltaX**2)*(sigma**2/2)*(matriz_J[i-1][j] - 2*matriz_J[i][j] + matriz_J[i+1][j])

    for i in range(0, n+1):
        x_i[i] = calcula_x(DeltaX, i, L)
        if x - x_i[i] < diferenca_x:
            diferenca_x = x - x_i[i]
            x_escolhido = x_i[i]
            i_escolhido = i
        for j in range(0, m+1):
            tau_j[j] = calcula_tau(DeltaT, j)
            V_i_j[i][j] = calcula_V_i_j(matriz_J[i][j], r, tau_j[j])
            print(matriz_J[i][j])
            S_i_j[i][j] = K + V_i_j[i][j]
            if tau_t - tau_j[j] < diferenca_tau:
                diferenca_tau = tau_j[j] - tau_t
                tau_escolhido = tau_j[j]
                j_escolhido = j
            if S_i_j[i][j] > maxvalor:
                maxvalor = S_i_j[i][j]
            elif S_i_j[i][j] < minvalor:
                minvalor = S_i_j[i][j]
            if i == 0:
                x_i[j] = calcula_x(DeltaX, L, j)

    #Condicionais do S0:
    if S0 > maxvalor:
        maxvalor = S0
    elif S0 < minvalor:
        minvalor = S0

    intervalo = (maxvalor - minvalor)/nintervalos

    return V_i_j[i_escolhido][j_escolhido]




#Elementos calculados a cada iteração
def calcula_x(DeltaX, i, L):
    return i * DeltaX - L

#Calculo de \tau(j)
def calcula_tau(DeltaTau, j):
    return DeltaTau * j

#Calculo iterativo de S
def calcula_S(K, x_i, r, sigma, tau_j):
    exponencial = np.exp(x_i - (r - sigma**2/2)*tau_j)
    S = K * exponencial
    return S

#Calculo iterativo de V
def calcula_V_i_j(u_i_j, r, tau_j):
    return u_i_j * np.exp(- 1 * (r * tau_j))

#Cálculo dos Deltas
def calcula_DeltaX(L, N):
    return 2*L/N

def calcula_DeltaTau(T, M):
    return T/M

#Calculo de \tau(t)
def tau_t(T, t):
    return T - t

#Calculo de x(S,t)
def x_ideal(S, K, r, sigma, tau_t):
    return np.log(S/K) + (r - sigma**2/2)*tau_t

#Calculo de V(S,t)
def V_S_t(x_i, x_proximo, V_i_j, x_ideal, V_proximo_j):
    return ((x_proximo - x_ideal)*V_i_j - (x_i - x_ideal)*V_proximo_j)/(x_proximo - x_i)


# calculo de uij com veteorizacao
def calcula_u_i_j_vetorizado(tau, x, sigma, N, M, deltaTau, deltaX, K, L, tauJ):
    # u_j+1 = A*u_j + u_j
    # A = Δtau/Δx^2 * sigma^2 / 2 * [[-2, 1 ...],
    #                               [1, -2, 1, ...]
    #                               ...
    #                               [... 1, -2, 1]
    #                               [...     1, -2]]

    # calculo da matriz A
    const = (deltaTau / (deltaX) ** 2) * (sigma ** 2 / 2)
    A = np.zeros(shape=(N, N))

    for i in range(0, N):
        A[i][i] = -2 * const
        if (i < N - 1):
            A[i + 1][i] = 1 * const
            A[i][i + 1] = 1 * const

    # definicao da matriz de u
    u = np.zeros(shape=(N, M))

    # calculo da primeira iteracao u0
    u_inicial = np.zeros(shape=(N, 1))
    for i in range(0, N):
        xi = i * deltaX - L
        u_inicial[i][0] = K * maximo(np.exp(xi) - 1, 0)

    u[:, [0]] = u_inicial

    # calcula demais iteracoes ate N-1
    u_atual = u_inicial
    for j in range(1, M - 1):
        u_prox = np.matmul(A, u_atual) + u_atual
        u_prox[j][0] = 0
        u_prox[N-1][0] = K * np.exp(L + sigma ** 2 * tauJ / 2)

        # salva resultado na matriz
        u[:, [j]] = u_prox

        # atual recebe proximo para iteracao seguinte
        u_atual = u_prox

    return u



#Sobre o cálculo de S_t
#Calcula V_i_j pela equação do calor
#A iteração X permite obter S_t
#Constrói-se um intervalo que abarque S_t, tal como S_0, que deve possuir espaçamentos equidistantes
#A partir dos valores do intervalo de S, calculam-se valores para S_t, de modo que seja possível construir um gráfico que os relacione

#Cálculo do Prêmio (By Pedro Bacic):
#Calcula-se V(S0, t = 0) e multiplica-se pela Quantidade de Ativos = Prêmio

def uIterativo(sigma, n, m, K, L, T):
    #Declaração de variáveis
    DeltaT = T/m
    DeltaX = 2*L/n
    u_i_j = np.zeros(shape=(n+1, m+1))
    const = (DeltaT / (DeltaX) ** 2) * (sigma ** 2 / 2)

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


def calculaVij(u, T, M, N, r):
    V = np.zeros(shape=(N + 1, N + 1))

    for j in range(M + 1):
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

def calculaOpcao(sigma, S0, K, N, M, L, T, r, t):
    # Declaração de variáveis

    u = uIterativo(sigma, N, M, K, L, T)
    V = calculaVij(u, T, M, N, r)

    # decide qual é o melhor Vij a se retornar
    i_xProximo, j_tauProximo = escolheMelhorVij(M, N, L, S0, K, r, sigma, T, t)

    return V[i_xProximo][j_tauProximo]



#Testes

#Teste Exercicio 4.1----------------------------------------------------------------------------------------------------------------------------------------

#Parâmetros:
#T = 1, r = 0.01, sigma = 0.01, K = 1, S0 = 1, Quantidade de ativos: 1000

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#4.1.1) Precificar a opção de compra
M = 50
N = 10000
L = 10
T = 1
sigma = 0.01
r = 0.01
qAtivos = 1000
K = 1
S0 = 1
tau = tau_t(T, 0)
x_t = x_ideal(S0, K, r, sigma, tau)
#Calcular os deltas:
Delta_tau = calcula_DeltaTau(T, M)
Delta_X = calcula_DeltaX(L, N)
#Calculando o Vij:
opcao = calculaOpcao(sigma, S0, K, N, M, L, T, r, 0) + K
print("A opção é precificada em R$" + str(opcao))
