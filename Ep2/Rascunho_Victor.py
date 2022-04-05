#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
import numpy as np

def itera(n, m, K, x, sigma, DeltaT, DeltaX, L, r, tau_j):
    V_i_j = np.zero(shape=(2, n))
    V_s_t = np.zero(shape=(2, n))
    matriz_J = np.zero(shape=(2, n))
    # para i = 0, u = 0
    matriz_J[0][0] = 0
    matriz_J[1][0] = 0

    max = maximo(np.exp(x) - 1, 0)

    # para j = 0, u = K*max(e^{x_i} - 1,0)
    for i in range(1, n): matriz_J[0][i] = K * max

    # para i = n, u = Ke^{L + \sigma^2\Tau_j/2}

    for j in range(1, m):
        for i in range(1, n-1):
            tau_j = calcula_tau(DeltaT, j)
            matriz_J[j][n] = K * np.exp(L + sigma ** 2 * tau_j/2)
            matriz_J[j+1][i] = matriz_J[j][i] + (DeltaT/DeltaX**2)*(sigma**2/2)*(matriz_J[j][i-1] - 2*matriz_J[j][i] + matriz_J[j][i+1])

    for i in range(0, m):
        for j in range(0, n):
            V_i_j[i][j] = calcula_V_i_j(matriz_J[i][j], r, tau_j)
            V_s_t[i][j] = V_S_t((calcula_x(DeltaX, i, L)),  calcula_x(DeltaX, i + 1, L), V_i_j[i][j], x, calcula_V_i_j(matriz_J[i+1][j], r, tau_j) )


    return matriz_J

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

def maximo(x, y):
    if(x >= y):
        return x
    else:
        return y

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
        u_prox[0][j] = 0
        u_prox[N][0] = K * np.exp(L + sigma ** 2 * tauJ / 2)

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