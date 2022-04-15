#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
import numpy as np

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

def calculaOpcaoInterpolacao(sigma, S0, K, N, M, L, T, r, t, S):
    # Declaração de variáveis
    DeltaX = 2*L/N
    DeltaT = T/M
    u = uIterativo(sigma, N, M, K, L, T)
    V = calculaVij(u, T, M, N, r)

    tauAnalitico = T - t
    x_Analitico = np.log(S / K) + (r - sigma ** 2 / 2) * tauAnalitico
    # decide qual é o melhor Vij a se retornar
    i_xProximo, j_tauProximo = escolheMelhorVij(M, N, L, S0, K, r, sigma, T, t)
    xi = i_xProximo * DeltaX - L
    xi_proximo = (i_xProximo + 1) * DeltaX - L
    V_interpolacao = V_S_t(xi, xi_proximo, V[i_xProximo][j_tauProximo], x_Analitico, V[i_xProximo][j_tauProximo + 1])

    return V_interpolacao


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
opcao_inter = calculaOpcaoInterpolacao(sigma, S0, K, N, M, L, T, r, 0, S0) + K
print("A opção é precificada em R$" + str(opcao_inter))