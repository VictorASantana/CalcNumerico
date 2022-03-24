#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
import numpy as np

def itera(n, m, K, x, sigma, DeltaT, DeltaX, L):
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
