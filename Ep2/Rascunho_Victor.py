#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
import numpy as np

def itera(n, K, x, sigma, DeltaT, DeltaX, L):
    matriz_J = np.zero(shape=(2, n))
    # para i = 0, u = 0
    matriz_J[0][0] = 0
    matriz_J[1][0] = 0

    max = maximo(np.exp(x) - 1, 0)

    # para j = 0, u = K*max(e^{x_i} - 1,0)
    for i in range(0, n): matriz_J[0][i] = K * max

    # para i = n, u = Ke^{L + \sigma^2\Tau_j/2}
    matriz_J[1][n] = K * np.exp(L + sigma**2 * DeltaT)

    return matriz_J

#Elementos calculados a cada iteração
def calcula_x(DeltaX, i, L):
    return i * DeltaX - L

def calcula_tau(DeltaTau, j):
    return DeltaTau * j

def calcula_S(K, x_i, r, sigma, tau_j):
    exponencial = np.exp(x_i - (r - sigma**2/2)*tau_j)
    S = K * exponencial
    return S

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