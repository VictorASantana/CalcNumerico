#Condições de contorno:
#u_{i}^{0} = Kmax(e^{x_i} - 1,0)
#u_{0}^{j} = 0
#u_{N}^{j} = Ke^{L + \sigma^2\Tau_j/2}
import numpy as np

def itera(n, K, x):
    matriz_J = np.zero(shape=(2, n))
    #inicialização:
    i = 0
    matriz_J[0][0] = 0
    matriz_J[1][0] = 0
    maximo = maxi(x)
    while(i < n):
        matriz_J[0][i] = K * maximo
        i = i + 1
    i = 0

    return matriz_J

#Elementos calculados a cada iteração
def calcula_x(DeltaX, i, L):
    x_i = 0
    return x_i

def calcula_S(K, x_i, r, sigma, tau_j):
    S = 0
    return S

def calcula_V_i_j(u_i_j, r, tau_j):
    return u_i_j * np.exp(- 1 * (r * tau_j))


def maxi(x):
    return x

