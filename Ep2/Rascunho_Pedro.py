import numpy as np

def deltaT(M, T):
    return T/M

def calcula_u_vetorizado(deltaT, sigma, taxaJuros, L, deltaX, K, N, M):

    # definicao da matriz A
    A = np.zeros(shape=(N-1, N-1))

    alpha = sigma**2 * deltaT
    beta = taxaJuros * deltaT

    A[0][0] = 1 - alpha - beta

    for i in range(1, N-1):
        v = 0.5 * (alpha * ((i+2)-1)**2 + beta * ((i+2)-1))
        d = 1 - alpha * (i+1)**2 - beta
        l = 0.5 * (alpha * (i-1)**2 - beta * (i-1))
        
        A[i][i] = d
        A[i-1][i] = v
        A[i][i-1] = l

    u = np.zeros(shape=(N, M))

    uAtual = np.zeros(shape=(N-1, 1))

    # =================================================================================================================
    # ATENCAO: o metodo precisa das condicoes de contorno que ainda nao foram implementadas pra funcionar
    # obs: eu nao sei quais sao as condicoes de contorno pro caso do documento bolado la
    # =================================================================================================================

    for j in range(M):
        z = np.zeros(shape=(N-1, 1))
        z[0][0] = 0.5 * (alpha - beta) * uAtual[0][0]
        z[N-2][0] = 0.5 * (alpha * (N-1)**2 + beta * (N-1)) * uAtual[N][0] # condicao de contorno que eu nao sei qual é

        uProximo = np.matmul(A, uAtual) + z

        u[:, [j]] = uProximo
    
    # calculo da primeira iteracao u0 (nao sei como faz)
    u_inicial = np.zeros(shape=(N, 1))
    for i in range(0, N):
        xi = i * deltaX - L

        u_inicial[i][0] = K * np.maximum(np.exp(xi) - 1, 0)