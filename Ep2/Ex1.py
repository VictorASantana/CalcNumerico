#Main
def processamento(N, L, sigma, K, T, r, t, S_0, vetorizado):
    M = calcula_M(N, sigma, T, L)
    return 0

#Retorna o M
def calcula_M(N, sigma, T, L):
    return (T*(sigma**2))/(2*L/N)**2

#Menu
N = 10000
L = 10
print("""Menu: 
1. Opção vetorizada.
2. Opção iterativa.""")
operacao = int(input("Opcao escolhida: "))
while(True):
    print("""Cenário: 
    1. Cenário fictício.
    2. Cenário de câmbio.
    3. Cenário real.
    0. Finalizar operação.""")
    cenario = int(input("Cenario escolhido: "))
    if(cenario == 1):
        #Executa o cenario ficticio
        if(operacao == 1):
            processamento(N, L, 0.01, 1, 1, 0.01, 0.5, 1, 1)
        else:
            processamento(N, L, 0.01, 1, 1, 0.01, 0.5, 1, 0)
    elif(cenario == 2):
        #Executa o cenario de cambio
        if (operacao == 1):
            processamento(N, L, 0.01, 1, 1, 0.01, 0.5, 1, 1)
        else:
            processamento(N, L, 0.01, 1, 1, 0.01, 0.5, 1, 0)
    elif(cenario == 3):
        #Executa o cenario real
        if (operacao == 1):
            processamento(N, L, 0.01, 1, 1, 0.01, 0.5, 1, 1)
        else:
            processamento(N, L, 0.01, 1, 1, 0.01, 0.5, 1, 0)
    elif(cenario == 0):
        break
    else:
        print("Opção inválida.")




