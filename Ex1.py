#Arquivo referente ao Exercicio 1
#Resolução de sistema linear a partir do método SOR.
#Implementar critério das linhas ou o critério de Sassenfeld, para verificar se o método SOR poderá ser aplicado
#Deve ser satisfeito que 1 <= \omega < 2 (\omega == 1 => Método de Gauss-Seidel)
#Objetivo: implementar o método da potência reverso com um critério de parada, estudar a ordem de convergência e o erro assintótico numéricos

#Usar o critério de convergência do autovetor como parada

#Definições importantes:
#x_k: vetor x na k-ésima iteração.
#x_converge: autovetor associado a \lambda_1

import numpy as np

#Verifica se a condição de parada foi atingida
def calcula_parada(x_k, x_converge, epsilon):
    if np.linalg.norm(x_k - x_converge) < epsilon:
        return True
    else:
        return False


#Calcula os autovalores de uma matriz
def calcula_autovalores(matriz_A):
    autovalores_ref  =  np.linalg.eig(matriz_A)


#Calcula os autovetores de uma matriz
def calcula_autovetores(matriz_A):
    autovetores_ref = np.linalg.eig(matriz_A)







