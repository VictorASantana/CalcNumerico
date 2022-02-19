# Arquivo referente ao exercicio 4
# ***************************************************************************
# * IMPORTANTE: para funcionamento, o arquivo utiliza a funcao              *
# * metodo_das_potencias presente no arquivo Ex1.py e este deve estar       *
# * no mesmo diretorio que o arquivo Ex4.py                                 *
# ***************************************************************************
# %%
import numpy as np
import matplotlib.pyplot as plt
from Ex1 import metodo_das_potencias

# ===========================================================================
# Exercicio 4.1
# Definicao e plotagem das redes G1 e G2
# ===========================================================================

#         Rede G1
#
#   0  1   2   3  4   5
# 5 a
# 4    a   b   c
# 3       abc     d
# 2    c   b   a  d
# 1 c      b      ad
# 0 c      b      d   d
#

n = 6

def plotagemRede1():
    # cada linha das matrizes representa
    # uma linha da rede desenvolvida
    x = np.array([[1, 2, 3, 4, 5],
                  [3, 3, 3, 3, 3],
                  [1, 1, 2, 3, 4],
                  [5, 5, 5, 5, 6]])
    y = np.array([[5, 4, 3, 2, 1],
                  [4, 3, 2, 1, 0],
                  [0, 1, 2, 3, 4],
                  [3, 2, 1, 0, 0]])
    # labels de cada estacao
    v = [[r"$v_{0}$", r"$v_{1}$", r"$v_{4}$", r"$v_{8}$", r"$v_{12}$"],
         [r"$v_{2}$", '', r"$v_{7}$", r"$v_{11}$", r"$v_{14}$"],
         [r"$v_{13}$", r"$v_{10}$", r"$v_{6}$", '', r"$v_{3}$"],
         [r"$v_{5}$", r"$v_{9}$", '', r"$v_{15}$", r"$v_{16}$"]]

    for linha in range(0, 4):
        plt.scatter(x[linha], y[linha], marker='o')
        plt.plot(x[linha], y[linha])

        # atribui as labels
        for i in range(0, 5):
            plt.text(x[linha][i] + 0.05, y[linha][i] + 0.2, v[linha][i])

    # especificacoes do grafico
    plt.title(r"Rede $G_{1}$")
    plt.grid(visible=True, axis='both', alpha=0.2)
    plt.figure(figsize=(10, 10))

    plt.savefig("rede1.png")

    plt.show()


#         Rede G2
#
#   0  1   2   3   4  5
# 5 b  b   b
# 4        bc  ab
# 3    c       a
# 2    c       a
# 1    c           a
# 0  d cd  d   d   d  a
#

def plotagemRede2():
    # cada linha das matrizes representa
    # uma linha da rede desenvolvida
    x = np.array([[4, 4, 4, 5, 6],
                  [1, 2, 3, 3, 4],
                  [3, 2, 2, 2, 2],
                  [1, 2, 3, 4, 5]])
    y = np.array([[4, 3, 2, 1, 0],
                  [5, 5, 5, 4, 4],
                  [4, 3, 2, 1, 0],
                  [0, 0, 0, 0, 0]])
    # labels de cada estacao
    v = [['', r"$v_{6}$", r"$v_{8}$", r"$v_{10}$", r"$v_{16}$"],
         [r"$v_{0}$", r"$v_{1}$", r"$v_{2}$", r"$v_{3}$", r"$v_{4}$"],
         ['', r"$v_{5}$", r"$v_{7}$", r"$v_{9}$", r"$v_{12}$"],
         [r"$v_{11}$", '', r"$v_{13}$", r"$v_{14}$", r"$v_{15}$"]]

    for linha in range(0, 4):
        plt.scatter(x[linha], y[linha], marker='o')
        plt.plot(x[linha], y[linha])

        # atribui as labels
        for i in range(0, 5):
            plt.text(x[linha][i] + 0.05, y[linha][i] + 0.2, v[linha][i])

    # especificacoes do grafico
    plt.title(r"Rede $G_{2}$")
    plt.grid(visible=True, axis='both', alpha=0.2)
    plt.figure(figsize=(10, 10))

    plt.savefig("rede2.png")

    plt.show()


# ===========================================================================
# Exercicio 4.2
# Obtencao das matrizes de adjacencia
# Definicao das matrizes de arestas eps1 e eps2
# ===========================================================================

# matriz de arestas que representa a rede de metro 1
eps1 = np.array(
    [[0, 1],
     [1, 4],
     [4, 8],
     [8, 12],
     [2, 4],
     [4, 7],
     [7, 11],
     [7, 14],
     [3, 4],
     [4, 6],
     [6, 10],
     [10, 13],
     [5, 9],
     [9, 12],
     [12, 15],
     [15, 16]]
)

# matriz de arestas que representa
# a rede de metro 2
eps2 = np.array(
    [[4, 6],
     [6, 8],
     [8, 10],
     [10, 16],
     [0, 1],
     [1, 2],
     [2, 3],
     [3, 4],
     [3, 5],
     [5, 7],
     [7, 9],
     [9, 12],
     [11, 12],
     [12, 13],
     [13, 14],
     [14, 15]]
)

# nVertices = nEstacoes/linha * linha - nConexoes
nVertices = 17

def calculaMatrizAdjacencia(matrizArestas):
    matrizAdjacencia = np.zeros(shape=(nVertices, nVertices), dtype=int)

    for aresta in matrizArestas:
        estacao1 = aresta[0]
        estacao2 = aresta[1]

        matrizAdjacencia[estacao1][estacao2] = 1
        matrizAdjacencia[estacao2][estacao1] = 1

    return matrizAdjacencia


# ===========================================================================
# Exercicio 4.3
# Calculo dos indices λ1 de cada rede
# Aplicação do método das potências
# ===========================================================================s

def comparaAutovaloresDominantesRedes(mArest1, mArest2):
    mAdjRede1 = calculaMatrizAdjacencia(mArest1)
    mAdjRede2 = calculaMatrizAdjacencia(mArest2)

    # vetor inicial para o metodo das potencias
    vetorIni = np.random.rand(nVertices, 1)

    autovalDomRede1 = metodo_das_potencias(mAdjRede1, vetorIni)[0]
    autovalDomRede2 = metodo_das_potencias(mAdjRede2, vetorIni)[0]

    print("λ(G1) = " + str(autovalDomRede1))
    print("λ(G2) = " + str(autovalDomRede2))

    if(autovalDomRede1 > autovalDomRede2):
        print("λ(G1) " + "maior que " + "λ(G2)")
    elif(autovalDomRede2 > autovalDomRede1):
        print("λ(G2) " + "maior que " + "λ(G1)")
    else:
        print("λ(G1) " + "igual a " + "λ(G2)")


# ===========================================================================
# Exercicio 4.4
# Calculo dos autovetores associados aos λ1 de cada rede
# Determinacao do vertice com maior centralidade
# ===========================================================================

def calculaVerticeMaiorCentralidade(matrizArestas):
    mAdjRede = calculaMatrizAdjacencia(matrizArestas)

    # vetor inicial para o metodo das potencias
    vetorIni = np.random.rand(nVertices, 1)

    autovalDom, autovetDom = metodo_das_potencias(mAdjRede, vetorIni)

    # vetor com centralidades de autovetor para cada vertice
    vetorCentralidades = np.zeros(shape=(nVertices, 1))

    # aplicacao da formula para centralidades de cada vertice
    for i in range(0, nVertices):
        somatoriaLinha = 0
        for j in range(0, nVertices):
            somatoriaLinha += mAdjRede[i][j] * autovetDom[j][0]
        vetorCentralidades[i][0] = somatoriaLinha/autovalDom

    # obtencao do vertice com maior centralidade
    maiorCentralidade = vetorCentralidades[0][0]
    verticeMaiorCentralidade = 1
    for i in range(1, nVertices):
        centralidadeAtual = vetorCentralidades[i][0]
        if(centralidadeAtual > maiorCentralidade):
            maiorCentralidade = centralidadeAtual
            verticeMaiorCentralidade = i+1

    return verticeMaiorCentralidade


# ===========================================================================
# Funcao de menu para chamada das funcoes via terminal
# ===========================================================================

def menu():
    operacao = 1

    print("Escolha o modo de operacao do algoritmo")
    print("1) Plotagem da Rede 1")
    print("2) Plotagem da Rede 2")
    print("3) Imprime matriz de arestas da Rede 1")
    print("4) Imprime matriz de arestas da Rede 2")
    print("5) Comparacao de autovalores dominantes em cada rede")
    print("6) Determina vertice de maior centralidade da Rede 1")
    print("7) Determina vertice de maior centralidade da Rede 2")

    operacao = int(input("Modo de operacao: "))

    if operacao == 1:
        plotagemRede1()
    elif operacao == 2:
        plotagemRede2()
    elif operacao == 3:
        print(eps1)
    elif operacao == 4:
        print(eps2)
    elif operacao == 5:
        comparaAutovaloresDominantesRedes(eps1, eps2)
    elif operacao == 6:
        vertice = calculaVerticeMaiorCentralidade(eps1)
        print("Vertice de maior centralidade de G1: " + str(vertice))
    elif operacao == 7:
        vertice = calculaVerticeMaiorCentralidade(eps2)
        print("Vertice de maior centralidade de G1: " + str(vertice))
    else:
        print("Modo de operacao invalido")

# Descomentar para testagem do algoritmo
menu()

# %%
