# Arquivo referente ao exercicio 4
# %%
import numpy as np
import matplotlib.pyplot as plt
import Ex1 as mdp

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

def comparaAutovaloresDominantesRedes():
    mAdjRede1 = calculaMatrizAdjacencia(eps1)
    mAdjRede2 = calculaMatrizAdjacencia(eps2)

    # vetor inicial para o metodo das potencias
    vetorIni = np.random.rand(nVertices, 1)

    autovalDomRede1 = mdp.metodo_das_potencias(mAdjRede1, vetorIni)[0]
    autovalDomRede2 = mdp.metodo_das_potencias(mAdjRede2, vetorIni)[0]

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

    autovalDom, autovetDom = mdp.metodo_das_potencias(mAdjRede, vetorIni)

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


# %%
