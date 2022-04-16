import metodos as mt

def imprimeMenu():

    print("""
    ********** EP2 de MAP3122 **********
    *** Opcoes no Mercado Financeiro ***
        
    Cenarios disponiveis
    1. Cenario ficticio
    2. Cenario de cambio
    3. Cenario real        
    """)

    cenario = int(input("Escolha o cenario desejado: "))
    print("""
    Métodos disponíveis:
    1. Método vetorizado
    2. Método iterativo
    """)
    metodo = int(input("Escolha o metodo desejado: "))

    M = 100
    N = 10000
    L = 10

    if cenario == 1:
        print("""
        1) Precificacao da opcao de compra
        2) Analise do lucro/prejuizo da opcao em 6 meses
        3) Analise com volatilidade a 2 por cento 
        4) Analise com volatilidade e taxa de juros a 10 por cento
        """)
        
        questao = int(input("Escolha o metodo: "))

        volatilidade = 0.01
        taxaJuros = 0.01
        valorAtual = 1
        quantidadeOpcoes = 1000
        T = 1
        precoExecucao = 1
        t = 0

        if questao == 1:
            with open("saida_cenario1_subitem1.txt", "a") as external_file:
                opcao = mt.calculaOpcao(volatilidade, valorAtual, precoExecucao, N, M, L, T, taxaJuros, t, metodo) + 1
                print("A opção é precificada em R$" + str(opcao), file=external_file)
                opcaoInterpolar = mt.calculaOpcaoInterpolacao(volatilidade, valorAtual, precoExecucao, N, M, L, T, taxaJuros, t, metodo) + precoExecucao
                print("A opção é precificada utilizando interpolação em R$" + str(opcaoInterpolar), file=external_file)

        elif questao == 2:
            t = 0.5
            mt.analiseLucroPrejuizo(volatilidade, N, M, L, precoExecucao, T, taxaJuros, t, valorAtual, quantidadeOpcoes, metodo)

        elif questao == 3:
            # analise do lucro com parametros diferentes
            volatilidade = 0.02
            
            opcao = mt.calculaOpcao(volatilidade, valorAtual, precoExecucao, N, M, L, T, taxaJuros, t, metodo) + precoExecucao
            print("A opção é precificada em R$" + str(opcao))
            
            opcaoInterpolar = mt.calculaOpcaoInterpolacao(volatilidade, valorAtual, precoExecucao, N, M, L, T, taxaJuros, t, metodo) + precoExecucao
            print("A opção é precificada utilizando interpolação em R$" + str(opcaoInterpolar))

            t = 0.5
            mt.analiseLucroPrejuizo(volatilidade, N, M, L, precoExecucao, T, taxaJuros, t, valorAtual, quantidadeOpcoes, metodo)

        elif questao == 4:
            # analise do lucro com parametros diferentes 2
            volatilidade = 0.1
            taxaJuros = 0.1

            opcao = mt.calculaOpcao(volatilidade, valorAtual, precoExecucao, N, M, L, T, taxaJuros, t, metodo) + precoExecucao
            print("A opção é precificada em R$" + str(opcao))

            opcaoInterpolar = mt.calculaOpcaoInterpolacao(volatilidade, valorAtual, precoExecucao, N, M, L, T, taxaJuros, t, metodo) + precoExecucao
            print("A opção é precificada utilizando interpolação em R$" + str(opcaoInterpolar))

            t = 0.5
            mt.analiseLucroPrejuizo(volatilidade, N, M, L, precoExecucao, T, taxaJuros, t, valorAtual, quantidadeOpcoes, metodo)
        else:
            print("Metodo invalido!")

    elif cenario == 2:

        volatilidade = 0.1692
        taxaJuros = 0.1075
        valorAtual = 5.6376
        quantidadeOpcoes = 100000
        T = 3/12
        precoExecucao = 5.7
        t = 0

        premio = mt.calculaPremio(volatilidade, valorAtual, precoExecucao, N, M, L, T, taxaJuros, t, quantidadeOpcoes, metodo)
        print(premio)

        t = T
        mt.analiseLucroPrejuizo(volatilidade, N, M, L, precoExecucao, T, taxaJuros, t, valorAtual, quantidadeOpcoes, metodo)

    elif cenario == 3:
        # mudar dados para realidade 
        volatilidade = 0.1692
        taxaJuros = 0.1075
        valorAtual = 5.6376
        quantidadeOpcoes = 100000
        T = 3/12
        precoExecucao = 5.7


    else:   
        print("Cenario invalido!")

imprimeMenu()