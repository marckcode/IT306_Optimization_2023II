from pyomo.environ import *
import pandas as pd


def antenas_opt():
    """
    Uma companhia de serviços públicos está implantando uma nova tecnologia que permite
    diminuir os custo de leitura dos consumos de seus clientes. A tecnologia utiliza um
    sistema que transmite via radiofrequencia os consumos desde os medidores de cada do-
    micilio para antenas centrais que logo retransmitem a leitura para a oficina central
    da companhia para gerar as faturas. Para assegurar uma leitura confíavel é necessário
    que a antena central de leitura esteja no máximo a 520 metros do medidor.
    
    - antenas.dat: 
        MED:  Conjunto de medidores
        ANT:  Conjunto de antenas
        BAI:  Conjunto de bairros
        DIST: Matriz de distancias
        DMAX: Distancia máxima de leitura
        CIO:  Custo de instalação e operação
        BB:   Bairro

    Returns:
        LOC:  Localização de antenas
        COVE: Matriz de covertura
    """
    
    m = AbstractModel()

    # declaração de conjuntos e parâmetros
    m.MED = Set()  # Conjunto de medidores
    m.ANT = Set()  # Conjunto de antenas
    m.BAI = Set()  # Conjunto de bairros

    m.DIST = Param(m.MED, m.ANT, domain=NonNegativeReals)  # Matriz de distâncias
    m.DMAX = Param(domain=NonNegativeReals)                # Distância máxima de leitura
    m.CIO  = Param(m.ANT, domain=NonNegativeReals)         # Custo de instalação e operação
    m.BB   = Param(m.ANT, domain=NonNegativeReals)         # Bairro

    # declaração das variáveis
    m.LOCA = Var(m.ANT, domain=Binary)  # localizacao de antenas
    m.COVE = Var(m.MED, m.ANT, domain=Binary)  # matriz de covertura  

    # definição da função objetivo
    # custo total de transporte
    def obj_func(m):
        return sum(m.CIO[a] * m.LOCA[a] for a in m.ANT)
    m.numero_antenas = Objective(rule=obj_func, sense=minimize)

    # definição das restrições
    def aloca_antena(m, med, ant):
        return m.COVE[med, ant] <= m.LOCA[ant]
    m.aloca_antena = Constraint(m.MED, m.ANT, rule=aloca_antena)

    def distancia_maxima(m, med, ant):
        return m.DIST[med, ant] * m.COVE[med, ant] <= m.DMAX
    m.distancia_maxima = Constraint(m.MED, m.ANT, rule=distancia_maxima)

    def covertura_minima(m, med):
        return sum(m.COVE[med, ant] for ant in m.ANT) >= 2
    m.covertura_minima = Constraint(m.MED, rule=covertura_minima)

    def posicao_bairro(m, b):
        return sum(m.LOCA[ant] for ant in m.ANT if m.BB[ant] == b) >= 1
    m.posicao_bairro = Constraint(m.BAI, rule=posicao_bairro)


    file = 'antenas.dat'
    instance = m.create_instance(file)
    solver = SolverFactory('cplex', executable="/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe")  # mudar locação
    results = solver.solve(instance, tee=False) # Mudar: tee=True -> show solver log

    results_list = []

     # turn data into a dataframe for better visualization of output
    loc_values = [instance.LOCA[l].value for l in instance.LOCA]
    df_loc = pd.DataFrame(loc_values, index=instance.ANT, columns=['LOC'])
    df_loc = df_loc.astype(int)
    print('\n')
    print('LOCA [*] := ')
    print(df_loc)

    # TRANS: 8x10 | ANT: 8 | MED: 10
    ant_values = [[instance.COVE[f, s].value for s in instance.ANT] for f in instance.MED]
    df_ant = pd.DataFrame(ant_values, index=instance.MED, columns=instance.ANT)
    df_ant = df_ant.astype(int)
    print('\n')
    print('COVE [*,*]')
    print(df_ant)
    print(f"\nCusto total de instalação e operação: {instance.numero_antenas()}")
    print('\n')


antenas_opt()