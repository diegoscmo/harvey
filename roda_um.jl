################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

 using Printf
 using SparseArrays
 using LinearAlgebra
 using IterativeEigensolvers
 using ProgressMeter

# Carrega rotinas
include("fem/rotinas_fem.jl")      # Rotinas de Elementos Finitos
include("opt/top_opt.jl")          # Rotinas de LA e Topológica

# Carrega Fobjs
include("fobj/11_norma_global.jl")
include("fobj/12_norma_local.jl")
include("fobj/13_norma_global_r-stress.jl")
include("fobj/14_norma_local_r-stress.jl")

#
# Entrada de dados + execução da rotina
#
function Roda_Um()

    #  Tipo de analise (1=NG_; 2=NL_; 3= NGS; 4=NLS )
    tipoan      = 1

    # Parâmetros da Função Objetivo, esses todos vao no nome
    freq        = 200.0     # Frequência de excitação
    beta        = "W"       # Amortecimento Proporcional - Rigidez "W" ou valor
    A           = 0.8      # Peso da primeira Fobj (negativo para resonant)
    P           = 2.0       # Parâmetros da Norma
    q           = 2.0       #
    Sy          = 4E20      # Tensão limite
    dini        = 0.50      # Volume inicial R=1 #0.52631277505 600@0.1w 0.532239322199 @ 600@1E-8b
    dmax        = 0.50      # Restrição de volume  F=1 #0.41345415 600@1E-8b
    QP          = 1.50      # Relaxação da tensão
    csi0        = 1.0       # Primeiro valor do Heaviside
    csimul      = 1.2      # Multiplicador do heaviside

    # Penalização inicial (usando sempre 3 rhos, vai no casefile, mas nao no nome)
    rho1        = 0.5      # Volume
    rho2        = 0.0      # R 5 (obsoleto, mas fica aqui de reserva)
    rho3        = 0.0      # Tensão (rho inicial) 0.05

    # Parâmetros do Lagrangiano Aumentado (esses vao no casefile, mas nao vao no nome)
    alfa        = 0.0       # Amortecimento Proporcional - Massa
    raiof       = 0.031     # Tamanho do filtro [m]
    max_fil     = 10        # Máximo de iteracoes externas
    max_hev     = 10        # Número de iterações com Heaviside
    max_int     = 10       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho_max     = 100.0     # Valor maximo de rho

    # Parâmetros da topológica fixos (tudo daqui pra baixo nao muda, nao vai do casefile)
    SP          = 3.00      # Parametro p do SIMP
    vmin        = 1E-3      # Densidade mínima

    # Parametro obsoleto (apenas para funções antigas)
    R_b         = .999      # Porcentagem da restrição (Est ou R)

    # Parâmetros do problema de FEM, 60x30 = 1800 // 70x35 = 2450 // 80x40 = 3200
    NX          = 100       # Nr. de elementos em X
    NY          = 50        # Nr. de elementos em Y
    LX          = 1.0       # Comprimento em X
    LY          = 0.5       # Comprimento em Y
    young       = 210E9     # Módulo de Young
    poisson     = 0.0       # Coeficiente de Poisson
    esp         = 1.0       # Espessura do retângulo
    p_dens      = 7860.0    # Densidade

    # Restrições de deslocamento (apoios):
    #        [ ponto_X0 ponto_Y0    ponto_XF    ponto_YF    dir (X=1 Y=2)]
    presos = [  0.0     0.0         0.0         LY          1 ;
                0.0     0.0         0.0         LY          2 ]

    # Regiões com elementos travados em 1.0 (quadrado)
    #        [ ponto_X0     ponto_Y0    ponto_X1     ponto_Y1
    travas = [  LX-0.025      0.2         LX           0.3  ]

    # Carregamentos:
    #       [  x0   y0  xf  yf  carreg      dir (X=1 Y=2)
    forcas = [ LX   0.2 LX  0.3 -10000.0    2   ]

    # Retorna o nome da analise e já cria o casefile
    dts,sub = Name_Game(tipoan,freq,alfa,beta,A,P,q,Sy,dini,dmax,QP,csi0,csimul,rho1,rho2,rho3,
                       raiof,max_fil,max_hev,max_int,tol_ext,tol_int,rho_max)

    # Corrige os rhos para um vetor
    rho = [rho1;rho2;rho3]

    # Executa.
    Top_Opt(tipoan, dts, sub, Sy, freq, alfa, beta, R_b, A, dini, max_fil, max_int, tol_ext,
              tol_int, rho, rho_max, SP, raiof, vmin, NX, NY, LX, LY,
            young, poisson, esp, p_dens, presos, forcas, travas, QP, csi0, csimul, dmax, P, q, false)

    println("\n\n\n  Fim da analise com filtro! \n\n")

    Top_Opt(tipoan, dts, sub, Sy, freq, alfa, beta, R_b, A, dini, max_hev, max_int, tol_ext,
              tol_int, rho, rho_max, SP, raiof, vmin, NX, NY, LX, LY,
            young, poisson, esp, p_dens, presos, forcas, travas, QP, csi0, csimul, dmax, P, q, true)

    println("\n\n\n  Fim da analise com heaviside! \n\n")

end
Roda_Um()
