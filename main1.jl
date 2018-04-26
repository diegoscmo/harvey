################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

using Printf; using SparseArrays; using LinearAlgebra; using IterativeEigensolvers; using ProgressMeter
# cd(ENV["HARVEY"]); include("main.jl"); main();

# Carrega rotinas
include("fem/rotinas_fem.jl")      # Rotinas de Elementos Finitos
include("opt/top_opt.jl")          # Rotinas de LA e Topológica

# Escolher função objetivo
#include("fobj/1_estatico.jl")
#include("fobj/2_dinamico.jl")
#include("fobj/3_din_r-est.jl")
#include("fobj/4_din_a-est.jl")
#include("fobj/5_potencia.jl")
#include("fobj/6_pot_r-est.jl")
#include("fobj/7_pot_a-est.jl")
#include("fobj/8_pot_a-est_r-r.jl")
#include("fobj/9_vol_r-stressdin.jl")
#include("fobj/10_pot_r-r_r-stressdin.jl")
#include("fobj/11_norma.jl")
include("fobj/12_norma_r-stress.jl")

#
# Entrada de dados + execução da rotina
#
function main()

    # Nome do Arquivo
    dts = "MIN_750_B_A099_P22_Q200_S2E5_V10_5000"

    # Parâmetros da Função Objetivo
    freq        = 750.0     # Frequência de excitação
    alfa        = 0.0       # Amortecimento Proporcional - Massa
    beta        = 1E-8     # Amortecimento Proporcional - Rigidez
    #beta        = 0.1/(2.0*pi*freq)
    A           = 0.99      # Peso da primeira Fobj (negativo para resonant)
    P           = 2.0       # Parâmetros da Norma
    q           = 2.0       # FIXME
    R_b         = .999      # Porcentagem da restrição (Est ou R)
    Sy          = 2E5       # Tensão limite
    dini        = 1.00      # Volume inicial R=1 #0.52631277505 600@0.1w 0.532239322199 @ 600@1E-8b
    dmax        = 0.50      # Restrição de volume  F=1 #0.41345415 600@1E-8b

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 30        # Máximo de iteracoes externas
    max_int     = 500       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho_max     = 100.0     # Valor maximo de rho

    # Penalização inicial (usando sempre 3 rhos)
    rho         = [0.5      # Volume
                   0.0      # R 5
                   0.05]    # Tensão (rho inicial) 0.05

    # Parâmetros da topológica
    SP          = 3.00      # Parametro p do SIMP
    QP          = 2.00       # Relaxação da tensão
    raiof       = 0.031     # Tamanho do filtro [m]
    vmin        = 1E-3      # Densidade mínima
    csi         = 0.5       # Valor inicial do parâmetro Heaviside

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

    # Executa.
    Top_Opt(dts, Sy, freq, alfa, beta, R_b, A, dini, max_ext, max_int, tol_ext,
              tol_int, rho, rho_max, SP, raiof, vmin, NX, NY, LX, LY,
            young, poisson, esp, p_dens, presos, forcas, travas, QP, csi, dmax, P, q)

    # Dá um tempinho e fecha
    println("  Fim! ", dts,"\n")
    sleep(60)
    exit()

end
main()
