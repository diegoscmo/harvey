################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

# using Printf; using SparseArrays; using LinearAlgebra; using IterativeEigensolvers
# cd(ENV["HARVEY"]); include("main.jl"); main();

# Carrega rotinas
include("fem\\Rotinas_Fem.jl")      # Rotinas de Elementos Finitos
include("opt\\Top_Opt.jl")          # Rotinas de LA e Topológica
include("opt\\Dif_Fin.jl")

# Escolher função objetivo
#include("fobj\\1_Estatico.jl")
#include("fobj\\2_Dinamico.jl")
#include("fobj\\3_Din_R-Est.jl")
#include("fobj\\4_Din_A-Est.jl")
#include("fobj\\5_Potencia.jl")
#include("fobj\\6_Pot_R-Est.jl")
#include("fobj\\7_Pot_A-Est.jl")
#include("fobj\\8_Pot_A-Est_R-R.jl")
#include("fobj\\11_Est_R-V_R-Stress.jl")
#include("fobj\\12_Vol_R-Stress.jl")
#include("fobj\\15_Pot_R-R_R-Stress.jl")
include("fobj\\16_Pot_R-R_R-StressDin.jl")

#
# Entrada de dados + execução da rotina
#
function main()

    # Nome do Arquivo
    dts = "16_101_S2E5_v049-5_SEMFSeVOL"
    #caso 92 1350

    # Parâmetros da Função Objetivo
    freq        = 700.0    # Frequência de excitação
    alfa        = 0.0       # Amortecimento Proporcional - Massa
    #beta        = 1E-8     # Amortecimento Proporcional - Rigidez
    beta        = 0.1/(2.0*pi*freq)
    A           = -1.0       # Peso da primeira Fobj (negativo para resonant)
    Ye          = 0.999     # Porcentagem da restrição (Est ou R)
    Sy          = 2E5       # Tensão limite
    dini        = 1.0       # Volume inicial
    dmax        = 0.99      # Restrição de volume

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 30        # Máximo de iteracoes externas
    max_int     = 500       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho         = 0.5       # Valor inicial de rho
    rho_max     = 10.00     # Valor maximo de rho
    mult_max    = 10.0      # Valor maximo dos multiplicadores

    # Parâmetros da topológica
    SP          = 3.00      # Parametro p do SIMP
    QP          = 2.6       # Relaxação da tensão
    raiof       = 0.031     # Tamanho do filtro [m]
    vmin        = 1E-3      # Densidade mínima
    csi         = 0.5       # Valor inicial do parâmetro Heaviside

    # Parâmetros do problema de FEM, 60x30 = 1800 // 140x70 = 9800
    NX          = 60        # Nr. de elementos em X
    NY          = 30        # Nr. de elementos em Y
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
    #            0.0     LY/2.0-2.0*dy         0.0         LY/2.0+2.0*dy          1 ;
    #            0.0     LY/2.0-2.0*dy         0.0         LY/2.0+2.0*dy          2 ]
    #presos = [  LX      0.0         LX          LY          1 ;
    #            LX      0.0         LX          LY          2 ]

    # Carregamentos:
    #        [ ponto_X  ponto_Y     força       dir (X=1 Y=2)]
        dy = LY/NY
    forcas = [ LX       LY/2.0-2.0*dy      -1125.0      2
               LX       LY/2.0-1.0*dy      -2250.0     2
               LX       LY/2.0             -2250.0     2
               LX       LY/2.0+1.0*dy      -2250.0     2
               LX       LY/2.0+2.0*dy      -1125.0      2  ]
    #forcas = [ LX     LY/2.0      -9000.0     2 ]

    # Executa.
    Top_Opt(dts, Sy, freq, alfa, beta, Ye, A, dini, max_ext, max_int, tol_ext,
               tol_int, rho, rho_max, mult_max, SP, raiof, vmin, NX, NY, LX, LY,
                           young, poisson, esp, p_dens, presos, forcas, QP, csi, dmax)
end
