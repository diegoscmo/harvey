################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

# using Printf
# cd(ENV["HARVEY"]); include("main.jl"); main();

# Carrega rotinas
include("fem\\Rotinas_Fem.jl")      # Rotinas de Elementos Finitos
include("opt\\Top_Opt.jl")          # Rotinas de LA e Topológica
include("fobj\\0_Todas_Fobj.jl")    # Funções Objetivo

#
# Entrada de dados + execução da rotina
#
function main()

    # Nome do Arquivo
    dts = "Teste"

    # Parâmetros da Função Objetivo
    caso        = 2         # Seleciona o caso (ver Fobjs)
    freq        = 100.0     # Frequência de excitação
    alfa        = 0.0       # Amortecimento Proporcional - Massa
    beta        = 1E-8      # Amortecimento Proporcional - Rigidez
    A           = 1.0       # Peso da primeira Fobj
    Ye          = 1.0       # Porcentagem da restrição (Est ou R)
    dini        = 0.49      # Volume inicial

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 30        # Máximo de iteracoes externas
    max_int     = 500       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho         = 1.00      # Valor inicial de rho
    rho_max     = 2.50      # Valor maximo de rho
    mult_max    = 10.0      # Valor maximo dos multiplicadores

    # Parâmetros da topológica
    SP          = 3.00      # Parametro p do SIMP
    raiof       = 0.031     # Tamanho do filtro [m]
    vmin        = 1E-9      # Densidade mínima

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
    #presos = [  LX      0.0         LX          LY          1 ;
    #            LX      0.0         LX          LY          2 ]

    # Carregamentos:
    #        [ ponto_X  ponto_Y     força       dir (X=1 Y=2)]
    forcas = [ LX       LY/2.0      -9000.0     2 ]
    #forcas = [ LX/2.0     0.0      -9000.0     2 ]

    # Executa.
    Top_Opt(dts, caso, freq, alfa, beta, Ye, A, dini, max_ext, max_int, tol_ext,
               tol_int, rho, rho_max, mult_max, SP, raiof, vmin, NX, NY, LX, LY,
                                    young, poisson, esp, p_dens, presos, forcas)
end
