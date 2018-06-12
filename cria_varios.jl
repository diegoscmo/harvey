################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

# Carrega rotinas
using Printf
using IterTools
using ProgressMeter
include("opt/top_opt.jl")          # Rotinas de LA e Topológica

#
# Entrada de dados + execução da rotina
#
function Cria_Varios()

    #  Tipo de analise (1=NG_; 2=NL_; 3= NGS; 4=NLS )
    tipoan      = [1]

    # Parâmetros da Função Objetivo, esses todos vao no nome
    freq        = [180.0;270.0;600.0;750.0]   # Frequência de excitação
    beta        = ["W";1E-8]       # Amortecimento Proporcional - Rigidez
    A           = [0.99;0.999;0.95]      # Peso da primeira Fobj (negativo para resonant)
    P           = [2.0; 4.0]       # Parâmetros da Norma
    q           = [2.0]       #
    Sy          = [4E20]      # Tensão limite
    dini        = [0.50;1.0]      # Volume inicial R=1 #0.52631277505 600@0.1w 0.532239322199 @ 600@1E-8b
    dmax        = [0.50]      # Restrição de volume  F=1 #0.41345415 600@1E-8b
    QP          = [1.50]      # Relaxação da tensão
    csi0        = [1.0; 5.0; 10.0]       # Primeiro valor do Heaviside
    csim        = [1.2 ; 1.3; 1.4; 1.5]      # Multiplicador do heaviside

    # Se nao quiser criar casos heaviside, deixe false!
    heavi = false

    V = collect(product(tipoan,freq,beta,A,P,q,Sy,dini,dmax,QP,csi0,csim))

    # Daqui pra baixo nao muda! #

    # Penalização inicial (usando sempre 3 rhos, vai no casefile, mas nao no nome)
    rho1        = 0.5      # Volume
    rho2        = 0.0      # R 5 (obsoleto, mas fica aqui de reserva)
    rho3        = 0.0      # Tensão (rho inicial) 0.05

    # Parâmetros do Lagrangiano Aumentado (esses vao no casefile, mas nao vao no nome)
    alfa        = 0.0       # Amortecimento Proporcional - Massa
    raiof       = 0.031     # Tamanho do filtro [m]
    max_fil     = 30        # Máximo de iteracoes externas
    max_hev     = 30        # Número de iterações com Heaviside
    max_int     = 500       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho_max     = 100.0     # Valor maximo de rho

    # Retorna o nome da analise e já cria o casefile
    for i=1:size(V,1)

        dts = Name_Game(V[i][1],V[i][2],alfa,V[i][3],V[i][4],V[i][5],V[i][6],V[i][7],V[i][8],
                            V[i][9],V[i][10],V[i][11],V[i][12], rho1,rho2,rho3,raiof,max_fil,
                                               max_hev,max_int,tol_ext,tol_int,rho_max,heavi)

        sleep(0.1)
    end

    println("  ",size(V,1)," casos criados!")

end
Cria_Varios()
sleep(10)
