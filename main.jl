################################################################################
# Lagrangiano Aumentado // Otimização Topológica
################################################################################

# Copiar para o console para executar:
# cd(ENV["HARVEY"]);
# include("main.jl"); main();

# Carrega os arquivos com as rotinas de elementos finitos
include("fem\\Gera_Malha.jl")       # GeraMalha, gl_livres_elemento
include("fem\\Quad4.jl")            # Kquad4, Kquad4_I
include("fem\\Monta_Global.jl")     # Global, Expande_Vetor
include("fem\\Gmsh.jl")             # Funções GMSH
include("fem\\Harmonica.jl")         # Cálculo Harmônico

# Carrega os arquivos com métodos as rotinas para otimização
include("opt\\Top_Opt.jl")
include("opt\\Steepest.jl")         # Método de Descida, Steepest
include("opt\\Wall_Search.jl")      # Procura em linha, Wall_Search
include("opt\\Filtros.jl")          # Filtros de densidades e sensibilidades
include("opt\\Saida.jl")            # Impressão das saídas em arquivo e console
include("opt\\Dif_Fin.jl")

# Carrega cálculo da Fobj, Lagrangiana e derivadas, selecionar um
#include("fobj\\1_Estatico.jl")
#include("fobj\\2_Dinamico.jl")
#include("fobj\\3_Din_R-Est.jl")
#include("fobj\\4_Din_A-Est.jl")
include("fobj\\5_Potencia.jl")
#include("fobj\\6_Pot_R-Est.jl")
#include("fobj\\7_Pot_A-Est.jl")
#include("fobj\\8_Pot_A-Est_R-R.jl")
#include("fobj\\0_Todas_Fobj.jl")

function main()
    # Parâmetros Harmônica / Fobj
    num         = 6
    caso        = 0
    freq        = 180.0
    alfa        = 0.0
    beta        = 0.1/(2.0*pi*freq)
    A           = 0.00
    Ye          = 1.5
    dini        = 0.49

    if alfa == 666.6
        alfa = 0.8*(2.0*pi*freq)
    end
    if beta == 666.6
        beta = 0.1/(2.0*pi*freq)
    end

    Top_Opt(num, caso, freq, alfa, beta, Ye, A, dini)
end
