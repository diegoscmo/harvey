################################################################################
#####                             Rotinas FEM                             ######
################################################################################

# Carrega os arquivos com as rotinas de elementos finitos
include("gera_malha.jl")       # GeraMalha, gl_livres_elemento
include("quad4.jl")          # Kquad4_I
include("monta_global.jl")     # Global, Expande_Vetor
include("gmsh.jl")             # Funções GMSH
include("harmonica.jl")        # Cálculo Harmônico
include("condicionamento.jl")
