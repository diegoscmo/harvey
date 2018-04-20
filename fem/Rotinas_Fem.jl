################################################################################
#####                             Rotinas FEM                             ######
################################################################################

# Carrega os arquivos com as rotinas de elementos finitos
include("Gera_Malha.jl")       # GeraMalha, gl_livres_elemento
include("Quad4_I.jl")          # Kquad4_I
include("Monta_Global.jl")     # Global, Expande_Vetor
include("Gmsh.jl")             # Funções GMSH
include("Harmonica.jl")        # Cálculo Harmônico
include("Condicionamento.jl")
