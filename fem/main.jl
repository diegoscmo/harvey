
################################################################################
# Elementos Finitos - Calculo Harmonico
################################################################################
# cd(ENV["HARVEY"]); cd("fem");

using Plots

# Carrega os arquivos com as rotinas auxiliares
include("gera_malha.jl")        # GeraMalha, gl_livres_elemento, Testa_Malha
include("quad4.jl")             # Kquad4, Kquad4_I
include("global.jl")            # Global, Expande_Vetor
include("harmonic.jl")          # Harmonica, Solve_FEM

################################################################################
# Dados do problema:
################################################################################

# Definição geométrica do problema (retângulo):
LX = 1.0       # Comprimento em X
LY = 0.5       # Comprimento em Y

# Malha:
NX = 20        # Nr. de elementos em X
NY = 10        # Nr. de elementos em Y

# Material:
young   = 210E9     # Módulo de Young
poisson = 0.0     # Coeficiente de Poisson (0.0 > viga)
esp     = 1.0       # Espessura do retângulo
rho     = 7860.0    # Densidade

# Restrições de deslocamento (apoios):
#        [ ponto_X_inicial ponto_Y_inicial ponto_X_final ponto_Y_final direção (X=1 Y=2)]
presos = [ 0.0             0.0             0.0           LY           1       ;
           0.0             0.0             0.0           LY           2       ]

# Carregamentos:
#        [ ponto_X         ponto_Y         força         direção (X=1 Y=2)]
forcas = [ LX              LY/2.0          -9000.0         2 ]


# Dados calculados e declaração de variáveis:
npresos = size(presos,1)       # Nr. de apoios
nforcas = size(forcas,1)       # Nr. de carregamentos
nelems  = NX*NY                   # Nr. de elementos
nnos    = (NX+1)*(NY+1)           # Nr. de nós

################################################################################
# Rotina principal
################################################################################

# Gera a malha
coord,ijk,nos_forcas,ID,nr_gl_livres = GeraMalha(LX,LY,NX,NY,npresos,presos,nforcas,forcas)
    # Onde:
    # coord = matriz com a posição de cada nó [x,y]
    # ijk = matriz conectividade
    # nos_forcas = matriz onde cada linha é uma força e as colunas são: [no, direção (X=1 Y=2), valor]
    # ID = matriz que informa a posição na matriz esparsa para cada grau de liberdade de cada nó
    # nr_gl_livres = número de graus de liberdades livres (tamanho da matriz esparsa)

# Gera a matrizes de rigidez e massa de um elemento finito - malha toda igual
(K0,) = Kquad4_I(1,coord,ijk,young,poisson,esp)
M0 = Mquad4(1,coord,ijk,esp,rho)

# Parametro p do SIMP
simp = 3.0
# Pseudo-densidades para a montagem global com o SIMP
dens = 0.49*ones(nelems)

# Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP
KG,MG,F = Global(nelems,ijk,ID,K0,M0,dens,simp,nforcas,nos_forcas,nr_gl_livres)

# Deslocamento estático
U,Fe = FEM_Solve(KG,F,nnos,ID)

# Amortecimento proporcional de Rayleigh
alfa = 0.0
beta = 0.0

# Dados para análise harmônica
maxfreq  = 400
stepfreq = 5

# Calcula e Plota a Flexibilidade Dinâmica
a,b = Harmonica(KG,MG,alfa,beta,F,nnos,ID,stepfreq,maxfreq)
z = plot(a,log.(b))

# Cálculo rápido de autovalores
autov = eigs(KG,MG,nev=3,which=:SM)[1]
print(sqrt.(autov)/(2.0*pi))
