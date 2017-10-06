
################################################################################
# Elementos Finitos - Calculo Harmonico
################################################################################
# cd(ENV["HARVEY"]); cd("fem");

using Plots

# Carrega os arquivos com as rotinas auxiliares
include("gera_malha.jl")        # GeraMalha, gl_livres_elemento, Testa_Malha
include("quad4.jl")             # Kquad4
include("global.jl")            # Global, Expande_Vetor

################################################################################
# Dados do problema:
################################################################################

# Definição geométrica do problema (retângulo):
LX = 1.0       # Comprimento em X
LY = 0.1       # Comprimento em Y

# Malha:
NX = 100        # Nr. de elementos em X
NY = 10        # Nr. de elementos em Y

# Material:
young   = 210E9     # Módulo de Young
poisson = 0.0     # Coeficiente de Poisson (0.0 > viga)
esp     = 1.0       # Espessura do retângulo
rho     = 7860.0    # Densidade

#w = 4.730041^2*sqrt((young*LY^2)/(12*rho*LX^4))/2/pi

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

  # Gera a matriz de rigidez de um elemento finito - malha toda igual
  K0 = 4.0*Kquad4(1,coord,ijk,young,poisson,esp)
  M0 = 4.0*Mquad4(1,coord,ijk,esp,rho)

  # Parametro p do SIMP
  simp = 3.0
  # Pseudo-densidades para a montagem global com o SIMP
  dens = 0.49*ones(nelems)

  # Monta matriz de rigidez Global
  KG,F = Global(nelems,ijk,ID,K0,dens,simp,nforcas,nos_forcas,nr_gl_livres)

  # Corrige os fatores densidade para a massa (Olhoff e Du)
  C1 =  6.0E5
  C2 = -5.0E6
  for i=1:nelems
      if dens[i] < 0.1
          dens[i] = C1*dens[i]^6 + C2*dens[i]^7
      end #if dens[i]
  end #for i

  # Monta matriz de massa Global
  MG, = Global(nelems,ijk,ID,M0,dens,simp,nforcas,nos_forcas,nr_gl_livres)

  # Amortecimento proporcional de Rayleigh
  alf = 0.0
  bet = 0.0
  CG = alf*MG + bet*KG

  maxfreq  = 300
  stepfreq = 5
  a = []
  b = []

  # Análise Harmônica, varredura
  #while false
  for f = 0:stepfreq:maxfreq
  KD = KG + 2*pi*f*im*CG - (2*pi*f)^2*MG

  # Soluciona o sistema de equações
  C = lufact(KD)
  U_sem_zeros = vec(C\F)

  # Expande o vetor de deslocamentos, adicionando os graus de liberdade travados:
  U  = Expande_Vetor(U_sem_zeros,nnos,ID)
  Fe = Expande_Vetor(F,nnos,ID)

  # Flexibilidade dinamica
  Y = abs(U'*Fe)^2

  push!(a,f)
  push!(b,Y)

  end

  # Plota Flexibilidade Dinamica
  z = plot(a,log.(b))

  # Cálculo rápido de autovalores
  autov = eigs(KG,MG,nev=3,which=:SM)[1]
  print(sqrt.(autov)/(2.0*pi))
