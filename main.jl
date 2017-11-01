################################################################################
# Lagrangiano Aumentado // Otimização Topológica
################################################################################
# cd(ENV["HARVEY"]); >> diretório da variável de sistema

# Horario de Execução
dtf = Dates.now()
dts = Dates.format(dtf,"YYYY-mm-dd-HH-MM-SS")
tic()

# Carrega os arquivos com as rotinas para otimização
include("opt\\F_Lagrangiana.jl")    # F_Lagrangiana
include("opt\\Descent.jl")          # Métodos de Descida
include("opt\\LS_Methods.jl")       # Métodos de busca em linha

# Carrega os arquivos com as rotinas de elementos finitos
include("fem\\Gera_Malha.jl")       # GeraMalha, gl_livres_elemento
include("fem\\Quad4.jl")            # Kquad4, Kquad4_I
include("fem\\Monta_Global.jl")     # Global, Expande_Vetor
include("fem\\Gmsh.jl")             # Funções GMSH

# Carrega cálculo das derivadas
include("Fobj_Estatico.jl")         # Fobj e Sensibilidades
include("Filtros.jl")               # Filtros de densidades e sensibilidades
include("Saidas.jl")                # Impressão das saídas em arquivo e console

# Parâmetros do Lagrangiano Aumentado
max_ext     = 10        # Máximo de iteracoes externas
tol_ext     = 1E-4      # Tolerância do laço externo
max_int     = 20        # Máximo de 'iterações internas
tol_int     = 1E-4      # Tolerância do laço interno
rho_max     = 3.0       # Valor maximo de rho
mult_max    = 1E2       # Valor maximo dos multiplicadores
descent     = "Steep"   # Metodo de descent (Steep, FletRe, DFP, BFGS)
search      = "Line"    # Metodo de busca (Line, Golden ou Back)

# Parâmetros da topológica
dens_ini    = 0.60      # Volume/Pseudo-densidades iniciais
simp        = 3.0       # Parametro p do SIMP
raiof       = 0.04      # Tamanho do filtro [m]

# Parâmetros do problema Harmônico
f    = 5.0      # Frequencia
alfa = 0.0      # Amortecimento proporcional de Rayleigh
beta = 1E-8

# Definição geométrica do problema (retângulo):
LX = 1.0       # Comprimento em X
LY = 0.5       # Comprimento em Y

# Malha:
NX = 20#120        # Nr. de elementos em X
NY = 10#60        # Nr. de elementos em Y

# Material:
young   = 210E9     # Módulo de Young
poisson = 0.0       # Coeficiente de Poisson (0.0 > viga)
esp     = 1.0       # Espessura do retângulo
rho     = 7860.0    # Densidade

# Restrições de deslocamento (apoios):
#        [ ponto_X_inicial ponto_Y_inicial ponto_X_final ponto_Y_final direção (X=1 Y=2)]
presos = [ 0.0             0.0             0.0           LY           1       ;
           0.0             0.0             0.0           LY           2       ]

# Carregamentos:
#        [ ponto_X         ponto_Y         força         direção (X=1 Y=2)]
forcas = [ LX              LY/2.0          -9000.0         2 ]


############################## PROGRAMA ##################################

# Dados calculados e declaração de variáveis:
npresos = size(presos,1)            # Nr. de apoios
nforcas = size(forcas,1)            # Nr. de carregamentos
nelems  = NX*NY                     # Nr. de elementos
nnos    = (NX+1)*(NY+1)             # Nr. de nós

# Gera a malha
coord,conect,nos_forcas,ID,nr_gl_livres = GeraMalha(LX,LY,NX,NY,npresos,presos,nforcas,forcas)

# Gera a matrizes de rigidez e massa de um elemento finito - malha toda igual
(K0,) = Kquad4_I(1,coord,conect,young,poisson,esp)
M0 = Mquad4(1,coord,conect,esp,rho)

# Inicializa matrizes
KG = zeros(Float64,nr_gl_livres,nr_gl_livres)
MG = zeros(Float64,nr_gl_livres,nr_gl_livres)
CG = zeros(Float64,nr_gl_livres,nr_gl_livres)
KD = zeros(Float64,nr_gl_livres,nr_gl_livres)
F  = zeros(Float64,nr_gl_livres)

# Pseudo-densidades para a montagem global com o SIMP
x  = dens_ini*ones(nelems)
xl = 0.00*ones(nelems)
xu = 1.00*ones(nelems)

# Correção de Hz para rad/s
w = 2.0*pi*f

# Prepara o H para filtros 63
H,sumH = Calc_H(coord,conect,nelems,raiof)

# Prepara o plot e primeira saída
Inicializa_Malha_Gmsh(string("z",dts,".pos"),nnos,nelems,conect,coord,2)
Adiciona_Vista_Escalar_Gmsh(string("z",dts,".pos"),"x",nelems,x,0.0)

# Obtem os valores de f e g em x0 para o calculo de c0
valor_fun,valor_res = F_Obj(x)

# Inicializa multiplicadores de Lagrange (u)
numres   = size(valor_res,1)
mult_res = zeros(numres)

# Define fator de penalização inicial (c)
rho = max.(1E-6,min(rho_max,(2.*abs(valor_fun)/norm(max.(valor_res,0.))^2.)))

# E calcula o criterio de atualizacao do c
crho_ant = max.(0.,norm(norm(max.(valor_res, -mult_res/rho))))

# Inicializa o contador de avaliacoes da F_Obj
count = 1
viewcount = 0.0

# Dá o display e inicializa o laço externo
Imprime(0,x,rho,mult_res,valor_fun,valor_res)
for i_ext=1:max_ext

    # Soluciona o problema interno (e salva a derivada)
    x,dL,count = Descent(x, mult_res, rho, xl, xu, max_int, tol_int, count, search, descent)

    # Atualiza o valor o multiplicador das restricoes (u)
    valor_fun,valor_res = F_Obj(x)
    mult_res = min.(mult_max,max.(rho*valor_res + mult_res, 0.))

    # Atualiza o multiplicador de penalizacao (c)
    crho_nov = max.(0.,norm(norm(max.(valor_res, -mult_res/rho))))
    if crho_nov >= .9*crho_ant
        rho = min.(1.1*rho, rho_max)
    end # if
    crho_ant = copy(crho_nov)

    # Imprime resultado atual e plota saida para o gmsh
    Adiciona_Vista_Escalar_Gmsh(string("z",dts,".pos"),string("x",string(mean(x))),nelems,x,Float64(i_ext))
    Imprime(i_ext,x,rho,mult_res,valor_fun,valor_res)

    # Verifica os criterios:
    if norm(dL) <= tol_ext &&                 # Condicao de gradiente
        norm(max.(valor_res,0.)) <= tol_ext &&   # Condicao de viabilidade
        norm(valor_res'*mult_res) < tol_ext #&&  # Condicao de complementariedade
        break
    end # if criterios

end # for i_ext

# Display Final
Imprime(-1,x,rho,mult_res,valor_fun,valor_res)
