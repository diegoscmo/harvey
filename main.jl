################################################################################
# Lagrangiano Aumentado // Otimização Topológica
################################################################################
# cd(ENV["HARVEY"]); >> diretório da variável de sistema

# Carrega os arquivos com as rotinas para otimização
include("opt\\F_Lagrangiana.jl")
include("opt\\Descent.jl")
include("opt\\LS_Methods.jl")

# Carrega os arquivos com as rotinas de elementos finitos
include("fem\\gera_malha.jl")        # GeraMalha, gl_livres_elemento, Testa_Malha
include("fem\\quad4.jl")             # Kquad4, Kquad4_I
include("fem\\global.jl")            # Global, Expande_Vetor
include("fem\\gmsh.jl")              # Funções GMSH

# Carrega cálculo das derivadas
include("sens.jl")                  # Análise de sensibilidade, função Dif_Fin
#include("opt\\Dif_Fin.jl")

# Parâmetros do Lagrangiano Aumentado
max_ext = 100         # Máximo de iteracoes externas
tol_ext = 1E-4      # Tolerância do laço externo
max_int = 50        # Máximo de 'iterações internas
tol_int = 1E-4      # Tolerância do laço interno
rho_max = 10.0       # Valor maximo de rho
mult_max = 1E2      # Valor maximo dos multiplicadores
descent = "Steep"   # Metodo de descent (Steep, FletRe, DFP, BFGS)
search = "Line"     # Metodo de busca (Line, Golden ou Back)

# Parâmetros do problema Harmônico
f    = 5.0      # Frequencia
alfa = 0.0      # Amortecimento proporcional de Rayleigh
beta = 1E-8

# Definição geométrica do problema (retângulo):
LX = 1.0       # Comprimento em X
LY = 0.5       # Comprimento em Y

# Malha:
NX = 120#120        # Nr. de elementos em X
NY = 60#60        # Nr. de elementos em Y

# Material:
young   = 210E9     # Módulo de Young
poisson = 0.0       # Coeficiente de Poisson (0.0 > viga)
esp     = 1.0       # Espessura do retângulo
rho     = 7860.0    # Densidade

# Volume/Pseudo-densidades iniciais
dens_ini = 0.49

# Restrições de deslocamento (apoios):
#        [ ponto_X_inicial ponto_Y_inicial ponto_X_final ponto_Y_final direção (X=1 Y=2)]
presos = [ 0.0             0.0             0.0           LY           1       ;
           0.0             0.0             0.0           LY           2       ]

# Carregamentos:
#        [ ponto_X         ponto_Y         força         direção (X=1 Y=2)]
forcas = [ LX              LY/2.0          -9000.0         2 ]

########### PROGRAMA ###########

# Dados calculados e declaração de variáveis:
npresos = size(presos,1)            # Nr. de apoios
nforcas = size(forcas,1)            # Nr. de carregamentos
nelems  = NX*NY                     # Nr. de elementos
nnos    = (NX+1)*(NY+1)             # Nr. de nós

# Gera a malha
coord,ijk,nos_forcas,ID,nr_gl_livres = GeraMalha(LX,LY,NX,NY,npresos,presos,nforcas,forcas)

# Gera a matrizes de rigidez e massa de um elemento finito - malha toda igual
(K0,) = Kquad4_I(1,coord,ijk,young,poisson,esp)
M0 = Mquad4(1,coord,ijk,esp,rho)

# Parametro p do SIMP
simp = 3.0

# Correção de Hz para rad/s
w = 2.0*pi*f

# Pseudo-densidades para a montagem global com o SIMP
x  = dens_ini*ones(nelems)
xl = 0.00*ones(nelems)
xu = 1.00*ones(nelems)

# Calcula valores para aplicação de filtros
# Mapeia coordenadas centrais dos elementos para filtros
centros = Coord_Centros(coord,ijk,nelems)

# Filtro de densidades4.39 pg65
rmin = 0.031 #m

# Calcula Hi
H = zeros(Float64,nelems,nelems)
for j=1:nelems
    for k=1:nelems
        dkj = sqrt((centros[k,1]-centros[j,1])^2.0+(centros[k,2]-centros[j,2])^2.0)
        if dkj <= rmin
            H[j,k]= rmin - dkj
        end
    end
    #else Hj=0
end


# Define o Função Objetivo
function F_Obj(x)

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    KG,MG,F = Global(nelems,ijk,ID,K0,M0,x,simp,nforcas,nos_forcas,nr_gl_livres)

    # Monta a matriz equivalente para análise harmônica
    #CG = alfa*MG + beta*KG
    #KD = KG + w*im*CG - w^2.0*MG
    KD = KG
       # Resolve o sistema
          U = vec(lufact(KD)\F);

    # Função objetivo

    #fun = real((0.5*w^2.0)*conj(U)'*CG*U)      # Potencia Ativa
    #fun = (abs(U'*F))^2.0                       # Flexibilidade dinamica
    fun = U'*KG*U

    # Funções de restrição
    res = [
            # Restrição de Volume
            -0.49+mean(x)
            # Restrição do R
        #   -1.0+(100.0+100.0*log(w^2.0*(real(U'*MG*U)/real(U'*KG*U))))
            ]
    return fun,res
end

####### Inicio da rotina de Otm #######

# Prepara o plot
Inicializa_Malha_Gmsh("zgmsh",nnos,nelems,ijk,coord,2)
# Plota saída para o gmsh
Adiciona_Vista_Escalar_Gmsh("zgmsh","dens",nelems,x,0.0)

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
global viewcount
# Dá o display e inicializa o laço externo
println("\n LagAug:: ",descent," // ",search)
for i_ext=1:max_ext

    # Display do valor atual
    @printf("\n  iter: %d \tfitness:\t%.3e\tvolume:\t%.3f ",i_ext, valor_fun,mean(x))
    @printf("\n \t\tmult. lagrange:\t")
    show(mult_res)

    # Soluciona o problema interno (e salva a derivada)
    x,dL,count = Descent(x, mult_res, rho, xl, xu, max_int, tol_int, count, search, descent)

    # Plota saída para o gmsh
    Adiciona_Vista_Escalar_Gmsh("zgmsh",string("dens ",string(mean(x))),nelems,x,Float64(i_ext))

    # Atualiza o valor o multiplicador das restricoes (u)
    valor_fun,valor_res = F_Obj(x)
    mult_res = min.(mult_max,max.(rho*valor_res + mult_res, 0.))

    # Atualiza o multiplicador de penalizacao (c)
    crho_nov = max.(0.,norm(norm(max.(valor_res, -mult_res/rho))))
    if crho_nov >= .9*crho_ant
        rho = min.(1.1*rho, rho_max)
    end # if
    crho_ant = copy(crho_nov)

    # Verifica os criterios:
    if norm(dL) <= tol_ext &&                 # Condicao de gradiente
        norm(max.(valor_res,0.)) <= tol_ext &&   # Condicao de viabilidade
        norm(valor_res'*mult_res) < tol_ext #&&  # Condicao de complementariedade
        break
    end # if criterios

end # for i_ext

# Imprime a saida dos resultados
@printf("\n  finalizado! \tfitness:\t%.3e\tvolume:\t%.3f ",valor_fun,mean(x))
@printf("\n \t\tmult. lagrange:\t")
show(mult_res)

# Aplica o filtro dedensidades
y = zeros(Float64,nelems)
for j=1:nelems
y[j] = (H[j,:]'*x)/sum(H[j,:])
end

Adiciona_Vista_Escalar_Gmsh("zgmsh",string("filt ",string(mean(y))),nelems,y,Float64(count+1))
