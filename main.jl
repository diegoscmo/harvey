################################################################################
# Lagrangiano Aumentado // Otimização Topológica
################################################################################
# cd(ENV["HARVEY"]); using StaticArrays;
# include("Main.jl");main()

# Limpa tudo
#workspace();

# Carrega os arquivos com as rotinas de elementos finitos
include("fem\\Gera_Malha.jl")       # GeraMalha, gl_livres_elemento
include("fem\\Quad4.jl")            # Kquad4, Kquad4_I
include("fem\\Monta_Global.jl")     # Global, Expande_Vetor
include("fem\\Gmsh.jl")             # Funções GMSH

# Carrega os arquivos com as rotinas para otimização
include("opt\\Opt_Methods.jl")          # Método de Descida

# Carrega cálculo das derivadas
include("Fobj_Estatico.jl")         # Fobj e Sensibilidades

# Etc
include("Filtros.jl")               # Filtros de densidades e sensibilidades
include("Saida.jl")                 # Impressão das saídas em arquivo e console

mutable struct finitos_var
    KG::SparseMatrixCSC{Float64,Int64}
    CG::SparseMatrixCSC{Float64,Int64}
    MG::SparseMatrixCSC{Float64,Int64}
    KD::SparseMatrixCSC{Complex{Float64},Int64}
    UE::Array{Float64,1}
    UD::Array{Float64,1}
end

struct finitos_fix
     F::Array{Float64,1}
    K0::StaticArrays.SArray{Tuple{8,8},Float64,2,64}
    M0::StaticArrays.SArray{Tuple{8,8},Float64,2,64}
    simp::Float64
    nelems::Int64
    conect::Array{Int64,2}
    NX::Int64
    NY::Int64
    ID::Array{Int64,2}
    w::Float64
    alfa::Float64
    beta::Float64
end

struct filtros
      raiof::Float64;
       vizi::Array{Int64,2}
       nviz::Array{Int64,1}
       dviz::Array{Float64,2}
     filtro::String
end

function main()

    # Horario de Execução
    dtf = Dates.now()
    dts = Dates.format(dtf,"YYYY-mm-dd-HH-MM-SS")
    tic();

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 100       # Máximo de iteracoes externas
    max_int     = 200       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho_ini     = 0.25
    rho_max     = 0.5       # Valor maximo de rho
    mult_max    = 10.0      # Valor maximo dos multiplicadores

    # Métodos de busca
    descent     = "BFGS"   # Steep, BFGS, FR, DFP
    lsearch     = "Equal"   # Equal, Golden, Back

    # Parâmetros da topológica
    dens_ini    = 0.49      # Volume/Pseudo-densidades iniciais
    simp        = 3.0       # Parametro p do SIMP
    raiof       = 0.031     # Tamanho do filtro [m]
    filtro      = "Dens"    # Filtros: (Off ou Dens)

    # Parâmetros do problema Harmônico
    f    = 0.0      # Frequencia
    alfa = 0.0      # Amortecimento proporcional de Rayleigh
    beta = 1E-8

    # Parâmetros do problema de FEM
    NX = 100   #80        # Nr. de elementos em X
    NY = 50   #40        # Nr. de elementos em Y

    # Definição geométrica do problema (retângulo):
    LX = 1.0       # Comprimento em X
    LY = 0.5       # Comprimento em Y

    # Material:
    young   = 210E9     # Módulo de Young
    poisson = 0.0       # Coeficiente de Poisson (0.0 > viga)
    esp     = 1.0       # Espessura do retângulo
    rho     = 7860.0    # Densidade

    # Restrições de deslocamento (apoios):
    #        [ ponto_X_inicial ponto_Y_inicial ponto_X_final ponto_Y_final direção (X=1 Y=2)]
    presos = [  0.0             0.0             0.0           LY           1       ;
                0.0             0.0             0.0           LY           2       ]

    # Carregamentos:
    #        [ ponto_X         ponto_Y         força         direção (X=1 Y=2)]
    forcas = [ LX              LY/2.0          -9000.0         2 ]

    # Nós e elementos
    nnos    = (NX+1)*(NY+1)             # Nr. de nós
    nelems  = NX*NY                     # Nr. de elementos

    # Pseudo-densidades para a montagem global com o SIMP
    x  = dens_ini*ones(nelems)
    xl = zeros(nelems)
    xu =  ones(nelems)

    # Gera a malha
    coord, conect, nos_forcas, ID, gdl_livres =
                         GeraMalha(nnos, nelems, LX, LY, NX, NY, presos, forcas)

    # Gera a matrizes de rigidez e massa de um elemento finito - malha toda igual
    (K0n,) = Kquad4_I(1, coord, conect, young, poisson, esp)
    M0n    = Mquad4(1, coord, conect, esp, rho)

    # Transforma para StaticArrays
    K0 = SMatrix{8,8}(K0n)
    M0 = SMatrix{8,8}(M0n)

    # Prepara os vizinhos para filtros 63
    vizi,nviz,dviz = Proc_Vizinhos(nelems, coord, conect, raiof)
    filt           = filtros(raiof, vizi, nviz, dviz, filtro)

    # Correção de Hz para rad/s
    w = 2.0*pi*f

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    KG,MG = Global_KM(x, nelems, conect, ID, K0, M0, simp)
    F     =  Global_F(ID, nos_forcas, gdl_livres)

    # Monta CG e KD (Harmônica)
    CG = alfa*MG + beta*KG
    KD = KG + w*im*CG - w^2.0*MG

    # Resolve o sistema pela primeira vez
    UE = vec(lufact(KG)\F);
    UD = vec(lufact(KD)\F);

    # Agrupa nos tipos
    fem_v = finitos_var(KG, CG, MG, KD, UE, UD)
    fem_f = finitos_fix(F, K0, M0, simp, nelems, conect, NX, NY, ID, w, alfa, beta)

    # Obtem os valores de f e g em x0 para o calculo de c0
    valor_fun, valor_res = F_Obj(x, fem_v, fem_f)#, 0.0)

    # Salva o valor da primeira função para
    #valor_zero = copy(valor_fun)
    #valor_fun  = 1.0

    # Inicializa multiplicadores de Lagrange (u)
    numres   = size(valor_res, 1)
    mult_res = zeros(numres)

    # Define fator de penalização inicial (c)
    #rho = max.(1E-6,min(rho_max,(2.0*abs(valor_fun)/norm(max.(valor_res,0.))^2.)))
    rho = rho_ini

    # E calcula o criterio de atualizacao do c
    crho_ant = max.(0.,norm(max.(valor_res, -mult_res/rho)))

    # Inicializa o contador de avaliacoes da F_Obj e loops internos
    count = 1
    n_int = 0

    # Prepara o plot e primeira saída
    Imprime_0(x, dts, valor_fun, rho, mult_res, valor_res, nelems, nnos,
                      conect, coord, max_ext, max_int, tol_ext, tol_int,
                      filtro, raiof, simp, f, descent, lsearch)

    # Inicia o laço externo do lagrangiano aumentado
    for i_ext=1:max_ext

        # Soluciona o problema interno (e salva a derivada)
        x, dL, count, n_int = Descent(x, valor_res, mult_res, rho, xl, xu,
                                max_int, tol_int, count, fem_v, fem_f, filt, descent, lsearch)#, valor_zero)

        # Verifica novos valores da função e restrições
        valor_fun, valor_res = F_Obj(x, fem_v, fem_f)#, valor_zero)

        # Atualiza o valor o multiplicador das restricoes (u)
        mult_res = min.(mult_max, max.(rho*valor_res + mult_res, 0.0))

        # Atualiza o multiplicador de penalizacao (c)
        crho_nov = max.(0.0, norm(max.(valor_res, -mult_res/rho)))
        if crho_nov >= 0.9*crho_ant
            rho = min.(1.1*rho, rho_max)
        end # if

        # Imprime resultado atual e plota saida para o gmsh
        Imprime_Atual(x, i_ext, n_int, count, dts, valor_fun, rho, mult_res, valor_res, nelems)

        # Verifica os criterios:
        if  norm(dL)                   <= tol_ext &&  # Condicao de gradiente
            norm(max.(valor_res, 0.0)) <= tol_ext &&  # Condicao de viabilidade
            norm(valor_res'*mult_res)  <  tol_ext     # Cond. de complementariedade
            break
        end # if criterios

    end # for i_ext

    # Display Final
    Imprime_F(dts, n_int, count)
end
