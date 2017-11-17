################################################################################
# Lagrangiano Aumentado // Otimização Topológica
################################################################################
# cd(ENV["HARVEY"]);include("Main.jl");main() >> diretório da variável de sistema

# Limpa tudo
workspace();

# Carrega os arquivos com as rotinas de elementos finitos
include("fem\\Gera_Malha.jl")       # GeraMalha, gl_livres_elemento
include("fem\\Quad4.jl")            # Kquad4, Kquad4_I
include("fem\\Monta_Global.jl")     # Global, Expande_Vetor
include("fem\\Gmsh.jl")             # Funções GMSH

# Carrega os arquivos com as rotinas para otimização
include("Steepest.jl")          # Método de Descida
include("Line_Search.jl")       # Método de busca em linha

# Carrega cálculo das derivadas
include("Fobj_Dinamico_Est.jl")         # Fobj e Sensibilidades
include("Dif_Fin.jl")         # Fobj e Sensibilidades

# Etc
include("Filtros.jl")               # Filtros de densidades e sensibilidades
include("Saida.jl")                 # Impressão das saídas em arquivo e console

type finitos
    KG; CG; MG; KD; F; UE; UD;
    K0; M0; simp;
    nelems; conect; NX; NY
    ID; nforcas; nos_forcas; nr_gl_livres;
    w; alfa; beta
end

type filtros
    raiof; vizi; nviz; dviz; filtro
end

function main()

    # Horario de Execução
    dtf = Dates.now()
    dts = Dates.format(dtf,"YYYY-mm-dd-HH-MM-SS")
    fname = string("z",dts,".pos")
    tic()

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 100       # Máximo de iteracoes externas
    max_int     = 200       # Máximo de iterações internas
    tol_ext     = 1E-5      # Tolerância do laço externo
    tol_int     = 1E-5      # Tolerância do laço interno
    rho_ini     = 0.2
    rho_max     = 5.0       # Valor maximo de rho
    mult_max    = 10.0      # Valor maximo dos multiplicadores

    # Parâmetros da topológica
    dens_ini    = 0.49      # Volume/Pseudo-densidades iniciais
    simp        = 3.0       # Parametro p do SIMP
    raiof       = 0.031     # Tamanho do filtro [m]
    filtro      = "Dens"    # Filtros: (Off ou Dens)

    # Parâmetros do problema Harmônico
    f    = 180.0      # Frequencia
    alfa = 0.0      # Amortecimento proporcional de Rayleigh
    beta = 1E-8

    # Malha:
    NX = 60   #80        # Nr. de elementos em X
    NY = 30   #40        # Nr. de elementos em Y

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

    # Dados calculados e declaração de variáveis:
    npresos = size(presos,1)            # Nr. de apoios
    nforcas = size(forcas,1)            # Nr. de carregamentos
    nelems  = NX*NY                     # Nr. de elementos
    nnos    = (NX+1)*(NY+1)             # Nr. de nós

    # Pseudo-densidades para a montagem global com o SIMP
    x  = dens_ini*ones(nelems)
    xl = 0.00*ones(nelems)
    xu = 1.00*ones(nelems)

    # Gera a malha
    coord, conect, nos_forcas, ID, nr_gl_livres =
        GeraMalha(nnos, nelems, LX, LY, NX, NY, npresos, presos, nforcas, forcas)

    # Gera a matrizes de rigidez e massa de um elemento finito - malha toda igual
    (K0,) = Kquad4_I(1, coord, conect, young, poisson, esp)
    M0    = Mquad4(1, coord, conect, esp, rho)

    # Prepara os vizinhos para filtros 63
    vizi,nviz,dviz = Proc_Vizinhos(nelems, coord, conect, raiof)
    filt           = filtros(raiof, vizi, nviz, dviz, filtro)

    # Correção de Hz para rad/s
    w = 2.0*pi*f

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    KG,MG,F = Global(nelems, conect, ID, K0, M0, x, simp, nforcas, nos_forcas, nr_gl_livres)

    # Monta CG e KD (Harmônica)
    CG = alfa*MG + beta*KG
    KD = KG + w*im*CG - w^2.0*MG

    # Resolve o sistema pela primeira vez
    UE = vec(lufact(KG)\F);
    UD = vec(lufact(KD)\F);
    fem = finitos(KG, CG, MG, KD, F, UE, UD, K0, M0, simp, nelems, conect, NX, NY,
                            ID, nforcas, nos_forcas, nr_gl_livres, w, alfa, beta)

    # Obtem os valores de f e g em x0 para o calculo de c0
    valor_fun, valor_res, fem = F_Obj(x, fem, 0.0)

    # Salva o valor da primeira função para
    valor_zero = copy(valor_fun)
    #valor_fun  = 1.0

    # Inicializa multiplicadores de Lagrange (u)
    numres   = size(valor_res, 1)
    mult_res = zeros(numres)

    # Define fator de penalização inicial (c)
    rho = max.(1E-6,min(rho_max,(2.0*abs(valor_fun)/norm(max.(valor_res,0.))^2.)))
    rho = rho_ini

    # E calcula o criterio de atualizacao do c
    crho_ant = max.(0.,norm(max.(valor_res, -mult_res/rho)))

    # Inicializa o contador de avaliacoes da F_Obj
    count = 1
    n_int = 0
    viewcount = 0.0

    # Prepara o plot e primeira saída
    Inicializa_Malha_Gmsh(fname, nnos, nelems, conect, coord, 2)
    Adiciona_Vista_Escalar_Gmsh(fname, "x", nelems, x, 0.0)
    Imprime(0, x, rho, mult_res, valor_fun, valor_res, dts, nelems,
            max_ext, max_int, tol_ext, tol_int, filtro, raiof, simp, f, n_int,count)

    # Inicia o laço externo do lagrangiano aumentado
    for i_ext=1:max_ext

        # Soluciona o problema interno (e salva a derivada)
        x, dL, count, fem, n_int = Steepest(x, valor_res, mult_res, rho, xl, xu,
                                max_int, tol_int, count, fem, filt, valor_zero)

        # Verifica novos valores da função e restrições
        valor_fun, valor_res, fem = F_Obj(x, fem, valor_zero)

        # Atualiza o valor o multiplicador das restricoes (u)
        mult_res = min.(mult_max, max.(rho*valor_res + mult_res, 0.0))

        # Atualiza o multiplicador de penalizacao (c)
        crho_nov = max.(0., norm(max.(valor_res, -mult_res/rho)))
        if crho_nov >= .9*crho_ant
            rho = min.(1.1*rho, rho_max)
        end # if

        # Imprime resultado atual e plota saida para o gmsh
        Adiciona_Vista_Escalar_Gmsh(fname, "xf", nelems, x, Float64(i_ext))
        Imprime(i_ext, x, rho, mult_res, valor_fun, valor_res, dts, nelems,
                max_ext, max_int, tol_ext, tol_int, filtro, raiof, simp, f,n_int,count)

        # Verifica os criterios:
        if  norm(dL)                   <= tol_ext &&  # Condicao de gradiente
            norm(max.(valor_res, 0.0)) <= tol_ext &&  # Condicao de viabilidade
            norm(valor_res'*mult_res)  <  tol_ext     # Cond. de complementariedade
            break
        end # if criterios

    end # for i_ext

    # Display Final
    Imprime(-1, x, rho, mult_res, valor_fun, valor_res, dts, nelems,
            max_ext, max_int, tol_ext, tol_int, filtro, raiof, simp, f,n_int,count)
end
