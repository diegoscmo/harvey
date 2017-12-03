################################################################################
# Lagrangiano Aumentado // Otimização Topológica
################################################################################

# Copiar para o console para executar:
# cd(ENV["HARVEY"]);
# include("main_all.jl"); main_all();

# Carrega os arquivos com as rotinas de elementos finitos
include("fem\\Gera_Malha.jl")       # GeraMalha, gl_livres_elemento
include("fem\\Quad4.jl")            # Kquad4, Kquad4_I
include("fem\\Monta_Global.jl")     # Global, Expande_Vetor
include("fem\\Gmsh.jl")             # Funções GMSH
include("fem\\Harmonica.jl")         # Cálculo Harmônico

# Carrega os arquivos com métodos as rotinas para otimização
include("opt\\Steepest_All.jl")         # Método de Descida, Steepest
include("opt\\Wall_Search_All.jl")      # Procura em linha, Wall_Search
include("opt\\Filtros.jl")          # Filtros de densidades e sensibilidades
include("opt\\Saida.jl")            # Impressão das saídas em arquivo e console

# Carrega cálculo da Fobj, Lagrangiana e derivadas, selecionar um
include("fobj\\0_Todas_Fobj.jl")

# Rotina principal
function main_all(num::Int64, caso::Int64, freq::Float64, alfa::Float64, beta::Float64, Ye::Float64, A::Float64)

    # Nome do Arquivo ou data de execução("OFF desliga")
    dts = string(num,"_",caso,"_f",freq,"_Y",Ye,"_A",A)

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 30       # Máximo de iteracoes externas
    max_int     = 500       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho_ini     = 1.00      # Valor inicial de rho
    rho_max     = 2.00      # Valor maximo de rho
    mult_max    = 10.0      # Valor maximo dos multiplicadores

    # Parâmetros da topológica
    dens_ini    = 0.49      # Volume/Pseudo-densidades iniciais
    SP          = 3.00      # Parametro p do SIMP
    raiof       = 0.031     # Tamanho do filtro [m]
    vmin        = 0.001     # 1%

    # Parâmetros do problema de FEM, 60x30 = 1800
    NX = 60              # Nr. de elementos em X
    NY = 30              # Nr. de elementos em Y

    # Definição geométrica do problema (retângulo):
    LX = 1.0       # Comprimento em X
    LY = 0.5       # Comprimento em Y

    # Material:
    young   = 210E9     # Módulo de Young
    poisson = 0.0       # Coeficiente de Poisson (0.0 > viga)
    esp     = 1.0       # Espessura do retângulo
    p_dens  = 7860.0    # Densidade

    # Restrições de deslocamento (apoios):
    #        [ ponto_X_inicial ponto_Y_inicial ponto_X_final ponto_Y_final direção (X=1 Y=2)]
    presos = [  0.0             0.0             0.0           LY           1       ;
    0.0             0.0             0.0           LY           2       ]

    # Carregamentos:
    #               [ ponto_X        ponto_Y         força           dir (X=1 Y=2)]
    forcas = [ LX              LY/2.0          -9000.0         2 ]

    # Nós e elementos
    nnos = (NX+1)*(NY+1)             # Nr. de nós
    nel  = NX*NY                     # Nr. de elementos

    # Gera a malha
    coord, ijk, nos_f, ID, gdll = GeraMalha(nnos, nel, LX, LY, NX, NY, presos, forcas)

    # Gera array com as forças
    F =  Global_F(ID, nos_f, gdll)

    # Gera a matrizes de rigidez e massa de um elemento e transforma para StaticArrays
    K0 = Kquad4_I(1, coord, ijk, young, poisson, esp)[1]
    M0 = Mquad4(1, coord, ijk, esp, p_dens)

    # Prepara os vizinhos para filtros 63
    vizi, nviz, dviz = Proc_Vizinhos(nel, coord, ijk, raiof)

    # Pseudo-densidades para a montagem global com o SIMP
    x = dens_ini*ones(nel)

    # Adquire valores iniciais se aplicável
    valor_0 = F_Obj(x, 0.0, [0.0], 0, nel, ijk, ID, K0, M0, SP, vmin,
    F, NX, NY, vizi, nviz, dviz, raiof, [0.0],
    caso, freq, alfa, beta, Ye, A)

    # Inicializa vetores e calcula para primeiro display
    valor_fun, valor_res = F_Obj(x, 0.0, [0.0], 1, nel, ijk, ID, K0, M0, SP, vmin,
    F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
    caso, freq, alfa, beta, Ye, A)
    numres = size(valor_res,1)
    mult_res  = zeros(Float64,numres)

    # Define fator de penalização inicial (c)
    rho = max.(1E-2,min(rho_max,(2.0*abs(valor_fun)/norm(max.(valor_res,0.))^2.)))

    # Usa rho_min se a função começar não restrita
    if rho == rho_max
        rho = rho_ini
    end

    # E calcula o criterio de atualizacao do c
    crho_ant = max.(0.,norm(max.(valor_res, -mult_res/rho)))

    # Prepara o plot e primeira saída
    Imprime_0(x, rho, mult_res, valor_fun, valor_res, max_ext, max_int, tol_ext,
    tol_int, dts, nnos, nel, ijk, coord)

    # Inicia o laço externo do lagrangiano aumentado
    for i_ext=1:max_ext

        # Soluciona o problema interno (e salva a derivada)
        x, dL = Steepest(x, rho, mult_res, max_int, tol_int, nel, ijk, ID, K0, M0,
        SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof, dts, valor_0,
        caso, freq, alfa, beta, Ye, A)

        # Verifica novos valores da função e restrições
        valor_fun, valor_res = F_Obj(x, 0.0, [0.0], 1, nel, ijk, ID, K0, M0, SP,
        vmin, F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
        caso, freq, alfa, beta, Ye, A)

        # Atualiza o valor o multiplicador das restricoes (u)
        mult_res = min.(mult_max, max.(rho*valor_res + mult_res, 0.0))

        # Atualiza o multiplicador de penalizacao (c)
        crho_nov = max.(0.0, norm(max.(valor_res, -mult_res/rho)))
        if crho_nov >= 0.9*crho_ant
            rho = min.(1.1*rho, rho_max)
        end # if

        # Imprime resultado atual e plota saida para o gmsh
        Imprime_Ext(x, rho, mult_res, valor_fun, valor_res, i_ext, dts, nel)

        # Verifica os criterios:
        if  norm(dL)                   <= tol_ext &&  # Condicao de gradiente
            norm(max.(valor_res, 0.0)) <= tol_ext &&  # Condicao de viabilidade
            norm(valor_res'*mult_res)  <  tol_ext     # Cond. de complementariedade
            break
        end # if criterios

    end # for i_ext

    @printf("\n\tAnálise Harmônica...")
    alfa = 0.0
    beta = 1E-8
    Harmonica(dts, 0.0, 5.0, 1000.0, alfa, beta, nel, ijk, ID, K0, M0, SP, vmin, F,
    vizi, nviz, dviz, raiof)

    @printf(" OK!\n")
end

function main()
    num = 1
    caso = 1
    freq = 0.0
    alfa = 0.0
    beta = 1E-8
    Ye = 1.0
    A = 1.0

    for num = 63:70 #1:86 
        if num == 1
            caso = 1
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 2
            caso = 2
            freq = 5.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 3
            caso = 2
            freq = 100.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 4
            caso = 2
            freq = 170.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 5
            caso = 2
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 6
            caso = 3
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            Ye = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 7
            caso = 3
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            Ye = 0.75
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 8
            caso = 3
            freq = 600.0
            alfa = 0.0
            beta = 1E-8
            Ye = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 9
            caso = 3
            freq = 600.0
            alfa = 0.0
            beta = 1E-8
            Ye = 0.5
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 10
            caso = 3
            freq = 600.0
            alfa = 0.0
            beta = 1E-8
            Ye = 0.3
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 11
            caso = 3
            freq = 600.0
            alfa = 0.0
            beta = 1E-8
            Ye = 0.2
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 12
            caso = 3
            freq = 750.0
            alfa = 0.0
            beta = 1E-8
            Ye = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 13
            caso = 3
            freq = 750.0
            alfa = 0.0
            beta = 1E-8
            Ye = 0.5
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 14
            caso = 3
            freq = 750.0
            alfa = 0.0
            beta = 1E-8
            Ye = 0.3
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 15
            caso = 3
            freq = 750.0
            alfa = 0.0
            beta = 1E-8
            Ye = 0.2
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 16
            caso = 4
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            A    = 0.99
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 17
            caso = 4
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            A    = 0.90
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 18
            caso = 4
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            A    = 0.50
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 19
            caso = 4
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            A    = 0.25
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 20
            caso = 4
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            A    = 0.05
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 21
            caso = 4
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            A    = 0.01
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 22
            caso = 2
            freq = 180.0
            alfa = 0.0
            beta = 1E-2
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 23
            caso = 2
            freq = 750.0
            alfa = 0.0
            beta = 1E-2
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 24
            caso = 2
            freq = 180.0
            beta = 0.0
            alfa = 1000.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 25
            caso = 2
            freq = 750.0
            beta = 0.0
            alfa = 1000.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 26
            caso = 2
            freq = 750.0
            beta = 1E-3
            alfa = 0.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 27
            caso = 2
            freq = 180.0
            beta = 1E-8
            alfa = 0.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 28
            caso = 2
            freq = 180.0
            beta = 1E-4
            alfa = 0.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 29
            caso = 2
            freq = 750.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 30
            caso = 2
            freq = 750.0
            alfa = 0.0
            beta = 1E-4
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 31
            caso = 5
            freq = 005.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 32
            caso = 5
            freq = 100.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 33
            caso = 5
            freq = 170.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 34
            caso = 5
            freq = 180.0
            alfa = 0.0
            beta = 1E-8
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 35
            caso = 5
            freq = 180.0
            w    = 2*pi*freq
            beta = 0.1/w
            alfa = 0.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 36
            caso = 5
            freq = 600.0
            w    = 2*pi*freq
            beta = 0.1/w
            alfa = 0.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 37
            caso = 5
            freq = 750.0
            w    = 2*pi*freq
            beta = 0.1/w
            alfa = 0.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 38
            caso = 5
            freq = 005.0
            beta = 0.0
            alfa = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 39
            caso = 5
            freq = 100.0
            beta = 0.0
            alfa = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 40
            caso = 5
            freq = 170.0
            beta = 0.0
            alfa = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 41
            caso = 5
            freq = 180.0
            beta = 0.0
            alfa = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 42
            caso = 5
            freq = 170.0
            w    = 2.0*pi*freq
            beta = 0.0
            alfa = 0.8*w
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 43
            caso = 5
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.0
            alfa = 0.8*w
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 44
            caso = 5
            freq = 1350.0
            w    = 2.0*pi*freq
            beta = 0.0
            alfa = 0.8*w
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 45
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 46
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.0
            alfa = 0.8*w
            Ye   = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 47
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 1.5
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 48
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 2.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 49
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 5.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 50
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 10.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 51
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 20.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 52
            caso = 6
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 0.5
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 53
            caso = 6
            freq = 600.0
            beta = 1E-8
            alfa = 0.0
            Ye   = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 54
            caso = 6
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 55
            caso = 6
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 0.3
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 56
            caso = 6
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 0.5
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 57
            caso = 6
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 2.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 58
            caso = 6
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 5.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 59
            caso = 6
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 10.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 60
            caso = 6
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 15.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 61
            caso = 6
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 1.0
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 62
            caso = 6
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 0.5
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 63
            caso = 6
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            Ye   = 0.3
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 64
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.99
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 65
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.95
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 66
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.90
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 67
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.80
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 68
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.70
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 69
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.50
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 70
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.30
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 71
            caso = 7
            freq = 750.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.10
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 72
            caso = 6
            freq = 1450.0
            beta = 1E-8
            alfa = 0.0
            Ye   = 1.00
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 73
            caso = 7
            freq = 1450.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.95
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 74
            caso = 7
            freq = 1450.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.90
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 75
            caso = 7
            freq = 1450.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.50
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 76
            caso = 7
            freq = 1450.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.30
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 77
            caso = 7
            freq = 1450.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.10
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 78
            caso = 7
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.999
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 79
            caso = 7
            freq = 180.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.99
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 80
            caso = 7
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.999
            main_all(num, caso, freq, alfa, beta, Ye, A)
        elseif num == 81
            caso = 7
            freq = 600.0
            w    = 2.0*pi*freq
            beta = 0.1/w
            alfa = 0.0
            A    = 0.99
            main_all(num, caso, freq, alfa, beta, Ye, A)
        end
    end
    quit()
end
