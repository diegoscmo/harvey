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
include("opt\\Steepest.jl")         # Método de Descida, Steepest
include("opt\\Wall_Search.jl")      # Procura em linha, Wall_Search
include("opt\\Filtros.jl")          # Filtros de densidades e sensibilidades
include("opt\\Saida.jl")            # Impressão das saídas em arquivo e console

# Carrega cálculo da Fobj, Lagrangiana e derivadas, selecionar um
#include("fobj\\1_Estatico.jl")
#include("fobj\\2_Dinamico.jl")
#include("fobj\\3_Din_R-Est.jl")
#include("fobj\\4_Din_A-Est.jl")
#include("fobj\\5_Potencia.jl")
#include("fobj\\6_Pot_R-Est.jl")
include("fobj\\7_Pot_A-Est.jl")


# Rotina principal
function main(caso::Int64)

    # Nome do Arquivo ou data de execução("OFF desliga")
    #dts = "1_Est_1800"
    #dts = "3_Din_R-Est_1800_f180-R100"
    #dts = "4_Din_A-Est_1800_f180-A099"
    #dts = "5_Pot_1800_f180_alfa1"
    #dts = "6_Pot_R-Est_1800_f600_Y100"
    #dts  = "7_Pot_A-Est_1800_f600_A0999"
    dts = "teste"

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 50       # Máximo de iteracoes externas
    max_int     = 300       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho_ini     = 1.00      # Valor inicial de rho
    rho_max     = 3.00      # Valor maximo de rho
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
                                                F, NX, NY, vizi, nviz, dviz, raiof, [0.0])

    # Inicializa vetores e calcula para primeiro display
    valor_fun, valor_res = F_Obj(x, 0.0, [0.0], 1, nel, ijk, ID, K0, M0, SP, vmin,
                                                F, NX, NY, vizi, nviz, dviz, raiof, valor_0)
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
                                 SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof, dts, valor_0)

        # Verifica novos valores da função e restrições
        valor_fun, valor_res = F_Obj(x, 0.0, [0.0], 1, nel, ijk, ID, K0, M0, SP,
                                        vmin, F, NX, NY, vizi, nviz, dviz, raiof, valor_0)

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
