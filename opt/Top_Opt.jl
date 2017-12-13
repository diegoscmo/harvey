################################################################################
# Lagrangiano Aumentado // Otimização Topológica
################################################################################

# Rotina principal
function Top_Opt(num::Int64, caso::Int64, freq::Float64, alfa::Float64,
                        beta::Float64, Ye::Float64, A::Float64, dini::Float64)

    # Nome do Arquivo ou data de execução("OFF desliga")
    dts = ""
    if caso == 0
        dts = "Teste"
    elseif caso == 1
        dts = string(num,"_",caso)
    elseif caso == 2
        dts = string(num,"_",caso,"_f",freq)
    elseif caso == 3
        dts = string(num,"_",caso,"_f",freq,"_Y",Ye)
    elseif caso == 4
        dts = string(num,"_",caso,"_f",freq,"_A",A)
    elseif caso == 5
        dts = string(num,"_",caso,"_f",freq)
    elseif caso == 6
        dts = string(num,"_",caso,"_f",freq,"_Y",Ye)
    elseif caso == 7
        dts = string(num,"_",caso,"_f",freq,"_A",A)
    elseif caso == 8
        dts = string(num,"_",caso,"_f",freq,"_A",A,"_R",Ye)
    else
        error("ERRO NO CASO")
    end

    # Parâmetros do Lagrangiano Aumentado
    max_ext     = 30       # Máximo de iteracoes externas
    max_int     = 500       # Máximo de iterações internas
    tol_ext     = 1E-6      # Tolerância do laço externo
    tol_int     = 1E-6      # Tolerância do laço interno
    rho         = 1.00      # Valor inicial de rho
    rho_max     = 2.50      # Valor maximo de rho
    mult_max    = 10.0      # Valor maximo dos multiplicadores

    # Parâmetros da topológica
    dens_ini    = dini      # Volume/Pseudo-densidades iniciais
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
#                LX              0.0             LX            LY           1       ;
#                LX              0.0             LX            LY           2       ]

    # Carregamentos:
    #               [ ponto_X        ponto_Y         força           dir (X=1 Y=2)]
    forcas = [ LX              LY/2.0          -9000.0         2 ]
    #forcas = [ LX/2.0              0.0          -9000.0         2 ]

    # Nós e elementos
    nnos = (NX+1)*(NY+1)             # Nr. de nós
    nel  = NX*NY                     # Nr. de elementos

    # Gera a malha
    coord, ijk, nos_f, ID, gdll = GeraMalha(nnos, nel, LX, LY, NX, NY, presos, forcas)

    # Gera array com as forças
    F =  Global_F(ID, nos_f, gdll)

    # Gera a matrizes de rigidez e massa de um elemento e transforma para StaticArrays
    K0 = Kquad4_I(1, coord, ijk, young, poisson, esp)
    M0 = Mquad4(1, coord, ijk, esp, p_dens)

    # Prepara os vizinhos para filtros 63
    vizi, nviz, dviz = Proc_Vizinhos(nel, coord, ijk, raiof)

    # Pseudo-densidades para a montagem global com o SIMP
    x = dens_ini*ones(nel)

    # Adquire valores iniciais se aplicável
    valor_0 = F_Obj(x, 0.0, [0.0], 0, nel, ijk, ID, K0, M0, SP, vmin, F, NX, NY,
                  vizi, nviz, dviz, raiof, [0.0], caso, freq, alfa, beta, A, Ye)

    # Inicializa vetores e calcula para primeiro display
    valor_fun, valor_res, to_plot = F_Obj(x, 0.0, [0.0], 1, nel, ijk, ID, K0, M0,
                                    SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                          valor_0, caso, freq, alfa, beta, A, Ye)
    numres = size(valor_res,1)
    mult_res  = zeros(Float64,numres)

    # Calcula o criterio de atualizacao do c
    crho_ant = max.(0.,norm(max.(valor_res, -mult_res/rho)))

    # Prepara o plot e primeira saída
    Imprime_0(x, rho, mult_res, valor_fun, valor_res, max_ext, max_int, tol_ext,
                                   tol_int, dts, nnos, nel, ijk, coord, to_plot)

    # Inicia o laço externo do lagrangiano aumentado
    for i_ext=1:max_ext

        # Soluciona o problema interno (e salva a derivada)
        x, dL = Steepest(x, rho, mult_res, max_int, tol_int, nel, ijk, ID, K0, M0,
                                     SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                      dts, valor_0, caso, freq, alfa, beta, A, Ye)

        # Verifica novos valores da função e restrições
        valor_fun, valor_res, to_plot = F_Obj(x, 0.0, [0.0], 1, nel, ijk, ID, K0, M0,
                                        SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                               valor_0, caso, freq, alfa, beta, A, Ye)

        # Atualiza o valor o multiplicador das restricoes (u)
        mult_res = min.(mult_max, max.(rho*valor_res + mult_res, 0.0))

        # Atualiza o multiplicador de penalizacao (c)
        crho_nov = max.(0.0, norm(max.(valor_res, -mult_res/rho)))
        if crho_nov >= 0.9*crho_ant
            rho = min.(1.1*rho, rho_max)
        end # if

        # Imprime resultado atual e plota saida para o gmsh
        Imprime_Ext(x, rho, mult_res, valor_fun, valor_res, i_ext, dts, nel, to_plot)

        # Verifica os criterios:
        if  norm(dL)                   <= tol_ext &&  # Condicao de gradiente
            norm(max.(valor_res, 0.0)) <= tol_ext &&  # Condicao de viabilidade
            norm(valor_res'*mult_res)  <  tol_ext     # Cond. de complementariedade
            break
        end # if criterios

    end # for i_ext

    @printf("\n\tAnálise Harmônica...")
    Harmonica(dts, 0.0, 2.0, 1000.0, alfa, beta, nel, ijk, ID, K0, M0, SP, vmin,
                                                     F, vizi, nviz, dviz, raiof)

    @printf(" OK!\n")
end
