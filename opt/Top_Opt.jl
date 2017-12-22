################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

# Carrega os arquivos com métodos as rotinas para otimização
include("Steepest.jl")         # Método de Descida, Steepest
include("Wall_Search.jl")      # Procura em linha, Wall_Search
include("Filtros.jl")          # Filtros de densidades e sensibilidades
include("Saida.jl")            # Impressão das saídas em arquivo e console
include("Dif_Fin.jl")          # Diferenças finitas, caso precise fazer validação

#
# Rotina principal, recebe parâmetros e executa a otimização topológica com LA
#
function Top_Opt(dts::AbstractString, caso::Int64, freq::Float64, alfa::Float64,
      beta::Float64, Ye::Float64, A::Float64, dini::Float64, max_ext::Int64, max_int::Int64,
      tol_ext::Float64, tol_int::Float64, rho::Float64, rho_max::Float64, mult_max::Float64,
              SP::Float64, raiof::Float64, vmin::Float64, NX::Int64, NY::Int64, LX::Float64,
                                LY::Float64, young::Float64, poisson::Float64, esp::Float64,
                        p_dens::Float64, presos::Array{Float64,2}, forcas::Array{Float64,2})

    # Cálcula número de nós e elementos
    nnos = (NX+1)*(NY+1)
    nel  = NX*NY

    # Gera a malha
    coord, ijk, nos_f, ID, gdll = GeraMalha(nnos, nel, LX, LY, NX, NY, presos, forcas)

    # Gera array com as forças
    F =  Global_F(ID, nos_f, gdll)

    # Gera a matrizes de rigidez e massa de um elemento, junto com a CBA para tensões
    K0, CBA = Kquad4_I(1, coord, ijk, young, poisson, esp)
    M0 = Mquad4(1, coord, ijk, esp, p_dens)

    # Prepara os vizinhos para filtros
    vizi, nviz, dviz = Proc_Vizinhos(nel, coord, ijk, raiof)

    # Pseudo-densidades para a montagem global com o SIMP
    x = dini*ones(nel)

    # Adquire valores iniciais se aplicável
    valor_0, = F_Obj(x, 0.0, [0.0], 0, nnos, nel, ijk, ID, K0, M0, SP, vmin, F, NX, NY,
                  vizi, nviz, dviz, raiof, [0.0], caso, freq, alfa, beta, A, Ye, CBA)

    # Inicializa vetores e calcula para primeiro display
    valor_fun, valor_res, to_plot, TS = F_Obj(x, 0.0, [0.0], 1, nnos, nel, ijk, ID, K0, M0,
                                    SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                          valor_0, caso, freq, alfa, beta, A, Ye, CBA)
    numres = size(valor_res,1)
    mult_res  = zeros(numres)

    # Calcula o criterio de atualizacao do c
    crho_ant = max.(0.,norm(max.(valor_res, -mult_res/rho)))

    # Prepara o plot e primeira saída
    Imprime_0(x, rho, mult_res, valor_fun, valor_res, max_ext, max_int, tol_ext,
             tol_int, dts, nnos, nel, ijk, coord, to_plot, TS, vizi, nviz, dviz, raiof)


    # Inicia o laço externo do lagrangiano aumentado
    for i_ext=1:max_ext

        # Soluciona o problema interno e salva a derivada
        x, dL = Steepest(x, rho, mult_res, max_int, tol_int, nnos, nel, ijk, ID, K0, M0,
                                     SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                      dts, valor_0, caso, freq, alfa, beta, A, Ye, CBA)

        # Verifica novos valores da função e restrições
        valor_fun, valor_res, to_plot, TS = F_Obj(x, 0.0, [0.0], 1, nnos, nel, ijk, ID, K0, M0,
                                        SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                               valor_0, caso, freq, alfa, beta, A, Ye, CBA)

        # Atualiza o valor o multiplicador das restricoes (u)
        mult_res = min.(mult_max, max.(rho*valor_res + mult_res, 0.0))

        # Atualiza o multiplicador de penalizacao (c)
        crho_nov = max.(0.0, norm(max.(valor_res, -mult_res/rho)))
        if crho_nov >= 0.9*crho_ant
            rho = min.(1.1*rho, rho_max)
        end # if

        # Imprime resultado atual e plota saida para o gmsh
        Imprime_Ext(x, rho, mult_res, valor_fun, valor_res, i_ext, dts,
                                    nel, to_plot, TS, vizi, nviz, dviz, raiof)

        # Verifica os criterios de parada:
        if  norm(dL)                   <= tol_ext &&  # Condicao de gradiente
            norm(max.(valor_res, 0.0)) <= tol_ext &&  # Condicao de viabilidade
            norm(valor_res'*mult_res)  <  tol_ext     # Cond. de complementariedade
            break
        end # if criterios

    end # for i_ext

    # Executa análise harmônica da configuração final
    println("\n\tAnálise Harmônica...")
    Harmonica(dts, 0.0, 2.5, 1000.0, alfa, beta, nel, ijk, ID, K0, M0, SP, vmin,
                                                     F, vizi, nviz, dviz, raiof)
    println(" OK!\n")
end
