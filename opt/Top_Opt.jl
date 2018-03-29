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
function Top_Opt(dts::AbstractString, Sy::Float64, freq::Float64, alfa::Float64,
      beta::Float64, Ye::Float64, A::Float64, dini::Float64, max_ext::Int64, max_int::Int64,
      tol_ext::Float64, tol_int::Float64, rho::Float64, rho_max::Float64, mult_max::Float64,
              SP::Float64, raiof::Float64, vmin::Float64, NX::Int64, NY::Int64, LX::Float64,
               LY::Float64, young::Float64, poisson::Float64, esp::Float64, p_dens::Float64,
              presos::Array{Float64,2}, forcas::Array{Float64,2}, QP::Float64, csi::Float64, dmax::Float64)

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
    valor_0, n_rho = F_Obj(x, [1.0], [0.0], 0, nnos, nel, ijk, ID, K0, M0, SP, vmin, F, NX, NY,
                  vizi, nviz, dviz, raiof, [1.0], Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax)

    # Inicializa os rhos
    rho = rho*ones(n_rho)

    # Inicializa vetores e calcula para primeiro display
    v_fun, v_res, to_plot, sigma = F_Obj(x, rho, [0.0], 1, nnos, nel, ijk, ID, K0, M0,
                                    SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                          valor_0, Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax)

    # Inicializa multiplicadores de penalização
    n_res   = size(v_res,1)
    mu_res  = zeros(n_res)

    # Calcula o criterio de atualizacao do c, começando pelos individuais
    crho_ant = zeros(n_rho)
    for j = 1:n_rho-1
        crho_ant[j] = max.(0.,norm(max.(v_res[j], -mu_res[j]/rho[j])))
    end
    # O ultimo critério engloba os demais (tensão)
    crho_ant[n_rho] = max.(0.,norm(max.(v_res[n_rho:end], -mu_res[n_rho:end]/rho[n_rho])))

    # Carrega execução anterior se houver
    itex = 1
    if isfile(string(dts,"_save.txt"))

        @printf("Utilizar execução anterior? Y/N")
        conf = input()
        if conf == "Y" || conf ==  "y"
            # Lê arquivo com densidades x
            x = Le_Densidades(string(dts,"_dens_o.pos"), nel)

            # Lê arquivo com execução anterior
            itex,rho,mu_res = Le_Ultima(dts,n_rho)

            println("AVISO: Carregando execução anterior, cuidado com as condições de contorno!")
        else
            Imprime_0(x, rho, mu_res, v_fun, v_res, max_ext, max_int, tol_ext,
                      tol_int, dts, nnos, nel, ijk, coord, to_plot, real(sigma), vizi, nviz, dviz, raiof)
        end

    else
        # Prepara o plot e primeira saída
        Imprime_0(x, rho, mu_res, v_fun, v_res, max_ext, max_int, tol_ext,
                  tol_int, dts, nnos, nel, ijk, coord, to_plot, real(sigma), vizi, nviz, dviz, raiof)
    end

    # Inicia o laço externo do lagrangiano aumentado
    for i_ext=itex:max_ext

        # Soluciona o problema interno e salva a derivada
        x, dL = Steepest(x, rho, mu_res, max_int, tol_int, nnos, nel, ijk, ID, K0, M0,
                                     SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                      dts, valor_0, Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax)

        # Verifica novos valores da função e restrições
        v_fun, v_res, to_plot, sigma = F_Obj(x, rho, mu_res, 1, nnos, nel, ijk, ID, K0, M0,
                                        SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof,
                                               valor_0, Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax)

        # Atualiza o valor o multiplicador das restricoes (u) e penalização (c)
        # Começando pelos individuais
        crho_nov = zeros(n_rho)
        for j=1:(n_rho-1)
            mu_res[j] = min.(mult_max, max.(rho[j]*v_res[j] + mu_res[j], 0.0))
            crho_nov[j] = max.(0.0, norm(max.(v_res[j], -mu_res[j]/rho[j])))

            if crho_nov[j] >= 0.9*crho_ant[j]
                rho[j] = min.(1.1*rho[j], rho_max)
            end # if
        end #for j

        # Depois os restantes
        mu_res[n_rho:end] = min.(mult_max, max.(rho[n_rho]*v_res[n_rho:end] + mu_res[n_rho:end], 0.0))
        crho_nov[n_rho] = max.(0.0, norm(max.(v_res[n_rho:end], -mu_res[n_rho:end]/rho[n_rho])))

        if crho_nov[n_rho] >= 0.9*crho_ant[n_rho]
            rho[n_rho] = min.(1.1*rho[n_rho], rho_max)
        end # if

        # Imprime resultado atual e plota saida para o gmsh
        Imprime_Ext(x, rho, mu_res, v_fun, v_res, i_ext, dts,
                                    nel, to_plot, real(sigma), vizi, nviz, dviz, raiof)

        # Verifica os criterios de parada:
        if  norm(dL)               <= tol_ext &&  # Condicao de gradiente
            norm(max.(v_res, 0.0)) <= tol_ext &&  # Condicao de viabilidade
            norm(v_res'*mu_res)    <  tol_ext     # Cond. de complementariedade
            break
        end # if criterios

    end # for i_ext

    # Análise na frequência
    fvarredura    = [0.5, 2.5, 1500.00]
    nos_monitor   = [nel 2]
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    arquivo_saida = string(dts,"_modal.pos")
    nmodos = 6

    Analise_Modal(nmodos,K0,M0,nel,nnos,ijk,ID,coord,arquivo_saida,xc,SP)
    #Analise_Harmonica(freq,arquivo_saida,nel,nnos,ijk,coord,ID,K0,M0,nos_f,
    #                       gdll,alfa,beta,xc,SP)
    arquivo_saida = string(dts,"_varred.pos")
    Varredura(arquivo_saida,fvarredura,nos_monitor,nel,nnos,
                            ijk,coord,ID,K0,M0,nos_f,gdll,alfa,beta,xc,SP)

    println(" OK!\n")
end
