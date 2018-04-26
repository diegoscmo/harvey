################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

# Carrega os arquivos com métodos as rotinas para otimização
include("steepest.jl")         # Método de Descida, Steepest
include("wall_search.jl")      # Procura em linha, Wall_Search
include("filtros.jl")          # Filtros de densidades e sensibilidades
include("saida.jl")            # Impressão das saídas em arquivo e console
include("dif_fin.jl")          # Diferenças finitas, caso precise fazer validação

#
# Rotina principal, recebe parâmetros e executa a otimização topológica com LA
#
function Top_Opt(dts::AbstractString, Sy::Float64, freq::Float64, alfa::Float64, beta::Float64, R_bar::Float64,
                 A::Float64, dini::Float64, max_ext::Int64, max_int::Int64, tol_ext::Float64, tol_int::Float64,
                 rho::Array{Float64,1}, rho_max::Float64, SP::Float64, raiof::Float64,
                 vmin::Float64, NX::Int64, NY::Int64, LX::Float64, LY::Float64, young::Float64, poisson::Float64,
                 esp::Float64, p_dens::Float64, presos::Array{Float64,2}, forcas::Array{Float64,2},
                 travas::Array{Float64,2}, QP::Float64, csi::Float64, dmax::Float64, P, q)

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

    # Prepara os vizinhos de nós para normas
    nos_viz = Proc_Vizinhos_Nos(ijk, nnos, nel)

    # Prepara elementos a travar
    trava_els = Proc_Travas(nel, coord, ijk, travas)

    # Pseudo-densidades para a montagem global com o SIMP, já trava elementos
    x = dini*ones(nel)
    x = Trava_Els(x, trava_els)

    # Abre diretório pras saídas
    if isdir(string("results/",dts)) == false
        mkdir(string("results/",dts))
    end

    # Adquire valores iniciais se aplicável
    valor_0, smu = F_Obj(x, [1.0], [0.0], 0, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin, F,
                              NX, NY, vizi, nviz, dviz, raiof, [1.0], Sy, freq, alfa, beta, A,
                                                      R_bar, CBA, QP, csi, dmax, nos_viz, dts,P,q)

    # Inicializa multiplicadores de penalização
    mu_res  = zeros(smu)
    n_rho   = size(rho,1)
    itex    = 1

    # Inicializa vetores e calcula para primeiro display
    v_fun, v_res, stat, sigma = F_Obj(x, rho, mu_res, 1, nnos, nel, ijk, coord, ID, K0, M0,
                                   SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof, valor_0, Sy,
                                  freq, alfa, beta, A, R_bar, CBA, QP, csi, dmax,nos_viz,dts,P,q)

    # Calcula o criterio de atualizacao do c, começando pelos individuais
    crit_rho = Atualiza_Rho(rho, zeros(n_rho), rho_max, v_res, mu_res, true)

    # Verifica se tem save_file
    load_flag = Load_Game(dts, nel, max_ext, max_int)

    #  Se tiver, ou continua ou faz o primeiro display pra começar do zero
    if load_flag == true

        # Carrega tudo
        x, itex, rho, mu_res = Le_Ultima(dts, nel, n_rho)

    else
        # Faz primeira impressão
        Imprime_0(x, rho, mu_res, v_fun, v_res, max_ext, max_int, tol_ext, tol_int, dts,
                   nnos, nel, ijk, coord, [0.0;0.0;stat], real(sigma), vizi, nviz, dviz, raiof)

        # Registra as variaveis em um arquivo
        Imprime_Caso(dts, Sy, freq, alfa, beta, R_bar, A, dini, max_ext, max_int, tol_ext,
                         tol_int, rho, rho_max, SP, raiof, vmin, NX, NY, LX, LY,
                               young, poisson, esp, p_dens, presos, forcas, QP, csi, dmax,P,q)
    end


    # Inicia o laço externo do lagrangiano aumentado
    for i_ext=itex:max_ext

        # Soluciona o problema interno e salva a derivada
        x, dL = Steepest(x, rho, mu_res, max_int, tol_int, nnos, nel, ijk, coord, ID, K0, M0,
                         SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof, dts, valor_0, Sy, freq,
                         alfa, beta, A, R_bar, CBA, QP, csi, dmax, trava_els, nos_viz, i_ext,P,q)

        # Verifica novos valores da função e restrições
        v_fun,v_res,stat,sigma = F_Obj(x, rho, mu_res, 1, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                          F, NX, NY, vizi, nviz, dviz, raiof, valor_0, Sy, freq, alfa,
                                          beta, A, R_bar, CBA, QP, csi, dmax, nos_viz, dts,P,q)

        # Atualiza o valor o multiplicador das restricoes (u) e penalização (c)
        rho, crit_rho, mu_res = Atualiza_Rho(rho, crit_rho, rho_max, v_res, mu_res)

        # Mini-exibição a cada 5 iterações
        if i_ext%5 == 0
            println("\n  LagAug:: ",dts)
        end

        # Imprime resultado atual e plota saida para o gmsh
        Imprime_Ext(x, rho, mu_res, v_fun, v_res, i_ext, dts, nel,
                    [i_ext;norm(dL);stat], real(sigma), vizi, nviz, dviz, raiof)

        # Verifica os criterios de parada:
        if  norm(dL)               <= tol_ext &&  # Condicao de gradiente
            norm(max.(v_res, 0.0)) <= tol_ext &&  # Condicao de viabilidade
            norm(v_res'*mu_res)    <  tol_ext     # Cond. de complementariedade
            break
        end # if criterios

    end # for i_ext

    # Agora que acabou a otimização, fazemos barba, cabelo e bigode
    fvarredura = 0.0:2.0:1500.0

    # Modal, Potência Ativa, Norma P, R....
    Varredura_Graph(dts, fvarredura, nel, nnos, ijk, coord,vizi, nviz, dviz, raiof,ID, K0, M0, F, alfa, beta, x, P, q, nos_viz, SP, vmin)

    # Diferenças finitas, só pra se sobrar tempo...
    #F_Obj_DFC(x, rho, mu_res, 3, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin, F,
    #                NX, NY, vizi, nviz, dviz, raiof, valor_0, Sy, freq, alfa, beta, A,
    #                                            R_bar, CBA, QP, csi, dmax,nos_viz,dts,P,q)

end

#
# Atualiza rhos, se flag==true, faz só a inicialização
#
function Atualiza_Rho(rho, crit_old, rho_max, v_res, mu_res, flag=false)

    # Copia o critério antigo
    crit_rho = copy(crit_old)

    # Inicializa
    n_rho    = size(rho,1)
    crit_rho = zeros(n_rho)

    # Varre os primeiros rhos para atualizar o critério
    for j = 1:(n_rho-1)
        criterio    = max.(v_res[j], -mu_res[j]/rho[j])
        crit_rho[j] = max(0.0, criterio)
    end

    # O ultimo critério engloba os demais multiplicadores (eg tensão)
    criterio        = norm(max.(v_res[n_rho:end], -mu_res[n_rho:end]/rho[n_rho]))
    crit_rho[n_rho] = max.(0.0, criterio)

    # Se for a inicialização, para por aqui
    if flag == true
        return crit_rho
    end

    # Se continuou, agora verifica se o rho precisa ser atualizado
    for j = 1:n_rho

        if crit_rho[j] >= 0.9*crit_old[j]
            rho[j] = min.(1.1*rho[j], rho_max)
        end # if

    end

    # Atualiza também os multiplicadores das restrições
    for j=1:(n_rho-1)
        mu_res[j] = max.(rho[j]*v_res[j] + mu_res[j], 0.0)
    end #for j

    # Do mesmo modo, o últimos multiplicadores são diferentes...
    mu_res[n_rho:end] = max.(norm(rho[n_rho]*v_res[n_rho:end] + mu_res[n_rho:end]), 0.0)
    return rho, crit_rho, mu_res

end
