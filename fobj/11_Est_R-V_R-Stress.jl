################################################################################
#####                       Flebixilidade Estática                        ######
################################################################################

#
# Define o Função Objetivo de Flexibilidade Estática
#
function F_Est_S(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
     nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
           M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
              NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
            raiof::Float64, Y0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64,
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64)

    # Filtra o x antes de qualquer coisa
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

    # Corrige x com vmin
    xc = broadcast(+,vmin,(1.0-vmin)*xf)

    # Remonta matriz de rigidez global, aqui é aplicado o SIMP
    KG, = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema novamente
    KGf = cholfact(Symmetric(KG))
    US = vec(KGf\F)

    # Flexibilidade Estática
    Est = 0.5*dot(F,US)

    # Retorna valor para primeira iteração
    if tipo == 0
        return [Est],0.0,0.0,0.0
    end

    # Função objetivo normalizada #FIXME
    valor_fun = Est #/ Y0[1]

    # Calcula as tensões sigma
    USx = Expande_Vetor(US, nnos, ID)
    sigma = Tquad4_I(xc, nel, SP, QP, ijk, CBA, USx)

    # Funções de restrição, volume normalizada
    valor_res = Array{Float64}(undef,1+nel*4)
    valor_res[1] = (mean(xf)-0.49)/0.51

    # Tensão limite (escoamento)
    Sy = 9E4

    # Tensor de von-Mises
    M = [1.0    -0.5    0.0
         -0.5   1.0     0.0
         0.0    0.0     3.0]

    # Restrições de tensão e von-Mises
    VM = zeros(nel,4)
    for j=1:nel
        for k=1:4
            # Localiza o sigma de j,k
            sigmax = sigma[j,(k*3-2):(k*3)]
            VM[j,k] = sqrt(transpose(sigmax)*M*sigmax)

            # Adiciona restrição no local correto [2:end]
            loc = j*4+k-3
            valor_res[loc] = VM[j,k]/Sy - 1.0
        end
    end

    # Se quiser a função obj normalizada
    if tipo == 1

        return valor_fun, valor_res, [valor_fun], sigma

    # Função Lagrangiana
    elseif tipo == 2

        L = valor_fun
        for j=1:size(valor_res,1)
            L += 0.5*rho*max(0.0, valor_res[j] + mult_res[j]/rho)^2.0
        end
        return L,0.0,0.0,0.0

    # Calculo da derivada
    elseif tipo == 3

        # Inicializa a derivada interna
        dL = Array{Float64}(undef,nel)

        # Valor para correção da derivada
        corr_min = 1.0 - vmin

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/NX/NY/0.51
        dL2 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

        # Resolve o adjunto da tensão
        ngl = size(US,1)
        aux_Sg = zeros(ngl,1)
        for j=1:nel

            # Determina os gls do elemento (vai ter que montar na reduzida)
            (glg,gll) = gl_livres_elemento(j,ijk,ID)

            for k=1:4
                # localiza j,k no vetor global de restrições
                loc = j*4+k-3
                # Calcula o vetor
                aux_S = rho*max(0.0, mult_res[loc]/rho + valor_res[loc])*(1.0/Sy)*
                        (transpose(sigma[j,k*3-2:k*3])*M / VM[j,k])*xc[j]^(SP-QP)*CBA[:,:,k]

                for m = 1:length(glg)
                    glgm =  glg[m]
                    gllm =  gll[m]
                    aux_Sg[glgm] += aux_S[gllm]
                end #m

            end #for k
        end #for j

        # Resolve o adjunto da tensão e expande
        adj_S = vec(KGf\(-0.5*F - aux_Sg))

        adj_Sx = Expande_Vetor(adj_S, nnos, ID)

        # Varre os elementos
        @inbounds for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Aloca vetor local e preenche
            USe = zeros(8)
            adj_Se = zeros(8)

            for k = 1:4
                nok        = nos_ele[k]
                USe[2*k-1] = USx[2*nok-1]
                USe[2*k]   = USx[2*nok]
                adj_Se[2*k-1] = adj_Sx[2*nok-1]
                adj_Se[2*k]   = adj_Sx[2*nok]
            end

            # Resolve o termo da tensão para o elemento
            dTdx = 0.0
            for k = 1:4
                # localiza j,k no vetor global de restrições
                loc = j*4+k-3
                dTdx = dTdx + rho*max(0.0,mult_res[loc]/rho + valor_res[loc])*(1.0/Sy)*
                        (transpose(sigma[j,k*3-2:k*3])*M / VM[j,k])*(SP-QP)*xc[j]^(SP-QP-1.0)*CBA[:,:,k]*USe
            end

            # derivada das matrizes de rigidez e massa 109, corrigida
            dKedx  = corr_min*SP*xf[j]^(SP-1.0)*K0

            # Derivada do LA
            dL[j] = dL2 + dTdx + transpose(adj_Se)*dKedx*USe

        end #j

        # Corrige aplicando a derivada do x filtrado em relação ao x original
        if filtra
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)
        end

        return dL,0.0,0.0,0.0

    end #tipo 3

end #fim da funcao
