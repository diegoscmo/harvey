################################################################################
#####          Potência + Flebibilidade Estática, restrição no R          ######
################################################################################

#
# Define o Função Objetivo de Potência + Flexibilidade Estática, restrição no R
#
function F_Obj(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1}, tipo::Int64,
     nnos::Int64, nel::Int64, ijk::Array{Int64,2}, coord::Array{Float64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
           M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
              NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
            raiof::Float64, Y0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64,
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64)

    # Número de rhos para penalização
    n_rho = 1

    # Filtra o x antes de qualquer coisa
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

    # Corrige x com vmin
    xc = @. vmin+(1.0-vmin)*xf

    # Volume
    vol = mean(xf)

    # Retorna valores zero
    if tipo == 0
        return [vol],n_rho,0.0,0.0
    end


    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(xc, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema dinamico
    w  = 2.0*pi*freq
    KD = KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG
    KDf = lufact(KD)
    UD = vec(KDf\F)

    # Função objetivo original
    valor_fun = vol

    # Funções de restrição, volume, R, stress
    valor_res = Array{Float64}(undef,nel*4)

    # Calcula tensões
    UDx = Expande_Vetor(UD, nnos, ID)
    sigma = Tquad4_I_Din(xc, nel, SP, QP, ijk, CBA, UDx, beta, w)

    # Tensor de von-Mises
    M = [ 1.0   -0.5    0.0
         -0.5    1.0    0.0
          0.0    0.0    3.0]

    # Calcula von Mises e adiciona restrição
    VM = Array{Float64}(undef,nel,4)
    for j=1:nel
        for k=1:4
            # Localiza o sigma de j,k
            sigmax = sigma[j,(k*3-2):(k*3)]
            VM[j,k] = real(sqrt(adjoint(sigmax)*M*sigmax))

            # Adiciona restrição no local correto [3:end]
            loc = j*4+k-4
            valor_res[loc] = VM[j,k]/Sy - 1.0
        end #k
    end #j

    # Se quiser só a F_Obj, retorna
    if tipo == 1
        maxS = maximum(VM)/Sy

        nmodos = 6
        freqs = Analise_Modal(nmodos,K0,M0,nel,nnos,ijk,ID,coord,xc,SP)
        @show(mean(xc)-dmax)
        @show(maxS)

        return valor_fun, valor_res, [valor_fun;mean(xf);maxS;freqs], sigma

    # Função Lagrangiana
    elseif tipo == 2
        L = valor_fun
        for j=1:size(valor_res,1)
            L += 0.5*rho[n_rho]*max(0.0, valor_res[j] + mult_res[j]/rho[n_rho])^2.0
        end

        return L,0.0,0.0,0.0

    # Calculo da derivada
    elseif tipo == 3

        # Inicializa a derivada interna
        dL = Array{Float64}(undef,nel)

        # Valor para correção da derivada
        crmin = 1.0 - vmin

        # Derivada do vol
        dVdx = 1.0/NX/NY

        #<.> de dTdx
        dTL = Array{Float64}(undef,nel*4)
        for j=1:nel*4
            dTL[j] = max(0.0, rho[n_rho]*valor_res[j] + mult_res[j])
        end

        # Problemas adjuntos, começando pela expressão da tensão
        ngl = size(UD,1)
        aux_Sg = zeros(Complex,ngl)

        # Algumas coisas que vão lá no meio
        PQ2 = 2.0*(SP-QP)
        b2w21 = (beta^2.0)*(w^2.0) + 1.0
        YY = zeros(8,8,4)
        for k=1:4
            YY[:,:,k] = transpose(CBA[:,:,k])*CBA[:,:,k]
        end

        #  Loop do adjunto
        for j=1:nel

            # Determina os gls do elemento (vai ter que montar na reduzida)
            (glg,gll) = gl_livres_elemento(j,ijk,ID)

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Monta matriz de deslocamentos locais
            UDe    = Array{Complex{Float64}}(undef,8)
            for k = 1:4
                nok        = nos_ele[k]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
            end

            for k=1:4
                # localiza j,k no vetor global de restrições
                loc = j*4+k-4
                # Calcula e salva os vetores
                aux_SL = dTL[loc]*crmin*(xc[j]^(PQ2))*b2w21*(-adjoint(UDe)*YY[:,:,k])/(Sy*VM[j,k])

                # Monta o aux_Sg sem os gdls presos
                for m = 1:length(glg)
                    glgm = glg[m]
                    gllm = gll[m]
                    aux_Sg[glgm] += aux_SL[gllm]
                end #m
            end #for k
        end #for j

        adj_P = vec(KDf\aux_Sg)

        # Expande o vetor de deslocamentos  (insere zeros)
        adj_Px = Expande_Vetor(adj_P, nnos, ID)
        UDx    = Expande_Vetor(UD, nnos, ID)

        # Varre os elementos
        for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Aloca vetor local e preenche
            UDe    = Array{Complex{Float64}}(undef,8)
            adj_Pe = Array{Complex{Float64}}(undef,8)
            for k = 1:4
                nok        = nos_ele[k]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
                adj_Pe[2*k-1] = adj_Px[2*nok-1]
                adj_Pe[2*k]   = adj_Px[2*nok]
            end # for k

            # derivada das matrizes de rigidez e massa 109, corrigida
            dKedx  = crmin*SP*xc[j]^(SP-1.0)*K0

            # Correção Olhoff & Du, 91 + Correção do valor minimo
            dMedx = crmin*M0
            if xc[j] < 0.1
                dMedx = (36E5*xc[j]^5.0 - 35E6*xc[j]^6.0)*dMedx
            end #if

            # Derivada da matriz dinamica
            dKDedx = dKedx*(1.0+im*w*beta) + dMedx*(-w^2.0+im*w*alfa)

            # Resolve o termo da tensão para o elemento
            dTdx = 0.0
            for k = 1:4
                # localiza j,k no vetor global de restrições
                loc = j*4+k-4
                # Acumula dTdx
                dTdx += dTL[loc]*(1.0/Sy)*(SP-QP)*(VM[j,k]/(crmin*xc[j]))

            end #for k

            # O que ficou pro adjunto
            d_adj = real(transpose(adj_Pe)*dKDedx*UDe)

            # Derivada do LA
            dL[j] = dVdx+real(dTdx + d_adj)

        end #for j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL,0.0,0.0,0.0

    end #tipo 3

end #fim da funcao
