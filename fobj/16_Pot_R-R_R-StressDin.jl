################################################################################
#####          Potência + Flebibilidade Estática, restrição no R          ######
################################################################################

#
# Define o Função Objetivo de Potência + Flexibilidade Estática, restrição no R
#
function F_Obj(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
     nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
           M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
              NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
            raiof::Float64, Y0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64,
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64)

    # Peso do problema de potência e do estático
    B   = 1.0 - abs(A)

    # Restrição do R
    R_bar = Ye

    # Filtra o x antes de qualquer coisa
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

    # Corrige x com vmin
    xc = @. vmin+(1.0-vmin)*xf

    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(xc, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema dinamico
    w  = 2.0*pi*freq
    KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
    #KDf = ldltfact(Hermitian(KD))
    KDf = lufact(KD)
    UD = vec(KDf\F)

    # Resolve o sistema estático
    KSf = cholfact(Symmetric(KG))
    US = vec(KSf\F)

    # Potencia ativa e correção
    Pa = 0.5*w*real(im*dot(F,UD))
    #Pa  = 100.0  + 10.0*log10(Pa1)

    # Flexibilidade estática
    FS  = 0.5*dot(F,US)

    # Retorna valores zero
    if tipo == 0
        # Caso A seja negativo, ainda tem que dar 1.0
        if A<0.0
            Pa = -Pa
        end
        return [Pa,FS],0.0,0.0,0.0
    end

    # Função objetivo original
    valor_fun = (A*Pa/Y0[1] + B*FS/Y0[2])

    # Restrição do R
    R_g = real((w^2.0)*(adjoint(UD)*MG*UD) - R_bar*(adjoint(UD)*KG*UD))
    #R_c = 1.0+1.0*log10(R)

    # Funções de restrição, volume, R, stress
    valor_res = Array{Float64}(undef,2+nel*4)
    valor_res[1] = (mean(xc)-dmax)/(1.0-dmax)
    valor_res[2] = R_g

    # Calcula tensões
    UDx = Expande_Vetor(UD, nnos, ID)
    sigma = Tquad4_I_Din(xc, nel, SP, QP, ijk, CBA, UDx, beta, w)

    # Tensor de von-Mises
    M = [ 1.0   -0.5    0.0
         -0.5    1.0    0.0
          0.0    0.0    3.0]

    # Calcula von Mises e adiciona restrição
    VM = zeros(nel,4)
    for j=1:nel
        for k=1:4
            # Localiza o sigma de j,k
            sigmax = sigma[j,(k*3-2):(k*3)]
            VM[j,k] = real(sqrt(adjoint(sigmax)*M*sigmax))

            # Adiciona restrição no local correto [3:end]
            loc = j*4+k-4
            valor_res[loc+2] = VM[j,k]/Sy - 1.0
        end #k
    end #j

    # Se quiser só a F_Obj, retorna
    if tipo == 1
        maxS = maximum(VM)/Sy

        return valor_fun, valor_res, [valor_fun;mean(xf);maxS], sigma

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
        dL = zeros(Float64,nel)

        # Valor para correção da derivada
        corr_min = 1.0 - vmin

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/(NX*NY*(1.0-dmax))

        # <.> dVdx
        dLV = rho*max(0.0, valor_res[1] + mult_res[1]/rho)*dVdx

        # Problemas adjuntos, começando pela expressão da tensão
        ngl = size(UD,1)
        aux_Sg = zeros(Complex,ngl)
        aux_S1 = zeros(Complex,8)

        for j=1:nel

            # Determina os gls do elemento (vai ter que montar na reduzida)
            (glg,gll) = gl_livres_elemento(j,ijk,ID)

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            UDe    = zeros(Complex,8)
            for k = 1:4
                nok        = nos_ele[k]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
            end

            for k=1:4
                # localiza j,k no vetor global de restrições
                loc = j*4+k-4
                # Calcula e salva os vetores
                aux_S1 = rho*max(0.0, valor_res[loc+2]+mult_res[loc+2]/rho)*((1.0+(beta^2.0)*(w^2.0))/(Sy*VM[j,k]))*
                         corr_min*xc[j]^(2.0*(SP-QP))*(-adjoint(UDe)*transpose(CBA[:,:,k])*CBA[:,:,k])

                # Monta o aux_Sg sem os gdls presos
                for m = 1:length(glg)
                    glgm = glg[m]
                    gllm = gll[m]
                    aux_Sg[glgm] += aux_S1[gllm]
                end #m
            end #for k
        end #for j

        # Parte do djunto de R
        aux_R = -max(0.0, mult_res[2] + rho*valor_res[2])*2.0*adjoint(UD)*(w^2.0*MG - R_bar*KG)

        # Parte do adjunto de P
        dPdx = -0.5*A*w*F/Y0[1]
        #dPdx = A*(10.0/(log(10.0)*Pa1)*dPdx1)/Y0[1]

        # Agrupa o adjunto e resolve
        aux_adj = aux_Sg + dPdx*im - transpose(aux_R)
        adj_R = vec(KDf\aux_adj)

        # Expande o vetor de deslocamentos  (insere zeros)
        adj_Rx = Expande_Vetor(adj_R, nnos, ID)
        UDx    = Expande_Vetor(UD, nnos, ID)
        USx    = Expande_Vetor(US, nnos, ID)

        # Varre os elementos
        @inbounds for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Aloca vetor local e preenche
            USe    = zeros(8)
            UDe    = zeros(Complex,8)
            adj_Re = zeros(Complex,8)
            for k = 1:4
                nok        = nos_ele[k]
                USe[2*k-1] = USx[2*nok-1]
                USe[2*k]   = USx[2*nok]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
                adj_Re[2*k-1] = adj_Rx[2*nok-1]
                adj_Re[2*k]   = adj_Rx[2*nok]
            end # for k

            # derivada das matrizes de rigidez e massa 109, corrigida
            dKedx  = corr_min*SP*xc[j]^(SP-1.0)*K0

            # Correção Olhoff & Du, 91 + Correção do valor minimo
            dMedx = corr_min*M0
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
                dTdx += rho*max(0.0,valor_res[loc+2]+mult_res[loc+2]/rho)*
                                        (1.0/Sy)*(SP-QP)*(VM[j,k]/(corr_min*xc[j]))

            end #for k

            # Parte da reativa
            dR = max(0.0, mult_res[2] + rho*valor_res[2])*adjoint(UDe)*(w^2.0*dMedx - R_bar*dKedx)*UDe

            # Parte da estática (adjunto)
            dS = -0.5*B*transpose(USe)*dKedx*USe/Y0[2]

            # adjunto Real
            d_adj = real(transpose(adj_Re)*dKDedx*UDe)

            # Derivada do LA
            dL[j] = real(dLV + dTdx + dR + dS + d_adj)
        end #for j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL,0.0,0.0,0.0

    end #tipo 3

end #fim da funcao
