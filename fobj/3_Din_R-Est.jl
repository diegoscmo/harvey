################################################################################
#####       Flebixilidade Dinâmica com restrição de tensão                ######
################################################################################

#
# Define o Função Objetivo de Flexibilidade Dinâmica com restrição de tensão
#
function F_Din_R_Est(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
               nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
               M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1},
               NX::Int64, NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
               dviz::Array{Float64,2}, raiof::Float64, Y0::Array{Float64,1},
               caso::Int64, freq::Float64, alfa::Float64, beta::Float64, A::Float64, Ye::Float64, filtra::Bool=true )

    # Filtra o x antes de qualquer coisa
    if filtra
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    else
        xf = copy(x)
    end

    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema estático
    US = vec(cholfact(Symmetric(KG))\F)

    # Resolve o sistema dinamico
    w  = 2.0*pi*freq
    KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
    UD = vec(lufact(KD)\F)

    # Flexibilidade Estática e Dinamica
    Din = abs(real(dot(F,UD)))
    Est = 0.5*dot(F,US)

    # Pega apenas os valores iniciais da restrição
    if tipo == 0
        Y0 = Array{Float64}(uninitialized,2)
        Y0[1] = Din
        Y0[2] = Est
        return Y0,0.0,0.0
    end

    # Função objetivo e restrição normalizadas
    valor_fun = Din/Y0[1]
    valor_est = Est/Y0[2]

    # Funções de restrição, volume normalizada
    valor_res = [ (mean(x)-0.49)/0.51
                valor_est-Ye ]

    # Se quiser só a F_Obj, retorna
    if tipo == 1
        return valor_fun, valor_res, [valor_fun;valor_est]

    # Função Lagrangiana
    elseif tipo == 2
        L = valor_fun + 0.5*rho*(max(0.0, valor_res[1] + mult_res[1]/rho)^2.0 +
                                 max(0.0, valor_res[2] + mult_res[2]/rho)^2.0 )
        return L,0.0,0.0

  # Calculo da derivada
    elseif tipo == 3

        # Inicializa a derivada interna
        dL = Array{Float64}(uninitialized,nel)

        # Valor para correção da derivada
        corr_min = 1.0 - vmin

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/NX/NY/0.51

        # Para flexibilidade dinamica, problema adjunto
        adj = dot(F,conj(UD)) / abs(dot(F,UD))

        # Parte da derivada referente a restrição de volume
        dL1 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

        # E a restrição da flexibilidade estática
        dL2 = max(0.0, mult_res[2] + rho*valor_res[2])

        # Expande o vetor de deslocamentos  (insere zeros)
        USx = Expande_Vetor(US, nnos, ID)
        UDx = Expande_Vetor(UD, nnos, ID, true)

        # Varre os elementos
        @inbounds for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Aloca vetor local e preenche
            USe = zeros(8)
            UDe = zeros(Complex,8)
            for k = 1:4
                nok        = nos_ele[k]
                USe[2*k-1] = USx[2*nok-1]
                USe[2*k]   = USx[2*nok]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
            end

            # derivada das matrizes de rigidez e massa 109, corrigida
            dKedx  = corr_min*SP*xf[j]^(SP-1.0)*K0

            # Correção Olhoff & Du, 91 + Correção do valor minimo
            dMedx = corr_min*M0
            if xf[j] < 0.1
                dMedx = (36E5*xf[j]^5.0 - 35E6*xf[j]^6.0)*dMedx
            end

            # Derivada da matriz dinamica
            dKDedx = dKedx*(1.0+im*w*beta) + dMedx*(-w^2.0+im*w*alfa)

            # Derivada da Flexibilidade dinamica
            dDdx = -real(adj*(transpose(UDe)*dKDedx*UDe))/Y0[1]

            # Derivada da Flexibilidade Estática
            dSdx = -0.5*transpose(USe)*dKedx*USe/Y0[2]

            # Derivada do LA
            dL[j] = dDdx + dL1 + dL2*dSdx

        end #j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL,0.0,0.0

    end #tipo 3

end #fim da funcao
