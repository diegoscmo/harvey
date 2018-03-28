################################################################################
#####                Potência + Flebibilidade Estática                    ######
################################################################################

#
# Define o Função Objetivo de Potência + Flexibilidade Estática
#
function F_Pot_A_Est(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
         nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
               M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
                  NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
                raiof::Float64, Y0::Array{Float64,1}, caso::Int64, freq::Float64, alfa::Float64,
               beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, filtra::Bool=true, QP::Float64=2.5)

    # Peso do problema de potência e do estático
    B   = 1.0 - A

    # Filtra o x antes de qualquer coisa
    if filtra
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    else
        xf = copy(x)
    end

    # Corrige x com vmin
    xc = broadcast(+,vmin,(1.0-vmin)*xf)

    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema dinamico
    w  = 2.0*pi*freq
    KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
    UD = vec(lufact(KD)\F)

    # Resolve o sistema estático
    US = vec(cholfact(Symmetric(KG))\F)

    # Potencia ativa e Pot Adaptada
    Pa   = 0.5*w*real(im*dot(F,UD))
    PadB = 100.0  + 10.0*log10(Pa)

    # Flexibilidade Estatica
    Est = 0.5*dot(F,US)

    # Retorna valores zero
    if tipo == 0
        Y0 = Array{Float64}(undef,2)
        Y0[1] = PadB
        Y0[2] = Est
       return Y0,0.0,0.0,0.0
    end

    # Função objetivo normalizada
    valor_fun1 = A*PadB/Y0[1]
    valor_fun2 = B*Est/Y0[2]
    valor_fun = valor_fun1 + valor_fun2

    # Funções de restrição, volume normalizada
    valor_res = [ (mean(xf)-0.49)/0.51 ]

    # Se quiser só a F_Obj, retorna
    if tipo == 1

        UDx = real(Expande_Vetor(UD, nnos, ID, true))
        TS = Tquad4_I(xf, nel, SP, QP, ijk, CBA, UDx)

        return valor_fun, valor_res, [valor_fun1;valor_fun2],TS

    # Função Lagrangiana
    elseif tipo == 2
        L = valor_fun + 0.5*rho*max(0.0, valor_res[1] + mult_res[1]/rho)^2.0
        return L,0.0,0.0,0.0

    # Calculo da derivada
    elseif tipo == 3

        # Inicializa a derivada interna
        dL = zeros(Float64,nel)

        # Valor para correção da derivada
        corr_min = 1.0 - vmin

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/NX/NY/0.51

        # Parte da derivada referente a restrição de volume
        dL1 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

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

            # Derivada da Potencia Ativa
            dPdx = -0.5*w*real(transpose(UDe)*dKDedx*UDe*im)
            dPBdx = (10.0/(log(10.0)*Pa))*dPdx/Y0[1]

            # Derivada da Flexibilidade Estática
            dSdx = -0.5*transpose(USe)*dKedx*USe/Y0[2]

            # Derivada do LA
            dL[j] = A*dPBdx + B*dSdx + dL1

        end #j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL,0.0,0.0,0.0

    end #tipo 3

end #fim da funcao
