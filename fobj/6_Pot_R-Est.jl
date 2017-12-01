# Define o Função Objetivo
function F_Obj(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
               nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
               M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1},
               NX::Int64, NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
               dviz::Array{Float64,2}, raiof::Float64, Y0::Array{Float64,1})

    # Parâmetros do problema dinãmico
    freq = 600.0
    w    = 2.0*pi*freq

    alfa = 0.0
    #beta = 0.0
    beta = 1E-8
    #alfa = 0.8*w
    #beta = 0.1/w


    # Porcentagem da flexibilidade estatica
    Ye = 1.0

    # Filtra o x antes de qualquer coisa
    x = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(x, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema dinamico
#    UD = Array{Complex{Float64}}(uninitialized,nel)
#    try
#        KD = Hermitian(sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG))
#        UD = vec(cholfact(KD)\F)
#    catch
        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)
#    end

    # Resolve o sistema estático
    US = vec(cholfact(Symmetric(KG))\F)

    # Retorna valores zero
    if tipo == 0
        Y0 = Array{Float64}(uninitialized,2)
        Y0[1] = 0.5*w*real(im*dot(F,UD))
        Y0[2] = 0.5*(abs(dot(F,US)))
       return Y0
    end

    # Função objetivo, flexibilidade estática
    valor_fun = 0.5*w*real(im*dot(F,UD)) /Y0[1]

    # Funções de restrição, volume normalizada
    valor_res = [ (mean(x)-0.49)/0.51
                0.5*(abs(dot(F,US)))/Y0[2]-Ye ]

    # Se quiser só a F_Obj, retorna
    if tipo == 1
        return valor_fun, valor_res

        # Função Lagrangiana
    elseif tipo == 2

        L = 0.0
        for j=1:size(valor_res,1)
            L = L + max(0.0, mult_res[j] + valor_res[j]/rho)^2.0
        end

        L = valor_fun + rho*L/2.0

        return L

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

        # Parte da derivada referente a restrição estática (falta dSdx)
        dL2 = max(0.0, mult_res[2] + rho*valor_res[2])

        # Varre os elementos
        @inbounds for j=1:nel

            # Localiza deslocamentos e lambdas
            nos = ijk[j,:]

            # Zera Ue, complexo para cálculos Harmônicos
            UDe = complex(zeros(Float64,8))
            USe = zeros(Float64,8)

            # Busca o U dos nós não restritos
            for k=1:4
                for q=1:2
                    loc = ID[nos[k],q]
                    if loc != 0
                        UDe[2*k-2+q] = UD[loc]
                        USe[2*k-2+q] = US[loc]
                    end #loc
                end #q
            end #k

            # derivada das matrizes de rigidez e massa 109, corrigida
            dKedx  = corr_min*SP*x[j]^(SP-1.0)*K0

            # Correção Olhoff & Du, 91 + Correção do valor minimo
            if x[j] >= 0.1
                dMedx = corr_min*M0
            else
                dMedx = corr_min*(6.0*(6.0E5)*x[j]^5.0+7.0*(-5.0E6)*x[j]^6.0)*M0
            end

            # Derivada da matriz dinamica
            dKDedx = dKedx*(1.0+im*w*beta) + dMedx*(-w^2.0+im*w*alfa)

            # Derivada da Potencia Ativa
            dPdx = -0.5*w*real(transpose(UDe)*dKDedx*UDe*im) / Y0[1]

            dSdx = real( -transpose(USe)*dKedx*USe ) / Y0[2]

            # Derivada do LA
            dL[j] = dPdx + dL1 + dL2*dSdx

        end #j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL

    end #tipo 3

end #fim da funcao
