# Define o Função Objetivo
function F_Obj(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
               nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
               M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1},
               NX::Int64, NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
               dviz::Array{Float64,2}, raiof::Float64, Y0::Array{Float64,1})

    # Filtra o x antes de qualquer coisa
    x = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

    # Remonta matriz de rigidez global, aqui é aplicado o SIMP
    KG, = Global_KM(x, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema novamente
    US = vec(cholfact(Symmetric(KG))\F)

    # Retorna valor para primeira iteração
    if tipo == 0
        Y0 = [0.5*dot(F,US)]
        return Y0
    end

    # Função objetivo, flexibilidade estática
    valor_fun = 0.5*dot(F,US) / Y0[1]

    # Funções de restrição, volume normalizada
    valor_res = [ (mean(x)-0.49)/0.51 ]

    # Se quiser a função obj normalizada
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

        dL2 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

        # Varre os elementos
        @inbounds for j=1:nel

            # Localiza deslocamentos e lambdas
            nos = ijk[j,:]

            # Zera Ue, complexo para cálculos Harmônicos
            USe = zeros(Float64,8)

            # Busca o U dos nós não restritos
            for k=1:4
                for q=1:2
                    loc = ID[nos[k],q]
                    if loc != 0
                        USe[2*k-2+q] = US[loc]
                    end #loc
                end #q
            end #k

            # derivada das matrizes de rigidez e massa 109, corrigida
            dKedx  = corr_min*SP*x[j]^(SP-1.0)*K0

            # Derivada da flexibilidade estática 57
            dSdx = real(-USe'*dKedx*USe) / Y0[1]

            # Derivada do LA
            dL[j] = dSdx + dL2

        end #j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL

    end #tipo 3

end #fim da funcao