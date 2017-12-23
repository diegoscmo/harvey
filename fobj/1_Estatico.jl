################################################################################
#####                       Flebixilidade Estática                        ######
################################################################################

#
# Define o Função Objetivo de Flexibilidade Estática
#
function F_Est(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
   nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
         M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
            NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
          raiof::Float64, Y0::Array{Float64,1}, caso::Int64, freq::Float64, alfa::Float64,
         beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, filtra::Bool=true)

    # Filtra o x antes de qualquer coisa
    if filtra
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    else
        xf = copy(x)
    end

    # Remonta matriz de rigidez global, aqui é aplicado o SIMP
    KG, = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema novamente
    US = vec(cholfact(Symmetric(KG))\F)

    # Flexibilidade Estática
    Est = 0.5*dot(F,US)

    # Retorna valor para primeira iteração
    if tipo == 0
        return [Est],0.0,0.0,0.0
    end

    # Função objetivo normalizada
    valor_fun = Est / Y0[1]

    # Funções de restrição, volume normalizada
    valor_res = [ (mean(xf)-0.49)/0.51 ]

    # Se quiser a função obj normalizada
    if tipo == 1

        USx = Expande_Vetor(US, nnos, ID)
        TS = Tquad4_I(xf, nel, SP, 2.8, ijk, CBA, USx)

        return valor_fun, valor_res, [valor_fun], TS

    # Função Lagrangiana
    elseif tipo == 2
        L = valor_fun + 0.5*rho*max(0.0, valor_res[1] + mult_res[1]/rho)^2.0
        return L,0.0,0.0,0.0

    # Calculo da derivada
    elseif tipo == 3

        # Inicializa a derivada interna
        dL = Array{Float64}(uninitialized,nel)

        # Valor para correção da derivada
        corr_min = 1.0 - vmin

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/NX/NY/0.51
        dL2 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

        # Expande o vetor de deslocamentos  (insere zeros)
        USx = Expande_Vetor(US, nnos, ID)

        # Varre os elementos
        @inbounds for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Aloca vetor local e preenche
            USe = zeros(8)
            for k = 1:4
                nok        = nos_ele[k]
                USe[2*k-1] = USx[2*nok-1]
                USe[2*k]   = USx[2*nok]
            end

            # derivada das matrizes de rigidez e massa 109, corrigida
            dKedx  = corr_min*SP*xf[j]^(SP-1.0)*K0

            # Derivada da flexibilidade estática 57
            dSdx = -0.5*real(transpose(USe)*dKedx*USe)/Y0[1]

            # Derivada do LA
            dL[j] = dSdx + dL2

        end #j

        # Corrige aplicando a derivada do x filtrado em relação ao x original
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL,0.0,0.0,0.0

    end #tipo 3

end #fim da funcao
