# Define o Função Objetivo
function F_Obj(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
    nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
    M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1},
    NX::Int64, NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
    dviz::Array{Float64,2}, raiof::Float64, Y0::Array{Float64,1},
    caso::Int64, freq::Float64, alfa::Float64, beta::Float64, A::Float64, Ye::Float64 )

    ###################### ESTÁTICO ######################
    if caso == 1

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Remonta matriz de rigidez global, aqui é aplicado o SIMP
        KG, = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

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
            return valor_fun, valor_res, [valor_fun]

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
            dL = Array{Float64}(uninitialized,nel)

            # Valor para correção da derivada
            corr_min = 1.0 - vmin

            # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
            dVdx = 1.0/NX/NY/0.51

            dL2 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

            # Varre os elementos
            @inbounds for j=1:nel

                # Localiza deslocamentos e lambdas
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                USe = zeros(8)

                # Busca o U dos nós não restritos
                for k=[1,2,3,4]
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        USe[2*k-1] = US[gdli]
                    end
                    if gdlj != 0
                        USe[2*k] = US[gdlj]
                    end
                end #k

                # derivada das matrizes de rigidez e massa 109, corrigida
                dKedx  = corr_min*SP*xf[j]^(SP-1.0)*K0

                # Derivada da flexibilidade estática 57
                dSdx = -0.5*transpose(USe)*dKedx*USe / Y0[1]

                # Derivada do LA
                dL[j] = dSdx + dL2

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3

        ###################### DINAMICO ######################
    elseif caso == 2

        # Parâmetros do problema dinãmico
        w    = 2.0*pi*freq

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
        KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

        # Resolve o sistema dinamico
        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)

        # Retorna valor para primeira iteração
        if tipo == 0
            Y0 = [abs(real(dot(F,UD)))]
            return Y0
        end

        # Função objetivo, flexibilidade estática
        valor_fun = abs(real(dot(F,UD)))/Y0[1]

        # Funções de restrição, volume normalizada
        valor_res = [ (mean(x)-0.49)/0.51 ]

        # Se quiser só a F_Obj, retorna
        if tipo == 1
            return valor_fun, valor_res, [valor_fun]

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
            dL = Array{Float64}(uninitialized,nel)

            # Valor para correção da derivada
            corr_min = 1.0 - vmin

            # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
            dVdx = 1.0/NX/NY/0.51

            # Para flexibilidade dinamica, problema adjunto
            adj = dot(F,conj(UD)) / abs(dot(F,UD))

            # Parte da derivada referente a restrição de volume
            dL2 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

            # Varre os elementos
            @inbounds for j=1:nel

                # Localiza deslocamentos e lambdas
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                UDe = zeros(Complex,8)

                # Busca o U dos nós não restritos
                for k=1:4
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        UDe[2*k-1] = UD[gdli]
                    end
                    if gdlj != 0
                        UDe[2*k] = UD[gdlj]
                    end
                end #k

                # derivada das matrizes de rigidez e massa 109, corrigida
                dKedx  = corr_min*SP*xf[j]^(SP-1.0)*K0

                # Correção Olhoff & Du, 91 + Correção do valor minimo
                dMedx = corr_min*M0

                if xf[j] < 0.1
                    dMedx = (36E5*xf[j]^5.0 - 35E6*xf[j]^6.0)*dMedx
                end

                # Derivada da matriz dinamica
                dKDedx = dKedx*(1.0+im*w*beta) + dMedx*(-(w^2.0)+im*w*alfa)

                # Derivada do LA - flexibilidade dinamica 57
                dfdx = real( -adj*(transpose(UDe)*dKDedx*UDe) ) /Y0[1]

                # Derivada do LA
                dL[j] = dfdx + dL2

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3

        ###################### DINAMICO C/REST EST ######################
    elseif caso == 3

        # Parâmetros do problema dinãmico
        w    = 2.0*pi*freq

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
        KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

        # Resolve o sistema estático
        US = vec(cholfact(Symmetric(KG))\F)

    #    # Resolve o sistema dinamico
        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)

        # Pega apenas os valores iniciais da restrição
        if tipo == 0
            Y0 = Array{Float64}(uninitialized,2)
            Y0[1] = abs(real(dot(F,UD)))
            Y0[2] = 0.5*(abs(dot(F,US)))
            return Y0
        end

        # Função objetivo, flexibilidade dinamica
        valor_fun = abs(real(dot(F,UD))) / Y0[1]

        valor_est = 0.5*(abs(dot(F,US)))/Y0[2]

        # Funções de restrição, volume normalizada
        valor_res = [ (mean(x)-0.49)/0.51
                    valor_est-Ye ]

        # Se quiser só a F_Obj, retorna
        if tipo == 1
            return valor_fun, valor_res, [valor_fun;valor_est]

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
            dL = Array{Float64}(uninitialized,nel)

            # Valor para correção da derivada
            corr_min = 1.0 - vmin

            # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
            dVdx = 1.0/NX/NY/0.51

            # Para flexibilidade dinamica, problema adjunto
            adj = dot(F,conj(UD)) / abs(dot(F,UD))

            # Parte da derivada referente a restrição de volume
            dL1 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

            dL2 = max(0.0, mult_res[2] + rho*valor_res[2])

            # Varre os elementos
            @inbounds for j=1:nel

                # Localiza deslocamentos e lambdas
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                USe = zeros(8)
                UDe = zeros(Complex,8)

                # Busca o U dos nós não restritos
                for k=1:4
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        USe[2*k-1] = US[gdli]
                        UDe[2*k-1] = UD[gdli]
                    end
                    if gdlj != 0
                        USe[2*k] = US[gdlj]
                        UDe[2*k] = UD[gdlj]
                    end
                end #k

                # derivada das matrizes de rigidez e massa 109, corrigida
                dKedx  = corr_min*SP*xf[j]^(SP-1.0)*K0

                # Correção Olhoff & Du, 91 + Correção do valor minimo
                dMedx = corr_min*M0

                if xf[j] < 0.1
                    dMedx = (36E5*xf[j]^5.0 - 35E6*xf[j]^6.0)*dMedx
                end

                # Derivada da matriz dinamica
                dKDedx = dKedx*(1.0+im*w*beta) + dMedx*(-w^2.0+im*w*alfa)

                # Derivada do LA - flexibilidade dinamica 57
                dDdx = -real(adj*(transpose(UDe)*dKDedx*UDe))/Y0[1]

                dSdx = -0.5*transpose(USe)*dKedx*USe/Y0[2]

                # Derivada do LA
                dL[j] = dDdx + dL1 + dL2*dSdx

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3

        ###################### DINAMICO + ESTATICO ######################
    elseif caso == 4

        # Parâmetros do problema dinãmico
        w    = 2.0*pi*freq

        # Peso do problema dinâmico e estático
        B   = 1.0 - A

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
        KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

        # Resolve o sistema estático
        US = vec(cholfact(Symmetric(KG))\F)

        # Resolve o sistema dinamico
        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)

        # Retorna valores zero
        if tipo == 0
            Y0 = Array{Float64}(uninitialized,2)
            Y0[1] = abs(dot(F,UD))
            Y0[2] = 0.5*abs(dot(F,US))
           return Y0
        end

        # Função objetivo, flexibilidade estática
        valor_fun1 = A*real(abs(dot(F,UD))/Y0[1])

        valor_fun2 = B*0.5*abs(dot(F,US))/Y0[2]

        valor_fun = valor_fun1 + valor_fun2

        # Funções de restrição, volume normalizada
        valor_res = [ (mean(x)-0.49)/0.51 ]

        # Se quiser só a F_Obj, retorna
        if tipo == 1
            return valor_fun, valor_res, [valor_fun1;valor_fun2]

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
            dL = Array{Float64}(uninitialized,nel)

            # Valor para correção da derivada
            corr_min = 1.0 - vmin

            # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
            dVdx = 1.0/NX/NY/0.51

            # Para flexibilidade dinamica, problema adjunto
            adj = dot(F,conj(UD)) / abs(dot(F,UD))

            # Parte da derivada referente a restrição de volume
            dL1 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

            # Varre os elementos
            @inbounds for j=1:nel

                # Localiza deslocamentos e lambdas
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                USe = zeros(8)
                UDe = zeros(Complex,8)

                # Busca o U dos nós não restritos
                for k=1:4
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        USe[2*k-1] = US[gdli]
                        UDe[2*k-1] = UD[gdli]
                    end
                    if gdlj != 0
                        USe[2*k] = US[gdlj]
                        UDe[2*k] = UD[gdlj]
                    end
                end #k

                # derivada das matrizes de rigidez e massa 109, corrigida
                dKedx  = corr_min*SP*xf[j]^(SP-1.0)*K0

                # Correção Olhoff & Du, 91 + Correção do valor minimo
                dMedx = corr_min*M0

                if xf[j] < 0.1
                    dMedx = (36E5*xf[j]^5.0 - 35E6*xf[j]^6.0)*dMedx
                end

                # Derivada da matriz dinamica
                dKDedx = dKedx*(1.0+im*w*beta) + dMedx*(-w^2.0+im*w*alfa)

                # Derivada do LA - flexibilidade dinamica 57
                dDdx = real( -adj*(transpose(UDe)*dKDedx*UDe) ) /Y0[1]

                dSdx = 0.5*real( -transpose(USe)*dKedx*USe ) /Y0[2]

                # Derivada do LA
                dL[j] = A*dDdx + B*dSdx + dL1

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3

        ###################### POTENCIA ######################
    elseif caso == 5

        # Parâmetros do problema dinãmico
        w    = 2.0*pi*freq

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
        KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

        # Resolve o sistema dinamico
        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)

        Pa   = 0.5*w*real(im*dot(F,UD))

        PadB = 100.0  + 10.0*log10(Pa)

        # Retorna valores zero
        if tipo == 0
           return [PadB]
        end

        # Função objetivo, flexibilidade estática
        valor_fun = PadB/Y0[1]

        # Funções de restrição, volume normalizada
        valor_res = [ (mean(x)-0.49)/0.51 ]

        # Se quiser só a F_Obj, retorna
        if tipo == 1

            return valor_fun, valor_res, [valor_fun]

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
            dL = Array{Float64}(uninitialized,nel)

            # Valor para correção da derivada
            corr_min = 1.0 - vmin

            # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
            dVdx = 1.0/NX/NY/0.51

            # Parte da derivada referente a restrição de volume
            dL1 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

            # Varre os elementos
            @inbounds for j=1:nel

                # Localiza deslocamentos e lambdas
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                UDe = zeros(Complex,8)

                # Busca o U dos nós não restritos
                for k=1:4
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        UDe[2*k-1] = UD[gdli]
                    end
                    if gdlj != 0
                        UDe[2*k] = UD[gdlj]
                    end
                end #k

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

                # Derivada do LA
                dL[j] = dPBdx + dL1

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3
        ###################### POTENCIA C/REST ESTATICA ######################
    elseif caso == 6

        # Parâmetros do problema dinãmico
        w    = 2.0*pi*freq

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
        KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)

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

        valor_est = 0.5*(abs(dot(F,US)))/Y0[2]

        # Funções de restrição, volume normalizada
        valor_res = [ (mean(x)-0.49)/0.51
                    valor_est-Ye ]

        # Se quiser só a F_Obj, retorna
        if tipo == 1
            return valor_fun, valor_res, [valor_fun;valor_est]

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
            dL = Array{Float64}(uninitialized,nel)

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
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                USe = zeros(8)
                UDe = zeros(Complex,8)

                # Busca o U dos nós não restritos
                for k=1:4
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        UDe[2*k-1] = UD[gdli]
                        USe[2*k-1] = US[gdli]
                    end
                    if gdlj != 0
                        UDe[2*k] = UD[gdlj]
                        USe[2*k] = US[gdlj]
                    end
                end #k

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
                dPdx = -0.5*w*real(transpose(UDe)*dKDedx*UDe*im) / Y0[1]

                dSdx = 0.5*real( -transpose(USe)*dKedx*USe ) / Y0[2]

                # Derivada do LA
                dL[j] = dPdx + dL1 + dL2*dSdx

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3
                ###################### POTENCIA + ESTATICO ######################
    elseif caso == 7

        # Parâmetros do problema dinãmico
        w    = 2.0*pi*freq

        # Peso do problema de potência e do estático
        B   = 1.0 - A

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
        KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

        # Resolve o sistema dinamico
        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)

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
        valor_fun1 = A*0.5*w*real(im*dot(F,UD))/Y0[1]
        valor_fun2 = B*0.5*abs(dot(F,US))/Y0[2]

        valor_fun = valor_fun1 + valor_fun2

        # Funções de restrição, volume normalizada
        valor_res = [ (mean(x)-0.49)/0.51 ]

        # Se quiser só a F_Obj, retorna
        if tipo == 1
            return valor_fun, valor_res, [valor_fun1;valor_fun2]

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

            # Varre os elementos
            @inbounds for j=1:nel

                # Localiza deslocamentos e lambdas
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                USe = zeros(8)
                UDe = zeros(Complex,8)

                # Busca o U dos nós não restritos
                for k=1:4
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        UDe[2*k-1] = UD[gdli]
                        USe[2*k-1] = US[gdli]
                    end
                    if gdlj != 0
                        UDe[2*k] = UD[gdlj]
                        USe[2*k] = US[gdlj]
                    end
                end #k

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
                dPdx = -0.5*w*real(transpose(UDe)*dKDedx*UDe*im)/Y0[1]

                dSdx = 0.5*real( -transpose(USe)*dKedx*USe )/Y0[2]

                # Derivada do LA
                dL[j] = A*dPdx + B*dSdx + dL1

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3

        ## Pot A Est R R
    elseif caso == 8

        # Parâmetros do problema dinãmico
        w    = 2.0*pi*freq

        # Peso do problema de potência e do estático
        B   = 1.0 - abs(A)

        # Restrição do R
        R_m = Ye

        # Filtra o x antes de qualquer coisa
        xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

        # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
        KG, MG = Global_KM(xf, nel, ijk, ID, K0, M0, SP, vmin)

        KD = sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG)
        UD = vec(lufact(KD)\F)

        # Resolve o sistema estático
        US = vec(cholfact(Symmetric(KG))\F)

        # Retorna valores zero
        if tipo == 0
            Y0 = Array{Float64}(uninitialized,2)
            Y0[1] = 0.5*w*real(im*dot(F,UD))
            if A<0.0
                Y0[1] = -Y0[1]
            end
            Y0[2] = 0.5*(abs(dot(F,US)))
           return Y0
        end

        # Função objetivo, flexibilidade estática
        valor_fun1 = A*0.5*w*real(im*dot(F,UD))/Y0[1]
        valor_fun2 = B*0.5*abs(dot(F,US))/Y0[2]

        valor_fun = valor_fun1 + valor_fun2

        # Energia cinética e potencial
        Ec = 0.25*(w^2.0)*(adjoint(UD)*MG*UD)
        Ep = 0.25*(adjoint(UD)*KG*UD)

        # R e R corrigido
        R  = real(Ec/Ep)
        R_c = 1.0+1.0*log10(R)

        # Funções de restrição, volume normalizada
        valor_res = [ (mean(x)-0.49)/0.51
                       real(R_c - R_m)]

        # Se quiser só a F_Obj, retorna
        if tipo == 1
            return valor_fun, valor_res, [valor_fun1;valor_fun2;real(Ep/Ec)]

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

            # Problema adjunto de R
            aux_R = 0.5*((Ec/(Ep^2.0))*KG*conj(UD) - ((w^2.0)/Ep)*MG*conj(UD))

            adj_R = vec(lufact(KD)\aux_R)

            # Parte da derivada referente a restrição de volume
            dL1 = max(0.0, mult_res[1] + rho*valor_res[1])*dVdx

            # Parte da derivada referente a restrição do R (falta dRdx)
            dL2 = max(0.0, mult_res[2] + rho*valor_res[2])

            # Varre os elementos
            @inbounds for j=1:nel

                # Localiza deslocamentos e lambdas
                nos_ele = ijk[j,:]

                # Zera Ue, complexo para cálculos Harmônicos
                USe = zeros(8)
                UDe = zeros(Complex,8)
                adj_Re = zeros(Complex,8)

                # Busca o U dos nós não restritos
                for k=1:4
                    no_k = nos_ele[k]
                    gdli = ID[no_k,1]
                    gdlj = ID[no_k,2]
                    if gdli != 0
                        UDe[2*k-1] = UD[gdli]
                        USe[2*k-1] = US[gdli]
                        adj_Re[2*k-1] = adj_R[gdli]
                    end
                    if gdlj != 0
                        UDe[2*k] = UD[gdlj]
                        USe[2*k] = US[gdlj]
                        adj_Re[2*k] = adj_R[gdlj]
                    end
                end #k

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
                dPdx = -0.5*w*real(transpose(UDe)*dKDedx*UDe*im)/Y0[1]

                # Derivada da flexibilidade estática
                dSdx = 0.5*real( -transpose(USe)*dKedx*USe )/Y0[2]

                # Derivada do R
                dRdx_1 = +(0.25*(w^2.0)/Ep)*(adjoint(UDe)*dMedx*UDe)
                dRdx_2 = -(0.25*Ec/(Ep^2.0))*(adjoint(UDe)*dKedx*UDe)
                dRdx_3 = +real(transpose(adj_Re)*(dKDedx*UDe))

                # Derivada do R corrigida
                dRdx = (1.0/(log(10.0)*R))*real((dRdx_1 + dRdx_2 + dRdx_3))

                # Derivada do LA
                dL[j] = A*dPdx + B*dSdx + dL1 +  dL2*dRdx

            end #j

            # Corrige aplicando a derivada do x filtrado em relação ao x original 65
            dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

            return dL

        end #tipo 3

    end # casos

    # Se chegou até aqui é pq deu problema
    error("ATENCAO AO VALOR DA VARIAVEL CASO")

end #fim da funcao
