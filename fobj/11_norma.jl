################################################################################
#####          Potência + Flebibilidade Estática, restrição no R          ######
################################################################################

#
# Define o Função Objetivo de Potência + Flexibilidade Estática, restrição no R
#
function F_Obj(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1}, tipo::Int64, nnos::Int64, nel::Int64,
               ijk::Array{Int64,2}, coord::Array{Float64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
               M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
               NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64,
               Y0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64, beta::Float64, A::Float64,
               R_bar::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64, nos_viz,
                dts::String,P::Float64,q::Float64)

    # Número de rhos para penalização
    n_rho = 3

    # Peso do problema de potência e do estático
    B   = 1.0 - abs(A)

    # Filtra o x antes de qualquer coisa
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)

    # Corrige x com vmin
    xc = @. vmin + (1.0 - vmin)*xf

    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(xc, nel, ijk, ID, K0, M0, SP, vmin)

    # Resolve o sistema dinamico
    w  = 2.0*pi*freq
    KD = KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG
    KDf = lufact(KD)
    UD = vec(KDf\F)

    # Resolve o sistema estático
    KSf = cholfact(Symmetric(KG))
    US = vec(KSf\F)

    # Gera um vetor com as "normas nodais de densidade"
    aj = zeros(nnos)
    for j=1:nnos

        # Recupera os elementos vizinhos deste nó
        viz = nos_viz[j]
        vaj = 0.0

        # Acumula os vizinhos
        for v in viz
            vaj += xc[v]^q
        end

        # Tira a norma q
        aj[j] = vaj^(1.0/q)
     end

     # Para montar o a_N
     soma = 0.0
     for j=1:nnos

       # Recupera aj para este nó
       vaj = aj[j]

       # Recupera os gls do nós
       gl1 = ID[j,1]
       gl2 = ID[j,2]

       # Primeiro gl
       if gl1 != 0
           ui = UD[gl1]
           soma += real(conj(ui)*vaj*ui)^(P/2.0)
       end
       # Segundo gl
       if gl2 != 0
           ui = UD[gl2]
           soma += real(conj(ui)*vaj*ui)^(P/2.0)
       end

     end

  # Aplica a norma
  Nor1 = soma^(1.0/P)

  # Corrigida por log
  Nor = 100.0 + 10.0*log10(Nor1)

   # Flexibilidade estática
    FS  = 0.5*dot(F,US)

    # Retorna valores zero
    if tipo == 0
        return [Nor,FS],0.0,0.0,0.0
    end

    # Função objetivo original
    valor_fun = (A*Nor/Y0[1] + B*FS/Y0[2])

    # Restrição do R
    Ec = w^2.0*real(adjoint(UD)*MG*UD)
    Ep = real(adjoint(UD)*KG*UD)

    R  = Ec/Ep

    # Funções de restrição, volume, R, stress
    valor_res = Array{Float64}(undef,3)
    valor_res[1] = mean(xf)/dmax - 1.0
    valor_res[2] = 0.0
    valor_res[3] = 0.0

    # Calcula tensões
    sigma = Array{Complex{Float64}}(undef,nel,12)

    # Calcula o Lagrangiano aumentado
    L = valor_fun + 0.5*rho[1]*max(0.0, valor_res[1] + mult_res[1]/rho[1])^2.0

    # Se quiser só a F_Obj, retorna
    if tipo == 1

        # Calcula a tensão para o display
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
                # Calcula o von-mises (forçando real pra não ficar resto)
                VM[j,k] = sqrt(real(adjoint(sigmax)*M*sigmax))
            end #k
        end #j
        maxS = maximum(VM)/Sy

        # Faz a análise modal e passa os valores para o display externo
        freqs = Analise_Modal(10,K0,M0,nel,nnos,ijk,ID,coord,vizi, nviz, dviz, raiof,x,SP,vmin,dts)
        Analise_Harmonica(freq,K0,M0,nel,nnos,ijk,ID,coord,vizi,nviz,dviz,raiof,alfa,beta,F,x,SP,vmin,dts,1)

        return valor_fun, valor_res, [L;Nor1;FS;mean(xf);R;maxS;freqs], sigma

        # Função Lagrangiana
    elseif tipo == 2

        return L,0.0,0.0,0.0

    # Calculo da derivada
    elseif tipo == 3

        # Inicializa a derivada interna
        dL = Array{Float64}(undef,nel)

        # Valor para correção da derivada
        crmin = 1.0 - vmin^SP

        # Derivada da restrição de volume (igual p/todos elementos)
        dVdx = max(0.0, rho[1]*valor_res[1] + mult_res[1])*(1.0/(NX*NY))*(1.0/dmax)

        ### Resolvendo os adjuntos, começando pela norma ###

        # Calcula a_N com o somatório
        a_N = soma^((1.0/P)-1.0)

        # Inicializa o auxiliar para resolver o adj reduzido e varre nós
        aux_N = zeros(Complex,size(UD,1))
        for j=1:nnos

            # Recupera aj para este nó
            vaj = aj[j]

            # Recupera os gls do nós
            gl1 = ID[j,1]
            gl2 = ID[j,2]

            # Primeiro gl
            if gl1 != 0
                ui = UD[gl1]
                bj  = real(conj(ui)*vaj*ui)^((P/2.0)-1.0)
                aux_N[gl1] = -a_N*bj*vaj*conj(ui)
            end

            # Segundo gl
            if gl2 != 0
                ui = UD[gl2]
                bj  = real(conj(ui)*vaj*ui)^((P/2.0)-1.0)
                aux_N[gl2] = -a_N*bj*vaj*conj(ui)
            end
        end

        # Resolve o adjunto da norma
        adj_N = vec(KDf\aux_N)

        # Expande os vetores de deslocamentos  (insere zeros). UDx já foi
        adj_Nx = Expande_Vetor(adj_N, nnos, ID)
        USx    = Expande_Vetor(US, nnos, ID)
        UDx    = Expande_Vetor(UD, nnos, ID)

        # Aloca vetores locais
        UDe    = Array{ComplexF64}(undef,8)
        adj_Ne = Array{ComplexF64}(undef,8)
        USe    = Array{Float64}(undef,8)

        # Varre os elementos
        for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Preenche vetores locais
            for k = 1:4
                nok        = nos_ele[k]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
                USe[2*k-1] = USx[2*nok-1]
                USe[2*k]   = USx[2*nok]
                adj_Ne[2*k-1] = adj_Nx[2*nok-1]
                adj_Ne[2*k]   = adj_Nx[2*nok]
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

            # Derivada da norma!
            dNdx0 = 0.0

            # Para cada nó do elemento j
            for k in nos_ele

                # Recupera os elementos vizinhos deste nó
                viz = nos_viz[k]

                # Acumula para calcular dajdx...
                vaj = 0.0
                for v in viz
                    vaj += xc[v]^(q)
                end

                # Faz a parte da derivada de aij
                dajdx = (vaj^(1.0/q - 1.0))*(crmin*xc[j]^(q-1.0))

                # Recupera vaj para este nó (usando vaj de novo, mas não tem problema)
                vaj = aj[k]

                # Recupera os gls do nó
                gl1 = ID[k,1]
                gl2 = ID[k,2]

                # Primeiro gl
                if gl1 != 0
                    ui = UD[gl1]
                    bj  = real(conj(ui)*vaj*ui)^((P/2.0)-1.0)
                    dNdx0 += (a_N/2.0)*bj*(real(ui)^2.0 + imag(ui)^2.0)*dajdx
                end

                # Segundo gl
                if gl2 != 0
                    ui = UD[gl2]
                    bj  = real(conj(ui)*vaj*ui)^((P/2.0)-1.0)
                    dNdx0 += (a_N/2.0)*bj*(real(ui)^2.0 + imag(ui)^2.0)*dajdx
                end

            end #k

            # Derivada da Norma
            dNdx1 = real(transpose(adj_Ne)*dKDedx*UDe)
            dNdx = (10.0/(log(10.0)*Nor1))*(dNdx0+dNdx1)/Y0[1]

            # Derivada da flexibilidade estática e correção
            dYdx = -0.5*transpose(USe)*dKedx*USe/Y0[2]

            # Derivada do LA
            dL[j] = A*dNdx + B*dYdx + dVdx

        end #for j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL,L,0.0,0.0

    end #tipo 3

end #fim da funcao
