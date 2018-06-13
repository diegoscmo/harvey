################################################################################
#####          Potência + Flebibilidade Estática, restrição no R          ######
################################################################################

#
# Define o Função Objetivo de Potência + Flexibilidade Estática, restrição no R
#
function F_Obj4(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1}, tipo::Int64, nnos::Int64, nel::Int64,
               ijk::Array{Int64,2}, coord::Array{Float64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
               M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
               NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64,
               Y0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64, beta::Float64, A::Float64,
               R_bar::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64, nos_viz,
               dts::String, P::Float64, q::Float64)

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
    KDf = lu(KD)
    UD = vec(KDf\F)

    # Resolve o sistema estático
    KSf = cholesky(Symmetric(KG))
    US = vec(KSf\F)

    # tolerancia para verificar se está na extremidade,
    # isso teria que vir de fora FIXME
    tol = 1.0/40/10.0

    # Gera um vetor com as "normas nodais de densidade"
    aj = zeros(nnos)
    for j=1:nnos

        # Verifica se está na extremidade... aí calcula aj FIXME
        xx = coord[j,1]
        yy = coord[j,2]
        if (xx>=1.0-tol) && (xx<=1.0+tol) && (yy>=0.2-tol) && (yy<=0.3+tol)

            # Recupera os elementos vizinhos deste nó
            viz = nos_viz[j]
            vaj = 0.0

            # Acumula os vizinhos
            for v in viz
                vaj += xc[v]^q
            end

            # Tira a norma q
            aj[j] = vaj^(1.0/q)

        # Se não estiver na tolerancia, zera o aj do nó FIXME
        else
            aj[j] = 0.0
        end

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
#       if gl1 != 0
#           ui = UD[gl1]
#           soma += real(conj(ui)*vaj*ui)^(P/2.0)
#       end
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
        return [Nor,FS],(2+nel*4),0.0,0.0
    end

    # Função objetivo original
    valor_fun = (A*Nor/Y0[1] + B*FS/Y0[2])

    # Restrição do R
    Ec = w^2.0*real(adjoint(UD)*MG*UD)
    Ep = real(adjoint(UD)*KG*UD)

    R  = Ec/Ep

    # Funções de restrição, volume, R, stress
    valor_res = Array{Float64}(undef,2+nel*4)
    valor_res[1] = mean(xf)/dmax - 1.0
    valor_res[2] = 0.0

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

            # Adiciona restrição no local correto [3:end]
            loc = j*4+k-4
            valor_res[loc+2] =  VM[j,k]/Sy - 1.0

        end #k
    end #j

    # Calcula o Lagrangiano aumentado
    L = valor_fun + 0.5*rho[1]*max(0.0, valor_res[1] + mult_res[1]/rho[1])^2.0

    # Parte da tensão
    for j = 3:size(valor_res,1)
        L += 0.5*rho[3]*max(0.0, valor_res[j] + mult_res[j]/rho[3])^2.0
    end

    # Se quiser só a F_Obj, retorna
    if tipo == 1

        maxS = maximum(VM)/Sy

        # Faz a análise modal e passa os valores para o display externo
        freqs = Analise_Modal(10,K0,M0,nel,csi,nnos,ijk,ID,coord,vizi, nviz, dviz, raiof,x,SP,vmin,dts)
        harm = Harmonica(dts,freq,x,csi,nnos,nel,coord,ijk,ID,vizi,nviz,dviz,raiof,alfa,beta,F,K0,M0,vmin)

        return valor_fun, valor_res, [L;Nor1;FS;mean(xf);R;maxS;freqs], sigma, harm

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

        #<.> de dTdx
        dSL = Array{Float64}(undef,nel*4)
        for j=1:nel*4
            dSL[j] = max(0.0, rho[3]*valor_res[j+2] + mult_res[j+2])
        end

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
#            if gl1 != 0
#                ui = UD[gl1]
#                bj  = real(conj(ui)*vaj*ui)^((P/2.0)-1.0)
#                aux_N[gl1] = -a_N*bj*vaj*conj(ui)
#            end

            # Segundo gl
            if gl2 != 0
                ui = UD[gl2]
                bj  = real(conj(ui)*vaj*ui)^((P/2.0)-1.0)
                aux_N[gl2] = -a_N*bj*vaj*conj(ui)
            end
        end

        # Resolve o adjunto da norma
        adj_N = vec(KDf\aux_N)

        # Termos em comum no assembly do adjunto da tensão
        PQ2 = 2.0*(SP-QP)
        b2w21Sy = ((beta^2.0)*(w^2.0) + 1.0)/Sy
        Psiek = Array{Float64}(undef,8,8,4)
        for k=1:4
            Psiek[:,:,k] = transpose(CBA[:,:,k])*M*CBA[:,:,k]
        end

        # Somatório do adjunto da tensão
        aux_Sg = zeros(Complex,size(UD,1))
        for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Monta o vetor de deslocamentos locais
            UDe    = Array{ComplexF64}(undef,8)
            for k = 1:4
                nok        = nos_ele[k]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
            end

            # Determina os gls do elemento (vai ter que montar na reduzida)
            (glg,gll) = gl_livres_elemento(j,ijk,ID)

            # Acumula o adjunto de um elemento
            for k=1:4
                # Calcula o vetor auxiliar
                aux_SL = dSL[j*4+k-4]*b2w21Sy*((crmin*xc[j]^(PQ2))/(VM[j,k]))*(-adjoint(UDe)*Psiek[:,:,k])

                # Acumula o auxiliar sem os gdls presos
                for m = 1:length(glg)
                    aux_Sg[glg[m]] += aux_SL[gll[m]]
                end #m
            end #for k
        end #for j

        # Resolve o adjunto da tensão
        adj_S = vec(KDf\aux_Sg)

        # Expande os vetores de deslocamentos  (insere zeros). UDx já foi
        USx    = Expande_Vetor(US, nnos, ID)
        UDx    = Expande_Vetor(UD, nnos, ID)
        adj_Nx = Expande_Vetor(adj_N, nnos, ID)
        adj_Sx = Expande_Vetor(adj_S, nnos, ID)

        # Aloca vetores locais
        USe    = Array{Float64}(undef,8)
        UDe    = Array{ComplexF64}(undef,8)
        adj_Ne = Array{ComplexF64}(undef,8)
        adj_Se = Array{ComplexF64}(undef,8)

        # Varre os elementos
        for j=1:nel

            # Identifica nos do elemento
            nos_ele = ijk[j,:]

            # Preenche vetores locais
            for k = 1:4
                nok        = nos_ele[k]
                USe[2*k-1] = USx[2*nok-1]
                USe[2*k]   = USx[2*nok]
                UDe[2*k-1] = UDx[2*nok-1]
                UDe[2*k]   = UDx[2*nok]
                adj_Se[2*k-1] = adj_Sx[2*nok-1]
                adj_Se[2*k]   = adj_Sx[2*nok]
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
#                if gl1 != 0
#                    ui = UD[gl1]
#                    bj  = real(conj(ui)*vaj*ui)^((P/2.0)-1.0)
#                    dNdx0 += (a_N/2.0)*bj*(real(ui)^2.0 + imag(ui)^2.0)*dajdx
#                end

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

            # Derivada da tensão e adjunto
            dSdx = real(transpose(adj_Se)*dKDedx*UDe)
            for k = 1:4
                # Acumula dTdx
                dSdx += dSL[j*4+k-4]*((SP-QP)/Sy)*(VM[j,k]/(crmin*xc[j]))
            end #for k

            # Derivada da flexibilidade estática e correção
            dYdx = -0.5*transpose(USe)*dKedx*USe/Y0[2]

            # Derivada do LA
            dL[j] = A*dNdx + B*dYdx + dVdx + dSdx

        end #for j

        # Corrige aplicando a derivada do x filtrado em relação ao x original 65
        dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

        return dL,L,0.0,0.0

    end #tipo 3

end #fim da funcao
