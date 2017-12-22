################################################################################
#####                                Filtros                              ######
################################################################################
#
# Identifica os vizinhos de cada elemento e a distancia entre eles
#
function Proc_Vizinhos(nel::Int64, coord::Array{Float64,2}, conect::Array{Int64,2}, raiof::Float64)

    # Inicializa a matriz de vizinhos, distancias e quantidade de vizinhos
    vizi = zeros(Int64,nel,1000)
    dviz = zeros(nel,1000) # Vai dar erro se o filtro for muito grande
    nviz = zeros(Int64,nel)

    # Inicializa centro dos elementos em estudo
    cj = zeros(Float64,2)
    ck = zeros(Float64,2)

    for j=1:nel

        # Define o centro do elemento utilizando nós 1 e 3
        n1 = conect[j,1]
        n3 = conect[j,3]
        cj[1] = (coord[n1,1]+coord[n3,1])/2.0
        cj[2] = (coord[n1,2]+coord[n3,2])/2.0

        for k=1:nel
            n1 = conect[k,1]
            n3 = conect[k,3]
            ck[1] = (coord[n1,1]+coord[n3,1])/2.0
            ck[2] = (coord[n1,2]+coord[n3,2])/2.0

            # Calcula a distancia entre os elementos
            dkj = norm(ck-cj)

            # Se for menor do que o raio, eles são vizinhos
            if dkj <= raiof

                # Incrementa o número de vizinhos e armazena vizi e dviz
                nviz[j] += 1
                vizi[j,nviz[j]] = k
                dviz[j,nviz[j]] = dkj

            end
        end
    end

    # Verifica o elemento com mais vizinho e reduz matrizes
    big = maximum(nviz)
    if big == 0
        big = 1
    end
    vizi = vizi[:,1:big]
    dviz = dviz[:,1:big]

    return vizi,nviz,dviz
end

#
# Aplica filtro de densidades, pág 65
#
function Filtro_Dens(xi::Array{Float64,1},nel::Int64,vizi::Array{Int64,2},
    nviz::Array{Int64,1},dviz::Array{Float64,2},raiof::Float64)

    # Inicializa vetor de saída
    xo = zeros(Float64,nel)

    # Varre os elementos
    @inbounds for ele = 1:nel

        # Zera acumuladores
        dividen = 0.0
        divisor = 0.0

        nvizk = nviz[ele]

        # Acumula somatórios para cada vizinho de k
        @inbounds for viz = 1:nvizk

            dist = dviz[ele,viz]
            #pesoH    = raiof - dist
            pesoH    = 1.0 - dist/raiof

            vz = vizi[ele,viz]

            dividen += pesoH*xi[vz]
            divisor += pesoH
        end #viz

        # Filtra elemento
        xo[ele] = dividen/divisor

    end #ele

    return xo
end
#
# Corrige derivada devido ao filtro de densidades (cadeia)
#
function dL_Dens(dL::Array{Float64,1},nel::Int64,vizi::Array{Int64,2},
    nviz::Array{Int64,1},dviz::Array{Float64,2},raiof::Float64)

    # Inicializa o vetor das derivadas corrigidas
    dLo = zeros(Float64,nel)

    # Varre m elementos
    for m=1:nel

        # Define vizinhos de m
        vizm = vizi[m,1:nviz[m]]

        # Varre os vizinhos de m
        @inbounds for k in vizm

            # Posição de m na tabela de vizinhanca
            pos = findfirst(equalto(m),vizi[k,1:nviz[k]])

            # Calcula o peso (H) entre m e k
            Hmk = 1.0 - dviz[k,pos]/raiof

            # Loop para denominador da filtragem, Hjk
            soma = 0.0
            @inbounds for j=1:nviz[k]
                soma += 1.0 - dviz[k,j]/raiof
            end #j

            # Acumula a correção de derivada
            dLo[m] += dL[k]*(Hmk/soma)

        end #k

    end #m

    return dLo
end



function Filtra(x,csi,nel)

    # Inicializa o x filtrado

    # Calcula filtro de densidades

    #se csi = 0.0, retorna filtro de densidades

    #se csi > 0.0, heaviside

    # Soma proveniente do filtro de densidades
    sum_xd = sum(xd)

    #resolve o filtro para csi
    a = 0.0
    b = 1.0
    c = 0.5

    tol = 1E-2

    x_a = Array{Float64}(uninitialized,nel)
    x_c = Array{Float64}(uninitialized,nel)
    xh  = Array{Float64}(uninitialized,nel)

    while true

        if (b-a)<=tol
            break
        end

        thcn_a = tanh(csi*a)
        thcn_c = tanh(csi*c)

        for j = 1:nel
            x_a[j] = (thcn_a + tanh(csi*(x[j]-a)))/(thcn_a+tanh(csi*(1.0-a)))
            x_c[j] = (thcn_c + tanh(csi*(x[j]-c)))/(thcn_c+tanh(csi*(1.0-c)))
        end

        fun_a = sum(x_a) - sum_xd
        fun_c = sum(x_c) - sum_xd

        if fun_a*fun_c < 0.0
            b = copy(c)
        else
            a = copy(c)
        end

        c = (a+b)/2.0

    end

    eta = c
    thcn = tanh(csi*eta)
    for j = 1:nel
        xh[j] = (thcn + tanh(csi*(x[j]-eta)))/(thcn+tanh(csi*(1.0-eta)))
    end

    return xh,eta
end
