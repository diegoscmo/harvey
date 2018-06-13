################################################################################
#####                                Filtros                              ######
################################################################################

#
# Verifica se elementos estão na área de travamento e devolve
#
function Proc_Travas(nel::Int64, coord::Array{Float64,2}, conect::Array{Int64,2}, travas::Array{Float64,2})

    # Número de regiões de trava
    ntravas = size(travas,1)

    # Inicializa o número de travados e o vetor para captura
    ntravados = 0
    trava_els = zeros(Int64,nel)

    # Varre pelo número de travas
    for i=1:ntravas

        # Renomeia as travas em questão
        xi = travas[i,1]
        yi = travas[i,2]
        xf = travas[i,3]
        yf = travas[i,4]

        # Varre os elementos
        for j=1:nel

            # Define o centro do elemento utilizando nós 1 e 3
            n1 = conect[j,1]
            n3 = conect[j,3]
            xcj = (coord[n1,1]+coord[n3,1])/2.0
            ycj = (coord[n1,2]+coord[n3,2])/2.0

            # Verifica se está dentro do quadrado
            if (xcj >= xi) && (xcj <= xf) && (ycj >= yi) && (ycj <= yf)

                # Se estiver, adiciona no vetor e incrementa o ntravados
                ntravados += 1
                trava_els[ntravados] = j

            end #if
        end #for.j
    end #for.i

    #  Remove o excesso de trava_els
    trava_els = trava_els[1:ntravados]

    return trava_els
end #Proc_Travas

#
# Bloqueia os elementos em trava_els
#
function Trava_Els(x::Array{Float64,1},trava_els::Array{Int64,1})

    # Varre e trava
    for j=1:size(trava_els,1)
        el = trava_els[j]
        x[el] = 1.0
    end

    return x

end #Trava_Els


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

function Proc_Vizinhos_Nos(ijk, nnos, nelems)

   # Declara uma lista que irá conter as listas
   lista = []

   # Loop por todos os nós
   for no=1:nnos

       # Lista de vizinhos deste nó
       lista_no = Int64[]

       # Loop pelos elementos
       for ele=1:nelems

           # Se o no estiver nas conectividades, coloca na lista
           if no in ijk[ele,:]
               push!(lista_no,ele)
           end

       end #ele

       # Adiciona a lista de vizinhos deste nó
       push!(lista,lista_no)

   end #no

   return lista

   end

#
# Aplica filtro de densidades, pág 65
#
function Filtro_Dens(xi::Array{Float64,1},nel::Int64,csi::Float64,vizi::Array{Int64,2},
    nviz::Array{Int64,1},dviz::Array{Float64,2},raiof::Float64)

    # Inicializa vetores de saída (ñ pode ser undef)
    xf = zeros(Float64,nel)

    # Heaviside também
    xh = zeros(Float64,nel)

    # Varre os elementos
    for ele = 1:nel

        # Zera acumuladores
        dividen = 0.0
        divisor = 0.0

        nvizk = nviz[ele]

        # Acumula somatórios para cada vizinho de k
        for viz = 1:nvizk

            dist = dviz[ele,viz]
            pesoH    = 1.0 - dist/raiof

            vz = vizi[ele,viz]

            dividen += pesoH*xi[vz]
            divisor += pesoH
        end #viz

        # Filtra elemento
        xfe     = dividen/divisor
        xf[ele] = xfe

        # Aplica Heaviside
        xh[ele] = 1.0 - exp(-csi*xfe)  + xfe*exp(-csi)

    end #ele

    return xh,xf
end

#
# Corrige derivada devido ao filtro de densidades e Heaviside (cadeia)
#

function dL_Dens(dL::Array{Float64,1},xf::Array{Float64,1},csi::Float64,nel::Int64,vizi::Array{Int64,2},
    nviz::Array{Int64,1},dviz::Array{Float64,2},raiof::Float64)

    # Inicializa o vetor das derivadas corrigidas (ñ pode ser undef)
    dLo = zeros(Float64,nel)

    # Varre m elementos
    for m=1:nel

        # Define vizinhos de m
        vizm = vizi[m,1:nviz[m]]

        # Varre os vizinhos de m
        for k in vizm

            # Posição de m na tabela de vizinhanca
            pos = findfirst(==(m),vizi[k,1:nviz[k]])

            # Calcula o peso (H) entre m e k
            Hmk = 1.0 - dviz[k,pos]/raiof

            # Loop para denominador da filtragem, Hjk
            soma = 0.0
            for j=1:nviz[k]
                soma += 1.0 - dviz[k,j]/raiof
            end #j

            # Acha o x filtrado do vizinho
            xfik  = xf[k]

            # Calcula o termo referente ao heaviside
            heavi = csi*exp(-csi*xfik) + exp(-csi)

            # Acumula a correção de derivada
            dLo[m] += dL[k]*(Hmk/soma)*heavi

        end #k

    end #m

    return dLo
end



##################### NÃO IMPLEMENTADO AINDA! #################################

#
#   Não implementado ainda!
#
function Filtro_Heavi_Tan(x::Array{Float64,1},nel::Int64,csi)

    # Inicializa o x filtrado

    # Calcula filtro de densidades
    #se csi = 0.0, retorna filtro de densidades
    #se csi > 0.0, heaviside

    # Soma proveniente do filtro de densidades
    sum_xd = sum(x)

    #resolve o filtro para csi
    a = 0.0
    b = 1.0
    c = 0.5

    tol = 1E-2

    x_a = Array{Float64}(undef,nel)
    x_c = Array{Float64}(undef,nel)
    xh  = Array{Float64}(undef,nel)

    # Procura o valor de eta
    while true

        # Criério de parada
        if (b-a)<=tol
            break
        end

        # Tanh dos valores inferiores e médios
        thcn_a = tanh(csi*a)
        thcn_c = tanh(csi*c)

        # Forma os vetores
        for j = 1:nel
            x_a[j] = (thcn_a + tanh(csi*(x[j]-a)))/(thcn_a+tanh(csi*(1.0-a)))
            x_c[j] = (thcn_c + tanh(csi*(x[j]-c)))/(thcn_c+tanh(csi*(1.0-c)))
        end

        # Verifica a diferença entre o valor original
        fun_a = sum(x_a) - sum_xd
        fun_c = sum(x_c) - sum_xd

        # Se houver uma raiz entre a e c, senão está em c e b, reitera
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

function dL_Heavi_Tan(dL::Array{Float64,1},nel::Int64,vizi::Array{Int64,2},
    nviz::Array{Int64,1},dviz::Array{Float64,2},raiof::Float64,xf::Array{Float64,1},csi::Float64,eta::Float64)

    # Inicializa o vetor das derivadas corrigidas
    dLo = zeros(Float64,nel)

    # Varre m elementos
    for m=1:nel

        # Define vizinhos de m
        vizm = vizi[m,1:nviz[m]]

        # Varre os vizinhos de m
        for k in vizm

            # Posição de m na tabela de vizinhanca
            pos = findfirst(equalto(m),vizi[k,1:nviz[k]])

            # Calcula o peso (H) entre m e k
            Hmk = 1.0 - dviz[k,pos]/raiof

            # Loop para denominador da filtragem, Hjk
            soma = 0.0
            for j=1:nviz[k]
                soma += 1.0 - dviz[k,j]/raiof
            end #j

            # Acumula a correção de derivada
            dH = (csi*sech(csi*(xf[k]-eta))^2.0)/(tanh(csi*eta)+tanh(csi*(1.0-eta)))
            dLo[m] += dL[k]*(Hmk/soma)*dH

        end #k

    end #m

    return dLo
end
