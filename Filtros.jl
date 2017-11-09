#
#   Identifica se há filtro em x e aplica
#
function Aplica_Filtro(x, filt)

    if filt.filtro == "Dens"
        x = Filtro_Dens(x, filt)
    elseif filt.filtro == "Off"
        x = copy(x)
    else
      error("Erro na declaracao do Filtro")
    end

    return x
end
#
#   Identifica se a derivada precisa de correção e a executa
#
function Derivada_Filtro(dLf, filt)

    # Declara a saída
    dLo = zeros(Float64,size(dLf,1))

    if filt.filtro == "Dens"
        dLo = dL_Dens(dLf, filt)
    elseif filt.filtro == "Off"
        dLo = copy(dLf)
    else
      error("Erro na declaracao do Filtro")
    end

    return dLo
end


#
# Identifica os vizinhos de cada elemento e a distancia entre eles
#
function Proc_Vizinhos(nelems,coord,conect,raiof)

    # Inicializa a matriz de vizinhos, distancias e quantidade de vizinhos
    vizi = zeros(Int64,nelems,nelems)
    dviz = zeros(Float64,nelems,nelems)
    nviz = zeros(Int64,nelems)

    # Inicializa centro dos elementos em estudo
    cj = zeros(Float64,2)
    ck = zeros(Float64,2)

    for j=1:nelems

        # Define o centro do elemento utilizando nós 1 e 3
        n1 = conect[j,1]
        n3 = conect[j,3]
        cj[1] = (coord[n1,1]+coord[n3,1])/2.0
        cj[2] = (coord[n1,2]+coord[n3,2])/2.0

        for k=1:nelems
            n1 = conect[k,1]
            n3 = conect[k,3]
            ck[1] = (coord[n1,1]+coord[n3,1])/2.0
            ck[2] = (coord[n1,2]+coord[n3,2])/2.0

            # Calcula a distancia entre os elementos
            dkj = sqrt((ck[1]-cj[1])^2.0+(ck[2]-cj[2])^2.0)

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
    vizi = vizi[:,1:big]
    dviz = dviz[:,1:big]

    return vizi,nviz,dviz
end

#
# Aplica filtro de densidades, pág 65
#
function Filtro_Dens(xi,filt)

    # Inicializa vetor de saída
    nelems = size(xi,1)
    xo = zeros(Float64,nelems)
    nviz  = filt.nviz
    vizi  = filt.vizi
    dviz  = filt.dviz
    raiof = filt.raiof

    # Varre os elementos
    for k = 1:nelems

        # Zera acumuladores
        dividen = 0.0
        divisor = 0.0

        # Acumula somatórios para cada vizinho de k
        for j = 1:nviz[k]
            #pesoH    = raiof - dviz[k,j]
            pesoH    = 1.0 - dviz[k,j]/raiof
            dividen += pesoH*xi[vizi[k,j]]
            divisor += pesoH
        end

        # Filtra elemento
        xo[k] = dividen/divisor

    end

    return xo
end
#
# Corrige derivada devido ao filtro de densidades (cadeia)
#
function dL_Dens(dLd, filt)

    # Inicializa o vetor das derivadas corrigidas
    dLo = zeros(Float64,size(dLd,1))
    nviz  = filt.nviz
    vizi  = filt.vizi
    dviz  = filt.dviz
    raiof = filt.raiof

    # Varre m elementos
    for m=1:size(dLd,1)

      # Define vizinhos de m
      vizm = vizi[m,1:nviz[m]]

      # Varre os vizinhos de m
      for k in vizm

        # Posição de m na tabela de vizinhanca
        pos = findfirst(vizi[k,1:nviz[k]],m)

        # Calcula o peso (H) entre m e k
        #Hmk = raiof - dviz[k,pos]
        Hmk = 1.0 - dviz[k,pos]/raiof

        # Loop para denominador da filtragem, Hjk
        soma = 0.0
        for j=1:nviz[k]
          #soma += raiof - dviz[k,j]
          soma += 1.0 - dviz[k,j]/raiof
        end #j

        # Acumula a correção de derivada
        dLo[m] += dLd[k]*(Hmk/soma)

      end #k

   end #m

    return dLo
end
