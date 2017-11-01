# Gera uma variável com as coordenadas do centro de cada elemento
# Depois converte isso para a matriz H
function Calc_H(coord,conect,nelems,rmin)

    # Inicializa
    centros = zeros(Float64,nelems,2)
    sumH    = zeros(Float64,nelems)
    H       = zeros(Float64,nelems,nelems)

    # Acumula as coordenadas referentes aos nós de cada elemento
    for j=1:nelems
        nos = conect[j,:]
        for k = 1:4
            centros[j,1] += coord[nos[k],1]
            centros[j,2] += coord[nos[k],2]
        end # for k

    end # for j

    # Faz a média, definindo as posições centrais
    centros = centros/4.0

    # Calcula Hj de k
    for j=1:nelems
        for k=1:nelems
            dkj = sqrt((centros[k,1]-centros[j,1])^2.0+(centros[k,2]-centros[j,2])^2.0)
            if dkj <= rmin
                H[j,k]= rmin - dkj
            end
        end
    end

    for j=1:nelems
        sumH[j] = sum(H[:,j])
    end

    return H,sumH
end

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

# Aplica filtro de densidades
function Filtra_Dens(xo)

    # Inicializa vetor
    xf = zeros(Float64,nelems)

    # Varre os elementos filtrando-os
    for k=1:nelems
        xf[k] = sum(H[k,:].*xo)/sumH[k]
    end

    return xf
end

# Corrige derivada devido ao filtro de densidades (derivada)
function Corrige_dL(x, dL)

    # Inicializa o vetor das derivadas corrigidas
    dLc = zeros(Float64,nelems)

    #FIXME? Jeito mais fácil de fazer? Aqui no filtro eu multiplico a matriz inteira,
    # talvez vale a pena listar os vizinhos...

    return dLc
end
