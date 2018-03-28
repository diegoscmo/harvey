################################################################################
#####         Elemento Quadrilateral Isoparamétrico Incompatível          ######
################################################################################
#
# Rotina que calcula a matriz de rigidez de um elemento da malha
#
function Kquad4_I(elem::Int64, coord::Array{Float64,2}, conect::Array{Int64,2},
    elast::Float64, poiss::Float64, espess::Float64)

    # Define os arrays locais
    coord_nos  = Array{Float64}(undef,2,4)
    K_I  = zeros(12,12)
    CBA  = zeros(3,8,4)

    # Determina as coordenadas dos nós
    for i=1:4
        coord_nos[1,i] = coord[conect[elem,i],1]
        coord_nos[2,i] = coord[conect[elem,i],2]
    end

    # tabela de pontos de Gauss
    valor = 1.0/sqrt(3.0)
    gpoint = [ -valor    valor     valor   -valor ;
               -valor   -valor     valor    valor ]

    # Agora define o tensor constitutivo EPT
    com = (elast/(1.0-poiss^2));
    De = [  com          com*poiss   0.0 ;
            com*poiss    com         0.0 ;
            0.0          0.0         com*(1.0-poiss)/2.0 ]

    # Integra numericamente com 4 pontos de Gauss
    for p=1:4
        (B,DJ) = Bquad4_I(gpoint[1,p],gpoint[2,p],coord_nos)
        K_I .= K_I + (transpose(B)*De*B)*DJ;
    end
    K_I = K_I*espess

    # Incompatível
    Kaa = view(K_I,1:8,1:8)
    #Kab = K_I[1:8,9:12]
    Kba = view(K_I,9:12,1:8)
    Kbb = view(K_I,9:12,9:12)

    invKbbKba = inv(Kbb)*Kba

    K = Kaa - transpose(Kba)*invKbbKba


    # Matriz A
    A = vcat(Matrix{Float64}(I,8,8),-invKbbKba)

    # Gera matrizes CBA
    for p = 1:4
        B, = Bquad4_I(gpoint[1,p],gpoint[2,p],coord_nos)
        CBA[:,:,p] = De*B*A
    end

    # Fim da funcao Kquad4
    return K,CBA

end #Kquad4

function Bquad4_I(r::Float64,s::Float64,coord_nos::Array{Float64,2})

    # Criando matrizes necessarias
    DN   = Array{Float64}(undef,2,6)
    B    = Array{Float64}(undef,3,12)

    #DN eh a matriz com as derivadas da funcao de forma em relacao a r e s
    DN[1,1] = (-1.0+s)/4.0
    DN[1,2] = (+1.0-s)/4.0
    DN[1,3] = (+1.0+s)/4.0
    DN[1,4] = (-1.0-s)/4.0
    DN[2,1] = (-1.0+r)/4.0
    DN[2,2] = (-1.0-r)/4.0
    DN[2,3] = (+1.0+r)/4.0
    DN[2,4] = (+1.0-r)/4.0
    # Derivadas incompatíveis
    DN[1,5] = -2.0*r
    DN[1,6] = 0.0
    DN[2,5] = 0.0
    DN[2,6] = -2.0*s

    #Calculo da matriz jacobiana
    J  = DN[:,1:4]*transpose(coord_nos)

    #Inverte a matriz jacobiana
    invJ = inv(J)

    #calcula o determinante da matriz Jacobiana
    DetJ = det(J)

    #Monta a matriz B
    Baux = invJ*DN
    for i=1:size(DN,2)
        B[1,2*i-1]  = Baux[1,i]
        B[1,2*i]    = 0.0
        B[2,2*i-1]  = 0.0
        B[2,2*i]    = Baux[2,i]
        B[3,2*i-1]  = Baux[2,i]
        B[3,2*i]    = Baux[1,i]
    end

    # Fim da rotina Bquad4
    return B, DetJ

end

function Tquad4_I(pdens::Array{Float64,1}, nel::Int64, SP::Float64, SQ::Float64,
                ijk::Array{Int64,2}, CBA::Array{Float64,3}, Ux::Array{Float64,1}) #U expandido

    # Define as dimensões e zera os arrays locais
    # [sigma_xx, sigma_yy, sigma_xy] * 4 pontos de Gauss
    sigma = Array{Float64}(undef,nel,12)
    SPQ   = SP-SQ

    # Loop pelos elementos
    for j=1:nel

        # Adquire os nós do elemento
        no1 = ijk[j,1]
        no2 = ijk[j,2]
        no3 = ijk[j,3]
        no4 = ijk[j,4]

        # Mapeia o graus de liberdade do elemento
        gdl = [2*no1-1,2*no1,2*no2-1,2*no2,
               2*no3-1,2*no3,2*no4-1,2*no4]

        # Deslocamentos do elemento
        Ue = Ux[gdl]

        # Tensões nos 4 pontos de Gauss
        S_aux = pdens[j]^SPQ

        for k=1:4
            sigma[j,(k*3-2):(k*3)]   = S_aux*(CBA[:,:,k]*Ue)
        end

    end # for j

    return sigma
end # function

function Tquad4_I_Din(pdens::Array{Float64,1}, nel::Int64, SP::Float64, SQ::Float64,
                ijk::Array{Int64,2}, CBA::Array{Float64,3}, Ux::Array{Complex{Float64},1}, beta::Float64, w::Float64) #U expandido complexo

    # Define as dimensões e zera os arrays locais
    # [sigma_xx, sigma_yy, sigma_xy] * 4 pontos de Gauss
    sigma = Array{Complex{Float64}}(undef,nel,12)
    SPQ   = SP-SQ

    # Loop pelos elementos
    for j=1:nel

        # Adquire os nós do elemento
        no1 = ijk[j,1]
        no2 = ijk[j,2]
        no3 = ijk[j,3]
        no4 = ijk[j,4]

        # Mapeia o graus de liberdade do elemento
        gdl = [2*no1-1,2*no1,2*no2-1,2*no2,
               2*no3-1,2*no3,2*no4-1,2*no4]

        # Deslocamentos do elemento
        Ue = Ux[gdl]

        # Tensões nos 4 pontos de Gauss
        S_aux = (1.0+beta*im*w)*pdens[j]^SPQ

        for k=1:4
            sigma[j,(k*3-2):(k*3)]   = S_aux*(CBA[:,:,k]*Ue)
        end

    end # for j

    return sigma
end # function


function Nquad4(r::Float64,s::Float64)

    # As funcoes de interpolacao
    N1 = 0.25*(1.0-r)*(1.0-s)
    N2 = 0.25*(1.0+r)*(1.0-s)
    N3 = 0.25*(1.0+r)*(1.0+s)
    N4 = 0.25*(1.0-r)*(1.0+s)

    # E o vetor com as funcoes
    #N = [N1 N2 N3 N4]

    N = [N1  0.0 N2  0.0 N3  0.0 N4  0.0
         0.0 N1  0.0 N2  0.0 N3  0.0 N4]

    return N

end

function Mquad4(elem::Int64, coordx::Array{Float64,2},conect::Array{Int64,2},
    espess::Float64,rho::Float64)

    # Define as dimensoes dos arrays locais
    M   = zeros(8,8)
    coord_nos = zeros(2,4)

    # determina as coordxenadas dos nos
    for i=1:4
        coord_nos[1,i] = coordx[conect[elem,i],1]
        coord_nos[2,i] = coordx[conect[elem,i],2]
    end

    # Tabela de pontos de Gauss
    valor = 1.0/sqrt(3.0)
    gpoint = [ -valor    valor     valor   -valor ;
               -valor   -valor     valor    valor ]

    # Integra nos quatro pontos de Gauss
    for p=1:4
        (B,DJ) = Bquad4_I(gpoint[1,p],gpoint[2,p],coord_nos)
        Nr     = Nquad4(gpoint[1,p],gpoint[2,p])
        M      = M + rho*(Nr'*Nr)*DJ
    end
    # Multiplica pela espura
    M = M * espess

end
