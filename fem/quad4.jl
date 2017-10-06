# Funções para quad4 (pesos de integração omitidos == 1.0)

# Rotina que calcula a matriz de rigidez de um elemento da malha
function Kquad4(elem::Int64, coordx::Array{Float64,2}, conect::Array{Int64,2},
    elast::Float64, poiss::Float64, espess::Float64)

    # Define os arrays locais
    B  = zeros(Float64,3,8)
    nos = zeros(Float64,2,4);
    K = zeros(Float64,8,8);
    DJ = 0.0

    # Determina as coordxenadas dos nós
    for i=1:4
        nos[1,i] = coordx[conect[elem,i],1]
        nos[2,i] = coordx[conect[elem,i],2]
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
        (B,DJ) = Bquad4(gpoint[1,p],gpoint[2,p],nos)
        K = K + (B'*De*B)*DJ;
    end

    # multiplica pela espess
    K = K*espess

    # Fim da funcao Kquad4
    return K

end #Kquad4

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

function Bquad4(r::Float64,s::Float64,nos::Array{Float64,2})

    # Criando matrizes necessarias
    J    = Array{Float64}(2,2)
    invJ = Array{Float64}(2,2)
    DN   = Array{Float64}(2,4)
    B    = Array{Float64}(3,8)

    #DN eh a matriz com as derivadas da funcao de forma em relacao a r e s
    DN[1,1] = (-1.0+s)/4.0
    DN[1,2] = (+1.0-s)/4.0
    DN[1,3] = (+1.0+s)/4.0
    DN[1,4] = (-1.0-s)/4.0
    DN[2,1] = (-1.0+r)/4.0
    DN[2,2] = (-1.0-r)/4.0
    DN[2,3] = (+1.0+r)/4.0
    DN[2,4] = (+1.0-r)/4.0

    #Calculo da matriz jacobiana
    J  = DN*nos'

    #Inverte a matriz jacobiana
    invJ = inv(J)

    #calcula o determinante da matriz Jacobiana
    DetJ = det(J)

    #Monta a matriz B
    Baux = invJ*DN
    for i=1:2:2*size(Baux,2)
        B[1,i]   = Baux[1,Int((i+1)/2)]
        B[2,i+1] = Baux[2,Int((i+1)/2)]
        B[3,i]   = Baux[2,Int((i+1)/2)]
        B[3,i+1] = Baux[1,Int((i+1)/2)]
    end

    # Fim da rotina Bquad4
    return B, DetJ

end


function Mquad4(elem::Int64, coordx::Array{Float64,2},conect::Array{Int64,2},
    espess::Float64,rho::Float64)

    # Define as dimensoes dos arrays locais
    M   = zeros(Float64,8,8)
    nos = zeros(2,4)

    # determina as coordxenadas dos nos
    for i=1:4
        nos[1,i] = coordx[conect[elem,i],1]
        nos[2,i] = coordx[conect[elem,i],2]
    end

    # Tabela de pontos de Gauss
    valor = 1.0/sqrt(3.0)
    gpoint = [ -valor    valor     valor   -valor ;
               -valor   -valor     valor    valor ]

    # Integra nos quatro pontos de Gauss
    for p=1:4
        (B,DJ) = Bquad4(gpoint[1,p],gpoint[2,p],nos)
        Nr     = Nquad4(gpoint[1,p],gpoint[2,p])
        M      = M + rho*(Nr'*Nr)*DJ
    end
    # Multiplica pela espura
    M = M * espess

end



#################################3




function Sigma_quad4(r::Float64,s::Float64,elem::Int64,U::Array{Float64,1},
    coordx::Array{Float64,2},conect::Array{Int64,2},
    elast::Float64,poiss::Float64)

    # Define as dimensões e zera os arrays locais
    sigma = zeros(3,1)
    eps   = zeros(3,1)
    ele_coordx = zeros(2,4)

    # Determina as coordxenadas dos nos
    for i=1:4
        pos = convert(Int64,conect[elem,i])
        ele_coordx[:,i] = coordx[pos,:]
    end

    # Agora define o tensor constitutivo EPT
    comum = (elast/(1.0-poiss^2))
    De = [   comum             comum*poiss          0.0;
    comum*poiss     comum                  0.0;
    0.0                 0.0           comum*(1.0-poiss)/2.0 ]

    # Extrai os deslocamentos do vetor global.
    no  = conect[elem,:]
    pos = [2*no[1]-1; 2*no[1]; 2*no[2]-1;  2*no[2];
    2*no[3]-1; 2*no[3]; 2*no[4]-1;  2*no[4]]

    # Recupera matriz B para o ponto de integracao.
    (B,) = Bquad4(r,s,ele_coordx)

    # Relaciona deformação e tensão.
    eps   = B*U[pos]
    sigma = De*eps

    return sigma
end


function Tensoes_quad4(nelems::Int64,coordx::Array{Float64,2},conect::Array{Int64,2},
    elast::Float64,poiss::Float64,U::Array{Float64,1})

    # Matriz com as tensoes (Sxx, Syy, Sxy).
    tensoes = zeros(nelems,3)

    # Tensão no centróide de cada elemento (r=0.0, s=0.0):
    for elem=1:nelems
        sigma = Sigma_quad4(0.0,0.0,elem,U,coordx,conect,elast,poiss)
        tensoes[elem,:] = sigma'
    end #fors
    return tensoes
end #Tensoes_quad4


function Fquad4(elem::Int64,coordx::Array{Float64,2},conect::Array{Int64,2},
    espess::Float64,sigma::Float64)

    # Sigma é a tensão calculada no elemento - Cte para o bilinear isoparamétrico

    # Define as dimensões e zera os arrays.
    nos_coordx = zeros(2,4)
    forcas  = zeros(4,1)

    # Determina as coordxenadas dos nos
    for i=1:4
        pos = convert(Int64,  conect[elem,i])
        nos_coordx[:,i] = coordx[pos,:]
    end

    # Tabela de pontos de Gauss
    valor = 1.0/sqrt(3.0)
    gpoint = [ -valor    valor     valor   -valor ;
    -valor   -valor     valor    valor ]

    # Integra nos quatro pontos de Gauss
    for p=1:4
        (B,DJ) = Bquad4(gpoint[1,p],gpoint[2,p],nos_coordx)
        Nr     = Nquad4(gpoint[1,p],gpoint[2,p])
        forcas = forcas + 1.0*(Nr'*sigma)*DJ # 1.0 eh a densidade
    end
    forcas = forcas*espess

    return forcas
end #Fquad4
