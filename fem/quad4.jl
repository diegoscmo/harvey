# Funções para quad4 (pesos de integração omitidos == 1.0)

# Rotina que calcula a matriz de rigidez de um elemento da malha
function Kquad4(ele::Int64, coord::Array{Float64,2},ijk::Array{Int64,2},
               young::Float64, poisson::Float64, esp::Float64)

# Define os arrays locais
  B  = zeros(Float64,3,8)
  no_coord = zeros(Float64,2,4);
  K = zeros(Float64,8,8);
  DJ = 0.0

# Nos do elemento
  no = ijk[ele,:]

# determina as coordenadas dos nos
  for i=1:4
     no_coord[1,i] = coord[no[i],1]
     no_coord[2,i] = coord[no[i],2]
  end

# tabela de pontos de Gauss
   valor = 1.0/sqrt(3.0)
   gpoint = [ -valor    valor     valor   -valor ;
              -valor   -valor     valor    valor ]

# Agora define o tensor constitutivo EPT
   com = (young/(1.0-poisson^2));
   De = [      com               com*poisson          0.0 ;
           com*poisson              com               0.0 ;
               0.0                  0.0           com*(1.-poisson)/2.0 ]

# Integra numericamente com 4 pontos de Gauss
   for p=1:4
       (B,DJ) = Bquad4(gpoint[1,p],gpoint[2,p],no_coord)
        K = K + (B'*De*B)*DJ;
   end

# multiplica pela espessura
   K = K*esp

# Fim da funcao Kquad4
  return K

end #Kquad4

function Nquad4(r::Float64,s::Float64)

    # As funcoes de interpolacao
          N1 = 1.0/4.0*(1.0-r)*(1.0-s)
          N2 = 1.0/4.0*(1.0+r)*(1.0-s)
          N3 = 1.0/4.0*(1.0+r)*(1.0+s)
          N4 = 1.0/4.0*(1.0-r)*(1.0+s)

    # E o vetor com as funcoes
          N = [N1 N2 N3 N4]

          return N

      end

function Bquad4(r::Float64,s::Float64,coord::Array{Float64,2})

# Define as dimensoes dos arrays locais
      J    = zeros(Float64,2,2)
      invJ = zeros(Float64,2,2)
      Ae   = zeros(Float64,4)
      Be   = zeros(Float64,4)

# Calcula a matriz jacobiana
       t1 = -1.0+s
       t4 = +1.0+s
      t13 = -1.0+r
      t15 = -1.0-r
       J  = 0.25* [(t1*coord[1,1]-t1*coord[1,2]+t4*coord[1,3]-t4*coord[1,4]) (t1*coord[2,1]-t1*coord[2,2]+t4*coord[2,3]-t4*coord[2,4]) ;
             (t13*coord[1,1]+t15*coord[1,2]-t15*coord[1,3]-t13*coord[1,4]) (t13*coord[2,1]+t15*coord[2,2]-t15*coord[2,3]-t13*coord[2,4])]

# Inverte a matriz jacobiana
   invJ = inv(J)

# Monta a matriz B
# dNi/dr
      Ae[1] = invJ[1,1]*t1+invJ[1,2]*t13
      Ae[2] = -invJ[1,1]*t1+invJ[1,2]*t15
      Ae[3] = invJ[1,1]*t4-invJ[1,2]*t15
      Ae[4] = -invJ[1,1]*t4-invJ[1,2]*t13

# dNi/ds
      Be[1] = invJ[2,1]*t1+invJ[2,2]*t13
      Be[2] = -invJ[2,1]*t1+invJ[2,2]*t15
      Be[3] = invJ[2,1]*t4-invJ[2,2]*t15
      Be[4] = -invJ[2,1]*t4-invJ[2,2]*t13
# B
      B = 0.25* [ Ae[1]    0.0  Ae[2]    0.0  Ae[3]    0.0  Ae[4]    0.0 ;
                    0.0  Be[1]    0.0  Be[2]    0.0  Be[3]    0.0  Be[4] ;
                  Be[1]  Ae[1]  Be[2]  Ae[2]  Be[3]  Ae[3]  Be[4]  Ae[4] ]

# calcula o determinante da matriz Jacobiana
      DJ = det(J)

# Fim da rotina Bquad4
  return B, DJ

end


function Mquad4(ele::Int64, coord::Array{Float64,2},ijk::Array{Int64,2},
               esp::Float64,rho::Float64)

# Define as dimensoes dos arrays locais
         M    = zeros(Float64,4,4)
    nos_coord = zeros(Float64,2,4)

# Determina as coordenadas dos nos
    for i=1:4
      pos = convert(Int64,ijk[ele,i])
      nos_coord[:,i] = coord[pos,:]
    end

# Tabela de pontos de Gauss
    valor = 1.0/sqrt(3.0)
    gpoint = [ -valor    valor     valor   -valor ;
               -valor   -valor     valor    valor ]

# Integra nos quatro pontos de Gauss
    for p=1:4
      (B,DJ) = Bquad4(gpoint[1,p],gpoint[2,p],nos_coord)
      Nr     = Nquad4(gpoint[1,p],gpoint[2,p])
      M = M + rho*(Nr'*Nr)*DJ
    end
# Multiplica pela espura
    M = M * esp

end


function Sigma_quad4(r::Float64,s::Float64,ele::Int64,U::Array{Float64,1},
                             coord::Array{Float64,2},ijk::Array{Int64,2},
                                          young::Float64,poisson::Float64)

# Define as dimensões e zera os arrays locais
        sigma = zeros(3,1)
        eps   = zeros(3,1)
    ele_coord = zeros(2,4)

# Determina as coordenadas dos nos
    for i=1:4
      pos = convert(Int64,ijk[ele,i])
      ele_coord[:,i] = coord[pos,:]
    end

# Agora define o tensor constitutivo EPT
   comum = (young/(1.0-poisson^2))
   De = [   comum             comum*poisson          0.0;
            comum*poisson     comum                  0.0;
            0.0                 0.0           comum*(1.0-poisson)/2.0 ]

# Extrai os deslocamentos do vetor global.
  no  = ijk[ele,:]
  pos = [2*no[1]-1; 2*no[1]; 2*no[2]-1;  2*no[2];
         2*no[3]-1; 2*no[3]; 2*no[4]-1;  2*no[4]]

# Recupera matriz B para o ponto de integracao.
  (B,) = Bquad4(r,s,ele_coord)

# Relaciona deformação e tensão.
    eps   = B*U[pos]
    sigma = De*eps

    return sigma
end


function Tensoes_quad4(nelems::Int64,coord::Array{Float64,2},ijk::Array{Int64,2},
                            young::Float64,poisson::Float64,U::Array{Float64,1})

# Matriz com as tensoes (Sxx, Syy, Sxy).
    tensoes = zeros(nelems,3)

# Tensão no centróide de cada elemento (r=0.0, s=0.0):
    for ele=1:nelems
      sigma = Sigma_quad4(0.0,0.0,ele,U,coord,ijk,young,poisson)
      tensoes[ele,:] = sigma'
    end #fors
    return tensoes
end #Tensoes_quad4


function Fquad4(ele::Int64,coord::Array{Float64,2},ijk::Array{Int64,2},
                               esp::Float64,sigma::Float64)

# Sigma é a tensão calculada no elemento - Cte para o bilinear isoparamétrico

# Define as dimensões e zera os arrays.
    nos_coord = zeros(2,4)
      forcas  = zeros(4,1)

# Determina as coordenadas dos nos
    for i=1:4
      pos = convert(Int64,ijk[ele,i])
      nos_coord[:,i] = coord[pos,:]
    end

# Tabela de pontos de Gauss
    valor = 1.0/sqrt(3.0)
    gpoint = [ -valor    valor     valor   -valor ;
               -valor   -valor     valor    valor ]

# Integra nos quatro pontos de Gauss
    for p=1:4
      (B,DJ) = Bquad4(gpoint[1,p],gpoint[2,p],nos_coord)
      Nr     = Nquad4(gpoint[1,p],gpoint[2,p])
      forcas = forcas + 1.0*(Nr'*sigma)*DJ # 1.0 eh a densidade
    end
    forcas = forcas*esp

    return forcas
end #Fquad4
