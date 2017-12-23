################################################################################
#####                     Montagem da Matriz Global                       ######
################################################################################
#
#
#
function Global_KM(densidades::Array{Float64,1}, nelems::Int64, ijk::Array{Int64,2},
                   ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                   simp::Float64, vminimo::Float64)

    # Corrige os fatores densidade para a massa (Olhoff e Du)
    massadens = copy(densidades)

    # Corrige massa de acordo com Olhoff & Du
    for j=1:nelems
        if massadens[j] < 0.1
            massadens[j] = 6E5*massadens[j]^6.0 - 5E6*massadens[j]^7.0
        end #if dens[i]
    end #for i

     # Monta a matriz Global de Rigidez e Massa
     contador = 1

     # Para correção das propriedades (E_min / rho_min)
     vsuper = 1.0 - vminimo

     # Inicializa vetores
     I = Array{Int64}(uninitialized,8*8*nelems)
     J = Array{Int64}(uninitialized,8*8*nelems)
     V = Array{Float64}(uninitialized,8*8*nelems)
     W = Array{Float64}(uninitialized,8*8*nelems)

     for el = 1:nelems
        # rotina para determinacao dos vetores para montagem da matriz esparsa
        # Modificação Olhoff & Du pag 94/95 eq 5.19-20
        kfator = densidades[el]^simp
        mfator = massadens[el]
        (glg,gll) = gl_livres_elemento(el,ijk,ID)
         @inbounds for i = 1:length(glg)
                        glgi =  glg[i]
                        glli =  gll[i]
                        @inbounds for j = 1: length(glg)
                                  glgj =  glg[j]
                                  gllj =  gll[j]
                                  # Monta só a parte superior
                                  #if glgj>=glgi
                                     I[contador] = glgi
                                     J[contador] = glgj

                                     # Corrige o valor mínimo na montagem da matriz
       	                             V[contador] = (vminimo + kfator*vsuper)*K0[glli,gllj]
                                     W[contador] = (vminimo + mfator*vsuper)*M0[glli,gllj]

                                     contador = contador + 1
                                #end # if U
                        end #j
       end #i
    end #e

    contador = contador - 1

    #KG = Symmetric(sparse(I[1:contador],J[1:contador],V[1:contador]))
    #MG = Symmetric(sparse(I[1:contador],J[1:contador],W[1:contador]))

    KG = sparse(I[1:contador],J[1:contador],V[1:contador])
    MG = sparse(I[1:contador],J[1:contador],W[1:contador])

    return KG,MG
end

function Global_F(ID::Array{Int64,2}, nos_forcas::Array{Float64,2}, ngl_efetivos::Int64)

    # Monta o vetor de forcas
    F = zeros(ngl_efetivos)
    nforcas = size(nos_forcas,1)

    for i=1:nforcas

          # No da forca
          nof = convert(Int64,nos_forcas[i,1])

          # Gl da forca
          glf = convert(Int64,nos_forcas[i,2])

          # Grau de liberdade efetivo
          gl = ID[nof,glf]

          # Se gl for valido, acrescenta ao vetor global de forcas
          if gl>0
             F[gl] = F[gl] + nos_forcas[i,3]
          end
    end
return F
end


# Rotina que expande um vetor global compactado (sem os gls restritos)
# para um vetor com a dimensão global original e com 0.0 nas
# posições restritas
function Expande_Vetor(entrada, nnos, ID, C::Bool=false)

   # Aloca o vetor de saida
   saida = zeros(2*nnos)

   # Se for complexo...
   if C == true
       saida = complex(saida)
   end

   # Loop pelas posicoes de ID, catando no vetor de entrada
   # a medida que ID permitir
   contador = 1
   for i=1:nnos
     for j=1:2
        gl = ID[i,j]
        if gl > 0
           saida[contador] = entrada[gl]
        end
        contador +=1
     end
  end
  return saida
end
