function Global(nelems,ijk,ID,K0,densidades,simp,nforcas,nos_forcas,ngl_efetivos)

     # Monta a matriz Global de Rigidez
     contador = 1

     I = zeros(Int64,8*8*nelems)
     J = zeros(Int64,8*8*nelems)
     V = zeros(8*8*nelems)
     for el = 1:nelems
        # rotina para determinacao dos vetores para montagem da matriz esparsa
        fator::Float64 = densidades[el]^simp
        (glg,gll) = gl_livres_elemento(el,ijk,ID)
         @inbounds  for i = 1:length(glg)
                        glgi =  glg[i]
                        glli =  gll[i]
                        @inbounds for j = 1: length(glg)
                                  glgj =  glg[j]
                                  gllj =  gll[j]
                                  I[contador] = glgi
                                  J[contador] = glgj
       	                          V[contador] = K0[glli,gllj]*fator
                                  contador = contador + 1
                        end #j
       end #i
    end #e

    contador = contador - 1

    KG = sparse(I[1:contador],J[1:contador],V[1:contador])

    # Monta o vetor de forcas
    F = zeros(ngl_efetivos)
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


    return KG,F
end


# Rotina que expande um vetor global compactado (sem os gls restritos)
# para um vetor com a dimensão global original e com 0.0 nas
# posições restritas
function Expande_Vetor(entrada, nnos, ID)

   # Aloca o vetor de saida
   saida = zeros(2*nnos)
   saida = complex(saida)
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
