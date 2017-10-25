function Global(nelems,ijk,ID,K0,M0,densidades,simp,nforcas,nos_forcas,ngl_efetivos)

    # Corrige os fatores densidade para a massa (Olhoff e Du)
    massadens = copy(densidades)
    C1 =  6.0E5
    C2 = -5.0E6
    # Corrige massa de acordo com Olhoff & Du
    for j=1:nelems
        if massadens[j] < 0.1
            massadens[j] = C1*massadens[j]^6 + C2*massadens[j]^7
        end #if dens[i]
    end #for i

     # Monta a matriz Global de Rigidez e Massa
     contador = 1

     I = zeros(Int64,8*8*nelems)
     J = zeros(Int64,8*8*nelems)
     V = zeros(8*8*nelems)
     W = zeros(8*8*nelems)
     for el = 1:nelems
        # rotina para determinacao dos vetores para montagem da matriz esparsa
        # Modificação Olhoff & Du pag 94/95 eq 5.19-20
        kfator::Float64 = (0.001+densidades[el]*(0.999))^simp
        mfator::Float64 = (0.001+massadens[el]*(0.999))
        (glg,gll) = gl_livres_elemento(el,ijk,ID)
         @inbounds  for i = 1:length(glg)
                        glgi =  glg[i]
                        glli =  gll[i]
                        @inbounds for j = 1: length(glg)
                                  glgj =  glg[j]
                                  gllj =  gll[j]
                                  I[contador] = glgi
                                  J[contador] = glgj
       	                          V[contador] = K0[glli,gllj]*kfator
                                  W[contador] = M0[glli,gllj]*mfator
                                  contador = contador + 1
                        end #j
       end #i
    end #e

    contador = contador - 1

    KG = sparse(I[1:contador],J[1:contador],V[1:contador])
    MG = sparse(I[1:contador],J[1:contador],W[1:contador])

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


    return KG,MG,F
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
