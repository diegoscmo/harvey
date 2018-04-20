################################################################################
#####                          Geração da Malha                           ######
################################################################################

#
#
#
"""
   Gera a malha,,,,
"""
function GeraMalha(nnos::Int64, nelems::Int64, LX::Float64, LY::Float64, NX::Int64, NY::Int64,
                   presos::Array{Float64,2}, forcas::Array{Float64,2})

   # Define o tamanho dos vetores das condições de contorno
   npresos = size(presos,1)            # Nr. de apoios
   nforcas = size(forcas,1)            # Nr. de carregamentos

   # Define a menor dimensao de um elemento, para fins de tol dimensional
   deltax = LX/NX
   deltay = LY/NY

   # Define a tol geometrica para as janelas
   tol= min(deltax,deltay)/10.0

   # Aloca matriz de coordenadas nodais
   elcoor = zeros(nnos,2)

   # Gera as coordenadas dos nos
   no::Int64 = 1
   coordy::Float64 = 0.0
   for coluna=1:(NY+1)
       coordx::Float64 = 0.0
       for linha=1:(NX+1)
           # Para cada coordenada y fixa,
           # varre todos os possiveis nos
           elcoor[no,1] = coordx
           elcoor[no,2] = coordy
           coordx = coordx + deltax
           no = no + 1
        end #linha
        coordy = coordy + deltay;
    end #coluna

   # Gera as conectividades
   ijk = zeros(Int64,nelems,4)

   noi::Int64 = 1
   ele::Int64 = 1
   for linha=1:NY
      for coluna=1:NX
         noj  = noi + 1
         nok  = noi + (NX + 2)
         nol  = noi + (NX + 1)
         ijk[ele,1] = noi
         ijk[ele,2] = noj
         ijk[ele,3] = nok
         ijk[ele,4] = nol
         ele = ele + 1
         noi = noi + 1
      end # coluna
      # Ja somamos 1 no final do loop anterior...
      noi = noi + 1
   end # linha

   # Para cada região com informacoes de condições de contorno,
   # varremos todos os nós que podem estar dentro do retângulo
   # Informacoes em cada linha de presos: xi, yi, xf, yf, direcao
   ngl_preso::Int64 = 0
   gls_presos = zeros(Int64,nnos,2)

   for i=1:npresos

       # Pega as informacoes do retangulo
       xi = presos[i,1]
       yi = presos[i,2]
       xf = presos[i,3]
       yf = presos[i,4]
       dir = convert(Int64,presos[i,5])

       # Loop por todos os nos da malha
       for no=1:nnos

           # Pega as coordenadas do nó
           x = elcoor[no,1]
           y = elcoor[no,2]

           if (x>=xi-tol) && (x<=xf+tol) && (y>=yi-tol) && (y<=yf+tol)
              ngl_preso +=1
              gls_presos[ngl_preso,1] = no
              gls_presos[ngl_preso,2] = dir
           end #if
       end #no
   end #i

   # Monta a matriz ID, que informa o grau de liberdade global
   # associado a um no/direcao
   ID = zeros(Int64,nnos,2)
   contador_global::Int64 = 1

   for no=1:nnos
     for gl=1:2
        # Procura por no,gl nas informacoes
        # de condicoes de contorno. Se existirem,
        # pula o gl e deixa com zero. Do contrario,
        # coloca o contador_global e incrementa
        flag_preso = false
        for i=1:ngl_preso
            if no==gls_presos[i,1] && gl==gls_presos[i,2]
               flag_preso = true
               break
            end #if
        end #i
         if !flag_preso
             ID[no,gl] = contador_global
             contador_global += 1
         end
    end #gl
   end #no
   contador_global = contador_global - 1


   ### MONTA AS FORCAS ###



   # Para cada região com informacoes de condições de contorno,
   # varremos todos os nós que podem estar dentro do retângulo
   # Informacoes em cada linha de presos: xi, yi, xf, yf, direcao

   nos_forcas  = zeros(nnos,3)

   nforcas_tot = 0

   for i=1:nforcas

       ngl_forcas::Int64 = 0
       gls_forcas = zeros(Int64,nnos,3)

       # Pega as informacoes do retangulo
       xi = forcas[i,1]
       yi = forcas[i,2]
       xf = forcas[i,3]
       yf = forcas[i,4]
       dir = convert(Int64,forcas[i,6])

       # Loop por todos os nos da malha pra ver onde estão as forças i
       for no=1:nnos

           # Pega as coordenadas do nó
           x = elcoor[no,1]
           y = elcoor[no,2]

           if (x>=xi-tol) && (x<=xf+tol) && (y>=yi-tol) && (y<=yf+tol)
              ngl_forcas +=1
              gls_forcas[ngl_forcas,1] = no
              gls_forcas[ngl_forcas,2] = dir
           end #if

       end #no

       # corta a gordura do vetor
      gls_forcas =  gls_forcas[1:ngl_forcas,:]

      #
      # Aplica forcas nodais "locais" (de i)
      # x,y,valor,direcao
      loc_forcas = zeros(ngl_forcas,3) # no, dir, valor

      # Insere as forças em um vetor local
      for j=1:ngl_forcas

          loc_forcas[j,1] = gls_forcas[j,1]
          loc_forcas[j,2] = gls_forcas[j,2]
          loc_forcas[j,3] = forcas[i,5]/(ngl_forcas-1)

          # Se for primeiro ou ultimo nó, divide a força por 2
          if j==1 || j==ngl_forcas
              loc_forcas[j,3] = loc_forcas[j,3]/2.0
          end

      end #loop das informacoes de forca


      # Acumula as forças no global
      nos_forcas[1+nforcas_tot:nforcas_tot+ngl_forcas,:] = loc_forcas[:,:]
      nforcas_tot += ngl_forcas

   end #i

   # trima o vetor de forcas
   nos_forcas = nos_forcas[1:nforcas_tot,:]

   # Retorna as informacoes da malha
   return elcoor, ijk, nos_forcas, ID, contador_global

end



function gl_livres_elemento(elemento::Int64,ijk::Array{Int64,2},ID::Array{Int64,2})

    # graus de liberdade locais
    gll = zeros(Int64,8)

    # graus de liberdade globais
    glg = zeros(Int64,8)

    contador_local::Int64 = 1
    contador_global::Int64 = 1
    for no in ijk[elemento,:]
       for j = 1:2
            gl = ID[no,j]
            if gl > 0
               glg[contador_global] = gl
               gll[contador_global] = contador_local
               contador_global +=1
            end
            contador_local +=1
       end
    end

    contador_global = contador_global - 1

    return (glg[1:contador_global],gll[1:contador_global])

end
