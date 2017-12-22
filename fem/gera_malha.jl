################################################################################
#####                          Geração da Malha                           ######
################################################################################

#
#
#
function GeraMalha(nnos::Int64, nelems::Int64, LX::Float64, LY::Float64, NX::Int64, NY::Int64,
                   presos::Array{Float64,2}, forcas::Array{Float64,2})

   # Define o tamanho dos vetores das condições de contorno
   npresos = size(presos,1)            # Nr. de apoios
   nforcas = size(forcas,1)            # Nr. de carregamentos

   # Define a menor dimensao de um elemento, para fins de tolerancia dimensional
   deltax = LX/NX
   deltay = LY/NY

   # Define a tolerancia geometrica para as janelas
   tolerancia= min(deltax,deltay)/10.0

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

           if (x>=xi) && (x<=xf) && (y>=yi) && (y<=yf)
              ngl_preso +=1
              gls_presos[ngl_preso,1] = no
              gls_presos[ngl_preso,2] = dir
           end #if

       end #no

   end #i


  #
  # Aplica forcas nodais
  # Informacoes em cada linha de forcas
  # x,y,valor,direcao
  nos_forcas = zeros(nforcas,3) # no, dir, valor

  for i=1:nforcas

      # Dados do nó
      x = forcas[i,1]
      y = forcas[i,2]
      valor = forcas[i,3]
      dir = forcas[i,4]

      # Loop pelos nós
      gravou = false
      for no=1:nnos

          xn = elcoor[no,1]
          yn = elcoor[no,2]

          if abs(x-xn)<=tolerancia && abs(y-yn)<=tolerancia
             nos_forcas[i,1] = no
             nos_forcas[i,2] = dir
             nos_forcas[i,3] = valor
             gravou = true
             break
           end #if
      end #loop dos nos

      # Testa se a forca achou o seu nó
      if !gravou
         error("A forca nao foi aplicada")
      end

  end #loop das informacoes de forca


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
    @inbounds for no in ijk[elemento,:]
       @inbounds  for j = 1:2
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
