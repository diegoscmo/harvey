function GeraMalha(LX::Float64,LY::Float64,NX::Int64,NY::Int64,npresos::Int64,presos::Array{Float64,2},
  nforcas::Int64, forcas::Array{Float64,2},lxmin,lymin,cube)

  # Define a menor dimensao de um elemento, para fins de tolerancia dimensional
  deltax = LX/NX
  deltay = LY/NY

  # Aloca matriz de coordenadas nodais
  nnos = (NX+1)*(NY+1)
  nelems = NX*NY
  elcoor = zeros(nnos,2)

  # Se a interpolacao cubica dos cantos for ativada
  if cube == true

    deltax = lxmin
    deltay = lymin

    # Gera um sistema de minimos quadrados para determinar os coefs
    a = zeros(4,2)
    n = 5.0
    for i=1:2
        # fit em x
      if i==1
        xfit = [0.0  1.0    NX/2.0  NX-1.0   NX]
        yfit = [0.0  lxmin  LX/2.0  LX-lxmin LX]
        # fit em y
      elseif i==2
        xfit = [0.0  1.0    NY/2.0  NY-1.0   NY]
        yfit = [0.0  lymin  LY/2.0  LY-lymin LY]
      end

      # Somas
      sx = sum(xfit)
      sx2 = sum(xfit.^2)
      sx3 = sum(xfit.^3)
      sx4 = sum(xfit.^4)
      sx5 = sum(xfit.^5)
      sx6 = sum(xfit.^6)
      sy = sum(yfit)
      sxy = sum(xfit.*yfit)
      sx2y = sum(xfit.^2.*yfit)
      sx3y = sum(xfit.^3.*yfit)

      # Determinacao dos coefs
      m1 = [ n   sx  sx2 sx3;
             sx  sx2 sx3 sx4;
             sx2 sx3 sx4 sx5;
             sx3 sx4 sx5 sx6 ]

      m2 = [ sy sxy sx2y sx3y ]'

      a[:,i] = m1\m2

    end #for
  end #if

  # Define a tolerancia geometrica para as janelas
  tolerancia= min(deltax,deltay)/10.0

  # Gera as coordenadas dos nos
  no= 1
  coordy = 0.0
  for coluna=1:(NY+1)
    coordx::Float64 = 0.0
    for linha=1:(NX+1)

      # Para cada coordenada y fixa,
      # varre todos os possiveis nos
      elcoor[no,1] = coordx
      elcoor[no,2] = coordy

      # O mesmo laco faz os dois tipos de posicionamento dos nos
      if cube == true
        coordx = a[1,1] + a[2,1]*linha + a[3,1]*linha^2 + a[4,1]*linha^3
      else
        coordx = coordx + deltax
      end
      no = no + 1
    end #linha

    if cube == true
      coordy = a[1,2] + a[2,2]*coluna + a[3,2]*coluna^2 + a[4,2]*coluna^3
    else
      coordy = coordy + deltay
    end
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
      x = round(elcoor[no,1],10)
      y = round(elcoor[no,2],10)

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

      xn = round(elcoor[no,1],10)
      yn = round(elcoor[no,2],10)

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
  const gll = zeros(Int64,8)

  # graus de liberdade globais
  const glg = zeros(Int64,8)


  const contador_local::Int64 = 1
  const contador_global::Int64 = 1
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
