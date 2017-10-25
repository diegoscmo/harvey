function DFP(x::Array{Float64,1},mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
   xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64, search::String)

  numvar = size(x,1)

  # Calcula o gradiente de L
  dL,count = Dif_Fin(x, mult_res, rho, count)

  # inicializa aproximacao da hessiana
  G = eye(Float64,numvar)

  # Calcula a norma para criterio de parada interno
  norma = norm(dL)

  # Direcao de minimizacao
  dir = (-G*dL)/norma

  for i=1:max_int

    # Verifica o criterio de parada interno
    if norma < tol_int
      break
    end

    # Bloqueia a direcao e zera o gradiente se bloqueado
    for j=1:numvar
      if dir[j] < 0.0 && x[j] <= xl[j]
        dir[j] = 0.0
        dL[j]  = 0.0
      end
      if dir[j] > 0.0 && x[j] >= xu[j]
        dir[j] = 0.0
        dL[j]  = 0.0
      end
    end

    # Line search nesta direcao
    alpha,count = LS_Methods(x, mult_res, rho, dir, xl, xu, dL, tol_int, count, search)

    # incrementa a estimativa do ponto
    x = x + alpha*dir

    # Regaculta o gradiente e a norma
    dL_ant = copy(dL)
    dL = Dif_Fin(x, mult_res, rho)
    count += 1 + numvar

    # Se houver atualizacao na derivada,
    y = dL - dL_ant
    if y == 0.0

      # Salva argumentos para atualizacaodo DFP
      v = alpha*dir
      A = (v*v') ./ (v'*y)
      B = ( -G*y*(G*y)' ) ./ ( y'*G*y )

      # Atualiza a aproximacao da hessiana
      G = G + A + B

    end

    # Calculo da nova direcao
    norma = norm(dL)
    dir = (-G*dL)/norma

  end #for i
  return x,dL,count
end #function
