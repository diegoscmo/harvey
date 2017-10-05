function Fletcher_Reeves(x::Array{Float64,1},mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
   xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64, search::String)

  numvar = size(x,1)
  dL = zeros(Float64,numvar)

  # Calcula o gradiente de L
  dL = Dif_Fin(x, mult_res, rho)
  count += 1 + numvar

  # Calcula a norma para criterio de parada interno
  norma = norm(dL)

  # Direcao de minimizacao
  dir = -dL/norma

  for i=1:max_int

    # Verifica o criterio de parada interno
    if norma < tol_int
      break
    end

    # Bloqueia a direcao
    for j=1:numvar
      if dir[j] < 0.0 && x[j] <= xl[j]
        dir[j] = 0.0
      end
      if dir[j] > 0.0 && x[j] >= xu[j]
        dir[j] = 0.0
      end
    end

    # Search nesta direcao
    alpha,count = LS_Methods(x, mult_res, rho, dir, xl, xu, dL, tol_int, count, search)

    # incrementa a estimativa do ponto
    x = x + alpha*dir

    # Salva os calculos da iteracao anterior
    dL_ant = copy(dL)
    norma_ant = copy(norma)
    dir_ant = copy(dir)

    # Regaculta o gradiente e a norma
    dL = Dif_Fin(x, mult_res, rho)
    count += 1 + numvar

    norma = norm(dL)

    # Atualiza betha para o criterio de fletcher reeves
    betha = ( norma / norma_ant )^2

    # Calculo da nova direcao
    dir = ( -dL + betha * dir_ant )/norma

  end #for i
  return x,dL,count
end #function
