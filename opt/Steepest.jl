function Steepest(x::Array{Float64,1},mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
   xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64, search::String)

  numvar = size(x,1)
  dL = zeros(Float64,numvar)

  for i=1:max_int

    # Calcula o gradiente de L
    dL = Dif_Fin(x, mult_res, rho)
    count += 1 + numvar

    # Calcula a norma para criterio de parada interno
    norma = norm(dL)
    if norma < tol_int
      break
    end

    # Direcao de minimizacao
    dir = -dL/norma

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

  end #for i
  return x,dL,count
end #function
