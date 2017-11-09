function Steepest(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
   xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64,
   fem,filt,valor_zero)

  numvar = size(x,1)
  dL = zeros(Float64,numvar)
  #contaestag = 0

  for i=1:max_int

    # Aplica filtro de densidade 65
    #x = Aplica_Filtro(x, filt) #FIXME, não funfa quando ligado!

    # Calcula o gradiente de L
    dL,count = Sensibilidade(x, valor_res, mult_res, rho, count, fem, filt, valor_zero)
    #dL,count = Dif_Fin(x, mult_res, rho, count, fem, valor_zero)

    # Calcula a norma para criterio de parada interno
    norma = norm(dL)
    if norma < tol_int
      break
    end

    # Direcao de minimizacao
    dir = -dL/norma

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

    # Search nesta direcao
    alpha,count = Line_Search(x, mult_res, rho, dir, xl, xu, tol_int, count, fem, valor_zero)

    # incrementa a estimativa do ponto
    x = x + alpha*dir

  end #for i

  return x,dL,count,fem
end #function
