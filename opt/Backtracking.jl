function Backtracking(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, dL::Array{Float64,1},
   dir::Array{Float64,1}, xl::Array{Float64,1}, xu::Array{Float64,1}, tol_int::Float64, count::Int64,
    fem_v, fem_f, filt, step_min::Float64)

  # Relaxação do alfa
    const tau = 0.5

  # Relaxação da inclinação inicial
    const cc = 0.1

  # Define um valor minimo de passo
    const minimo = step_min

  # Calcula o valor do custo no ponto atual
    xf = Aplica_Filtro(x, filt)
    const f0 =  F_Lagrangiana(xf, mult_res, rho, fem_v, fem_f)
    count += 1

    # Normaliza a direção de busca, só para garantir...
    dir = dir/norm(dir)

  # Fator de comparação do método
    const direita = -cc*dot(dL,dir)

  # Verifica se não temos um alfa limite. Do contrário, utilizamos o máximo
    const alfa = 0.1

  # Bloqueio alpha
    a = +Inf
    @inbounds for j=1:size(x,1)
        if dir[j]<0
            a = (xl[j]-x[j])/dir[j]
        end
        if dir[j]>0
            a = (xu[j]-x[j])/dir[j]
        end
        alfa = min(a, alfa)
    end


  #  println("\n ****************** ")
  # Loop do Método
    for i=1:1000

        xn = x + alfa*dir

        xf = Aplica_Filtro(xn, filt)
        const fu =  F_Lagrangiana(xf, mult_res, rho, fem_v, fem_f)
        count += 1

  #      println(f0," ",fu," ",f0-fu," ",direita," ",alfa)

        if f0-fu < alfa*direita
            alfa = alfa * tau
            if alfa < minimo
                break
            end
        else
            break
        end

    end #i
  #  println("\n =================== ")

    return alfa,count

end # Armijo
