function LS_Backtracing(x, mult_res, rho, dir, xl, xu, dL, tol_int, count)

  # Relaxação do alfa
    const tau = 0.5

  # Relaxação da inclinação inicial
    const cc = 0.1

  # Define um valor minimo de passo
    const minimo = 1E-12

  # Calcula o valor do custo no ponto atual
    const f0 =  F_Lagrangiana(x, mult_res, rho)
    count += 1

    # Normaliza a direção de busca, só para garantir...
    dir = dir/norm(dir)

  # Fator de comparação do método
    const direita = -cc*dot(dL,dir)

  # Verifica se não temos um alfa limite. Do contrário, utilizamos o máximo
    const alfa = 10.0

  # Bloqueio alpha
    a = +Inf
    for j=1:size(x,1)
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

        xf = x + alfa*dir
        fu = F_Lagrangiana(xf,mult_res, rho)
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
