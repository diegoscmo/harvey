function Golden_Search(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64,
   dir::Array{Float64,1}, xl::Array{Float64,1}, xu::Array{Float64,1}, tol_int::Float64, count::Int64,
    fem_v, fem_f, filt)

    # Define um valor minimo de passo
    const minimo = 1E-12
    # o passo inicial,
    const delta = 0.1
    # e o golden ratio
    const GR = (sqrt(5.)-1)/2.
    const GR2 = GR^2

    # Bloqueio do alpha
    alpha = 1.0
    a_check = alpha
    for j=1:size(x,1)
        if dir[j] < 0.0
            a_check = (xl[j]-x[j])/dir[j]
        end
        if dir[j] > 0.0
            a_check = (xu[j]-x[j])/dir[j]
        end
        alpha = min(a_check, alpha)
    end

    # Calcula o valor do custo no ponto atual
    a = x
    af = Aplica_Filtro(a, filt)
    L_a = F_Lagrangiana(af, mult_res, rho, fem_v, fem_f)
    count += 1

    # Bracketing superior
    while true
        b = x + delta*dir

        # Calcula funcao lagrangiana nesta nova posicao
        bf = Aplica_Filtro(b, filt)
        L_b = F_Lagrangiana(bf, mult_res, rho, fem_v, fem_f)
        count += 1

        # Se nao minimizar, diminui o passo delta
        if L_a < L_b
            delta = delta/10.
            if delta < minimo
                return minimo,count
            end
        # Se minimizar, definiu o bracketing superior
        else
            break
        end
    end #while


    # Define os pontos intermediarios do Golden Search
    L = delta*dir
    b = x + L
    l1 = a + GR2*L
    l2 = a + GR *L
    l1f = Aplica_Filtro(l1, filt)
    fl1 = F_Lagrangiana(l1f, mult_res, rho, fem_v, fem_f)
    l2f = Aplica_Filtro(l2, filt)
    fl2 = F_Lagrangiana(l2f, mult_res, rho, fem_v, fem_f)
    count += 2

    # Enquanto for menor que a tolerancia,
    while norm(L) > tol_int

      # Verifica qual secao remover
      if fl1 > fl2
        # Remove a parte inferior
        a = copy(l1)
        l1 = copy(l2)
        fl1 = copy(fl2)
        L = (b - a)
        l2 = a + GR*L
        l2f = Aplica_Filtro(l2, filt)
        fl2 = F_Lagrangiana(l2f, mult_res, rho, fem_v, fem_f)
        count += 1
      else
        # Remove a parte superior
        b = copy(l2)
        l2 = copy(l1)
        fl2 = copy(fl1)
        L = (b - a)
        l1  = a + GR2*L
        l1f = Aplica_Filtro(l1, filt)
        fl1 = F_Lagrangiana(l1f, mult_res, rho, fem_v, fem_f)
        count += 1
      end #if

      # Recalcula delta com o maior dir
      (dir,ind) = findmax(dir)
      delta = L[ind]/dir
      #println(delta)

      # Se delta ficar muito pequeno, define como o minimo
      if delta < minimo
        return minimo,count
      end
    end #while

    return delta,count

end
