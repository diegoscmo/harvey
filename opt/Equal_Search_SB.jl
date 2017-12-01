function Equal_Search_SB(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64,
   dir::Array{Float64,1}, xl::Array{Float64,1}, xu::Array{Float64,1}, tol_int::Float64, count::Int64,
    fem_v, fem_f, filt, step_min::Float64)

    # Define um valor minimo de passo
    const minimo = step_min
    # E o passo inicial
    const delta = 0.1

    # Aoca as variaveis da rotina
    const alfa = 0.0
    const al = 0.0
    const aa = 0.0
    const La = 0.0
    const au = 0.0
    const al = 0.0
    const Ll = 0.0

    # Bloqueio do alpha
    #alpha = 1e10
    #a_check = alpha
    #for j=1:size(x,1)
#           if dir[j]<0.0
#            a_check = (xl[j]-x[j])/dir[j]
#        end
#        if dir[j]>0.0
#            a_check = (xu[j]-x[j])/dir[j]
#        end
#        alpha = min(a_check, alpha)
#    end

    # O passo máximo não tem bloqueio
    alpha = delta

    # Verifica o número de variáveis para zerar quando chegar nas paredes
    numvar = size(x,1)

    # Calcula o valor do custo no ponto atual
    xf = Aplica_Filtro(x, filt)
    Ll = F_Lagrangiana(xf, mult_res, rho, fem_v, fem_f)#, valor_zero)
    count += 1

    # Bracketing
    while true
        aa = delta
        xn = x + aa*dir

        # Bloqueia as variáveis que passaram
        for j=1:numvar
            if x[j] >= xu[j]
                x[j] = xu[j]
            end
            if x[j] <= xl[j]
                x[j] = xl[j]
            end
        end


        # Calcula funcao lagrangiana nesta nova posicao
        xf = Aplica_Filtro(xn, filt)
        La = F_Lagrangiana(xf, mult_res, rho, fem_v, fem_f)#, valor_zero)
        count += 1

        # Se nao minimizar, diminui o passo delta
        if La>Ll
            delta = delta/10.
            if delta < minimo
                return minimo,count
            end
        # Se minimizar, definiu o bracketing superior
        else
            break
        end
    end #while

    # Refina o Bracketing
    while true
        au = aa + delta
        xn = x + au*dir

        # Se passar do limite de alpha, retorna alpha
        if au>alpha
            return alpha,count
        end

        # Bloqueia as variáveis que passaram
        for j=1:numvar
            if x[j] >= xu[j]
                x[j] = xu[j]
            end
            if x[j] <= xl[j]
                x[j] = xl[j]
            end
        end

        # Calcula funcao lagrangiana nesta nova posicao
        xf = Aplica_Filtro(xn, filt)
        Lu = F_Lagrangiana(xf, mult_res, rho, fem_v, fem_f)#, valor_zero)
        count += 1

        # Aproxima o search pelo brackeing inferior
        if La>Lu
            al = aa
            aa = au
            Ll = La
            La = Lu
        else
            break
        end
    end

    while true
        # Verifica o criterio de parada
        if (au-al) <= tol_int
            break
        end

        delta = delta/10.0
        aa = al
        La = Ll

        while true
            au = aa + delta
            xn = x + au*dir

            # Bloqueia as variáveis que passaram
            for j=1:numvar
                if x[j] >= xu[j]
                    x[j] = xu[j]
                end
                if x[j] <= xl[j]
                    x[j] = xl[j]
                end
            end

            # Sai se passar do bloqueio
            if au>alpha
                return alpha,count
            end

            xf = Aplica_Filtro(xn, filt)
            Lu = F_Lagrangiana(xf, mult_res, rho, fem_v, fem_f)#, valor_zero)
            count += 1
            if La>Lu
                al = aa
                aa = au
                Ll = La
                La = Lu
            else
                break
            end
        end
    end

    alpha = (au+al)/2.0

    return alpha,count

end
