function Line_Search(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64,
   dir::Array{Float64,1}, xl::Array{Float64,1}, xu::Array{Float64,1}, tol_int::Float64, count::Int64,
    fem,filt,valor_zero)

    # Define um valor minimo de passo
    const minimo = 1E-12
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
    alpha = 1e10
    a_check = alpha
    for j=1:size(x,1)
        if dir[j]<0.0
            a_check = (xl[j]-x[j])/dir[j]
        end
        if dir[j]>0.0
            a_check = (xu[j]-x[j])/dir[j]
        end
        alpha = min(a_check, alpha)
    end

    # Calcula o valor do custo no ponto atual
    xf = Aplica_Filtro(x, filt)
    Ll = F_Lagrangiana(xf, mult_res, rho, fem, valor_zero)
    count += 1

    # Bracketing
    while true
        aa = delta
        xn = x + aa*dir

        # Calcula funcao lagrangiana nesta nova posicao
        xf = Aplica_Filtro(xn, filt)
        La = F_Lagrangiana(xf, mult_res, rho, fem, valor_zero)
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

        # Calcula funcao lagrangiana nesta nova posicao
        xf = Aplica_Filtro(xn, filt)
        Lu = F_Lagrangiana(xf, mult_res, rho, fem, valor_zero)
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

            # Sai se passar do bloqueio
            if au>alpha
                return alpha,count
            end

            xf = Aplica_Filtro(xn, filt)
            Lu = F_Lagrangiana(xf, mult_res, rho, fem, valor_zero)
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
