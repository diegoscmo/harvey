function Equal_Search(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1},
                      dir::Array{Float64,1}, tol_int::Float64, minimo::Float64,
                      nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2},
                      K0::Array{Float64,2},
                      #K0::StaticArrays.SArray{Tuple{8,8},Float64,2,64},
                      M0::Array{Float64,2},
                      #M0::StaticArrays.SArray{Tuple{8,8},Float64,2,64},
                      SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
                      NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
                      dviz::Array{Float64,2}, raiof::Float64)

    # E o passo inicial
    delta = 0.1

    # Aoca as variaveis da rotina
    alfa = 0.0
    al = 0.0
    aa = 0.0
    La = 0.0
    au = 0.0
    al = 0.0
    Ll = 0.0

    # Bloqueio do alpha
    alpha = 1e10
    a_check = alpha
    for j=1:size(x,1)
        if dir[j]<0.0
            a_check = (0.0-x[j])/dir[j]
        end
        if dir[j]>0.0
            a_check = (1.0-x[j])/dir[j]
        end
        alpha = min(a_check, alpha)
    end

    # Calcula o valor do custo no ponto atual
    Ll = F_Obj(x, rho, mult_res, 2, nel, ijk, ID, K0, M0, SP, vmin,
                                 F, NX, NY, vizi, nviz, dviz, raiof)
    conta_line = 1

    # Bracketing
    while true
        aa = delta
        xn = x + aa*dir

        # Calcula funcao lagrangiana nesta nova posicao
        La = F_Obj(xn, rho, mult_res, 2, nel,  ijk, ID, K0, M0, SP, vmin,
                                       F, NX, NY, vizi, nviz, dviz, raiof)
        conta_line += 1

        # Se nao minimizar, diminui o passo delta
        if La>Ll
            delta = delta/10.
            if delta < minimo
                return minimo,conta_line
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
            return alpha,conta_line
        end

        # Calcula funcao lagrangiana nesta nova posicao
        Lu = F_Obj(xn, rho, mult_res, 2, nel, ijk, ID, K0, M0, SP, vmin,
                                      F, NX, NY, vizi, nviz, dviz, raiof)
        conta_line += 1

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
                return alpha,conta_line
            end

            Lu = F_Obj(xn, rho, mult_res, 2, nel, ijk, ID, K0, M0, SP, vmin,
                                          F, NX, NY, vizi, nviz, dviz, raiof)
            conta_line += 1

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

    return alpha, conta_line

end
