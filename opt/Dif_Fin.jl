function Dif_Fin(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, count::Int64, tipo=1)
    #
    # tipos, 1 = à frente, 2 = central
    #

    # Step para as diferenças finitas
    h = 3.0*sqrt(eps())

    # Só para não precisar pegar de fora...
    numvar = size(x,1)

    # Declara
    dL = zeros(Float64,numvar)

    # Diferenças finitas à frente
    if tipo == 1
        # Valor da funcao Lagrangiana no ponto x
        L0 = F_Lagrangiana(x, mult_res, rho)

        # Gradiente em cada direção
        dL = zeros(size(x,1))
        for i=1:numvar
            b = x[i]
            x[i] = b + h
            L = F_Lagrangiana(x, mult_res, rho)
            dL[i] = (L - L0)/h
            x[i] = b
        end
        count = count + 1 + numvar


    elseif tipo == 2
        # Gradiente em cada direção
        for i=1:numvar
            b = x[i]
            x[i] = b + h
            L   = F_Lagrangiana(x, mult_res, rho)
            x[i] = b - h
            L0 = F_Lagrangiana(x, mult_res, rho)
            dL[i] = (L - L0)/(2.0*h)
            x[i] = b
        end
        count = count + 2*numvar
    end
    return dL,count
end
