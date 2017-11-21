function Dif_Fin(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, count::Int64, fem, valor_zero)
    
    # Step para as diferenças finitas
    h = 3.0*sqrt(eps())

    # Só para não precisar pegar de fora...
    numvar = size(x,1)

    # Declara
    dL = zeros(Float64,numvar)

    # Diferenças finitas à frente

        # Valor da funcao Lagrangiana no ponto x
        L0 = F_Lagrangiana(x, mult_res, rho, fem, valor_zero)

        # Gradiente em cada direção
        dL = zeros(size(x,1))
        for i=1:numvar
            b = x[i]
            x[i] = b + h
            L = F_Lagrangiana(x, mult_res, rho, fem, valor_zero)
            dL[i] = (L - L0)/h
            x[i] = b
        end
        count = count + 1 + numvar

    return dL,count
end
