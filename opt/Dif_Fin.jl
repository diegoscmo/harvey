# Diferenças Finitas para validação!

function Sensibilidade(x::Array{Float64,1}, valor_res::Array{Float64,1}, mult_res::Array{Float64,1},
                       rho::Float64, fem_v, fem_f, filt)

    # Step para as diferenças finitas
    h = 3.0*sqrt(eps())

    # Só para não precisar pegar de fora...
    numvar = size(x,1)

    # Declara
    dL = zeros(Float64,numvar)

    # Diferenças finitas à frente

        # Valor da funcao Lagrangiana no ponto x
        L0 = F_Lagrangiana(x, mult_res, rho, fem_v, fem_f)

        # Gradiente em cada direção
        dL = zeros(size(x,1))
        for i=1:numvar
            b = x[i]
            x[i] = b + h
            L = F_Lagrangiana(x, mult_res, rho, fem_v, fem_f)
            #dL[i] = (L - L0)/h

            x[i] = b - h
            L0 = F_Lagrangiana(x, mult_res, rho, fem_v, fem_f)
            dL[i] = (L - L0)/(2.0*h)

            x[i] = b
        end

    dL = Derivada_Filtro(dL, filt)

    return dL
end
