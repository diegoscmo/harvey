function Dif_Fin(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64)

    # Valor da funcao Lagrangiana no ponto x
    #L0 = F_Lagrangiana(x, mult_res, rho)
    y = copy(x)
    # Step para as diferenças finitas
    h = 3*sqrt(eps())
    # Gradiente em cada direção
    dL = zeros(size(x,1))
    for i=1:numvar
        b = y[i]
        y[i] = b + h
        L = F_Lagrangiana(y, mult_res, rho)
        y[i] = b - h
        L0 = F_Lagrangiana(y, mult_res, rho)
        dL[i] = (L - L0)/(2*h)
        y[i] = b
    end
    return dL
end
