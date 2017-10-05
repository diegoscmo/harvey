function F_Lagrangiana(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64)

    # Obtem o valor do objetivo e das restricoes
    valor_fun,valor_res = F_Obj(x)

    # Dimensoes do problema
    numvar = size(x,1)
    numres = size(valor_res,1)

    # Define o valor de L(x) no ponto x
    L = 0.0
    for j=1:numres
        L = L + max(0.0, mult_res[j]+valor_res[j]/rho)^2.0
    end

    L = valor_fun + rho*L/2.0

    return L

end
