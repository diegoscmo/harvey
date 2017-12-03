function Wall_Search(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1},
                      dir::Array{Float64,1}, tol_int::Float64, minimo::Float64,
                      nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2},
                      K0::Array{Float64,2}, M0::Array{Float64,2},
                      SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
                      NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
                      dviz::Array{Float64,2}, raiof::Float64,passo0::Float64,L0::Float64, valor_0::Array{Float64,1},
                      caso::Int64, freq::Float64, alfa::Float64, beta::Float64, Ye::Float64, A::Float64)

    # Declaração de variáveis
    conta_line = 0      # Contador
    L1 = 0.0

    # Parâmetros
    p_zero = passo0       # Passo inicial
    p_min  = minimo       # Passo mínimo

    # Enquanto o paso for maior do que o mínimo
    alfa = p_zero
    while true

        # Define o novo x para verificação
        xn = x + alfa*dir

        # Bloqueia o x se passar
        @inbounds for j=1:nel
            xn[j] = min(1.0,max(xn[j],0.0))
        end

        # Cálcula a função na posição nova
        L1 = F_Obj(xn, rho, mult_res, 2, nel, ijk, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              caso, freq, alfa, beta, A, Ye)
        conta_line += 1

        if L1 < L0
            break
        else
            alfa = alfa/2.0
        end

        if alfa < p_min
            return p_min, conta_line, L1
        end

    end

    return alfa, conta_line, L1

end
