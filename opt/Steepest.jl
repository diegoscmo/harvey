function Steepest(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1},
                    max_int::Int64, tol_int::Float64, nel::Int64, ijk::Array{Int64,2},
                    ID::Array{Int64,2},
                    K0::Array{Float64,2},
                    M0::Array{Float64,2},
                    SP::Float64, vmin::Float64,
                    F::Array{Float64,1}, NX::Int64, NY::Int64, vizi::Array{Int64,2},
                    nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64 ,dts::String, valor_0::Array{Float64,1})

    # Parâmetros, passo mínimo e máximo de iterações estagnado
    minimo    = 1E-4
    max_break = 300

    # Inicializa os vetores para derivadas
    dL = zeros(Float64,nel)

    # Inicializa o contador para break estagnado e de avaliações da Fobj
    breaker = 0
    i_int   = 0
    norma   = 0.0

    passo0 = 10.0

    L0 = F_Obj(x, rho, mult_res, 2, nel, ijk, ID, K0, M0, SP, vmin,
                                 F, NX, NY, vizi, nviz, dviz, raiof, valor_0)
    contaev = 2

    tmp = @elapsed for i=1:max_int

        # Calcula sensibilidade
        dL = F_Obj(x, rho, mult_res, 3, nel, ijk, ID, K0, M0, SP, vmin,
                                     F, NX, NY, vizi, nviz, dviz, raiof, valor_0)
        contaev += 1

        # Bloqueia a direcao e zera o gradiente se bloqueado
        @inbounds for j=1:nel
            if dL[j] > 0.0 && x[j] <= 0.0
                dL[j]  = 0.0
            elseif dL[j] < 0.0 && x[j] >= 1.0
                dL[j]  = 0.0
            end
        end

        # Direcao de minimizacao
        norma = norm(dL)
        dir = -dL/norma

        if norma < tol_int && i_int > 1
            break
        end

        # Search nesta direcao
        (alpha, conta_line, L0) = Wall_Search(x, rho, mult_res, dir, tol_int,
                                minimo, nel, ijk, ID, K0, M0, SP, vmin, F, NX, NY, vizi, nviz, dviz, raiof, passo0, L0, valor_0)
        contaev += conta_line
        i_int += 1
        # incrementa a estimativa do ponto
        x = x + alpha*dir

        passo0 = alpha*2.0
        if i%5 == 0
            @printf(".")
        end

        # Bloqueia o x se passar
        @inbounds for j=1:nel
            x[j] = min(1.0,max(x[j],0.0))
        end

        # Critério adicional de saída
        if alpha <= minimo
            breaker += 1
            if breaker >= max_break
                break
            end
        end

    end #for i
    @printf("\n")
    # Da o display do laço interno
    Imprime_Int(i_int, contaev, norma, dts, tmp)

    return x, dL
end #function
