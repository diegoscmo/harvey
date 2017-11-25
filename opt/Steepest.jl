function Steepest(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
                  xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64,
                  fem_v, fem_f, filt, lsearch::String, step_min::Float64) #,valor_zero)

    # Lê número de variáveis
    numvar = size(x,1)

    # Inicializa os vetores para derivadas
    dL = zeros(Float64,numvar)

    # Inicializa o contador de iterações internas
    n_int = 0

    # Critério adicional de saída, vezes com passo mínimo
    minimo = step_min
    breaker = 0
    max_break = fem_f.nbreaker

    for i=1:max_int

        # Aplica filtro de densidade 65
        xf = Aplica_Filtro(x, filt)

        # Atualiza as matrizes de elementos finitos filtrados
        fun, res = F_Obj(xf, fem_v, fem_f)#, valor_zero)
        count += 1

        # Calcula o gradiente de L
        dL = Sensibilidade(xf, valor_res, mult_res, rho, fem_v, fem_f, filt)#, valor_zero)

        # Bloqueia a direcao e zera o gradiente se bloqueado
        @inbounds for j=1:numvar
            if dL[j] > 0.0 && x[j] <= (xl[j]+minimo)
                dL[j]  = 0.0
                x[j]   = xl[j]
            end
            if dL[j] < 0.0 && x[j] >= (xu[j]-minimo)
                dL[j]  = 0.0
                x[j]   = xu[j]
            end
        end

        # Direcao de minimizacao
        norma = norm(dL)
        dir = -dL/norma

        if norma < tol_int && i > 1
            break
        end

        # Search nesta direcao
        alpha,count = LineSearch(x, mult_res, rho, dL, dir, xl, xu, tol_int, count, fem_v, fem_f, filt, lsearch, step_min)

        # incrementa a estimativa do ponto
        x = x + alpha*dir

        # Número de passos internos
        n_int += 1

        # Extrai derivada para verificação
        if false
            fmesh = string("sens_pot_loo.txt")
            if isfile(fmesh)
                rm(fmesh)
            end
            saida  = open(fmesh,"a")
            for z = 1:numvar
                println(saida,dL[z])
            end
            close(saida)
            error("\n Derivada em ",fmesh)
        end

        # Critério adicional de saída
        if alpha <= minimo
            breaker += 1
            if breaker >= max_break
                break
            end
        end
        @printf(".")

    end #for i
    @printf("\n")



    return x,dL,count,n_int
end #function
