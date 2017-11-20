function Fletcher_Reeves(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
    xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64,
    fem_v, fem_f, filt, lsearch::String)

    # Lê número de variáveis
    numvar = size(x,1)

    # Inicializa os vetores para derivadas
    dL = zeros(Float64,numvar)
    norma = zeros(Float64,numvar)
    dir = (Float64,numvar)

    # Inicializa o contador de iterações internas
    n_int = 0

    for i=1:max_int

        # Salva os calculos da iteracao anterior
        dL_ant = copy(dL)
        norma_ant = copy(norma)
        dir_ant = copy(dir)

        # Aplica filtro de densidade 65
        xf = Aplica_Filtro(x, filt)

        # Atualiza as matrizes de elementos finitos filtrados
        fun, res = F_Obj(xf, fem_v, fem_f)#, valor_zero)
        count += 1

        # Calcula o gradiente de L
        dL = Sensibilidade(xf, valor_res, mult_res, rho, fem_v, fem_f, filt)

        # Calcula a nova norma
        norma = norm(dL)

        # Direcao de minimizacao
        if i == 1
            dir = -dL/norma
        else
            # Calculo da nova direcao
            dir = ( -dL + betha * dir_ant )/norma
        end

        # Bloqueia a direcao, posição e a derivada se necessário
        @inbounds for j=1:numvar
            if dir[j] < 0.0 && x[j] <= xl[j]
                dir[j] = 0.0
                dL[j]  = 0.0
                x[j]   = xl[j]
            end
            if dir[j] > 0.0 && x[j] >= xu[j]
                dir[j] = 0.0
                dL[j]  = 0.0
                x[j]   = xu[j]
            end
        end

        # Nova norma, agora com bloqueio
        norma = norm(dL)

        # Verifica o criterio de parada interno
        if norma < tol_int
            break
        end

        # Atualiza betha para o criterio de fletcher reeves
        betha = ( norma / norma_ant )^2

        # Line search nesta direcao
        alpha,count = LineSearch(x, mult_res, rho, dir, xl, xu, tol_int, count, fem_v, fem_f, filt, lsearch)

        # incrementa a estimativa do ponto
        x = x + alpha*dir

        # E o número de iterações internas
        n_int += 1

    end #for i

    return x, dL, count, n_int
end #function
