function BFGS(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
    xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64,
    fem_v, fem_f, filt, lsearch::String)

    # Lê número de variáveis
    numvar = size(x,1)

    # Inicializa os vetores para derivadas
    dL_ant = zeros(Float64,numvar)
    dL = zeros(Float64,numvar)

    # Inicializa o contador de iterações internas
    n_int = 0

    # inicializa aproximacao da hessiana
    G = eye(Float64,numvar)

    for i=1:max_int

        # Armazena a derivada anterior
        dL_ant = copy(dL)

        # Se houver atualizacao na derivada, (2a iteração adiante)
        y = norm(dL - dL_ant)
        if y > tol_int && i > 1

            # Salva argumentos para atualizacaodo BFGS
            v = alpha*dir
            A = (v*v') ./ (v'*y)
            C = 1.0 + ( ( y'*(G*y) ) ./ ( v'*y ) )
            D = (v*y'*G + G*y*v') ./ (v'*y)

            # Atualiza a aproximacao da hessiana
            G = G + C*A - D

        end # if y

        # Aplica filtro de densidade 65
        xf = Aplica_Filtro(x, filt)

        # Atualiza as matrizes de elementos finitos filtrados
        fun, res = F_Obj(xf, fem_v, fem_f)#, valor_zero)
        count += 1

        # Calcula o gradiente de L
        dL = Sensibilidade(xf, valor_res, mult_res, rho, fem_v, fem_f, filt)

        # Calcula a norma para criterio de parada interno
        norma = norm(dL)

        # Direcao de minimizacao
        dir = (-G*dL)/norma

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

        # Line search nesta direcao
        alpha,count = LineSearch(x, mult_res, rho, dL, dir, xl, xu, tol_int, count, fem_v, fem_f, filt, lsearch)

        # incrementa a estimativa do ponto
        x = x + alpha*dir

        # E o número de iterações internas
        n_int += 1

    end #for i

    return x, dL, count, n_int
end #function
