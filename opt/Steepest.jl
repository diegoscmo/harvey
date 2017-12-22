################################################################################
#####                       Steepest Descent                              ######
################################################################################
#
# Método de busca baseado no gradiente
#
function Steepest(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, max_int::Int64,
             tol_int::Float64, nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2},
    K0::Array{Float64,2}, M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1},
       NX::Int64, NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
             raiof::Float64 ,dts::String, valor_0::Array{Float64,1}, caso::Int64, freq::Float64,
                   alfa::Float64, beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3})

    # Parâmetros, passo mínimo e máximo de iterações estagnado
    minimo    = 1E-5*(nel^0.5)
    max_break = Int(round(0.7*max_int))

    # Inicializa os vetores para derivadas
    dL = Array{Float64}(uninitialized,nel)

    # Inicializa o contador para break estagnado e de avaliações da Fobj
    breaker = 0
    i_int   = 0
    contaev = 2
    norma   = 0.0

    # Passo Inicial
    passo0 = 0.25*(nel^0.5)

    # Lagrangiano no ponto de partida
    L0, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, ID, K0, M0, SP, vmin,
                                 F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                 caso, freq, alfa, beta, A, Ye, CBA)

    # Laço interno + cronômetro
    tmp = @elapsed for i=1:max_int

        # Calcula sensibilidade
        dL, = F_Obj(x, rho, mult_res, 3, nnos, nel, ijk, ID, K0, M0, SP, vmin,
                                     F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                     caso, freq, alfa, beta, A, Ye, CBA)
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

        # Verifica saída pela norma
        if norma < tol_int && i_int > 1
            break
        end

        # Search na direção do Steepest
        (alpha, conta_line, L0) = Wall_Search(x, rho, mult_res, dir, tol_int,
                                minimo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F, NX, NY, vizi,
                                nviz, dviz, raiof, passo0, L0, valor_0,
                                caso, freq, alfa, beta, A, Ye, CBA)
        contaev += conta_line
        i_int += 1

        # Incrementa a estimativa do ponto
        x = x + alpha*dir

        # Define o passo máximo da próxima iteração
        passo0 = alpha*2.0

        # Mini-exibição a cada 5 iterações
        if i%5 == 0
            print(".")
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
    println("")
    # Da o display do laço interno
    Imprime_Int(i_int, contaev, norma, dts, tmp)

    return x, dL
end #function
