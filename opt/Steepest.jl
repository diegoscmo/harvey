################################################################################
#####                       Steepest Descent                              ######
################################################################################
#
# Método de busca baseado no gradiente
#
function Steepest(T::Int64, x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1},
                max_int::Int64, tol_int::Float64, nnos::Int64, nel::Int64, ijk::Array{Int64,2},
                coord::Array{Float64,2}, ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64, NY::Int64, vizi::Array{Int64,2},
                nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64 ,dts::String,
                valor_0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64, beta::Float64,
                A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64,
                trava_els::Array{Int64,1}, nos_viz, i_ext::Int64, P::Float64, q::Float64, F_Obj)

    # Na primeira iteração externa, dobra o número de internas
    if i_ext == 1
        max_int = 2*max_int
    end

    # Parâmetros, passo mínimo e máximo de iterações estagnado
    minimo    = 1E-5*(nel^0.5)
    max_break = Int(round(0.7*max_int))

    # Inicializa os vetores para derivadas
    dL = Array{Float64}(undef,nel)

    # Inicializa o contador para break estagnado e de avaliações da Fobj
    breaker = 0
    i_int   = 0
    contaev = 2
    norma   = 0.0

    # Passo Inicial
    passo0 = 0.10*(nel^0.5)

    # Laço interno + cronômetro
    @showprogress "  iteracao $i_ext... " for i=1:max_int

        # Calcula sensibilidade #F_Obj_DF = diferenças finitas + impressão
        dL,L0, = F_Obj[T](x, rho, mult_res, 3, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                     F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                     Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)
        contaev += 1

        # Bloqueia a direcao
        for j=1:nel
            if (dL[j] > 0.0 && x[j] <= 0.0) || (dL[j] < 0.0 && x[j] >= 1.0)
                dL[j] = 0.0
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
        (alpha, conta_line) = Wall_Search(T, x, rho, mult_res, dir, tol_int, minimo, nnos, nel,
                                ijk, coord, ID, K0, M0, SP, vmin, F, NX, NY, vizi, nviz, dviz,
                                raiof, passo0, L0, valor_0, Sy, freq, alfa, beta, A, Ye, CBA,
                                QP, csi, dmax, trava_els, nos_viz, dts, P, q, F_Obj)
        contaev += conta_line
        i_int += 1

        # Incrementa a estimativa do ponto
        x = x + alpha*dir

        # Define o passo máximo da próxima iteração
        passo0 = alpha*2.0

        # Bloqueia o x se passar
        for j=1:nel
            x[j] = min(1.0,max(x[j],0.0))
        end

        # Trava elementos
        x = Trava_Els(x,trava_els)

        # Critério adicional de saída
        if alpha <= minimo
            breaker += 1
            if breaker >= max_break
                break
            end
        end

    end #for i

    # Da o display do laço interno
    Imprime_Int(i_int, contaev, norma)

    return x, dL
end #function
