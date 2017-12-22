################################################################################
#####                       Wall Line Search                              ######
################################################################################
#
# Busca em linha
#
function Wall_Search(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1},
                      dir::Array{Float64,1}, tol_int::Float64, minimo::Float64,
                      nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2},
                      K0::Array{Float64,2}, M0::Array{Float64,2},
                      SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
                      NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
                      dviz::Array{Float64,2}, raiof::Float64, passo0::Float64,L0::Float64, valor_0::Array{Float64,1},
                      caso::Int64, freq::Float64, alfa::Float64, beta::Float64, A::Float64, Ye::Float64, CBA )

    # Declaração de variáveis
    conta_line = 0      # Contador
    L1 = 0.0            # Valor da Função Objetivo LA
    alpha = passo0      # Passo inicial

    # Enquanto o paso for maior do que o mínimo
    while true

        # Define o novo x para verificação
        xn = x + alpha*dir

        # Bloqueia o x se passar
        @inbounds for j=1:nel
            xn[j] = min(1.0,max(xn[j],0.0))
        end #for j

        # Cálcula a função na posição nova
        L1, = F_Obj(xn, rho, mult_res, 2, nnos, nel, ijk, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              caso, freq, alfa, beta, A, Ye, CBA)
        conta_line += 1

        # Se diminuir sai do laço, do contrario diminui o passo
        if L1 < L0
            break
        else
            alpha = alpha/2.0
        end # if L1

        # Se passar do mínimo, retorna-o
        if alpha < minimo
            return minimo, conta_line, L1
        end # if alpha

    end #while

    return alpha, conta_line, L1

end
