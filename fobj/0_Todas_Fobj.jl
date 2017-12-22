################################################################################
#####                          Funções Objetivo                           ######
################################################################################

# Carrega as rotinas de cada função objetivo
include("1_Estatico.jl")
include("2_Dinamico.jl")
include("3_Din_R-Est.jl")
include("4_Din_A-Est.jl")
include("5_Potencia.jl")
include("6_Pot_R-Est.jl")
include("7_Pot_A-Est.jl")
include("8_Pot_A-Est_R-R.jl")
include("11_Estatico_Stress.jl")

#
# Função indice das funções objetivo, entra com o caso e escolhe a função.
#
function F_Obj(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
    nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
    M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1},
    NX::Int64, NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
    dviz::Array{Float64,2}, raiof::Float64, Y0::Array{Float64,1},
    caso::Int64, freq::Float64, alfa::Float64, beta::Float64, A::Float64, Ye::Float64, CBA, filtra::Bool=true )

    # Caso Estático
    if caso == 1
        return F_Est(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    # Dinâmico
    elseif caso == 2
         return F_Din(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
           NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    # Dinâmico com restrição estática
    elseif caso == 3
        return F_Din_R_Est(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    # Dinâmico + Estático
    elseif caso == 4
        return F_Din_A_Est(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    # Potência
    elseif caso == 5
        return F_Pot(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    # Potência com restrição estática
    elseif caso == 6
        return F_Pot_R_Est(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    # Potência + Estática
    elseif caso == 7
        return F_Pot_A_Est(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    # Potência + Estática com restrição R
    elseif caso == 8
        return F_Pot_A_Est_R_R(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye)

    elseif caso == 11
        return F_Est_S(x, rho, mult_res, tipo, nnos, nel, ijk, ID, K0, M0, SP, vmin, F,
            NX, NY, vizi, nviz, dviz, raiof, Y0, caso, freq, alfa, beta, A, Ye, CBA)

    else
    # Se chegou até aqui é pq deu problema
    error("ATENCAO AO VALOR DA VARIAVEL CASO")

end #fim da funcao

end
