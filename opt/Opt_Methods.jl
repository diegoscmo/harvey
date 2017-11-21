# Carrega os metodos de busca
include("Steepest.jl")
include("Fletcher_Reeves.jl")
include("DFP.jl")
include("BFGS.jl")

# Carrega os metodos de procura em linha
include("Equal_Search.jl")
include("Golden_Search.jl")
include("Backtracking.jl")

#
# Seleção do método de descida
#
function Descent(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1},rho::Float64,xl::Array{Float64,1},
                  xu::Array{Float64,1}, max_int::Int64, tol_int::Float64, count::Int64,
                  fem_v, fem_f, filt, descent::String, lsearch::String, step_min::Float64)

  # Search nesta direcao
  if descent == "Steep"
    x, dL, count, n_int = Steepest(x, valor_res, mult_res, rho, xl, xu,
                            max_int, tol_int, count, fem_v, fem_f, filt, lsearch, step_min)
  elseif descent == "FR"
    x, dL, count, n_int = Fletcher_Reeves(x, valor_res, mult_res, rho, xl, xu,
                            max_int, tol_int, count, fem_v, fem_f, filt, lsearch, step_min)
  elseif descent == "DFP"
    x, dL, count, n_int = DFP(x, valor_res, mult_res, rho, xl, xu,
                            max_int, tol_int, count, fem_v, fem_f, filt, lsearch, step_min)
  elseif descent == "BFGS"
    x, dL, count, n_int = BFGS(x, valor_res, mult_res, rho, xl, xu,
                            max_int, tol_int, count, fem_v, fem_f, filt, lsearch, step_min)
  else
    error("Erro na declaracao do método de Descent")
  end

  return x,dL,count,n_int

end

#
# Seleção do método de busca em linha
#
function LineSearch(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, dL::Array{Float64,1},
   dir::Array{Float64,1}, xl::Array{Float64,1}, xu::Array{Float64,1}, tol_int::Float64, count::Int64,
    fem_v, fem_f, filt, lsearch::String, step_min::Float64)

  # Search nesta direcao
  if lsearch == "Golden"
    alpha,count = Golden_Search(x, mult_res, rho, dir, xl, xu, tol_int, count, fem_v, fem_f, filt, step_min)
elseif lsearch == "Equal"
    alpha,count = Equal_Search(x, mult_res, rho, dir, xl, xu, tol_int, count, fem_v, fem_f, filt, step_min)
elseif lsearch == "Back"
    alpha,count = Backtracking(x, mult_res, rho, dL, dir, xl, xu, tol_int, count, fem_v, fem_f, filt, step_min)
  else
    error("Erro na declaracao do método de Line Search")
  end

  return alpha,count

end
