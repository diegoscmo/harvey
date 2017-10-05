
# Carrega os metodos de procura em linha
include("Line_Search.jl")
include("Golden_Search.jl")
include("LS_Backtracing.jl")

function LS_Methods(x, mult_res, rho, dir, xl, xu, dL, tol_int, count, search)

  # Search nesta direcao
  if search == "Golden"
    alpha,count = Golden_Search(x, mult_res, rho, dir, xl, xu, tol_int, count)
  elseif search == "Line"
    alpha,count = Line_Search(x, mult_res, rho, dir, xl, xu, tol_int, count)
  elseif search == "Back"
    alpha,count = LS_Backtracing(x, mult_res, rho, dir, xl, xu, dL, tol_int, count)
  else
    error("Erro na declaracao do metodo de Line Search")
  end

  return alpha,count

end
