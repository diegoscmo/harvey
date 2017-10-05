
# Carrega os metodos baseados em gradiente
include("Steepest.jl")
include("Fletcher_Reeves.jl")
include("DFP.jl")
include("BFGS.jl")


function Descent(x, mult_res, rho, xl, xu, max_int, tol_int, count, search, descent)

  # Search nesta direcao
  if descent == "Steep"
    x,dL,count = Steepest(x, mult_res, rho, xl, xu, max_int, tol_int, count, search)
  elseif descent == "FletRe"
    x,dL,count = Fletcher_Reeves(x, mult_res, rho, xl, xu, max_int, tol_int, count, search)
  elseif descent == "DFP"
    x,dL,count = DFP(x, mult_res, rho, xl, xu, max_int, tol_int, count, search)
  elseif descent == "BFGS"
    x,dL,count = BFGS(x, mult_res, rho, xl, xu, max_int, tol_int, count, search)
  else
    error("Erro na declaracao do metodo de Descent")
  end

  return x,dL,count

end
