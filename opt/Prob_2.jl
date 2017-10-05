
function F_Obj(x::Array{Float64,1})

  # Funcao objetivo
  fun = 3.0*x[1] + sqrt(3.0)*x[2]

  # Funcoes de restricao
  res = [ 18.0/x[1] + 6.0*sqrt(3.0)/x[2] - 3.0,
          5.73 - x[1],
          7.17 - x[2] ]

  return fun,res
end

# Número de variáveis e restrições
numvar = 2
numres = 3

# Chute Inicial (x0)
x=[5.0 ; 5.0]

# Limites laterais
xl=[-1E6  ; -1E6]
xu=[+1E6  ; +1E6]
