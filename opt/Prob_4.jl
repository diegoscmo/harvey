
function F_Obj(x::Array{Float64,1})

  # Funcao objetivo
  fun = x[1]^3.

  # Funcoes de restricao
  res = [ 0. ]

  return fun,res
end

# Número de variáveis
numvar = 1
# e restricoes (sempre >= 1)
numres = 1

# Chute Inicial (x0)
x=[5.]

# Limites laterais
xl=[-5.]
xu=[+5.]
