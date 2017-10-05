
function F_Obj(x::Array{Float64,1})

  # Função objetivo
  fun = (x[1]-3.)^2+(x[2]-1.)^2.

  # Funções de restrição
  res = [ 5.-x[1]-x[2] ]

  return fun,res
end

# Número de variáveis e restrições
numvar = 2
numres = 1

# Chute Inicial (x0)
x=[5. ; 5.]

# Limites laterais
xl=[0.  ; 0.]
xu=[10. ;10.]

# esperado 3. e 2.
