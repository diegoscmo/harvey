
function F_Obj(x::Array{Float64,1})

  # Funcao objetivo
  fun = -x[1]^3. - 2.*x[2]^2. + 10.*x[1] - 6. - 2.*x[2]^3.

  # Funcoes de restricao
  res = [ x[1]*x[2] - 10.,
          -x[1],
          x[2] - 10. ]

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
#xl=[0.0; 2.0]
#xu=[10.0;10.0]

# Parametros adicionais (problema do capeta!)
rho_max = 1E6
mult_max = 1E6
