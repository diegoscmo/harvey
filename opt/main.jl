################################################################################
# Lagrangiano Aumentado
################################################################################
# cd(ENV["TOPOPT"]); cd("opt"); >> diretório da variável de sistema

# Carrega a biblioteca de funções auxiliáres
include("F_Lagrangiana.jl")
include("Dif_Fin.jl")
include("Descent.jl")
include("LS_Methods.jl")

# Parâmetros padroes do Algoritmo
max_ext = 200    # Máximo de iteracoes externas
tol_ext = 1E-2    # Tolerância do laço externo
max_int = 500     # Máximo de 'iterações internas
tol_int = 1E-3    # Tolerância do laço interno
rho_max = 1E1     # Valor maximo de rho
mult_max = 1E2    # Valor maximo dos multiplicadores
descent = "BFGS"   # Metodo de descent (Steep, FletRe, DFP, BFGS)
search = "Line"   # Metodo de busca (Line, Golden ou Back)

# Define o problema ser resolvido
include("Prob_1.jl")

# Obtem os valores de f e g em x0 para o calculo de c0
valor_fun,valor_res = F_Obj(x)

# Inicializa multiplicadores de Lagrange (u)
mult_res = zeros(numres)

# Define fator de penalização inicial (c)
rho = max.(1E-6,min(rho_max,(2.*abs(valor_fun)/norm(max.(valor_res,0.))^2.)))

# E calcula o criterio de atualizacao do c
crho_ant = max.(0.,norm(norm(max.(valor_res, -mult_res/rho))))

# Inicializa o contador de avaliacoes da F_Obj
count = 1

# Dá o display e inicializa o laço externo
println("\n LagAug:: ",descent," // ",search)
for i_ext=1:max_ext

  # Display do valor atual
  @printf("\n  iter: %d \t fitness: %.3e \t \n x: \t",i_ext, valor_fun)
  print(x)
  @printf("\n lagrange: \t")
  print(mult_res)

  # Soluciona o problema interno (e salva a derivada)
  x,dL,count = Descent(x, mult_res, rho, xl, xu, max_int, tol_int, count, search, descent)

  # Atualiza o valor o multiplicador das restricoes (u)
  valor_fun,valor_res = F_Obj(x)
  mult_res = min.(mult_max,max.(rho*valor_res + mult_res, 0.))

  # Atualiza o multiplicador de penalizacao (c)
  crho_nov = max.(0.,norm(norm(max.(valor_res, -mult_res/rho))))
  if crho_nov >= .9*crho_ant
    rho = min.(1.1*rho, rho_max)
  end # if
  crho_ant = copy(crho_nov)

  # Verifica os criterios:
  if norm(dL) <= tol_ext &&                 # Condicao de gradiente
     norm(max.(valor_res,0.)) <= tol_ext &&   # Condicao de viabilidade
     norm(valor_res'*mult_res) < tol_ext #&&  # Condicao de complementariedade
     #sum(valor_res .>= 0.) == numres         # Verific. dos mult. de desigualdade  #!modificado 05/10
      break
  end # if criterios

end # for i_ext

# Imprime a saida dos resultados
println("\n Aval. da F_Obj = [$count]")
println("\n   Pt de minimo = $x")
println("\n     Valor de F = [$valor_fun]")
println("\n     Restricoes = $valor_res")
