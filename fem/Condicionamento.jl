# Teste de condicionamento de matrizes

#
# Checa condicionamento de uma matriz de Hilbert
#
function Cond_H(t,tol)

    println("Montando Hermitiana de $t por $t")
    H = sparse(zeros(t,t))
    for i=1:t
        for j=1:t
            H[i,j] = 1.0/(i+j-1.0)
        end
    end

    Checa_Cond(H,tol)

end

#
# Checa o condicionamento de uma matriz A qualquer
#
function Checa_Cond(A,tol::Float64,imax::Int64)

    #maior = Pot_Maior(A,tol,imax)
    #menor = Pot_Menor(A,tol,imax)

    maior = eigs(A,which=:LM)[1][1]
    menor = eigs(A,which=:SM)[1][1]
    cond = abs(maior)/abs(menor)

    return cond
end

#
# Busca o maior autovalor pelo método da potência
#
function Pot_Maior(A,tol,imax)

    # Primeira estimativa de x
    n = size(A,1)
    x = ones(n)

    # Primeira estimativa de lambda
    lambda = 0.0

    # Primeira estimativa do intervalo de convergência
    dif = 1.0

    iter = 0

    # Itera
    while dif > tol

        # Calcula a nova estimativa para o x
        xk = A*x

        # Calcula a estimativa do autovalor
        lambda_k = xk[1]/x[1]

        # Intervalo para calculo de tolerancia
        dif = abs(lambda_k - lambda)

        # Atualiza x e lambda
        x = copy(xk)/norm(xk)
        lambda = copy(lambda_k)

        # Incrementa iteração e verifica saída
        iter += 1
        if iter > imax
            break
        end

    end #while

    return lambda
end

#
# Busca o menor autovalor pelo método da potência
#
function Pot_Menor(A,tol,imax)

    # Primeira estimativa de x
    n = size(A,1)
    x = ones(n)

    # Primeira estimativa de lambda
    lambda = 0.0

    # Primeira estimativa do intervalo de convergência
    dif = 1.0

    invA = lufact(A)

    iter = 0

    # Itera
    while dif > tol

        # Calcula a nova estimativa para o x
        xk = A\x

        # Calcula a estimativa do autovalor
        lambda_k = x[1]/xk[1]

        # Intervalo para calculo de tolerancia
        dif = abs(lambda_k - lambda)

        # Atualiza x e lambda
        x = copy(xk)/norm(xk)
        lambda = copy(lambda_k)

        # Incrementa iteração e verifica saída
        iter += 1
        if iter > imax
            break
        end

    end #while

    return lambda
end

# Tennativa de melhoria de condicionamento com QR
function Resolve_QR(A,b)

n=size(b,1)

# Agora resolvendo com R*x = Q^H*b
QRF = qr(Array(A),full=true)
Q = QRF[1]
R = QRF[2]

# b2= Q^H*b
b2 = adjoint(Q)*b

# Executando a retrosubstituição de R*x = b2
x = zeros(Complex,n)
for i=n:(-1):1
     # Armazena o somatório
     somat = 0.0
     for j=(i+1):n
         somat = somat + R[i,j]*x[j]
     end #j
     x[i] = (b2[i] - somat)/R[i,i]
 end #i

return x

end

function Retro_QR(Q,R,b,n)

    # b2= Q^H*b
    b2 = adjoint(Q)*b

    # Executando a retrosubstituição de R*x = b2
    x = Array{Complex{Float64}}(undef,n)
    for i=n:(-1):1
         # Armazena o somatório
         somat = 0.0
         for j=(i+1):n
             somat = somat + R[i,j]*x[j]
         end #j
         x[i] = (b2[i] - somat)/R[i,i]
     end #i

    return x

end
