#
# Varredura da Flexibilidade Dinâmica pela frequência em Hertz
#
function Harmonica(K,M,alfa,beta,F,nnos,ID,fstep,fmax)

 fs = []
 ha = []

 C = alfa*M + beta*K

for f = 0:fstep:fmax

# Monta a matriz equivalente para análise harmônica
KD = K + 2*pi*f*im*C - (2*pi*f)^2*M

# Resolve o sistema
U,Fe = FEM_Solve(KD,F,nnos,ID)

# Flexibilidade dinamica
Y = abs(U'*Fe)^2

push!(fs,f)
push!(ha,Y)

end

return fs,ha

end

# Resolve o sistema de equações e expande os vetores
function FEM_Solve(KD,F,nnos,ID)

    # Soluciona o sistema de equações
    L = lufact(KD)
    U_sem_zeros = vec(L\F)

    # Expande o vetor de deslocamentos, adicionando os graus de liberdade travados:
    U  = Expande_Vetor(U_sem_zeros,nnos,ID)
    Fe = Expande_Vetor(F,nnos,ID)

    return U,Fe
end
