# Define o Função Objetivo
function F_Obj(x)

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    global KG,MG,CG,KD,F,KD,U
    KG,MG,F = Global(nelems,conect,ID,K0,M0,x,simp,nforcas,nos_forcas,nr_gl_livres)

    # Monta a matriz equivalente para análise harmônica
    #CG = alfa*MG + beta*KG
    #KD = KG + w*im*CG - w^2.0*MG
    KD = KG
       # Resolve o sistema
    U = vec(lufact(KD)\F);

    # Função objetivo

    #fun = real((0.5*w^2.0)*conj(U)'*CG*U)      # Potencia Ativa
    #fun = (abs(U'*F))^2.0                       # Flexibilidade dinamica
    fun = U'*KG*U                              # Flexibilidade Estática

    # Funções de restrição
    res = [
            # Restrição de Volume (Normalizada)
            (mean(x)-0.49)/0.51
            # Restrição do R
        #   -1.0+(100.0+100.0*log(w^2.0*(real(U'*MG*U)/real(U'*KG*U))))
            ]
    return fun,res
end

function Dif_Fin(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, count::Int64)

    # derivada de R 300/301
    #Ep = real(U'*MG*U)
    #Ec = real(U'*KG*U)

    # Resolve problema adjunto
    #adj = 0.5*((Ec/Ep^2.0)*KG*conj(U)-(w^2.0/Ep)*MG*conj(U))
    #lambda = KD\adj

    # Para flexibilidade dinamica
    #a = (F'*conj(U))/sqrt((F'*real(U))^2.0 + (F'*imag(U))^2.0)

    dL = zeros(Float64,nelems)
    for j=1:nelems

    # Localiza deslocamentos e lambdas
    nos = conect[j,:]

    # Inicializa
    Ue = zeros(Float64,8)
    Ue = complex(Ue)
    #lambdae = zeros(Float64,8)
    #lambdae = complex(lambdae)
    for k=1:4
      for q=1:2
        loc = ID[nos[k],q]
        if loc != 0
          Ue[2*k-2+q] = U[loc]
    #      lambdae[2*k-2+q] = lambda[loc]
        end
      end
    end

    # derivada das matrizes de rigidez e massa 109
    dMedx  = M0
    if x[j] < 0.1       # Correção de Olhoff & Du, 91
        dMedx = (6.0*6.0E5*x[j]^5.0+7.0*(-5.0E6)*x[j]^6.0)*M0
    end
    dKedx  = simp*x[j]^(simp-1.0)*K0
    #dKDedx = dKedx*(1.0+im*w*beta)+dMedx*(-w^2.0+im*w*alfa)

    # Derivada da restrição de Volume Normalizada, 58
    dVdx = 1.0/NX/NY/0.51

    # Renomeia as restrições e multiplicadores
    resV = valor_res[1]
    mulV =  mult_res[1]

    # Derivada do LA - flexibilidade estatica 57
    dL[j] = real(-Ue'*dKedx*Ue + rho*max(0.0,mulV+resV)*dVdx)

    # derivada da flexibilidade dinamica 109
    #dfdx = -real(a*(Ue'*dKDedx*Ue))

    #derivada da potencia ativa 109/181
    #dfdx = (0.5)*w*real(im*Ue'*dKDedx*Ue)

    # Derivada R
    #dRdx = ((w^2.0/Ep)*(conj(Ue)'*dMedx*Ue)-(Ec/(4.0*Ep^2.0))*(conj(Ue)'*dKedx*Ue))+
    #                                             real(lambdae'*(dKDedx*Ue)) #-dFdx))

    end

    # filtro de sensibilidade 63
    #dLf = zeros(Float64,nelems)
    # topans 75
    for j=1:nelems
    #    dLf[j] = (1.0/(x[j]*sum(H[j,:])))*sum(H[j,:].*x.*dL)
    end

    return dL,count
end
