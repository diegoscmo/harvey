
function Dif_Fin(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, count::Int64)

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP
    KG,MG,F = Global(nelems,ijk,ID,K0,M0,x,simp,nforcas,nos_forcas,nr_gl_livres)

    # Monta a matriz equivalente para análise harmônica
    CG = alfa*MG + beta*KG
    KD = KG + w*im*CG - w^2.0*MG

    # Resolve o sistema
    U = vec(lufact(KD)\F)
    count = count+1
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
    nos = ijk[j,:]

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
    if x[j] < 0.1   #olhoff & du
        dMedx = (6.0*6.0E5*x[j]^5.0+7.0*(-5.0E6)*x[j]^6.0)*M0
    end
    dKedx  = simp*x[j]^(simp-1.0)*K0
    dKDedx = dKedx*(1.0+im*w*beta)+dMedx*(-w^2.0+im*w*alfa)

    #derivada da potencia ativa 109/181
    #dfdx = (0.5)*w*real(im*Ue'*dKDedx*Ue)

    # derivada da flexibilidade dinamica 109
    #dfdx = -real(a*(Ue'*dKDedx*Ue))

    # derivada da flexibilidade estatica 57
    dfdx = -simp*x[j]^(simp-1.0)*Ue'*K0*Ue

    # Derivada R
    #dRdx = ((w^2.0/Ep)*(conj(Ue)'*dMedx*Ue)-(Ec/(4.0*Ep^2.0))*(conj(Ue)'*dKedx*Ue))+
    #                                             real(lambdae'*(dKDedx*Ue)) #-dFdx))

    # derivada de V 58
    dVdx = 1.0/NX/NY # volume de um elemento se o todo é =1.0

    # Derivada do Lagrangiano Aumentado
    #dHxdx = [dVdx; dRdx]
    #dL[j] = real(dPaedx + (valor_res.*dHxdx)'*(rho*valor_res+mult_res))

    dL[j] = real(dfdx + dVdx*(rho*valor_res[1]+mult_res[1]))

    end

    # filtro de sensibilidade 63
    dLf = zeros(Float64,nelems)

    for j=1:nelems
#        dLf[j] = (1.0/(x[j]*sum(H[j,:])))*sum(H[j,:].*x.*dL)
    end

    return dL,count
end
