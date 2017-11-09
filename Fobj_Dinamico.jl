# Define o Função Objetivo
function F_Obj(x)

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    global KG,MG,CG,KD,F,U
    KG,MG,F = Global(nelems,conect,ID,K0,M0,x,simp,nforcas,nos_forcas,nr_gl_livres)

    # Monta a matriz equivalente para análise harmônica
    CG = alfa*MG + beta*KG
    KD = KG + w*im*CG - w^2.0*MG

    # Resolve o sistema
    U = vec(lufact(KD)\F);

    # Função objetivo, flexibilidade dinamica
    fun = abs(F'*U)
    try
        fun = fun/valor_zero
    end

    # Funções de restrição, volume normalizada
    res = [ (mean(x)-0.49)/0.51 ]

    return fun,res
end

function Dif_Fin(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, count::Int64)

    # Inicializa a derivada
    dLi = zeros(Float64,nelems)

    # Para flexibilidade dinamica
    a = (F'*conj(U))/abs(U'*F)

    # Varre os elementos
    for j=1:nelems

        # Localiza deslocamentos e lambdas
        nos = conect[j,:]

        # Zera Ue, complexo para cálculos Harmônicos
        Ue = complex(zeros(Float64,8))

        # Busca o U dos nós não restritos
        for k=1:4
            for q=1:2
                loc = ID[nos[k],q]
                if loc != 0
                    Ue[2*k-2+q] = U[loc]
                end #loc
            end #q
        end #k

        # derivada das matrizes de rigidez e massa 109
        dMedx  = M0
        if x[j] < 0.1       # Correção de Olhoff & Du, 91
            dMedx = (6.0*6.0E5*x[j]^5.0+7.0*(-5.0E6)*x[j]^6.0)*M0
        end
        dKedx  = simp*x[j]^(simp-1.0)*K0
        dKDedx = dKedx*(1.0+im*w*beta)+dMedx*(-w^2.0+im*w*alfa)

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/NX/NY/0.51

        # Renomeia as restrições e multiplicadores
        resV = valor_res[1]
        mulV =  mult_res[1]

        #derivada da flexibilidade dinamica 109
        dfdx = real(-a*(Ue'*dKDedx*Ue))/valor_zero

        # Derivada do LA - flexibilidade estatica 57 #FIXME, ABRIR DERIVADAS
        dLi[j] = dfdx + rho*max(0.0,mulV+resV)*dVdx

    end #j

    # Corrige aplicando a derivada do x filtrado em relação ao x original 65
    dLi = Derivada_Filtro(dLi, raiof, vizi, nviz, dviz, filtro)

    count += 1

    return dLi,count
end
