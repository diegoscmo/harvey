# Define o Função Objetivo
function F_Lagrangiana(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, fem, valor_zero)

    # Obtem o valor do objetivo e das restricoes
    valor_fun,valor_res,fem = F_Obj(x,fem,valor_zero)

    # Dimensoes do problema
    numvar = size(x,1)
    numres = size(valor_res,1)

    # Define o valor de L(x) no ponto x
    L = 0.0
    for j=1:numres
        L = L + max(0.0, mult_res[j]+valor_res[j]/rho)^2.0
    end

    L = valor_fun + rho*L/2.0

    return L

end


function F_Obj(x, fem ,valor_zero)

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    fem.KG,fem.MG,fem.F = Global(fem.nelems,fem.conect,fem.ID,fem.K0,fem.M0,x,fem.simp,fem.nforcas,fem.nos_forcas,fem.nr_gl_livres)
    fem.CG = fem.alfa*fem.MG + fem.beta*fem.KG
    fem.KD = fem.KG + fem.w*im*fem.CG - fem.w^2.0*fem.MG

    # Resolve o sistema novamente
    fem.U = vec(lufact(fem.KD)\fem.F);

    # Função objetivo, flexibilidade estática
    fun = real((0.5*fem.w^2.0)*conj(fem.U)'*fem.CG*fem.U)

    if valor_zero != 0.0
        fun = fun/valor_zero
    end

    # Funções de restrição, volume normalizada
    res = [ (mean(x)-0.49)/0.51
            -1.0 + (1.0+1.0*log(fem.w^2.0*(real(fem.U'*fem.MG*fem.U)/
                                        real(fem.U'*fem.KG*fem.U))))
            ]

    return fun, res , fem
end

function Sensibilidade(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1}, rho::Float64, count::Int64, fem, filt, valor_zero)

    # Inicializa a derivada interna
    dLi = zeros(Float64,fem.nelems)

    # derivada de R 300/301
    Ep = real(fem.U'*fem.MG*fem.U)
    Ec = real(fem.U'*fem.KG*fem.U)

    # Resolve problema adjunto
    adj = 0.5*((Ec/Ep^2.0)*fem.KG*conj(fem.U)-(fem.w^2.0/Ep)*fem.MG*conj(fem.U))
    lambda = fem.KD\adj

    # Varre os elementos
    for j=1:fem.nelems

        # Localiza deslocamentos e lambdas
        nos = fem.conect[j,:]

        # Zera Ue, complexo para cálculos Harmônicos
        Ue = complex(zeros(Float64,8))
        Lbe = zeros(Float64,8)
        # Busca o U dos nós não restritos
        for k=1:4
            for q=1:2
                loc = fem.ID[nos[k],q]
                if loc != 0
                    Ue[2*k-2+q] = fem.U[loc]
                    Lbe[2*k-2+q] = lambda[loc]
                end #loc
            end #q
        end #k

        # derivada das matrizes de rigidez e massa 109
        dMedx  = fem.M0
        if x[j] < 0.1       # Correção de Olhoff & Du, 91
            dMedx = (6.0*6.0E5*x[j]^5.0+7.0*(-5.0E6)*x[j]^6.0)*fem.M0
        end
        dKedx  = fem.simp*x[j]^(fem.simp-1.0)*fem.K0
        dKDedx = dKedx*(1.0+im*fem.w*fem.beta) + dMedx*(-fem.w^2.0+im*fem.w*fem.alfa)

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/fem.NX/fem.NY/0.51

        # Derivada R
        dRdx = ((fem.w^2.0/Ep)*(conj(Ue)'*dMedx*Ue)-(Ec/(4.0*Ep^2.0))*(conj(Ue)'*dKedx*Ue))+real(Lbe*(dKDedx*Ue)) #-dFdx))


        # Renomeia as restrições e multiplicadores
        resV = valor_res[1]
        mulV =  mult_res[1]
        resR = valor_res[2]
        mulR =  mult_res[2]

        #derivada da potencia ativa 109/181
        dfdx = 0.5*fem.w*real(im*Ue'*dKDedx*Ue)/valor_zero

        # Derivada do LA #FIXME, ABRIR DERIVADAS NA MAO
        dLi[j] = real(dfdx + rho*(max(0.0,mulV/rho+resV)*dVdx + max(0.0,mulR/rho+resR)*dRdx))

    end #j

    # Corrige aplicando a derivada do x filtrado em relação ao x original 65
    dLi = Derivada_Filtro(dLi, filt)

    count += 1

    return dLi,count
end
