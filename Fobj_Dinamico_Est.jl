# Define o Função Objetivo
function F_Lagrangiana(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, fem, valor_zero)

    # Obtem o valor do objetivo e das restricoes
    valor_fun, valor_res, fem = F_Obj(x, fem, valor_zero)

    # Dimensoes do problema
    numvar = size(x,1)
    numres = size(valor_res,1)

    # Define o valor de L(x) no ponto x
    L = 0.0
    for j=1:numres
        L = L + max(0.0, mult_res[j] + valor_res[j]/rho)^2.0
    end

    L = valor_fun + rho*L/2.0

    return L

end


function F_Obj(x, fem ,valor_zero)

    # Acessa variáveis usadas mais de uma vez na rotina
    FL = fem.F
    wL = fem.w

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    KGL, MGL = Global_KM(fem.nelems,fem.conect,fem.ID,fem.K0,fem.M0,x,fem.simp)
    CGL = fem.alfa*MGL + fem.beta*KGL
    KDL = KGL + wL*im*CGL - (wL^2.0)*MGL

    # Resolve o sistema novamente
    UDL = vec(lufact(KDL)\FL);
    UEL = vec(lufact(KGL)\FL);

    # Pesos
    A   = 0.01
    B   = 1.0 - A

    # Função objetivo, flexibilidade dinamica
    fun = A*abs(FL'*UDL) + B*abs(FL'*UEL)

    #if valor_zero != 0.0
    #    fun = fun/valor_zero
    #end

    # Funções de restrição, volume normalizada
    res = [ (mean(x)-0.49)/0.51 ]

    # Salva em fem
    fem.KG = KGL
    fem.MG = MGL
    fem.CG = CGL
    fem.KD = KDL
    fem.UE = UEL
    fem.UD = UDL

    return fun, res, fem
end

function Sensibilidade(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1}, rho::Float64, count::Int64, fem, filt, valor_zero)

    # Inicializa a derivada interna
    dLi = zeros(Float64,fem.nelems)
    UDL = fem.UD
    UEL = fem.UE
    FL  = fem.F
    wL  = fem.w
    M0L = fem.M0
    spL = fem.simp

    # PESOS!
    A   = 0.01
    B   = 1.0 - A

    # Para flexibilidade dinamica
    a = (FL'*conj(UDL))/abs(FL'*UDL)

    # Varre os elementos
    for j=1:fem.nelems

        # Localiza deslocamentos e lambdas
        nos = fem.conect[j,:]

        # Zera Ue, complexo para cálculos Harmônicos
        UDe = complex(zeros(Float64,8))
        UEe = zeros(Float64,8)

        # Busca o U dos nós não restritos
        for k=1:4
            for q=1:2
                loc = fem.ID[nos[k],q]
                if loc != 0
                    UEe[2*k-2+q] = UEL[loc]
                    UDe[2*k-2+q] = UDL[loc]
                end #loc
            end #q
        end #k

        # derivada das matrizes de rigidez e massa 109
        dMedx  = 0.999*M0L      #correções adicionais da derivada! FIXME
        if x[j] < 0.1       # Correção de Olhoff & Du, 91
            dMedx = (6.0*(6.0E5)*x[j]^5.0+7.0*(-5.0E6)*x[j]^6.0)*M0L
        end
        dKedx  = 0.999*spL*x[j]^(spL-1.0)*fem.K0 #correções adicionais da derivada! FIXME
        dKDedx = dKedx*(1.0+im*wL*fem.beta) + dMedx*(-wL^2.0+im*wL*fem.alfa)

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/fem.NX/fem.NY/0.51

        # Renomeia as restrições e multiplicadores
        resV = valor_res[1]
        mulV =  mult_res[1]

        # Derivada do LA - flexibilidade dinamica + estatica 57
        df1dx  = real(-a*(UDe'*dKDedx*UDe))     #/valor_zero

        df2dx  = real(-UEe'*dKedx*UEe)

        dLi[j] = A*df1dx + B*df2dx + max(0.0, mulV + rho*resV)*dVdx

    end #j

    # Corrige aplicando a derivada do x filtrado em relação ao x original 65
    dLi = Derivada_Filtro(dLi, filt)

    return dLi,count
end
