# Define o Função Objetivo
function F_Lagrangiana(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64, fem_v, fem_f)#, valor_zero)

    # Obtem o valor do objetivo e das restricoes
    valor_fun, valor_res = F_Obj(x, fem_v, fem_f)

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


function F_Obj(x, fem_v, fem_f)

    # Acessa variáveis usadas mais de uma vez na rotina
    FL = fem_f.F
    wL = fem_f.w

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    KGL, MGL = Global_KM(x,fem_f.nelems,fem_f.conect,fem_f.ID,fem_f.K0,fem_f.M0,fem_f.simp,fem_f.vminimo)
    CGL = fem_f.alfa*MGL + fem_f.beta*KGL
    KDL = KGL + wL*im*CGL - (wL^2.0)*MGL

    # Resolve o sistema novamente
    UDL = vec(lufact(KDL)\FL);

    # Função objetivo, flexibilidade dinamica
    fun = abs(dot(FL,UDL))

    # Funções de restrição, volume normalizada
    res = [ (mean(x)-0.49)/0.51 ]

    # Salva em fem
    fem_v.KG = KGL
    fem_v.MG = MGL
    fem_v.CG = CGL
    fem_v.KD = KDL
    fem_v.UD = UDL

    return fun, res
end

function Sensibilidade(x::Array{Float64,1}, valor_res::Array{Float64,1}, mult_res::Array{Float64,1},
                       rho::Float64, fem_v, fem_f, filt)

    # Carrega variáveis locais
    nelems = fem_f.nelems
    conect = fem_f.conect
    IDL = fem_f.ID
    UDL = fem_v.UD
    FL  = fem_f.F
    wL  = fem_f.w
    M0L = fem_f.M0
    spL = fem_f.simp

    # Inicializa a derivada interna
    dLi = zeros(Float64,nelems)

    # Valor para correção da derivada
    corr_min = 1.0 - fem_f.vminimo

    # Para flexibilidade dinamica
    a = dot(FL,conj(UDL)) / abs(dot(FL,UDL))

    # Varre os elementos
    for j=1:nelems

        # Localiza deslocamentos e lambdas
        nos = conect[j,:]

        # Zera Ue, complexo para cálculos Harmônicos
        Ue = complex(zeros(Float64,8))

        # Busca o U dos nós não restritos
        for k=1:4
            for q=1:2
                loc = IDL[nos[k],q]
                if loc != 0
                    Ue[2*k-2+q] = UDL[loc]
                end #loc
            end #q
        end #k

        # derivada das matrizes de rigidez e massa 109, corrigida
        dKedx  = corr_min*spL*x[j]^(spL-1.0)*fem_f.K0

        # Correção Olhoff & Du, 91 + Correção do valor minimo
        if x[j] >= 0.1
            dMedx  = corr_min*M0L
        else
            dMedx = corr_min*(6.0*(6.0E5)*x[j]^5.0+7.0*(-5.0E6)*x[j]^6.0)*M0L
        end

        # Derivada da matriz dinamica
        dKDedx = dKedx*(1.0+im*wL*fem_f.beta) + dMedx*(-wL^2.0+im*wL*fem_f.alfa)

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/fem_f.NX/fem_f.NY/0.51

        # Renomeia as restrições e multiplicadores
        resV = valor_res[1]
        mulV =  mult_res[1]

        # Derivada do LA - flexibilidade dinamica 57
        dfdx = real( -a*(Ue'*dKDedx*Ue) )     #/valor_zero
        dLi[j] = dfdx + max(0.0, mulV + rho*resV)*dVdx

    end #j

    # Corrige aplicando a derivada do x filtrado em relação ao x original 65
    dLi = Derivada_Filtro(dLi, filt)

    return dLi
end
