# Define o Função Objetivo
function F_Lagrangiana(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64,
                       fem_v, fem_f) #, valor_zero)

    # Obtem o valor do objetivo e das restricoes
    valor_fun, valor_res = F_Obj(x, fem_v, fem_f) #, valor_zero)

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


function F_Obj(x::Array{Float64,1}, fem_v, fem_f)#, valor_zero)

    # Remonta matriz de rigidez global, aqui é aplicado o SIMP
    KGL,MGL = Global_KM(x,fem_f.nelems,fem_f.conect,fem_f.ID,fem_f.K0,fem_f.M0,fem_f.simp,fem_f.vminimo)
    FL  = fem_f.F

    # Resolve o sistema novamente
    UEL = vec(lufact(KGL)\FL)

    # Função objetivo, flexibilidade estática
    fun = dot(FL,UEL)

    # Funções de restrição, volume normalizada
    res = [ (mean(x)-0.49)/0.51 ]

    # Salva em fem
    fem_v.KG = KGL
    fem_v.UE = UEL

    return fun, res
end

function Sensibilidade(x::Array{Float64,1}, valor_res::Array{Float64,1}, mult_res::Array{Float64,1},
                       rho::Float64, fem_v, fem_f, filt)#, valor_zero)

    nelems = fem_f.nelems
    conect = fem_f.conect
    ID     = fem_f.ID
    spL   = fem_f.simp
    UEL    = fem_v.UE

    # Inicializa a derivada interna
    dLi = zeros(Float64,nelems)

    # Valor para correção da derivada
    corr_min = 1.0 - fem_f.vminimo

    # Varre os elementos
    for j=1:nelems

        # Localiza deslocamentos e lambdas
        nos = conect[j,:]

        # Zera Ue, complexo para cálculos Harmônicos
        UEe = zeros(Float64,8)

        # Busca o U dos nós não restritos
        for k=1:4
            for q=1:2
                loc = ID[nos[k],q]
                if loc != 0
                    UEe[2*k-2+q] = UEL[loc]
                end #loc
            end #q
        end #k

        # derivada das matrizes de rigidez e massa 109, corrigida
        dKedx  = corr_min*spL*x[j]^(spL-1.0)*fem_f.K0

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/fem_f.NX/fem_f.NY/0.51

        # Renomeia as restrições e multiplicadores
        resV = valor_res[1]
        mulV =  mult_res[1]

        # Derivada do LA - flexibilidade estatica 57
        dfdx = real(-UEe'*dKedx*UEe)
        dLi[j] = dfdx + max(0.0,mulV + rho*resV)*dVdx

    end #j

    # Corrige aplicando a derivada do x filtrado em relação ao x original 65
    dLi = Derivada_Filtro(dLi, filt)

    return dLi
end
