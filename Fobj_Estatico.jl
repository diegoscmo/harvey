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


function F_Obj(x, fem, valor_zero)

    # Remonta matriz de rigidez global, aqui é aplicado o SIMP
    KGL = Global_K(fem.nelems,fem.conect,fem.ID,fem.K0,x,fem.simp)
    FL  = fem.F

    # Resolve o sistema novamente
    UEL = vec(lufact(KGL)\fem.F);

    # Função objetivo, flexibilidade estática
    fun = FL'*UEL

    #if valor_zero != 0.0
    #    fun = fun/valor_zero
    #end

    # Funções de restrição, volume normalizada
    res = [ (mean(x)-0.49)/0.51 ]

    # Salva em fem
    fem.KG = KGL
    fem.UE = UEL

    return fun, res, fem
end

function Sensibilidade(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1}, rho::Float64, count::Int64, fem, filt, valor_zero)

    # Inicializa a derivada interna
    dLi = zeros(Float64,fem.nelems)

    # Varre os elementos
    for j=1:fem.nelems

        # Localiza deslocamentos e lambdas
        nos = fem.conect[j,:]

        # Zera Ue, complexo para cálculos Harmônicos
        Ue = complex(zeros(Float64,8))

        # Busca o U dos nós não restritos
        for k=1:4
            for q=1:2
                loc = fem.ID[nos[k],q]
                if loc != 0
                    Ue[2*k-2+q] = fem.UE[loc]
                end #loc
            end #q
        end #k

        # derivada das matrizes de rigidez 109
        dKedx  = 0.999*fem.simp*x[j]^(fem.simp-1.0)*fem.K0

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/fem.NX/fem.NY/0.51

        # Renomeia as restrições e multiplicadores
        resV = valor_res[1]
        mulV =  mult_res[1]

        # Derivada do LA - flexibilidade estatica 57
        dfdx = real(-Ue'*dKedx*Ue)/valor_zero
        dLi[j] = dfdx + rho*max(0.0,mulV + rho*resV)*dVdx

    end #j

    # Corrige aplicando a derivada do x filtrado em relação ao x original 65
    dLi = Derivada_Filtro(dLi, filt)

    return dLi,count
end
