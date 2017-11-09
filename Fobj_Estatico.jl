# Define o Função Objetivo
function F_Lagrangiana(x::Array{Float64,1}, mult_res::Array{Float64,1}, rho::Float64,
                    fem,valor_zero)

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


function F_Obj(x,fem,valor_zero)

    # Monta matrizes de rigidez e massa Globais, aqui é aplicado o SIMP e  CORREÇÔES OLHOF&DU!
    fem.KG,fem.MG,fem.F = Global(fem.nelems,fem.conect,fem.ID,fem.K0,fem.M0,x,fem.simp,fem.nforcas,fem.nos_forcas,fem.nr_gl_livres)
    # Resolve o sistema novamente
    fem.U = vec(lufact(fem.KG)\fem.F);

    # Função objetivo, flexibilidade estática
    fun = fem.U'*fem.KG*fem.U

    if valor_zero != 0.0
        fun = fun/valor_zero
    end

    # Funções de restrição, volume normalizada
    res = [ (mean(x)-0.49)/0.51 ]

    return fun, res , fem
end

function Sensibilidade(x::Array{Float64,1}, valor_res, mult_res::Array{Float64,1}, rho::Float64, count::Int64,fem,filt,valor_zero)

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
                    Ue[2*k-2+q] = fem.U[loc]
                end #loc
            end #q
        end #k

        # derivada das matrizes de rigidez 109
        dKedx  = fem.simp*x[j]^(fem.simp-1.0)*fem.K0

        # Derivada da restrição de Volume Normalizada (1.0-0.49), 58
        dVdx = 1.0/fem.NX/fem.NY/0.51

        # Renomeia as restrições e multiplicadores
        resV = valor_res[1]
        mulV =  mult_res[1]

        # Derivada do LA - flexibilidade estatica 57
        dfdx = real(-Ue'*dKedx*Ue)/valor_zero
        dLi[j] = dfdx + rho*max(0.0,mulV+resV)*dVdx

    end #j

    # Corrige aplicando a derivada do x filtrado em relação ao x original 65
    dLi = Derivada_Filtro(dLi, filt)

    count += 1

    return dLi,count
end
