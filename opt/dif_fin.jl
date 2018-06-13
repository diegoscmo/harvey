################################################################################
#####                       Diferenças Finitas                            ######
################################################################################
#
#
#
function F_Obj_DFC(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1}, tipo::Int64,
     nnos::Int64, nel::Int64, ijk::Array{Int64,2}, coord::Array{Float64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
           M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
              NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
            raiof::Float64, valor_0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64,
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64, nos_viz,dts,P,q)


    # Calcula sensibilidade
    dL_orig, = F_Obj(x, rho, mult_res, 3, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                 F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                 Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax, nos_viz,dts,P,q)


    # Step para as diferenças finitas
    #h = 3.0*sqrt(eps())
    h = 1E-6

    # Calcula o L0 caso derive só pra frente ou pra trás
    L0, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                      F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                      Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

    # Declara
    dL = Array{Float64}(undef,nel)

    # Gradiente em cada direção
    @showprogress "  verificando por dfc..." for i=1:nel

        b = x[i]

        # Se estiver no limite superior, derivada pra trás
        if b > (1.0-h)

            x[i] = b - h

            L, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                                  F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                                  Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

            dL[i] = (L0 - L)/h

        # Se estive no limite inferior, derivada pra frente
        elseif b < h

            x[i] = b + h

             L, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                                F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                            Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

            dL[i] = (L - L0)/h

        # Se no intervalo normal, derivadas centrais
        else

            x[i] = b + h

            L2, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

            x[i] = b - h

            L1, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)
            dL[i] = (L2 - L1)/(2.0*h)

        end

        x[i] = b

    end

    # NAO CORRIGE O FILTRO!!

    # media da direção da derivada
    media = median(dL_orig./dL)

    println()

    file = string("results/",dts,"/derivada_dfc_",h,".txt")
    if isfile(file);  rm(file);  end
    saida = open(file,"a")
    println("Média de Sensibilidade/DFC = $media")
    println(saida,"Média de Sensibilidade/DFC = $media")
    println(saida,"Sensibilidade       Diferenças Finitas Centrais  x")
    for z=1:nel;  println(saida,dL_orig[z],"\t\t",dL[z],"\t\t",x[z]) ; end
    close(saida);

    return dL,0.0,0.0,0.0
end

function F_Obj_DFF(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1}, tipo::Int64,
     nnos::Int64, nel::Int64, ijk::Array{Int64,2}, coord::Array{Float64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
           M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
              NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
            raiof::Float64, valor_0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64,
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64,nos_viz,dts,P,q)

    # Step para as diferenças finitas
    #h = 3.0*sqrt(eps())
    h = 1E-8

    L0, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                          F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                          Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

    # Declara
    dL = Array{Float64}(undef,nel)

    # Gradiente em cada direção
    for i=1:nel
        b = x[i]
        x[i] = b + h

        L, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)
        dL[i] = (L - L0)/h
        x[i] = b
    end

    return dL,L0,0.0,0.0
end
