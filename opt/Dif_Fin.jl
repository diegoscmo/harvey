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
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64, nos_viz,dts)


    dL_orig, = F_Obj(x, rho, mult_res, 3, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                 F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                 Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax, nos_viz,dts)

    # Step para as diferenças finitas
    #h = 3.0*sqrt(eps())
    h = 1E-8

    # Arruma as densidades
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(xc, nel, ijk, ID, K0, M0, SP, vmin)

    # Checa o condicionamento
    w  = 2.0*pi*freq
    KD = KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG
    condor = Checa_Cond(KD,1E-1,1000)

    println(" Calculando derivada por Diferenças Finitas...")

    # Declara
    dL = Array{Float64}(undef,nel)

    # Gradiente em cada direção
    for i=1:nel
        b = x[i]
        x[i] = b + h


        L, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts)
        #dL[i] = (L - L0)/h
        x[i] = b - h

        L0, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts)
        dL[i] = (L - L0)/(2.0*h)

        x[i] = b
    end

    dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

    # media da direção da derivada
    media = sum(dL_orig./dL)/nel

    println()

    file = string("results\\",dts,"\\",dts,"_diffin.txt")
    if isfile(file);  rm(file);  end
    saida = open(file,"a")
    println("Condicionamento de KD = $condor")
    println("Média de Sensibilidade/DFC = $media")
    println(saida,"Condicionamento de KD = $condor")
    println(saida,"Média de Sensibilidade/DFC = $media")
    println(saida,"Sensibilidade       Diferenças Finitas Centrais")
    for z=1:nel;  println(saida,dL_orig[z],"   ",dL[z]) ; end
    close(saida);  error("Diferenças finitas... OK")

    return dL,0.0,0.0,0.0
end

function F_Obj_DFF(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1}, tipo::Int64,
     nnos::Int64, nel::Int64, ijk::Array{Int64,2}, coord::Array{Float64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
           M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
              NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
            raiof::Float64, valor_0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64,
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64,nos_viz,dts)

    # Step para as diferenças finitas
    #h = 3.0*sqrt(eps())
    h = 1E-8

    L0, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                          F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                          Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts)

    # Declara
    dL = Array{Float64}(undef,nel)

    # Gradiente em cada direção
    for i=1:nel
        b = x[i]
        x[i] = b + h

        L, = F_Obj(x, rho, mult_res, 2, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts)
        dL[i] = (L - L0)/h
        x[i] = b
    end

    dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

    return dL,L0,0.0,0.0
end
