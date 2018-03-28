################################################################################
#####                       Diferenças Finitas                            ######
################################################################################
#
#
#
function F_Obj_DF(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
     nnos::Int64, nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
           M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1}, NX::Int64,
              NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2},
            raiof::Float64, valor_0::Array{Float64,1}, Sy::Float64, freq::Float64, alfa::Float64,
           beta::Float64, A::Float64, Ye::Float64, CBA::Array{Float64,3}, QP::Float64, csi::Float64, dmax::Float64)

    dL_orig, = F_Obj(x, rho, mult_res, 3, nnos, nel, ijk, ID, K0, M0, SP, vmin,
                        F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                        Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax)

    # Step para as diferenças finitas
    h = 3.0*sqrt(eps())

    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Declara
    dL = Array{Float64}(undef,nel)

    # Gradiente em cada direção
    @inbounds for i=1:nel
        b = xc[i]
        xc[i] = b + h


        L, = F_Obj(xc, rho, mult_res, 2, nnos, nel, ijk, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax)
        #dL[i] = (L - L0)/h
        xc[i] = b - h

        L0, = F_Obj(xc, rho, mult_res, 2, nnos, nel, ijk, ID, K0, M0, SP, vmin,
                                              F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                              Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax)
        dL[i] = (L - L0)/(2.0*h)

        xc[i] = b
    end

    dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)


    file = "diffin.txt"
    if isfile(file);  rm(file);  end
    saida = open(file,"a")
    for z=1:nel;  println(saida,dL_orig[z],"   ",dL[z]) ; end
    close(saida);  error("Diferenças finitas... OK")

    return dL,0.0
end
