################################################################################
#####                       Diferenças Finitas                            ######
################################################################################
#
#
#
function DF(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1}, tipo::Int64,
                                      nel::Int64, ijk::Array{Int64,2}, ID::Array{Int64,2}, K0::Array{Float64,2},
                                      M0::Array{Float64,2}, SP::Float64, vmin::Float64, F::Array{Float64,1},
                                      NX::Int64, NY::Int64, vizi::Array{Int64,2}, nviz::Array{Int64,1},
                                      dviz::Array{Float64,2}, raiof::Float64, Y0::Array{Float64,1},
                                      L0, caso, freq, alfa, beta, A, Ye)

    # Step para as diferenças finitas
    #h = 3.0*sqrt(eps())
    h = 1E-6

    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    # Declara
    dL = Array{Float64}(uninitialized,nel)

    # Gradiente em cada direção
    @inbounds for i=1:nel
        b = xf[i]
        xf[i] = b + h


        L = F_Obj(xf, rho, mult_res, 2, nel, ijk, ID, K0, M0, SP, vmin,
                                     F, NX, NY, vizi, nviz, dviz, raiof, Y0,
                                     caso, freq, alfa, beta, A, Ye,false)
        #dL[i] = (L - L0)/h
        xf[i] = b - h

        L0 = F_Obj(xf, rho, mult_res, 2, nel, ijk, ID, K0, M0, SP, vmin,
                                     F, NX, NY, vizi, nviz, dviz, raiof, Y0,
                                     caso, freq, alfa, beta, A, Ye,false)
        dL[i] = (L - L0)/(2.0*h)

        xf[i] = b
    end

    dL = dL_Dens(dL, nel, vizi, nviz, dviz, raiof)

    return dL
end
