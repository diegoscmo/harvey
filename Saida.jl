#
# Dá o display dos resultados no console e arquivo
#

function Imprime_Atual(x::Array{Float64,1}, i_ext::Int64, n_int::Int64, count::Int64, dts::String, valor_fun::Float64,
                      rho::Float64, mult_res::Array{Float64,1}, valor_res::Array{Float64,1}, nelems::Int64, ndL::Float64)

    if dts !=  "OFF"
    # Abre os arquivos
    fmesh = string("z",dts,".pos")
    file1 = string("z",dts,".txt")
    file2 = string("z",dts,"-2.txt")

    saida  = open(file1,"a")
    saida2 = open(file2,"a")
    end
    vol = mean(x)

    t = 0.0
    if i_ext > 0
        t = toq()/60.0;
        tic();

        @printf("  \t\tpassos internos: %d | T: %.3f mins | evals: %d | dL: %.3e\n",n_int,t,count,ndL)

        if dts !=  "OFF"

        @printf(saida,"  \t\tpassos internos: %d | T: %.3f mins | evals: %d | dL: %.3e\n",n_int,t,count,ndL)
        end
    end

    # Displau dos valores atuais
    @printf("\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ", i_ext, valor_fun, vol)
    @printf("\n  rho: %.4f", rho)
    @printf(" \tmults. lagrange: ")
    show(mult_res);
    @printf("\n\t\tval. restrições: ")
    show(valor_res);
    @printf("\n")

    if dts !=  "OFF"
    # Displau dos valores atuais
    @printf(saida,"\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ", i_ext, valor_fun, vol)
    @printf(saida,"\n  rho: %.4f", rho)
    @printf(saida," \tmults. lagrange: ")
    show(saida,mult_res);
    @printf(saida,"\n\t\tval. restrições: ")
    show(saida,valor_res);
    @printf(saida,"\n")

    #FObj e volume
    println(saida2,valor_fun,"  \t",vol," \t",t)

    Adiciona_Vista_Escalar_Gmsh(fmesh, "x", nelems, x, Float64(i_ext))

    end
    # Saida gmsh

end


function Imprime_0(x, dts, valor_fun, rho, mult_res, valor_res, nelems, nnos,
    conect, coord, max_ext, max_int, tol_ext, tol_int, filtro, raiof, simp, f, descent, lsearch)

    if dts !=  "OFF"
    # Abre os arquivos
    fmesh = string("z",dts,".pos")
    file1  = string("z",dts,".txt")

    saida  = open(file1,"a")
end

    println("\n  LagAug:: ",dts," | ",descent," | ",lsearch)
    println("\n  Nelems:: ",nelems," | Iter:: ",max_ext,"/",max_int," | Tol:: ",tol_ext,"/",tol_int)
    println("\n  Filtro:: ",filtro," | raiof:: ",raiof, " | pSIMP:: ",simp," | freq:: ", f)

    if dts !=  "OFF"

        # Imprime malha no gmsh
        Inicializa_Malha_Gmsh(fmesh, nnos, nelems, conect, coord, 2)

    println(saida,"\n  LagAug:: ",dts," | ",descent," | ",lsearch)
    println(saida,"\n  Nelems:: ",nelems," | Iter:: ",max_ext,"/",max_int," | Tol:: ",tol_ext,"/",tol_int)
    println(saida,"\n  Filtro:: ",filtro," | raiof:: ",raiof, " | pSIMP:: ",simp," | freq:: ", f)

    Imprime_Atual(x, 0, 0, 0, dts, valor_fun, rho, mult_res, valor_res, nelems, 0.0)
close(saida)
    end




end

function Imprime_F(dts, n_int, count)

    # Abre os arquivos
    if dts !=  "OFF"
    file1  = string("z",dts,".txt")

    saida  = open(file1,"a")

end
    t = toq()/60.0;
    tic();

    @printf("  \t\tpassos internos: %d | T: %.3f mins | evals: %d\n",n_int,t,count)
    @printf("\tFIM!\n")

    if dts != "OFF"
    @printf(saida,"  \t\tpassos internos: %d | T: %.3f mins | evals: %d\n",n_int,t,count)
    @printf(saida,"\tFIM!\n")
    end
    if dts !=  "OFF"
    close(saida)
end
end
