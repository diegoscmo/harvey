#
# Dá o display dos resultados no console e arquivo
#

function Imprime_Int(i_int::Int64, contaev::Int64, norma::Float64, dts, tmp)

    t = tmp/60.0
    @printf("  \t\tpassos internos: %d | T: %.3f mins | evals: %d | dL: %.3e\n",i_int,t,contaev,norma)

    if dts !=  "OFF"
        file1 = string(dts,"_B.txt")
        saida  = open(file1,"a")
        @printf(saida,"  \t\tpassos internos: %d | T: %.3f mins | evals: %d | dL: %.3e\n",i_int,t,contaev,norma)
        close(saida)
    end
end

function Imprime_Ext(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1},
                     valor_fun::Float64, valor_res::Array{Float64,1}, i_ext::Int64,
                     dts::String, nel::Int64,to_plot::Array{Float64,1})

    vol = mean(x)

    # Displau dos valores atuais
    @printf("\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ", i_ext, valor_fun, vol)
    @printf("\n  rho: %.3f", rho)
    @printf(" \tmults. lagrange: ")
    show(mult_res);
    @printf("\n\t\tval. restrições: ")
    show(valor_res);
    @printf("\n")

    if dts !=  "OFF"
    # Abre os arquivos
        fmesh = string(dts,"_A.pos")
#        fmesh2 = string(dts,"_F.pos")
        file1 = string(dts,"_B.txt")
        file2 = string(dts,"_C.txt")

        saida  = open(file1,"a")
        saida2 = open(file2,"a")

        @printf(saida,"\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ", i_ext, valor_fun, vol)
        @printf(saida,"\n  rho: %.3f", rho)
        @printf(saida," \tmults. lagrange: ")
        show(saida,mult_res);
        @printf(saida,"\n\t\tval. restrições: ")
        show(saida,valor_res);
        @printf(saida,"\n")

        #Fobj e volume

        s_fobj = size(to_plot,1)

        if s_fobj == 1
            println(saida2,to_plot[1]," ",vol)
        elseif s_fobj == 2
            println(saida2,to_plot[1]," ",to_plot[2]," ",vol)
        elseif s_fobj == 3
            println(saida2,to_plot[1]," ",to_plot[2]," ",to_plot[3]," ",vol)
        end

        close(saida)
        close(saida2)
        # Imprime GMSH
        Adiciona_Vista_Escalar_Gmsh(fmesh, "x", nel, x, Float64(i_ext))
#        Adiciona_Vista_Escalar_Gmsh(fmesh2, "x", nel, Filtro_Dens(x, nel, vizi, nviz, dviz, raiof), Float64(i_ext))

    end

end


function Imprime_0(x::Array{Float64,1}, rho::Float64, mult_res::Array{Float64,1},
                     valor_fun::Float64, valor_res::Array{Float64,1}, max_ext::Int64, max_int::Int64,
                     tol_ext::Float64, tol_int::Float64, dts::String, nnos::Int64,
                      nel::Int64, ijk::Array{Int64,2}, coord::Array{Float64,2},to_plot::Array{Float64})
                      #vizi::Array{Int64,2},nviz::Array{Int64,1},dviz::Array{Float64,2},raiof::Float64)

    println("\n  LagAug:: ",dts)
    println("\n  Nelems:: ",nel," | Iter:: ", max_ext,"/", max_int," | Tol:: ",tol_ext,"/",tol_int)

    if dts !=  "OFF"

        # Abre os arquivos
        fmesh = string(dts,"_A.pos")
#        fmesh2 = string(dts,"_F.pos")
        file1  = string(dts,"_B.txt")
        file2 = string(dts,"_C.txt")

        if isfile(fmesh)
            rm(fmesh)
        end
#        if isfile(fmesh2)
#            rm(fmesh2)
#        end
        if isfile(file1)
            rm(file1)
        end
        if isfile(file2)
            rm(file2)
        end

        saida  = open(file1,"a")

        # Imprime malha no gmsh
        Inicializa_Malha_Gmsh(fmesh, nnos, nel, ijk, coord, 2)
#        Inicializa_Malha_Gmsh(fmesh2, nnos, nel, ijk, coord, 2)
        println(saida,"\n  LagAug:: ",dts)
        println(saida,"\n  Nelems:: ",nel," | Iter:: ",max_ext,"/",max_int," | Tol:: ",tol_ext,"/",tol_int)
        close(saida)

    end



    Imprime_Ext(x, rho, mult_res, valor_fun, valor_res, 0, dts, nel,to_plot)#, vizi, nviz, dviz, raiof)

end

#
#
function Le_Densidades(arquivo::String, nel::Int64)

    # Le arquivo
   texto = readlines(arquivo)

   array = texto[(end-nel):(end-1)]

   # Inicializa o array
   dens = Array{Float64}(uninitialized,nel)

   for i=1:nel

       # Parse da linha
       dados  = split(array[i]," ")[2]

       # Extrai a ultima densidade
       dens[i] = parse(Float64,dados)

   end #i

   # retorna o vetor de densidades
   return dens

end
