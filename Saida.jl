#
# Dá o display dos resultados no console e arquivo
#
function Imprime(i_ext,x,rho,mult_res,valor_fun,valor_res,
                dts,nelems,max_ext,max_int,tol_ext,tol_int,
                filtro,raiof,simp,f,n_int,count)

    # Abre arquivo para armazenar iterações
    file  = string("z",dts,".txt")
    file2 = string("z",dts,"-2.txt")
    saida = open(file,"a")
    saida2 = open(file2,"a")

    # Imprime no console, primeiro print
    if     i_ext == 0
        println("\n  LagAug:: ",dts)
        println("\n  Nelems:: ",nelems," // Iter:: ",max_ext,"/",max_int," // Tol:: ",tol_ext,"/",tol_int)
        println("\n  Filtro:: ",filtro," // raiof:: ",raiof, " // pSIMP ",simp," // freq: ", f)

        println(saida,"\n  LagAug:: ",dts)
        println(saida,"\n  Nelems:: ",nelems," // Iter:: ",max_ext,"/",max_int," // Tol:: ",tol_ext,"/",tol_int)
        println(saida,"\n  Filtro:: ",filtro," // raiof:: ",raiof, " // pSIMP ",simp," // freq: ", f)

    end
    # Imprime a ultima saída com o tempo de execução
    if i_ext == -1
        tc = toq()
        @printf("  \t\tpassos internos: %d",n_int)
        @printf(saida,"  \t\tpassos internos: %d",n_int)
        @printf("\t\tevals: %d\n",count)
        @printf(saida,"\t\tevals: %d\n",count)

        @printf("\tFIM! %.2f minutos\n",tc/60.)
        @printf(saida,"\tFIM! %.2f minutos\n",tc/60.)
    else
        # numero de passos internos na iteração anterior
        if n_int > 0
            @printf("  \t\tpassos internos: %d",n_int)
            @printf(saida,"  \t\tpassos internos: %d",n_int)

            @printf("\t\tevals: %d\n",count)
        @printf(saida,"\t\tevals: %d\n",count)
        end

    # Displau dos valores atuais
        @printf("\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ",i_ext, valor_fun,mean(x))
        @printf("\n  rho: %.4f",rho)
        @printf(" \tmults. lagrange: ")
        show(mult_res);
        @printf("\n\t\tval. restrições: ")
        show(valor_res);
        @printf("\n")

        # Despeja o print no arquivo
        @printf(saida,"\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ",i_ext, valor_fun,mean(x))
        @printf(saida,"\n  rho: %.4f",rho)
        @printf(saida," \t\tmults. lagrange: ")
        show(saida,mult_res);
        @printf(saida,"\n\t\t\tval. restrições: ")
        show(saida,valor_res);
        @printf(saida,"\n")
    end

        println(saida2,valor_fun,"\t ",mean(x))

        close(saida)
        close(saida2)

    return
end
