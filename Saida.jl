#
# Dá o display dos resultados no console e arquivo
#
function Imprime(i_ext,x,rho,mult_res,valor_fun,valor_res,
                dts,nelems,max_ext,max_int,tol_ext,tol_int,
                filtro,raiof,simp,f)

    # Abre arquivo para armazenar iterações
    file  = string("z",dts,".txt")
    saida = open(file,"a")

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
        @printf("\tFIM! %.2f minutos\n",tc/60.)
        @printf(saida,"\tFIM! %.2f minutos\n",tc/60.)
    else
    # Displau dos valores atuais
        @printf("\n  iter: %d \t\tfitness: %.4f\tvolume:\t%.3f ",i_ext, valor_fun,mean(x))
        @printf("\n   rho: %.2f",rho)
        @printf(" \t mult. lagrange: ")
        show(mult_res);
        @printf("\n\t\tval. restrições: ")
        show(valor_res);
        @printf("\n")

        # Despeja o print no arquivo
        @printf(saida,"\n  iter: %d \t\tfitness: %.4f\tvolume:\t%.3f ",i_ext, valor_fun,mean(x))
        @printf(saida,"\n   rho: %.2f",rho)
        @printf(saida," \t mult. lagrange: ")
        show(saida,mult_res);
        @printf(saida,"\n\t\tval. restrições: ")
        show(saida,valor_res);
        @printf(saida,"\n")
    end
        close(saida)

    return
end
