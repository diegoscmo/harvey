

function Imprime(i_ext,x,rho,mult_res,valor_fun,valor_res)

    # Abre arquivo para armazenar iterações
    file  = string("z",dts,".txt")
    saida = open(file,"a")

    # Imprime no console, primeiro print
    if     i_ext == 0
        println("\n  LagAug:: ",descent," // ",search, " // ",dts)
        println(saida,"\n  LagAug:: ",descent," // ",search, " // ",dts)

    # Imprime a ultima saída com o tempo de execução
    elseif i_ext == -1
        tc = toq()
        @printf("\tFIM! %.2f minutos\n",tc/60.)
        @printf(saida,"\tFIM! %.2f minutos\n",tc/60.)

    # Displau dos valores atuais
    else
        @printf("\n  iter: %d \t\tfitness:\t%.3e\tvolume:\t%.3f ",i_ext, valor_fun,mean(x))
        @printf("\n   rho %.2f",rho)
        @printf(" \t mult. lagrange:\t")
        show(mult_res);
        @printf("\n\t\tval. restrições:\t")
        show(valor_res);
        @printf("\n")

        # Despeja o print no arquivo
        @printf(saida,"\n  iter: %d \t\tfitness:\t%.3e\tvolume:\t%.3f ",i_ext, valor_fun,mean(x))
        @printf(saida,"\n   rho %.2f",rho)
        @printf(saida," \t mult. lagrange:\t")
        show(saida,mult_res);
        @printf(saida,"\n\t\tval. restrições:\t")
        show(saida,valor_res);
        @printf(saida,"\n")
    end
        close(saida)

    return
end
