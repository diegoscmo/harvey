################################################################################
#####                Saída de Arquivos + Display no Console               ######
################################################################################

#
# Dá o display dos resultados no console e arquivo
#
function Imprime_Int(i_int::Int64, contaev::Int64, norma::Float64, dts::AbstractString, tmp::Float64)

    t = tmp/60.0
    @printf("  \t\tpassos internos: %d | T: %.3f mins | evals: %d | dL: %.3e\n",i_int,t,contaev,norma)

    if dts !=  "OFF"
        file1 = string(dts,"_log.txt")
        saida  = open(file1,"a")
        @printf(saida,"  \t\tpassos internos: %d | T: %.3f mins | evals: %d | dL: %.3e\n",i_int,t,contaev,norma)
        close(saida)
    end
end

function Imprime_Ext(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1},
                     valor_fun::Float64, valor_res::Array{Float64,1}, i_ext::Int64,
                     dts::String, nel::Int64, to_plot::Array{Float64,1}, TS::Array{Float64,2},
                     vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64)

    xf =  Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    xh, = Filtro_Heavi(x, nel, 0.5)
    vol = mean(xf)

    # Displau dos valores atuais
    @printf("\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ", i_ext, valor_fun, vol)
    @printf("\n  rho: %.3f", mean(rho))
    @printf(" \tmults. lagrange: ")
    show(mean(mult_res));
    @printf("\n\t\tval. restrições: ")
    show(mean(valor_res));
    @printf("\n")

    if dts !=  "OFF"
    # Abre os arquivos
        fmesh  = string(dts,"_dens_o.pos")
        fmesh2 = string(dts,"_dens_f.pos")
        fmesh3 = string(dts,"_stress.pos")
        fmesh4 = string(dts,"_dens_H.pos")
        file1  = string(dts,"_log.txt")
        file2  = string(dts,"_monitora.txt")
        file3  = string(dts,"_save.txt")

        saida  = open(file1,"a")
        saida2 = open(file2,"a")

        @printf(saida,"\n  iter: %d \torigin. fitness: %.6e\tvolume:\t%.5f ", i_ext, valor_fun, vol)
        @printf(saida,"\n  rho: %.3f", mean(rho))
        @printf(saida," \tmults. lagrange: ")
        show(saida,mean(mult_res));
        @printf(saida,"\n\t\tval. restrições: ")
        show(saida,mean(valor_res));
        @printf(saida,"\n")

        close(saida)

        # Plots de valores
        s_toplot = size(to_plot,1)
        for j = 1:s_toplot
            @printf(saida2, "%.3e ", to_plot[j])
        end
        @printf(saida2,"\n")
        close(saida2)

        # Imprime para retomar
        if isfile(file3); rm(file3); end
        if i_ext>0
            saida3 = open(file3,"a")
            println(saida3,i_ext)
            for j=1:size(rho,1)
                println(saida3,rho[j])
            end
        for j=1:size(mult_res,1)
                println(saida3,mult_res[j])
            end
            close(saida3)
        end

        Adiciona_Vista_Escalar_Gmsh(fmesh, "x", nel, x, Float64(i_ext))
        Adiciona_Vista_Escalar_Gmsh(fmesh2, "xf", nel, xf, Float64(i_ext))
        Adiciona_Vista_Nodal_Tensor_Gmsh(fmesh3, nel, "S", TS, Float64(i_ext))
        Adiciona_Vista_Escalar_Gmsh(fmesh4, "xh", nel, xh, Float64(i_ext))

    end

end


function Imprime_0(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1},
                     valor_fun::Float64, valor_res::Array{Float64,1}, max_ext::Int64, max_int::Int64,
                     tol_ext::Float64, tol_int::Float64, dts::String, nnos::Int64,
                      nel::Int64, ijk::Array{Int64,2}, coord::Array{Float64,2}, to_plot::Array{Float64,1}, TS::Array{Float64,2},
                      vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64)

    # Verifica se precisa imprimir arquivos
    if dts !=  "OFF"

        # Abre os arquivos
        fmesh  = string(dts,"_dens_o.pos")
        fmesh2 = string(dts,"_dens_f.pos")
        fmesh3 = string(dts,"_stress.pos")
        fmesh4 = string(dts,"_dens_h.pos")
        file1  = string(dts,"_log.txt")
        file2  = string(dts,"_monitora.txt")
        file3  = string(dts,"_save.txt")

        # Remove se já existir
        if isfile(fmesh);  rm(fmesh);  end
        if isfile(fmesh2); rm(fmesh2); end
        if isfile(fmesh3); rm(fmesh3); end
        if isfile(fmesh4); rm(fmesh4); end
        if isfile(file1);  rm(file1);  end
        if isfile(file2);  rm(file2);  end
        if isfile(file3);  rm(file3);  end

        # Imprime pro log
        saida  = open(file1,"a")
        println(saida,"\n  LagAug:: ",dts)
        println(saida,"\n  Nelems:: ",nel," | Iter:: ",max_ext,"/",max_int," | Tol:: ",tol_ext,"/",tol_int)
        close(saida)

        # Imprime malha no gmsh
        Inicializa_Malha_Gmsh(fmesh, nnos, nel, ijk, coord, 2)
        Inicializa_Malha_Gmsh(fmesh2, nnos, nel, ijk, coord, 2)
        Inicializa_Malha_Gmsh(fmesh3, nnos, nel, ijk, coord, 2)
        Inicializa_Malha_Gmsh(fmesh4, nnos, nel, ijk, coord, 2)
    end #if

    # Dá o print do console
    println("\n  LagAug:: ",dts)
    println("\n  Nelems:: ",nel," | Iter:: ", max_ext,"/", max_int," | Tol:: ",tol_ext,"/",tol_int)

    # Vai pra função que imprime a externa
    Imprime_Ext(x, rho, mult_res, valor_fun, valor_res, 0, dts, nel,
                                    to_plot, TS, vizi, nviz, dviz, raiof)

end

#
# Carrega o vetor de densidades de um arquivo
#
function Le_Densidades(arquivo::String, nel::Int64)

    # Le arquivo
   texto = readlines(arquivo)

   array = texto[(end-nel):(end-1)]

   # Inicializa o array
   dens = Array{Float64}(undef,nel)

   for i=1:nel

       # Parse da linha
       dados  = split(array[i]," ")[2]

       # Extrai a ultima densidade
       dens[i] = parse(Float64,dados)

   end #i

   # retorna o vetor de densidades
   return dens

end

#
# Carrega itex, rho e mult_res
#
function Le_Ultima(arquivo::String,n_rho::Int64)

    # Le arquivo
    arquivo  = string(arquivo,"_save.txt")
    texto = readlines(arquivo)

    # Número da próxima iteração
   i_ext = parse(Int64,texto[1])+1

   # Número de restrições (tudo menos rhos e i_ext)
   nrest = size(texto,1)-n_rho-1

   # Inicializa
   rho  = Array{Float64,1}(undef,n_rho)
   mult_res = Array{Float64,1}(undef,nrest)

   # Captura valores de rho
   for j=2:(n_rho+1)
       rho[j-1]   = parse(Float64,texto[j])
   end

   # Captura valores dos multplicadores
   for j=(n_rho+2):size(texto,1)
       mult_res[j-n_rho-1] = parse(Float64,texto[j])
   end

   return i_ext,rho,mult_res

end

# Leitura da Y/N
function input(prompt::String="")::String
   print(prompt)
   return chomp(readline())
end
