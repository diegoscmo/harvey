10################################################################################
#####                Saída de Arquivos + Display no Console               ######
################################################################################

#
# Dá o display dos resultados no console e arquivo
#
function Imprime_Int(i_int::Int64, contaev::Int64, norma::Float64)

    @printf("  interno: %d\t evals: %d \t dL: %.3e\n", i_int, contaev, norma)

end

function Imprime_Ext(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1},
                     valor_fun::Float64, valor_res::Array{Float64,1}, i_ext::Int64, csi::Float64,
                     dts::String, nel::Int64, to_plot::Array{Float64,1}, TS::Array{Float64,2},
                     vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64)

    # Recalcula o volume
    xf, =  Filtro_Dens(x, nel, csi, vizi, nviz, dviz, raiof)
    vol = mean(xf)

    # Displau dos valores atuais
    @printf("\n  lagrangian: \t %.4e\t volume:  %.5f ", valor_fun, vol)
    @printf("\n  rho: \t\t rests: \t mults: ")
    @printf("\n  %.4e\t %.4e\t %.4e", rho[1], valor_res[1], mult_res[1])
    @printf("\n  %.4e\t %.4e\t %.4e", rho[2], valor_res[2], mult_res[2])
    @printf("\n  %.4e\t %.4e\t %.4e", rho[3], norm(valor_res[3:end]), norm(mult_res[3:end]))
    @printf("\n")

    # Abre os arquivos de saída
    fmesh  = string("results/",dts,"/densidades_originais.pos")
    fmesh2 = string("results/",dts,"/densidades_filtradas.pos")
    fmesh3 = string("results/",dts,"/analise_tensoes_nos.pos")
    file2  = string("results/",dts,"/monitoramento.txt")
    file3  = string("results/",dts,"/save_file.txt")

    # Plots de valores do arquivo de monitoramento
    saida2 = open(file2,"a")

    # Imprime o cabeçadlho na primeira
    if i_ext == 0
        @printf(saida2, "i_ext     dL        L         Fobj1     Fobj2     Vol       R         maxS      freqs...\n")
    end

    # Imprime os dados dessa iteração
    for j = 1:size(to_plot,1)
        @printf(saida2, "%.3e ", to_plot[j])
    end

    # Fecha o arquivo
    @printf(saida2,"\n")
    close(saida2)

    # Imprime o save_file para retomar
    if isfile(file3); rm(file3); end

    # Imprime só se já houver alguma coisa calculada
    if i_ext > 0

        # Abre o arquivo
        saida3 = open(file3,"a")

        # Imprime o número da iteração
        println(saida3,i_ext)

        # Valor do CSI para Heaviside
        println(saida3,csi)

        # Imprime os valores de rho
        for j=1:size(rho,1)
            println(saida3,rho[j])
        end

        # Imprime os valores dos multiplicadores
        for j=1:size(mult_res,1)
            println(saida3,mult_res[j])
        end

        # Fecha arquivo
        close(saida3)
    end

    # Adiciona as leitura de densidades e tensões
    Adiciona_Vista_Escalar_Gmsh(fmesh, "x", nel, x, Float64(i_ext))
    Adiciona_Vista_Escalar_Gmsh(fmesh2, "xf", nel, xf, Float64(i_ext))
    Adiciona_Vista_Nodal_Tensor_Gmsh(fmesh3, nel, "S", TS, Float64(i_ext))

end #function


function Imprime_0(x::Array{Float64,1}, rho::Array{Float64,1}, mult_res::Array{Float64,1},
                     valor_fun::Float64, valor_res::Array{Float64,1}, max_ext::Int64, max_int::Int64,
                     tol_ext::Float64, tol_int::Float64, dts::String, nnos::Int64,
                      nel::Int64, ijk::Array{Int64,2}, coord::Array{Float64,2}, to_plot::Array{Float64,1}, TS::Array{Float64,2},
                      vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64)

    # Abre os arquivos
    fmesh  = string("results/",dts,"/densidades_originais.pos")
    fmesh2 = string("results/",dts,"/densidades_filtradas.pos")
    fmesh3 = string("results/",dts,"/analise_tensoes_nos.pos")
    file2  = string("results/",dts,"/monitoramento.txt")
    file3  = string("results/",dts,"/save_file.txt")

    # Remove se já existir
    if isfile(file2);  rm(file2);  end
    if isfile(file3);  rm(file3);  end

    # Imprime malha no gmsh
    Inicializa_Malha_Gmsh(fmesh, nnos, nel, ijk, coord, 2)
    Inicializa_Malha_Gmsh(fmesh2, nnos, nel, ijk, coord, 2)
    Inicializa_Malha_Gmsh(fmesh3, nnos, nel, ijk, coord, 2)

    # Dá o print do console
    printstyled("\n  LagAug:: ",dts,"\n",color=:10)
    printstyled("\n  Nelems:: ",nel," | Iter:: ", max_ext,"/", max_int,"\n",color=:10)

    # Vai pra função que imprime a externa
    Imprime_Ext(x, rho, mult_res, valor_fun, valor_res, 0, 0.0, dts, nel,
                                    to_plot, TS, vizi, nviz, dviz, raiof)

end

function Imprime_Caso(dts, Sy, freq, alfa, beta, R_b, A, dini, max_ext, max_int, tol_ext,
           tol_int, rho, rho_max, SP, raiof, vmin, NX, NY, LX, LY,
                       young, poisson, esp, p_dens, presos, forcas, QP, csi, dmax,P,q)

    if isdir(string("results/",dts)) == false
       mkdir(string("results/",dts))
    end

    file = string("results/",dts,"/variaveis_do_caso.txt")
    if isfile(file);  rm(file);  end
    saida = open(file,"a")

    matriz1 = [ "filename:" dts; " " " ";
                "# Parametros da Fobj" " "; "freq" freq; "alfa" alfa; "beta" beta;
                "A   " A; "P   " P; "q   " q; "Rb  " R_b; "Sy  " Sy; "dini" dini; "dmax" dmax; " " " ";
                "# Parametros do Lagrangiano Aumentado" " "; "max_ext" max_ext;
                "max_int" max_int; "tol_ext" tol_ext; "tol_int" tol_int;
                "rho_max" rho_max;  " " " "; "# Penalização inicial" "";]

    matriz2 = [ " " " ";"# Parâmetros da topológica" " "; "SP   " SP; "QP   " QP;
                "raiof" raiof; "vmin " vmin; "nel  " NX*NY]

    for j=1:size(matriz1,1)
        println(saida,matriz1[j,1], " ", matriz1[j,2])
    end
    for j=1:size(rho,1)
        println(saida,"rho$j", " ", rho[j])
    end
    for j=1:size(matriz2,1)
        println(saida,matriz2[j,1], " ", matriz2[j,2])
    end

    close(saida)
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
# Carrega x, itex, rho e mult_res
#
function Le_Ultima(dts, nel::Int64, n_rho::Int64)

   # Nome do save_file
   arquivo = string("results/",dts,"/save_file.txt")

   # Le arquivo
   texto = readlines(arquivo)

    # Número da próxima iteração
   i_ext = parse(Int64,texto[1])+1

   csi   = parse(Float64,texto[2])

   # Número de restrições (tudo menos rhos e i_ext)
   nrest = size(texto,1)-n_rho-2

   # Inicializa
   rho  = Array{Float64,1}(undef,n_rho)
   mult_res = Array{Float64,1}(undef,nrest)

   # Captura valores de rho
   for j=3:(n_rho+2)
       rho[j-2]   = parse(Float64,texto[j])
   end

   # Captura valores dos multplicadores
   for j=(n_rho+3):size(texto,1)
       mult_res[j-n_rho-2] = parse(Float64,texto[j])
   end

   # Arquivo com as densidades_originais
   arq_dens = string("results/",dts,"/densidades_originais.pos")

   # Lê arquivo com densidades x
   x = Le_Densidades(arq_dens, nel)

   return x, i_ext, rho, mult_res, csi

end


# Leitura da Y/N
function input(prompt::String="")::String
   print(prompt)
   return chomp(readline())
end


#
# Verifica se tem uma execução não finalizada
#
function Load_Game(dts, nel, max_ext, max_int)

    # Procura se o save file esta lá
    if isfile(string("results/",dts,"/save_file.txt"))

        # Se estive Vai ficar perguntando até receber Y ou N
        while true

            printstyled("\nUtilizar execução anterior? Y/N",bold=true,color=:yellow)
            conf = input()
            if conf == "Y" || conf ==  "y"

                # Imprime aviso
                printstyled("AVISO: Carregando execução anterior, cuidado com as condições de contorno!\n",color=:yellow)
                printstyled("\n  LagAug:: ",dts,"\n",color=:10)
                printstyled("\n  Nelems:: ",nel," | Iter:: ", max_ext,"/", max_int,"\n",color=:10)
                return true

            elseif conf == "N" || conf =="n"
                return false
            end #if

        end
    end

end #function
