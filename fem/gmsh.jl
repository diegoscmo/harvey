#
# Cria o cabecalho com informacoes da malha
# para posterior adicao de vistas com saidas
#
function Inicializa_Malha_Gmsh(nome_arquivo::String,nnos,nelems,conec,coord,dimensao=2)

    # Abre o arquivo para escrita
    saida = open(nome_arquivo,"w")


    # Cabecalho do gmsh
    println(saida,"\$MeshFormat")
    println(saida,"2.2 0 8")
    println(saida,"\$EndMeshFormat")

    # Nodes
    println(saida,"\$Nodes")
    println(saida,nnos)
    for i=1:nnos
        println(saida,i," ",coord[i,1]," ",coord[i,2]," 0.0 ")
    end
    println(saida,"\$EndNodes")

    # Conectividades
    tipo_elemento = 3
    if dimensao == 3
        tipo_elemento = 5
    end

    println(saida,"\$Elements")
    println(saida,nelems)
    for i=1:nelems
        con = string(i)*" "*string(tipo_elemento)*" 0 "*string(conec[i,1])
        for j=2:size(conec,2)
            con = con * " " * string(conec[i,j])
        end
        println(saida,con)
    end
    println(saida,"\$EndElements")

    # Fecha o arquivo ... por hora
    close(saida)


end # Gera_Malha_Gmsh






#
# Adiciona uma vista escalar nodal a um arquivo (que ja deve ter o cabecalho)
# O vetor que contem os valores nodais deve ter dimensao nnos
#
function Adiciona_Vista_Nodal_Escalar_Gmsh(nome_arquivo::String,nome_vista::String,nnos::Int64,
                                           escalares::Array{Float64,1},tempo::Float64)


    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end


    # Verifica se a dimensao esta correta
    if size(escalares,1)!=nnos
        error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: vetor com escalares deve ter dimensao nnos")
    end

    #
    #
    println(saida,"\$NodeData")
    println(saida,"1")
    println(saida,"\" $nome_vista $tempo\"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"1")
    println(saida,nnos)
    for i=1:nnos
        println(saida,i," ",escalares[i])
    end
    println(saida,"\$EndNodeData")

    # Fecha por hora
    close(saida)

end


#
# Adiciona uma vista escalar a um arquivo (que ja deve ter o cabecalho)
# O vetor que contem os valores centroidais deve ter dimensao nelems
#
function Adiciona_Vista_Escalar_Gmsh(nome_arquivo::String,nome_vista::String,nelems::Int64,
                                           escalares::Array,tempo::Float64)


    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end


    # Verifica se a dimensao esta correta
    if size(escalares,1)!=nelems
        error("ERROR::Adiciona_Vista_Nodal_Escalar_Gmsh:: vetor com escalares deve ter dimensao nnos")
    end

    #
    #
    println(saida,"\$ElementData")
    println(saida,"1")
    println(saida,"\" $nome_vista $tempo\"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"1")
    println(saida,nelems)
    for i=1:nelems
        println(saida,i," ",round(escalares[i],15))
    end
    println(saida,"\$EndElementData")

    # Fecha por hora
    close(saida)

end



#
# Adiciona uma vista vetorial nodal a um arquivo (que ja deve ter o cabecalho)
# O vetor que contem os valores nodais deve ser expandido, ou seja,
# deve ter sido criado por
#
function Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos::Int64,nome_arquivo::String,nome_vista::String,
                                            vetor::Array{Float64},tempo::Float64)


    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Nodal_Vetorial_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end


    # Verifica se a dimensao esta correta
    dim_total = 2*nnos
    if size(vetor,1)!= dim_total
        error("ERROR::Adiciona_Vista_Nodal_Vetorial_Gmsh:: vetor com escalares deve ter dimensao $dim_total")
    end

    #
    #
    println(saida,"\$NodeData")
    println(saida,"1")
    println(saida,"\" $nome_vista\"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"3")
    println(saida,nnos)
    for no=1:nnos
        pos1 = 2*(no-1)+1; val1 = vetor[pos1]
        pos2 = 2*(no-1)+2; val2 = vetor[pos2]
        val3 = 0.0
        println(saida,no," ",val1," ",val2," ",val3 )
    end
    println(saida,"\$EndNodeData")

    # Fecha por hora
    close(saida)

end


#
# Adiciona uma vista tensorial simétrica centroidal a um arquivo (que ja deve ter o cabecalho)
# valores é uma lista com as tres componentes [xx yy xy] e deve ter dimensão nelems
#
function Adiciona_Vista_Centroidal_Tensor2D_Gmsh(nome_arquivo::String,nelems,nome_vista::AbstractString,
                                                    valores,tempo::Float64)


    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Centroidal_Tensor2D_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end


    # Verifica se a dimensao esta correta
    if length(valores)!=3*nelems
        error("ERROR::Adiciona_Vista_Centroidal_Tensor2D_Gmsh:: lista com as componentes deve ter dimensao nelems x [xx yy xy]")
    end

    #
    #
    println(saida,"\$ElementData")
    println(saida,"1")
    println(saida,"\" $nome_vista\"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"9")
    println(saida,nelems)
    for i=1:nelems
        comp = valores[i,:]
        println(saida,i," ",comp[1], " ", comp[3], " ", 0.0, " ", comp[3], " ",comp[2]," ",0.0, " ",0.0, " ",0.0," ",0.0)
    end
    println(saida,"\$EndElementData")

    # Fecha por hora
    close(saida)
end





########################################################################################################################################

#
# Tensoes é uma matrizinha com 4 linhas (pto de gauss) por 3 colunas (tensão)
#

function Mapeia_Nos_ConsistenteQuad(tensoes)

    # Esta matriz foi obtida no maxima e os cálculos estão
    # na pasta de documentação.
    invA  = [1.866025403784438 -0.5 0.1339745962155612 -0.5;
                  -0.5 1.866025403784438 -0.5 0.1339745962155612;
                   0.1339745962155612 -0.5 1.866025403784438 -0.5;
                  -0.5 0.1339745962155612 -0.5 1.866025403784438]


    # e calcula os valores nodais
    S = invA*tensoes

    # Já que estamos nos ocupando da visualização,
    # vamos montar a string para o gmsh
    # Cada nó vai ter que conter uma sequência como a
    # comp[1], " ", comp[3], " ", 0.0, " ", comp[3], " ",comp[2]," ",0.0, " ",0.0, " ",0.0," ",0.0
    saida = " "
    for i=1:4
        saida = saida * " " * string(S[i,1]) * " " * string(S[i,3]) * " 0.0 " * string(S[i,3]) * " " * string(S[i,2]) *  " 0.0 " *  " 0.0 " * " 0.0 "  * " 0.0 "
    end

    return saida
end


# Se processado for true, então os valores já são os nodais. Do contrário,
# extrapolamos os valores para os nós intermente.
function Adiciona_Vista_Nodal_Tensor_Gmsh(nome_arquivo::AbstractString,nelems,nome_vista::AbstractString,
                                          valores,tempo::Float64,processado=false)



    # Tenta abrir o arquivo para append
    saida = try
                open(nome_arquivo,"a")
    catch
        error("ERROR::Adiciona_Vista_Nodal_Tensor_Gmsh:: Nao foi possivel acessar $nome_arquivo")
    end


    # Verifica se a dimensao esta correta

    if length(valores)!=4*3*nelems
        error("ERROR::Adiciona_Vista_Nodal_Tensor_Gmsh:: lista com as componentes deve ter dimensao nelems")
    end

    if processado
        error("ERROR::Adiciona_Vista_Nodal_Tensor_Gmsh::  Ainda não implementado")
    end

    #
    #
    println(saida,"\$ElementNodeData")
    println(saida,"1")
    println(saida,"\" $nome_vista\"")
    println(saida,"1")
    println(saida,tempo)
    println(saida,"3")
    println(saida,"0")
    println(saida,"9")
    println(saida,nelems)
    for i=1:nelems

            # extrai as tensoes deste elemento
            tensoes_elemento = valores[4*(i-1)+1:4*(i-1)+4,:]

            texto = Mapeia_Nos_ConsistenteQuad(tensoes_elemento)
            println(saida,i," 4 ", texto)

    end
    println(saida,"\$EndElementNodeData")

    # Fecha por hora
    close(saida)

end
