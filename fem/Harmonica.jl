################################################################################
#####                          Análise Harmônica                          ######
################################################################################

function Analise_Harmonica(freq::Float64, arquivo_saida, nelems::Int64, nnos::Int64,
                           ijk::Array{Int64,2},coord::Array{Float64,2},
                    ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                    nos_forcas, ngl_efetivos, alfa::Float64, beta::Float64,
                    densidades=[], simp=1.0, vminimo=1E-3)

    println("\n Iniciando a análise harmônica para $freq Hz ")

    # Se não foram informadas densiades...
    if densidades==[]
        densidades=ones(nelems)
    end

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(densidades,nelems,ijk,ID,K0,M0,simp,vminimo)

    # Monta o vetor de forças globais
    F = Global_F(ID,nos_forcas,ngl_efetivos)

    # Converte a frequência para radianos/s
    w = 2.0*pi*freq

    # Monta a matriz dinâmica
    KD = K - (w^2.0)*M + w*im*(alfa*M + beta*K)

    # Soluciona o sistema...como pode ser ou não posdef, vamos usar LU
    UD = vec(lufact(KD)\F)

    # Grava para visualização no gmsh
    Inicializa_Malha_Gmsh(arquivo_saida,nnos,nelems,ijk,coord)

    # Grava a parte real do vetor de deslocamentos
    desloc = Expande_Vetor(real.(UD), nnos, ID)

    # Grava no arquivo
    Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos, arquivo_saida," Freq $(freq) Hz",
                                      desloc,0.0)


    println("Análise Harmonica terminada com sucesso")
end


#
# Realiza uma varredura em frequência (Harmônica)
#
function Varredura(arquivo_saida, fvarredura,nos_monitor,
                   nelems::Int64, nnos::Int64,
                   ijk::Array{Int64,2},coord::Array{Float64,2},
                   ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                   nos_forcas, ngl_efetivos, alfa::Float64, beta::Float64,
                   densidades=[], simp=1.0, vminimo=1E-3)


    # Gera lista de frequências angulares para varredura
    lista_w = 2*pi.*collect(fvarredura)

    println("\n Iniciando varredura harmônica de $(minimum(lista_w)/(2*pi))
                            até $(maximum(lista_w)/(2*pi))Hz ")

    # Numero de termos na lista de varreduras
    numero_w = length(lista_w)

    println("\n Total de análises:  $(numero_w)  em ")

    # Abre um arquivo para monitoramento da resposta
    # de alguns gls
    arquivo_monitor = open("monitor.txt","w")


    # Se não foram informadas densiades...
    if densidades==[]
        densidades=ones(nelems)
    end

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(densidades,nelems,ijk,ID,K0,M0,simp,vminimo)

    # Monta o vetor de forças globais
    F = Global_F(ID,nos_forcas,ngl_efetivos)

    # Grava para visualização no gmsh
    Inicializa_Malha_Gmsh(arquivo_saida,nnos,nelems,ijk,coord)

    # Vamos fazer uma manha para inicialização dos fatores LDL
    # na matriz C que será reaproveitada em ldltfact!(C,KD) no
    # loop principal...
    w = lista_w[1]
    KD = Hermitian(K - (w^2.0)*M + w*im*(alfa*M + beta*K))
    C = ldltfact(KD)

    # Loop pela lista de frequências
    @time for w in lista_w

        # Monta a matriz dinâmica
        KD = Hermitian(K - (w^2.0)*M + w*im*(alfa*M + beta*K))

        # Soluciona o sistema...como pode ser ou não posdef, vamos usar LDL
        ldltfact!(C,KD)
        UD = vec(C\F)

        # Grava a resposta no arquivo de monitoramento
        for i=1:size(nos_monitor,1)

            no = nos_monitor[i,1]
            gl = nos_monitor[i,2]
            glg = ID[no,gl]
            ud = UD[glg]
            pr = real(ud)
            pc = imag(ud)
            pn = norm(ud)
            print(arquivo_monitor," $no $gl $pr $pc $pn ")

        end
        print(arquivo_monitor,"\n")

        # Grava a parte real do vetor de deslocamentos para visualização
        desloc = Expande_Vetor(real.(UD), nnos, ID)

        # Grava no arquivo
        f_Hz = w/(2.0*pi)
        Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos, arquivo_saida," Freq $(f_Hz) Hz",
                                           desloc,w)


    end # w in lista_w

    close(arquivo_monitor)
    println(" Varredura Harmonica terminada com sucesso")
end

#
# Calcula os <numero_modos> menores modos e frequências naturais da estrutura
#
function Analise_Modal(numero_modos,K0,M0,nelems,nnos,ijk,ID,coord,arquivo_saida,
               densidades=[],simp=3.0,vminimo=1E-3)

    println("\n Iniciando a análise modal para $(numero_modos) modos ")

    # Se não foram informadas densiades...
    if densidades==[]
        densidades=ones(nelems)
    end

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(densidades,nelems,ijk,ID,K0,M0,simp,vminimo)

    # Utilizamos o ARPACK
    @time AVL = IterativeEigensolvers.eigs(K,M,nev=numero_modos,which=:SM,ritzvec=true)

    # Converte as frequências para Hertz
    auto_vals = sqrt.(real(AVL[1]))/2.0/pi
    auto_vets = convert.(Float64,AVL[2])

    # Vamos gravar os modos para um arquivo de visualização no gmsh
    Inicializa_Malha_Gmsh(arquivo_saida,nnos,nelems,ijk,coord)

    # Para cada modo, convertemos para FULL e geramos a visualização
    for modo=1:numero_modos

           # Expande o modo
           desloc = Expande_Vetor(auto_vets[:,modo], nnos, ID)

           # Grava no arquivo
           Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos,arquivo_saida," Freq $(auto_vals[modo]) Hz",
                                             desloc,convert(Float64,modo))

    end #modo

    println("Problema modal terminado com sucesso")

end
