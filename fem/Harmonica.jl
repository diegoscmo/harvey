################################################################################
#####                          Análise Harmônica                          ######
################################################################################

#
# Calcula os <numero_modos> menores modos e frequências naturais da estrutura
#
function Analise_Modal(numero_modos,K0,M0,nelems,nnos,ijk,ID,coord,vizi, nviz, dviz, raiof,
               x,simp=3.0,vmin=1E-3,dts="",i_ext=0)

    # Corrige as densidades
    xf = Filtro_Dens(x, nelems, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nelems,ijk,ID,K0,M0,simp,vmin)

    # Utilizamos o ARPACK
    AVL = IterativeEigensolvers.eigs(K,M,nev=numero_modos,which=:SM,ritzvec=true)

    # Converte as frequências para Hertz
    auto_vals = sqrt.(real(AVL[1]))/2.0/pi
    auto_vets = convert.(Float64,AVL[2])

    # Vamos gravar os modos para um arquivo de visualização no gmsh
    if dts != ""
        arquivo_saida = string("results\\",dts,"\\",dts,"_modal",i_ext,".pos")

        if isfile(arquivo_saida);  rm(arquivo_saida);  end

        Inicializa_Malha_Gmsh(arquivo_saida,nnos,nelems,ijk,coord)

        # Para cada modo, convertemos para FULL e geramos a visualização
        for modo=1:numero_modos

           # Expande o modo
           desloc = Expande_Vetor(auto_vets[:,modo], nnos, ID)

           # Grava no arquivo
           Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos,arquivo_saida," Freq $(auto_vals[modo]) Hz",
                                             desloc,convert(Float64,modo))
        end #modo
    end

    return auto_vals

end


#
# Calcula os <numero_modos> menores modos e frequências naturais da estrutura
#
function Analise_Harmonica(freq,K0,M0,nelems,nnos,ijk,ID,coord,vizi, nviz, dviz, raiof, alfa,beta,F,
               x,simp=3.0,vmin=1E-3,dts="",i_ext=1)

    # Corrige as densidades
    xf = Filtro_Dens(x, nelems, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nelems,ijk,ID,K0,M0,simp,vmin)

    w = 2.0*pi*freq

    # Monta a matriz dinâmica
    KD = sparse(K - (w^2.0)*M + w*im*(alfa*M + beta*K))

    # Soluciona o sistema... como não é Hermitiana, vamos com LU
    UD = vec(lufact(KD)\F)

    # Vamos gravar os modos para um arquivo de visualização no gmsh
    if dts != ""

        arquivo_saida = string("results\\",dts,"\\",dts,"_harmonica.pos")

        if i_ext == 1

            if isfile(arquivo_saida);  rm(arquivo_saida);  end
            Inicializa_Malha_Gmsh(arquivo_saida,nnos,nelems,ijk,coord)

        end

       # Expande o modo
       desloc = real(Expande_Vetor(UD, nnos, ID))

       # Grava no arquivo
       Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos,arquivo_saida," $i_ext - $freq Hz", desloc,convert(Float64,freq))
    end

    return

end













#
# Realiza uma varredura em frequência (Harmônica)
#
function Varredura_Graph(dts, fvarredura,
                   nelems::Int64, nnos::Int64,
                   ijk::Array{Int64,2},coord::Array{Float64,2},vizi, nviz, dviz, raiof,
                   ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                   F, alfa::Float64, beta::Float64,
                   x, simp=3.0, vminimo=1E-3)

    # Corrige as densidades
    xf = Filtro_Dens(x, nelems, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Gera lista de frequências angulares para varredura
    lista_w = 2.0*pi.*collect(fvarredura)

    # Abre o arquivo de saída
    file  = string("results\\",dts,"\\",dts,"_graph.txt")
    if isfile(file);  rm(file);  end
    saida = open(file,"w")

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nelems,ijk,ID,K0,M0,simp,vminimo)

    # loop principal...
    w = lista_w[1]

    # Loop pela lista de frequências
    @time for w in lista_w

        # Monta a matriz dinâmica
        KD = sparse(K - (w^2.0)*M + w*im*(alfa*M + beta*K))

        # Soluciona o sistema... como não é Hermitiana, vamos com LU
        UD = vec(lufact(KD)\F)

        # Calcula o valor da Pot At
        #FD = abs(real(dot(F,UD)))
        FD = 0.5*w*real(im*dot(F,UD))

        # Frequencia em Hz
        f_Hz = w/(2.0*pi)

        # Grava a resposta no arquivo de monitoramento
        println(saida,"$f_Hz $FD")

    end # w in lista_w

    close(saida)
    println(" FRF atualizada ")
end

function Varredura_R(dts, fvarredura,
                   nelems::Int64, nnos::Int64,
                   ijk::Array{Int64,2},coord::Array{Float64,2},vizi, nviz, dviz, raiof,
                   ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                   F, alfa::Float64, beta::Float64,
                   x, simp=3.0, vminimo=1E-3)

    # Corrige as densidades
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Gera lista de frequências angulares para varredura
    lista_w = 2.0*pi.*collect(fvarredura)

    # Abre o arquivo de saída
    file  = string("results\\",dts,"\\",dts,"_R.txt")
    if isfile(file);  rm(file);  end
    saida = open(file,"w")

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nelems,ijk,ID,K0,M0,simp,vminimo)

    # loop principal...
    w = lista_w[1]

    # Loop pela lista de frequências
    @time for w in lista_w

        # Monta a matriz dinâmica
        KD = sparse(K - (w^2.0)*M + w*im*(alfa*M + beta*K))

        # Soluciona o sistema... como não é Hermitiana, vamos com LU
        UD = vec(lufact(KD)\F)

        # Monta pro R
        Ec = w^2.0*real(adjoint(UD)*M*UD)
        Ep = real(adjoint(UD)*K*UD)

        R  = Ec/Ep

        # Frequencia em Hz
        f_Hz = w/(2.0*pi)

        # Grava a resposta no arquivo de monitoramento
        println(saida,"$f_Hz $R")

    end # w in lista_w

    close(saida)
    println(" R atualizado ")
end


#
# Realiza uma varredura em frequência (Harmônica)
#
function Varredura_Phi(dts, fvarredura,
                   nelems::Int64, nnos::Int64,
                   ijk::Array{Int64,2},coord::Array{Float64,2},vizi, nviz, dviz, raiof,
                   ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                   F, alfa::Float64, beta::Float64, x, P=12.0,
                   simp=3.0, vminimo=1E-3)

    # Corrige as densidades
    xf = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Gera lista de frequências angulares para varredura
    lista_w = 2.0*pi.*collect(fvarredura)

    # Abre o arquivo de saída
    file  = string("results\\",dts,"\\",dts,"_phi_",P,".txt")
    if isfile(file);  rm(file);  end
    saida = open(file,"w")

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nelems,ijk,ID,K0,M0,simp,vminimo)

    # loop principal...
    w = lista_w[1]

    # Loop pela lista de frequências
    @time for w in lista_w

        # Monta a matriz dinâmica
        KD = sparse(K - (w^2.0)*M + w*im*(alfa*M + beta*K))

        # Soluciona o sistema... como não é Hermitiana, vamos com LU
        UD = vec(lufact(KD)\F)

        phi = fzinha(UD,F,P)

        # Frequencia em Hz
        f_Hz = w/(2.0*pi)

        # Grava a resposta no arquivo de monitoramento
        println(saida,"$f_Hz $phi")

    end # w in lista_w

    close(saida)
    println(" PHI atualizada ")
end

function fzinha(UD,F,P)

    phi = 0.0

    #norma_F =  sqrt(dot(F,F))

    #proj_UF = ((dot(UD,F))/norma_F)*F

    for u in UD

        phi += (conj(u)*u)^(P/2.0)

    end #for u

    phi = real(phi^(1.0/P))

    return phi

end
