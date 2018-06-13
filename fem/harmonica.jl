################################################################################
#####                          Análise Harmônica                          ######
################################################################################

#
# Calcula os <numero_modos> menores modos e frequências naturais da estrutura
#
function Analise_Modal(numero_modos,K0,M0,nel, csi, nnos,ijk,ID,coord,vizi, nviz, dviz, raiof,
               x,simp=3.0,vmin=1E-3,dts="")

    # Corrige as densidades
    xf, = Filtro_Dens(x, nel, csi, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nel,ijk,ID,K0,M0,simp,vmin)

    # Utilizamos o ARPACK
    AVL = IterativeEigensolvers.eigs(K,M,nev=numero_modos,which=:SM,ritzvec=true)

    # Converte as frequências para Hertz
    auto_vals = sqrt.(real(AVL[1]))/2.0/pi
    auto_vets = convert.(Float64,AVL[2])

    # Vamos gravar os modos para um arquivo de visualização no gmsh
    if dts != ""
        if csi == 0.0
            arquivo_saida = string("results/",dts,"/f_modalans.pos")
        else
            arquivo_saida = string("results/",dts,"/h_modalans.pos")
        end

        Inicializa_Malha_Gmsh(arquivo_saida,nnos,nel,ijk,coord)

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
# Calcula os <numero_modos> menores modos e frequências naturais da
# Pra rodar externo!
#
function Harmonica(dts,freq,x,csi,nnos,nel,coord,ijk,ID,vizi,nviz,dviz,raiof,alfa,beta,F,K0,M0,vmin,simp=3.0)

    # Corrige as densidades
    xf, = Filtro_Dens(x, nel, csi, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nel,ijk,ID,K0,M0,simp,vmin)

    w = 2.0*pi*freq

    # Monta a matriz dinâmica
    KD = sparse(K - (w^2.0)*M + w*im*(alfa*M + beta*K))

    # Soluciona o sistema... como não é Hermitiana, vamos com LU
    UD = vec(lu(KD)\F)

   # Expande o modo
   desloc = real(Expande_Vetor(UD, nnos, ID))

   return desloc

end



#
# Calcula os <numero_modos> menores modos e frequências naturais da estrutura
#
function Analise_Harmonica(freq,K0,M0,nel, csi, nnos,ijk,ID,coord,vizi, nviz, dviz, raiof, alfa,beta,F,
               x,simp=3.0,vmin=1E-3,dts="",i_ext=1)

    # Corrige as densidades
    xf, = Filtro_Dens(x, nel, csi, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf

    # Monta as matrizes globais de massa e de rigidez
    K,M = Global_KM(xc,nel,ijk,ID,K0,M0,simp,vmin)

    w = 2.0*pi*freq

    # Monta a matriz dinâmica
    KD = sparse(K - (w^2.0)*M + w*im*(alfa*M + beta*K))

    # Soluciona o sistema... como não é Hermitiana, vamos com LU
    UD = vec(lu(KD)\F)

    # Vamos gravar os modos para um arquivo de visualização no gmsh
    if dts != ""

       if csi == 0.0
           arquivo_saida = string("results/",dts,"/f_topology.pos")
       else
           arquivo_saida = string("results/",dts,"/h_topology.pos")
       end
       # Expande o modo
       desloc = real(Expande_Vetor(UD, nnos, ID))

       # Grava no arquivo
       Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos,arquivo_saida," $i_ext - $freq Hz", desloc,convert(Float64,freq))
    end

    return

end

#
# Calcula os <numero_modos> menores modos e frequências naturais da estrutura
#
function Varredura_dL(fvarredura,x, rho, mult_res, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                         F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                         Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

    lista_w = collect(fvarredura)

    arquivo_saida = string("results/",dts,"/analise_derivada.pos")

    Inicializa_Malha_Gmsh(arquivo_saida,nnos,nel,ijk,coord)

    @showprogress "  gerando derivadas... " for freq in lista_w

        dL, = F_Obj(x, rho, mult_res, 3, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                 F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                 Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

        norma = norm(dL)
        dir = -dL/norma

                 # Grava no arquivo
         Adiciona_Vista_Escalar_Gmsh(arquivo_saida, "$freq", nel, dir, Float64(freq))

    end

    return

end

#
# Calcula os <numero_modos> menores modos e frequências naturais da estrutura
#
function Varredura_dL_ele(ele,fvarredura,x, rho, mult_res, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                         F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                         Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

    lista_w = collect(fvarredura)

    file = string("results/",dts,"/varredura_derivada.txt")

    if isfile(file);  rm(file);  end
    saida = open(file,"w")

    @showprogress "  varrendo derivadas... " for freq in lista_w

        dL, = F_Obj(x, rho, mult_res, 3, nnos, nel, ijk, coord, ID, K0, M0, SP, vmin,
                                 F, NX, NY, vizi, nviz, dviz, raiof, valor_0,
                                 Sy, freq, alfa, beta, A, Ye, CBA, QP, csi, dmax,nos_viz,dts,P,q)

        norma = norm(dL)
        dir = -dL/norma

        dele = dir[ele]

        # Grava no arquivo
        println(saida,"$freq   $dele")

    end
    close(saida)

    return

end





#
# Realiza uma varredura em frequência (Harmônica)
#
function Varredura_Graph(dts, fvarredura,
                   nel::Int64, csi, nnos::Int64,
                   ijk::Array{Int64,2},coord::Array{Float64,2},vizi, nviz, dviz, raiof,
                   ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                   F, alfa::Float64, beta::Float64,
                   x,P,q,nos_viz, heavi, simp=3.0, vmin=1E-3)

    # Monta as matrizes globais de massa e de rigidez
    xf, = Filtro_Dens(x, nel, csi, vizi, nviz, dviz, raiof)
    xc = @. vmin+(1.0-vmin)*xf
    K,M = Global_KM(xc,nel,ijk,ID,K0,M0,simp,vmin)
    C   = (alfa*M + beta*K)

    # Gera lista de frequências angulares para varredura
    lista_w = 2.0*pi.*collect(fvarredura)

    # Diferencia heaviside de filtro
    if !heavi
        pfix = "f"
    else
        pfix = "h"
    end

    # Abre os arquivos de saída
    file  = string("results/",dts,"/",pfix,"_zomegas.txt")
    if isfile(file);  rm(file);  end
    saida_w = open(file,"w")

    file  = string("results/",dts,"/",pfix,"_zpotenat.txt")
    if isfile(file);  rm(file);  end
    saida_p = open(file,"w")

    file  = string("results/",dts,"/",pfix,"_znormas.txt")
    if isfile(file);  rm(file);  end
    saida_n = open(file,"w")


    # Loop pela lista de frequências
    @showprogress "  gerando graficos... " for w in lista_w

        # Monta a matriz dinâmica
        KD = sparse(K - (w^2.0)*M + w*im*C)

        # Soluciona o sistema... como não é Hermitiana, vamos com LU
        UD = vec(lu(KD)\F)

        # Verifica o condicionamento na frequência
    #    cond = Checa_Cond(KD,1E-4,1000)

        # Calcula o valor da Pot At
        Pa = 0.5*w*real(im*dot(F,UD))

        # Calcula o valor de r
    #    Ec = w^2.0*real(adjoint(UD)*M*UD)
    #    Ep = real(adjoint(UD)*K*UD)
    #    R  = Ec/Ep

        # Calcula o valor da Norma
        phi = fzinha_MOD(UD,nnos,nos_viz,xc,ID,P,q)

        # Frequencia em Hz
        f_Hz = w/(2.0*pi)

        # Grava a resposta no arquivo de monitoramento
        println(saida_w,f_Hz)
        println(saida_p,Pa)
    #    println(saida_r,R)
        println(saida_n,phi)
    #    println(saida_c,cond)

    end # w in lista_w

    close(saida_w)
    close(saida_p)
    #close(saida_r)
    close(saida_n)
    #close(saida_c)
end


function fzinha_MOD(UD,nnos,nos_viz,xc,ID,P,q)

    # Gera um vetor com as "normas nodais de densidade"
    aj = zeros(nnos)
    for j=1:nnos

        # Recupera os elementos vizinhos deste nó
        viz = nos_viz[j]
        vaj = 0.0
        for v in viz
            vaj += xc[v]^q
        end
        aj[j] = vaj^(1.0/q)
    end

    # Para montar o a_N
    a_N = 0.0
    for j=1:nnos

       # Recupera aj para este nó
       vaj = aj[j]

       # Recupera os gls do nós
       gl1 = ID[j,1]
       gl2 = ID[j,2]

       # Primeiro gl
       if gl1 != 0
           ui = UD[gl1]
           a_N += real(conj(ui)*vaj*ui)^(P/2.0)
       end
       # Segundo gl
       if gl2 != 0
           ui = UD[gl2]
           a_N += real(conj(ui)*vaj*ui)^(P/2.0)
       end

    end #if

    # Calcula a norma
    Nor1 = a_N^(1.0/P)

    return Nor1

end


function Varredura_Cond(dts, fvarredura,
                   nel::Int64,csi, nnos::Int64,
                   ijk::Array{Int64,2},coord::Array{Float64,2},vizi, nviz, dviz, raiof,
                   ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                   F, alfa::Float64, beta::Float64,
                   x, tol,miter,simp=3.0, vmin=1E-3)

   # Monta as matrizes globais de massa e de rigidez
   xf, = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)
   xc = @. vmin+(1.0-vmin)*xf
   K,M = Global_KM(xc,nel,ijk,ID,K0,M0,simp,vmin)

   # Gera lista de frequências angulares para varredura
   lista_w = 2.0*pi.*collect(fvarredura)

   # Abre os arquivos de saída
   file  = string("results/",dts,"/varredura_omegac.txt")
   if isfile(file);  rm(file);  end
   saida_w = open(file,"w")

   file  = string("results/",dts,"/varredura_condic.txt")
   if isfile(file);  rm(file);  end
   saida_c = open(file,"w")
   # Loop pela lista de frequências
   @showprogress "  gerando condicionamento... " for w in lista_w

       # Monta a matriz dinâmica
       KD = sparse(K - (w^2.0)*M)# + w*im*(alfa*M + beta*K))

       #KDM = Matrix(KD)
       # Verifica o condicionamento na frequência
       #cdn = Checa_Cond(KD,tol,miter)
       #cdn = cond(Matrix(KD))
       cdn = det(KD)

       # Frequencia em Hz
       f_Hz = w/(2.0*pi)

       # Grava a resposta no arquivo de monitoramento
       println(saida_w,f_Hz)
       println(saida_c,cdn)

   end # w in lista_w
   close(saida_w)
   close(saida_c)
end
