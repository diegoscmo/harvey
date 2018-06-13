################################################################################
#####                Saída de Arquivos + Display no Console               ######
################################################################################

#
# Dá o display dos resultados no console e arquivo
#
function Imprime_Ext(x::Array{Float64,1}, rho::Array{Float64,1}, mu_res::Array{Float64,1},
                     valor_fun::Float64, valor_res::Array{Float64,1}, i_ext::Int64, csi::Float64,
                     dts::String, nel::Int64, to_plot::Array{Float64,1}, TS::Array{Float64,2}, harm,
                     vizi::Array{Int64,2}, nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64,
                     heavi::Bool, nnos::Int64, ijk, coord, loadf=false)

     # Atualiza a trava
     Lock_Game(dts)

     # Escolhendo nome de filtro ou heaviside
     if !heavi
         pfix = "f"
         disp = "Filtro"
     else
         pfix = "h"
         disp = "Heaviside"
     end

     # Nomeia os arquivos
     dens_file = string("results/",dts,"/",pfix,"_densities.pos")
     topo_file = string("results/",dts,"/",pfix,"_topology.pos")
     moni_file = string("results/",dts,"/",pfix,"_monitora.txt")

     # Se nao for carregar uma execucao anterior, limpa os arquivos se existirem
     if i_ext == 0

            # Remove se já existir
            if isfile(moni_file);  rm(moni_file);  end

            # Inicializa a malha para as densidades
            Inicializa_Malha_Gmsh(dens_file, nnos, nel, ijk, coord, 2)

     end

    # Mini-exibição a cada 5 iterações, ou na primeira depois do load
    if i_ext%5 == 0 || to_plot[1] == 0.0
        printstyled("\n  LagAug:: ",dts," :: ",disp,"\n",color=:10)
    end

    # Recalcula o volume
    xf, = Filtro_Dens(x, nel, csi, vizi, nviz, dviz, raiof)
    vol = mean(xf)

    # Displau dos valores atuais
    @printf("\n  lagrangian: \t %.4e\t volume:  %.5f\tcsi:  %.1f ", valor_fun, vol, csi)
    @printf("\n  rho: \t\t rests: \t mults: ")
    @printf("\n  %.4e\t %.4e\t %.4e", rho[1], valor_res[1], mu_res[1])
    @printf("\n  %.4e\t %.4e\t %.4e", rho[2], valor_res[2], mu_res[2])
    @printf("\n  %.4e\t %.4e\t %.4e", rho[3], norm(valor_res[3:end]), norm(mu_res[3:end]))
    @printf("\n")

    # Se nao tiver carregando um ja iniciado, printa pros arquivos
    if !loadf

        # Plots de valores do arquivo de monitoramento
        open_moni = open(moni_file,"a")

        # Imprime o cabeçadlho na primeira
        if i_ext == 0
            @printf(open_moni, "i_ext     dL        L         Fobj1     Fobj2     Vol       R         maxS      freqs...\n")
        end

        # Imprime os dados dessa iteração
        for j = 1:size(to_plot,1)
            @printf(open_moni, "%.3e ", to_plot[j])
        end

        # Fecha o arquivo
        @printf(open_moni,"\n"); close(open_moni)

        # Adiciona a densidade filtrada
        Adiciona_Vista_Escalar_Gmsh(dens_file, "xf", nel, xf, Float64(i_ext))

        # Deleta o topo_file antigo e adiciona densidade, harmonica e tensao
        Inicializa_Malha_Gmsh(topo_file, nnos, nel, ijk, coord, 2)
        Adiciona_Vista_Escalar_Gmsh(topo_file, "xf", nel, xf, Float64(i_ext))
        Adiciona_Vista_Nodal_Vetorial_Gmsh(nnos, topo_file,"H", harm, Float64(i_ext))
        Adiciona_Vista_Nodal_Tensor_Gmsh(topo_file, nel, "S", TS, Float64(i_ext))

        # Salva se nao for a primeira iteracao
        if i_ext > 0
            Save_Game(x, i_ext, rho, mu_res, csi, dts, heavi)
        end
    end

end #function

#
#   Para iteracoes internas..
#
function Imprime_Int(i_int::Int64, contaev::Int64, norma::Float64)

    @printf("  interno: %d\t evals: %d \t dL: %.3e\n", i_int, contaev, norma)

end

#
#   Grava variaveis para proximas execucoes
#
function Save_Game(x, i_ext, rho, mu_res, csi, dts, heavi)

    # Determina o arquivo para o save
    if heavi == false
        save_file = string("results/",dts,"/f_savefile.txt")
    else
        save_file = string("results/",dts,"/h_savefile.txt")
    end

    # Remove o ultimo save_file
    if isfile(save_file); rm(save_file); end

    # Abre o arquivo
    open_save = open(save_file,"a")

    # Imprime as variaveis de interesse
    println(open_save, i_ext)
    println(open_save, csi)
    println(open_save, size(rho,1))
    println(open_save, size(mu_res,1))
    println(open_save, size(x,1))

    # Imprime os valores de rho, mu_res e x
    for j=1:size(rho,1);    println(open_save,rho[j]);     end
    for j=1:size(mu_res,1); println(open_save,mu_res[j]);  end
    for j=1:size(x,1);      println(open_save,x[j]);       end

    # Fecha arquivo
    close(open_save)

end


#
#   Verifica se tem um save_file e determina se precisa executar filtro ou heaviside
#
function Check_Game(dts, sub, max_ext, heavi)

    # Nome dos arquivos para tentar carregar
    fil_save = string("results/",dts,"/f_savefile.txt")
    hev_save = string("results/",dts,"/",sub,"/h_savefile.txt")

    # Se nao for heaviside, verifica o arquivo e retorna loadf caso tenha algo
    if !heavi && isfile(fil_save)

        loadf = true

    # Verifica se tem um save heaviside
    elseif heavi && isfile(hev_save)

        loadf = true

    # Se nao achou nada, entao nao precisa carregar
    else
        loadf =  false

    end

    # Se vai ter qualquer tipo de load, o usuario recebe um tempo pra cancelar
    if loadf
        printstyled("\n Sobrescrevendo valores anteriores, CTRL+C para cancelar (5s)",bold=true,color=:yellow)
        printstyled("\n Se deseja iniciar do zero, delete manualmente o save_file\n",bold=true,color=:yellow)

        sleep(5.0)
    end

    # SE for heaviside, ja atualiza o dts
    if heavi
        dts = string(dts,"/",sub)
    end

    return loadf, dts

end #function

#
#   Carrega x, itex, csi, rho e mu_res de uma execucao anterior
#
function Load_Game(dts0, sub, heavi)

   # Nome do save_file
   arquivo = string("results/",dts0,"/f_savefile.txt")

   # Se for pra heaviside, aí muda o arquivo
   if heavi == true
       arquivo = string("results/",dts0,"/",sub,"/h_savefile.txt")
   end

   # Le arquivo
   texto = readlines(arquivo)

    # Número da ultima iteração, csi, penalizacoes, nrest, nel
   i_ext = parse(Int64,texto[1])
   csi   = parse(Float64,texto[2])
   nrho  = parse(Int64,texto[3])
   nrest = parse(Int64,texto[4])
   nel   = parse(Int64,texto[5])

   # Captura valores de rho
   rho = Array{Float64,1}(undef,nrho)
   ir  = 5
   for j = 1:nrho
       rho[j]    = parse(Float64,texto[j+ir])
   end

   # Captura valores dos multplicadores
   mu_res = Array{Float64,1}(undef,nrest)
   im     = ir + nrho
   for j = 1:nrest
       mu_res[j] = parse(Float64,texto[j+im])
   end

   # Captura valores de x_original
   x      = Array{Float64,1}(undef,nel)
   ix = im + nrest
   for j=1:nel
       x[j]      = parse(Float64,texto[j+ix])
   end

   return x, i_ext, rho, mu_res, csi

end

#
#   Recebe os parametros, cria a pasta e o casefile, retorna o nome do arquivo
#
function Name_Game(tipoan,freq,alfa,beta,A,P,q,Sy,dini,dmax,QP,csi0,csim,rho1,rho2,rho3,
                       raiof,max_fil,max_hev,max_int,tol_ext,tol_int,rho_max,heavi=true)

        # Verifica o tipo de analise pro nome do arquivo
        if tipoan == 1;  tiponome = "NG_"; end
        if tipoan == 2;  tiponome = "NL_"; end
        if tipoan == 3;  tiponome = "NGS"; end
        if tipoan == 4;  tiponome = "NLS"; end

        # Verifica se é minimizacao ou maximizacao
        if A < 0.0
            minmax = "_MAX_"
        else
            minmax = "_MIN_"
        end

        # Valor de beta, verifica se é proporcional
        if beta == 0.1/(2.0*pi*freq)
            beta = "W"
        end

        # Converte a frequencia pra omitir o zero
        ifreq = Int64(freq)

        # COloca dados de tensao so em 3 e 4
        if tipoan != 3 && tipoan != 4
            #  Agrupa para criar algo do tipo "NG__MIN_f180a08p20q20s4E20vi05vm05q15bWc10s12" # b1E-8
            dts = string(tiponome,minmax,"f",ifreq,"a",A,"p",P,"q",q,"vi",dini,"vm",dmax,"b",beta)
        else
            dts = string(tiponome,minmax,"f",ifreq,"a",A,"p",P,"q",q,"vi",dini,"vm",dmax,"b",beta,"sy",Sy,"q",QP)
        end

        # Nome da pasta interna, heaviside
        sub = string("c",csi0,"cs",csim)

        # Tira os pontos
        dts = replace(dts,"."=>"")
        sub = replace(sub,"."=>"")

        if !heavi
            sub = ""
        end

        # Cria diretorio e escreve o case_file
        New_Game(dts,sub,tipoan,freq,alfa,beta,A,P,q,Sy,dini,dmax,QP,csi0,csim,rho1,rho2,rho3,
                               raiof,max_fil,max_hev,max_int,tol_ext,tol_int,rho_max)

       return dts,sub

end

function New_Game(dts,sub,tipoan,freq,alfa,beta,A,P,q,Sy,dini,dmax,QP,csi0,csim,rho1,rho2,rho3,
                       raiof,max_fil,max_hev,max_int,tol_ext,tol_int,rho_max)

    # Cria o diretorio se nao tiver ainda
    if !isdir(string("results/",dts))
        mkdir(string("results/",dts))
    end

    # Escreve a casefile la dentro
    case_file = string("results/",dts,"/0_casefile.txt")

    # Remove se ja tem
    if isfile(case_file); rm(case_file); end

    # Abre
    open_case = open(case_file,"a")

    # Salva todas as variaveis no arquivo
    println(open_case,"tipoan:  ",tipoan)
    println(open_case,"freq:    ",freq)
    println(open_case,"alfa:    ",alfa)
    println(open_case,"beta:    ",beta)
    println(open_case,"A:       ",A)
    println(open_case,"P:       ",P)
    println(open_case,"q:       ",q)
    println(open_case,"dini:    ",dini)
    println(open_case,"dmax:    ",dmax)
    println(open_case,"Sy:      ",Sy)
    println(open_case,"QP:      ",QP)
    println(open_case,"rho1:    ",rho1)
    println(open_case,"rho2:    ",rho2)
    println(open_case,"rho3:    ",rho3)
    println(open_case,"raiof:   ",raiof)
    println(open_case,"max_fil: ",max_fil)
    println(open_case,"max_int: ",max_int)
    println(open_case,"tol_ext: ",tol_ext)
    println(open_case,"tol_int: ",tol_int)
    println(open_case,"rho_max: ",rho_max)

    # Fecha
    close(open_case)

    # Cria tambem o diretorio interno caso seja de interesse
    if sub !=""
        if !isdir(string("results/",dts,"/",sub))
            mkdir(string("results/",dts,"/",sub))
        end
        # Escreve o case do heavi tambem
        heav_file = string("results/",dts,"/",sub,"/0_heavifile.txt")
        if isfile(heav_file); rm(heav_file); end

        heav_case = open(heav_file,"a")

        # Salva as do heavi tambem
        println(heav_case,"max_hev: ",max_hev)
        println(heav_case,"csi0:    ",csi0)
        println(heav_case,"csim:    ",csim)

        close(heav_case)
    end

end

#
#   Busca o case_file e devolve as variaveis
#
function Read_Game(dts,sub="")

    # Busca o case_file
    case_file = string("results/",dts,"/0_casefile.txt")

    texto = readlines(case_file)

     # Pega uma variavel de cada linha
    tipoan  = parse(Int64,split(texto[1],':')[2])
    freq    = parse(Float64,split(texto[2],':')[2])
    alfa    = parse(Float64,split(texto[3],':')[2])
    beta    = split(split(texto[4],':')[2],' ')[end]
    A       = parse(Float64,split(texto[5],':')[2])
    P       = parse(Float64,split(texto[6],':')[2])
    q       = parse(Float64,split(texto[7],':')[2])
    dini    = parse(Float64,split(texto[8],':')[2])
    dmax    = parse(Float64,split(texto[9],':')[2])
    Sy      = parse(Float64,split(texto[10],':')[2])
    QP      = parse(Float64,split(texto[11],':')[2])
    rho1    = parse(Float64,split(texto[12],':')[2])
    rho2    = parse(Float64,split(texto[13],':')[2])
    rho3    = parse(Float64,split(texto[14],':')[2])
    raiof   = parse(Float64,split(texto[15],':')[2])
    max_fil =   parse(Int64,split(texto[16],':')[2])
    max_int =   parse(Int64,split(texto[17],':')[2])
    tol_ext = parse(Float64,split(texto[18],':')[2])
    tol_int = parse(Float64,split(texto[19],':')[2])
    rho_max = parse(Float64,split(texto[20],':')[2])

    # Pega o heaviside tambem se tiver
    if sub != ""
        # Busca o heav_file tambem
        heav_file = string("results/",dts,"/",sub,"/0_heavifile.txt")

        texto = readlines(heav_file)
        # Do heav tambem
        max_hev =   parse(Int64,split(texto[1],':')[2])
        csi0    = parse(Float64,split(texto[2],':')[2])
        csim    = parse(Float64,split(texto[3],':')[2])

    else
        max_hev = 0
        csi0 = 0.0
        csim = 0.0
    end

    # Caso especial do beta, pode ser definido como string
    if beta != "W"
        beta = parse(Float64,beta)
    end

    return tipoan,freq,alfa,beta,A,P,q,Sy,dini,dmax,QP,csi0,csim,rho1,rho2,rho3,
                           raiof,max_fil,max_hev,max_int,tol_ext,tol_int,rho_max

end

#
#   Verifica disponibilidade da tarefa, 0=OK, 1=Ocupado, 2=Completo
#
function Lock_Game(dts)

    lock_file = string("results/",dts,"/.lock")
    done_file = string("results/",dts,"/.done")

    # Se estiver pronto
    if isfile(done_file)

        return 2

    # Se estiver travado, atualiza a trava
    elseif isfile(lock_file)
        open_lock = open(lock_file,"a")
        println(open_lock,time())
        close(open_lock)
        return 1

    # Se nao tiver nenhum dos dois cria o lock
    else

        close(open(lock_file,"a"))
        return 0

    end

end

function End_Game(dts)

    lock_file = string("results/",dts,"/.lock")
    done_file = string("results/",dts,"/.done")

    if isfile(lock_file); rm(lock_file); end

    close(open(done_file,"a"))

end
