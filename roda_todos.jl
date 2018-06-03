################################################################################
###           Lagrangiano Aumentado // Otimização Topológica                 ###
################################################################################

 using Printf
 using SparseArrays
 using LinearAlgebra
 using IterativeEigensolvers
 using ProgressMeter

# Carrega rotinas
include("fem/rotinas_fem.jl")      # Rotinas de Elementos Finitos
include("opt/top_opt.jl")          # Rotinas de LA e Topológica

# Carrega Fobjs
include("fobj/11_norma_global.jl")
include("fobj/12_norma_local.jl")
include("fobj/13_norma_global_r-stress.jl")
include("fobj/14_norma_local_r-stress.jl")

#
# Entrada de dados + execução da rotina
#
function Roda_Todos()

    # Iniciamos deixando todos os parametros fixos aqui...
    # Parâmetros da topológica fixos (tudo daqui pra baixo nao muda, nao vai do casefile)
    SP          = 3.00      # Parametro p do SIMP
    vmin        = 1E-3      # Densidade mínima

    # Parametro obsoleto (apenas para funções antigas)
    R_b         = .999      # Porcentagem da restrição (Est ou R)

    # Parâmetros do problema de FEM, 60x30 = 1800 // 70x35 = 2450 // 80x40 = 3200
    NX          = 100       # Nr. de elementos em X
    NY          = 50        # Nr. de elementos em Y
    LX          = 1.0       # Comprimento em X
    LY          = 0.5       # Comprimento em Y
    young       = 210E9     # Módulo de Young
    poisson     = 0.0       # Coeficiente de Poisson
    esp         = 1.0       # Espessura do retângulo
    p_dens      = 7860.0    # Densidade

    # Restrições de deslocamento (apoios):
    #        [ ponto_X0 ponto_Y0    ponto_XF    ponto_YF    dir (X=1 Y=2)]
    presos = [  0.0     0.0         0.0         LY          1 ;
                0.0     0.0         0.0         LY          2 ]

    # Regiões com elementos travados em 1.0 (quadrado)
    #        [ ponto_X0     ponto_Y0    ponto_X1     ponto_Y1
    travas = [  LX-0.025      0.2         LX           0.3  ]

    # Carregamentos:
    #       [  x0   y0  xf  yf  carreg      dir (X=1 Y=2)
    forcas = [ LX   0.2 LX  0.3 -10000.0    2   ]

    # Agora vamos ao que interessa, loop nas pastas.
    # Enquanto houver pastas para verificar...
    lista  = 1
    standby = 0
    while sum(lista) > 0

        println("\n   Varrendo diretorios...")

        # Vai para a pasta de resultados
        cd("results")

        # Acumula as pastas, arquivos etc
        diretorios = [];
        for (root,dirs,files) in walkdir(".")
                push!(diretorios,dirs)
        end

        # E retorna para o nivel anterior
        cd("..")

        # Separa só os diretorios
        diretorios = diretorios[1]

        # Numero de diretorios
        nd = length(diretorios)

        println("   ",nd," diretorios encontrados...")

        # Lista dos diretorios numerados
        lista = collect(1:nd)

        # Primeiro loop interno
        # Remove diretorios completos ou travados
        for j=1:nd

            # Renomeia o diretorio de trabalho
            dts = diretorios[j]

            # Verifica se a trava esta obsoleta...
            Limpa_Trava(dts)

            # Nome dos arquivos de trava/done
            lock_file = string("results/",dts,"/.lock")
            done_file = string("results/",dts,"/.done")

            # Se houver um dos dois arquivos, zera posição na lista
            if isfile(lock_file) || isfile(done_file)
                lista[j] = 0
            end

        end

        nz = count(!iszero, lista)

        println("   ",nz," diretorios disponiveis...")

        # Agora podemos calcular nos diretorios livres
        for j in lista
            if j != 0

                # Pega o nome da pasta para rodar
                dts = diretorios[j]

                # Carrega as variaveis
                tipoan,freq,alfa,beta,A,P,q,Sy,dini,dmax,QP,csi0,csim,rho1,rho2,rho3,
                                       raiof,max_fil,max_hev,max_int,tol_ext,tol_int,rho_max = Read_Game(dts)
                # Arruma rho
                rho = [rho1;rho2;rho3]

                println("   Rodando diretorio - [",j,"/",nd,"]")

                # Executa.
                Top_Opt(tipoan, dts, Sy, freq, alfa, beta, R_b, A, dini, max_fil, max_int, tol_ext,
                          tol_int, rho, rho_max, SP, raiof, vmin, NX, NY, LX, LY,
                        young, poisson, esp, p_dens, presos, forcas, travas, QP, csi0, csim, dmax, P, q, max_hev)
            end
        end

        # Toda vez que iterar vai rever a lista, só vai parar quando todos estiverem lock ou done.

        # Caso nao houver nada disponivel
        if nz == 0

            # Itera o numero de standby, se passou de 5, fecha
            standby += 1
            if standby >= 5
                exit()
            end

            # Caso contrario espera 20min
            println("   Mais nada para calcular, StandBy por 20min... [",standby,"/5]\n")
            sleep(20.0*60.0)
        end
    end # while

end # function

#
#   Compara a ultima alteracao no lock file e no monitora para excluir travas obsoletas
#   As travas sao liberadas quando o arquivo ficou parado por mais de X h
#
function Limpa_Trava(dts)

    # Nomeia os arquivos de interesse
    lock_file = string("results/",dts,"/.lock")

    # Se tiver um lock file...
    if isfile(lock_file)

        # Pega a data da ultima alteracao do lock file
        time_lock = mtime(lock_file)

        # Cria um arquivo .now so pra pegar o tempo atual sem usar Dates
        now_file = string("results/",dts,"/.now")
        close(open(now_file,"a"))
        time_now = mtime(now_file)
        rm(now_file)

        # Computa a diferenca
        time_diff = time_now - time_lock

        # Se a ultima alteracao for maior do que X h (X*60.0*60.0 segundos), limpa o lock
        X = 1.
        if time_diff > X*3600.0
            rm(lock_file)
            println("   trava removida - ",dts)
        end
    end

end

# Já roda.
Roda_Todos()
