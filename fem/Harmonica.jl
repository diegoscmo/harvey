    function Harmonica(dts::String, f_ini::Float64 ,f_delta::Float64 ,f_fim::Float64,
                    alfa::Float64, beta::Float64, nel::Int64, ijk::Array{Int64,2},
                    ID::Array{Int64,2}, K0::Array{Float64,2}, M0::Array{Float64,2},
                    SP::Float64, vmin::Float64, F::Array{Float64,1}, vizi::Array{Int64,2},
                    nviz::Array{Int64,1}, dviz::Array{Float64,2}, raiof::Float64)

    # Abre saída
    file  = string(dts,"_H.txt")
    if isfile(file)
        rm(file)
    end
    saida = open(file,"a")

    # Carrega o x do arquivo
    x = Le_Densidades(string(dts,"_A.pos"),nel)

    # Filtra
    x = Filtro_Dens(x, nel, vizi, nviz, dviz, raiof)


    # Monta matriz de rigidez e massa global, aqui é aplicado o SIMP
    KG, MG = Global_KM(x, nel, ijk, ID, K0, M0, SP, vmin)

    # Inicializa frequencia e vetores
    freq = collect(f_ini:f_delta:f_fim)
    tamf = size(freq,1)

    for i = 1:tamf

        # Determina a frequencia em rad/s
        f = freq[i]
        w = 2.0*pi*f

        # Resolve o sistema dinamico
        UD = Array{Complex{Float64}}(uninitialized,nel)
        try
            KD = Hermitian(sparse(KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG))
            UD = vec(cholfact(KD;shift=-1E2)\F)
        catch
            KD = KG + w*im*(alfa*MG + beta*KG) - (w^2.0)*MG
            UD = vec(lufact(KD)\F)
        end

        # Salva frequencia e flexibilidade
        Y = log10(abs(real(dot(F,UD))))

        println(saida,f," ",Y)

    end #i

    close(saida)

end #function
