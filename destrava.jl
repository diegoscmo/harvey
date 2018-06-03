#
# Rotina boba para destravar os locks
#
function Destrava()

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

    # Remove travas de cada diretorios
    for j=1:nd

            # Renomeia o diretorio de trabalho
            dts = diretorios[j]

            # Nome dos arquivos de trava/done
            lock_file = string("results/",dts,"/.lock")

            # Se houver um dos dois arquivos, zera posição na lista
            if isfile(lock_file)
                rm(lock_file)
            end

        end

end
Destrava()
