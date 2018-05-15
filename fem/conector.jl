function box_patch_elemento(i::Int64, j::Int64, nx::Int64, ny::Int64, Nx::Int64, Ny::Int64)

    # Testa para ver se é um par válido
    @assert ((i>0 && i<=Nx) && (j>0 && j<=Ny)) "vizinhos_elemento: i ou j inválido"

    # Limites do Box
    mx = max(1,  i - nx)
    Mx = min(Nx, i + nx)

    my = max(1,  j - ny)
    My = min(Ny, j + ny)

    # Iteradores
    Ii = mx:Mx
    Ij = my:My

    return Ii, Ij

end

# Dado o elemento ele, retorna o par (i,j)
function ele2par(ele::Int64,Nx::Int64,Ny::Int64)

    # A divisão de ele por Nx (dois inteiros) terá
    # uma parte inteira e uma parte fracionária.
    # No Julia, temos duas funções interessantes:
    # div(ele,m.Nx) -> devolve a parte inteira
    # rem(ele,m.Nx) -> devolve um inteiro que permite
    #                  obter a parte fracionária.
    # Por exemplo, 12/7=1.714285.....
    # div(12,7) = 1
    # rem(12,5) = 5, tal que 5/7 = 0.714285....

    # Testa para ver se é um elemento válido
    @assert (ele>0 && ele<=(Nx*Ny)) "ele2par:: ele inválido"

    # Relações descritas acima
    j = div(ele,Nx) + 1
    i = rem(ele,Nx)

    # Temos um caso a ser tratado ( caso o elemento esteja
    # exatamente na ultima posição de cada camada ou
    # na primeira camada). O mesmo vale para a primeira e
    # a última camada..
    if i==0
        i = Nx
        j -= 1
    end

    # Retorna o par
    return i,j

end

# Dado o par (i,j), retorna o elemento
function par2ele(i::Int64,j::Int64,Nx::Int64, Ny::Int64)

     # Testa para ver se é um par válido
     @assert ((i>0 && i<=Nx) && (j>0 && j<=Ny)) "par2ele: i ou j inválido"

     return (j-1)*(Nx) + i

end

function densidades_box(x,ele,Nx,Ny)

    # Adquire o part i,j do elemento
    i,j = ele2par(ele,Nx,Ny)

    # Devolve os contadores com os elementos da vizinhança (1,1)
    Ii,Ij = box_patch_elemento(i, j, 1, 1, Nx, Ny)

    # Aloca x local
    xl = Float64[]

    for k in Ii
        for l in Ij

            # Verifica o elemento k,l
            viz = par2ele(k,l,Nx,Ny)

            # Dá push sem pegar o central
            if viz != l
                push!(xl,x[viz])
            end

        end #for l
    end # for k

    return xl

end

function Heavi(a,b)

    h = 1.0/(1.0 + exp(-2.0*b*(a - 0.2 )))

end

function conector(ele,x,Nx,Ny)

    xcentral = Heavi(x[ele],20.0)

    xviz     = densidades_box(x,ele,Nx,Ny)

    phi = 1.0

    for i in xviz
        phi = phi*abs(xcentral - Heavi(i,10.0))
    end

    phi = xcentral*phi

end

function plota_con(x,nnos,ijk,coord,nel,dts,Nx,Ny)

    fmesh  = string("results/",dts,"/densidades_conects.pos")

    Inicializa_Malha_Gmsh(fmesh, nnos, nel, ijk, coord, 2)

    xcon = zeros(nel)

    for i = 1:nel

        xcon[i] = conector(i, x, Nx, Ny)

    end

    Adiciona_Vista_Escalar_Gmsh(fmesh, "xc", nel, xcon, Float64(0.0))

end
