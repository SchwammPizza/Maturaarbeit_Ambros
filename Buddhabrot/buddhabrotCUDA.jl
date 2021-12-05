# Dies ist die mit CUDA optimierte Fassung, somit die 2.

using Images, Colors, CUDA

@time begin
    #Varierende variabeln
    const n = Int(2668)
    const m = Int(floor(n/2*3))

    const iteration = 150
    const anzahlThreads = 256

    const zoomPoint = -.5 + 0im
    const zoom = 1

    #Berechnete variabeln
    const zoomPointAsMatrixPoint = ((-imag(zoomPoint) + 1)*n*zoom/2 + 1, (real(zoomPoint) + 2)*m*zoom/3 + 1)
    println(zoomPointAsMatrixPoint)

    # verschiebung des Bildes im gesamt array
    const horizontal = Int(floor(zoomPointAsMatrixPoint[1] - n/2))        # zuerst die auf der Komplexenebene rechtere
    const vertical = Int(floor(zoomPointAsMatrixPoint[2] - m/2))      # zuerst die auf der Komplexenebene hoechere
    println(string(horizontal) * " " * string(vertical))

    # erstellen der array
    mandelbrot = CUDA.zeros(Int8, round(Int, n*zoom), round(Int, m*zoom))

    # erstellen der Funktionen

    # bearbeitet das mandelbrot array so das nunroch 1 und 0 gibt, 1 fuer drausen und 0 fuer in der Menge
    function mandelbrotmenge!(nzoom::Int64, mzoom::Int64, n::Int64, m::Int64, iteration::Int64, mandelbrot, f::CuDeviceVector{ComplexF64, 1}) 
        indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        strideX = blockDim().x * gridDim().x
        indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
        strideY = blockDim().y * gridDim().y
        
        for i = indexX:strideX:nzoom
            for j = indexY:strideY:mzoom
                y = (2*i - n)/n * 1im # definitionbereich = [-1im, 1im]
                x = (3*j - 2*m)/m # definitionbereich = [-2, 1]
                z = x + y
                c = x + y
                
                for r = 1:iteration
                    if z in f
                        break
                    elseif abs(z) >= 2
                        mandelbrot[i, j] += 1
                        break
                    end
                    f[r] = z
                    z = z^2 + c
                end
            end
        end
        return nothing
    end

    function bench_mandel!(nzoom::Int64, mzoom::Int64, n::Int64, m::Int64, iteration::Int64, mandelbrot)
        f = CUDA.zeros(ComplexF64, iteration)
        f .= 3
        numblocks = ceil(Int, length(mandelbrot)/anzahlThreads)
        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks mandelbrotmenge!(nzoom, mzoom, n, m, iteration, mandelbrot, f)
        end
    end

    # schaut was die Maximale Iterationzahl war, und speichert in maxValues wie oft dort ein Punkt ankam
    function berechnungBuddhaBrot!(nzoom::Int64, mzoom::Int64, n::Int64, m::Int64, iteration::Int64, maxValues, mandelbrot, f::CuDeviceVector{ComplexF64, 1}) 
        indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        strideX = blockDim().x * gridDim().x
        indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
        strideY = blockDim().y * gridDim().y
        
        for i = indexX:strideX:nzoom
            for j = indexY:strideY:mzoom
                if mandelbrot[i, j] != 0
                    y = (2*i - n)/n * 1im # definitionbereich = [-1im, 1im]
                    x = (3*j - 2*m)/m # definitionbereich = [-2, 1]
                    z = x + y
                    if !((abs(z) > 1) & (x >= 0))
                        c = x + y
                        for r = 1:iteration
                            y = CUDA.floor(Int64, (imag(z)+1)*nzoom/2)
                            x = CUDA.floor(Int64, (real(z)+2)*mzoom/3)

                            if z in f
                                break
                            elseif abs(z) >= 2
                                break
                            end
                            if (1 < x < mzoom) & (1 < y < nzoom)
                                maxValues[y, x] += 1
                            end
                            
                            f[r] = z
                            z = z^2 + c
                        end
                    end
                end
            end
        end
        return nothing
    end

    function bench_buddhi!(nzoom::Int64, mzoom::Int64, n::Int64, m::Int64, iteration::Int64, maxValues, mandelbrot)
        f = CUDA.zeros(ComplexF64, iteration)
        f .= 3
        numblocks = ceil(Int, length(mandelbrot)/anzahlThreads)
        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks berechnungBuddhaBrot!(nzoom, mzoom, n, m, iteration, maxValues, mandelbrot, f)
        end
    end

    function zeichnen(n::Int64, m::Int64, horizontal::Int64, vertical::Int64, maxValues, img, maxLanding::Int128)
        indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        strideX = blockDim().x * gridDim().x
        indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
        strideY = blockDim().y * gridDim().y

        for i = indexX:strideX:n
            for j = indexY:strideY:m
                img[i, j] = RGB{Float64}(maxValues[i+horizontal-1, j+vertical-1]/maxLanding, maxValues[i+horizontal-1, j+vertical-1]/maxLanding, maxValues[i+horizontal-1, j+vertical-1]/maxLanding)
            end
        end
    end

    function bench_zeich!(n::Int64, m::Int64, horizontal::Int64, vertical::Int64, maxValues, img, maxLanding::Int128)
        numblocks = ceil(Int, length(maxValues)/anzahlThreads)
        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks zeichnen(n, m, horizontal, vertical, maxValues, img, maxLanding)
        end
    end

    # Hauptprogramm
    if (real(zoomPoint) - 3/(zoom*2) < -2) || (real(zoomPoint) + 3/(zoom*2) > 1) || (imag(zoomPoint) - 1/(zoom) < -1) || (imag(zoomPoint) + 1/(zoom) > 1)
        println("Zoom auserhalb der Vordefinierte Bildreichweite")
        exit()
    end

    bench_mandel!(round(Int64, n*zoom), round(Int64, m*zoom), n, m, iteration, mandelbrot)
    maxValues = CUDA.zeros(Int128, round(Int,n*zoom), round(Int,m*zoom))
    bench_buddhi!(round(Int64, n*zoom), round(Int64, m*zoom), n, m, iteration, maxValues, mandelbrot)
    mandelbrot = nothing
    
    img = CUDA.zeros(RGB{Float64}, n, m)
    println(maximum(maxValues))
    bench_zeich!(n, m, horizontal, vertical, maxValues, img, maximum(maxValues))
    
    # Bildstellung
    img_cpu = zeros(RGB{Float64}, n, m)
    img_cpu .= img
    save(string(@__DIR__) * "/Pictures/gpu/BuddhabrotmengeWithZoomGPU$(zoom)ToPoint$(zoomPoint)WithIteration$(iteration)withResolution$(m)x$(n)high.png", img_cpu)
end
