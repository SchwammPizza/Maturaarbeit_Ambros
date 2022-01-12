# Dies ist die letzte und somit Optimierte Fassung der Arbeit

using Images, Colors, CUDA

@time begin
    #Varierende variabeln
    n = Int(2668) # muss eine Gerade Zahl sein
    m = Int(floor(n/2*3))

    iteration = 1000

    anzahlThreads = 256

    zoomPoint = -1.25 + 0im
    zoom = 3. # zoom >= 1

    #Berechnete variabeln
    zoomPointAsMatrixPoint = ((-imag(zoomPoint) + 1)*n*zoom/2 + 1, (real(zoomPoint) + 2)*m*zoom/3 + 1)

    # verschiebung des Bildes im gesamt array
    s = 0
    horizontal = floor(Int, zoomPointAsMatrixPoint[1] - n/2)        # zuerst die auf der Komplexenebene hoechere
    if (horizontal < 1)
        s = horizontal
        horizontal = 1
    end
    horizontal2 = floor(Int, zoomPointAsMatrixPoint[1] + n/2) - s
    if (horizontal2 > n*zoom)
        horizontal -= (horizontal2 - n*zoom - 1)
        horizontal2 = n*zoom
    end
    s = 0
    vertical = floor(Int, zoomPointAsMatrixPoint[2] - m/2)      # zuerst die auf der Komplexenebene rechtere
    if (vertical < 1)
        s = vertical
        vertical = 1
    end
    vertical2 = floor(Int, zoomPointAsMatrixPoint[2] + m/2) - s
    if (vertical2 > m*zoom)
        vertical -= (vertical2 - m*zoom - 1)
        vertical2 = m*zoom
    end   
    zoomPointAsMatrixPoint = nothing
    # erstellen der Funktionen

    # bearbeitet das mandelbrot array so das nunroch 1 und 0 gibt, 1 fuer drausen und 0 fuer in der Menge
    function mandelbrotberechnung!(nn::Int64, mm::Int64, iteration::Int64, mandelbrotPart, f::CuDeviceVector{ComplexF64, 1}, r::Int8) 
        indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        strideX = blockDim().x * gridDim().x
        indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
        strideY = blockDim().y * gridDim().y

        for i = indexX:strideX:nn
            for j = indexY:strideY:mm
                y = floor((r-1)/2)*-1 + i/nn # definitionbereich = [-1im, 1im]
                x = 2*abs(r-2.5)-3 + 2*j/mm  # definitionbereich = [-2, 1]
                z = x + y*1im
                if !((abs(z) > 1) & (x >= 0))
                    c = x + y*1im
                        
                    for q = 1:iteration
                        if z in f
                            break
                        elseif abs(z) >= 2
                            mandelbrotPart[i, j] += 1
                            break
                        end
                        f[q] = z
                        z = z^2 + c
                    end
                end
            end
        end
        return nothing
    end

    function bench_mandel!(iteration::Int64, mandelbrotPart, r::Int8)
        f = CUDA.zeros(ComplexF64, iteration)
        f .= 3

        nn = CUDA.size(mandelbrotPart, 1)
        mm = CUDA.size(mandelbrotPart, 2)

        numblocks = ceil(Int, length(mandelbrotPart)/anzahlThreads)
        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks mandelbrotberechnung!(nn, mm, iteration, mandelbrotPart, f, r)
        end
    end

    # schaut was die Maximale Iterationzahl war, und speichert in maxValues wie oft dort ein Punkt ankam
    function berechnungBuddhaBrot!(nzoom::Int64, mzoom::Int64, nn::Int64, mm::Int64, iteration::Int64, maxValues, mandelbrot, f::CuDeviceVector{ComplexF64, 1}, r::Int8) 
        indexY = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        strideY = blockDim().x * gridDim().x
        indexX = (blockIdx().y - 1) * blockDim().y + threadIdx().y
        strideX = blockDim().y * gridDim().y
        
        for i = indexY:strideY:nn
            for j = indexX:strideX:mm
                if mandelbrot[i, j] != 0
                    y = floor((r-1)/2)*-1 + i/nn # definitionbereich = [-1im, 1im]
                    x = 2*abs(r-2.5)-3 + 2*j/mm  # definitionbereich = [-2, 1]
                    z = x + y*1im
                    c = x + y*1im
                    for q = 1:iteration
                        y = CUDA.floor(Int64, (imag(z)+1)*(nzoom-1)/2+1)
                        x = CUDA.floor(Int64, (real(z)+2)*(mzoom-1)/3+1)

                        if z in f # Kontrolle
                            break
                        elseif abs(z) >= 2
                            break
                        end
                        if (1 <= x <= mzoom) & (1 <= y <= nzoom)
                            maxValues[y, x] += 1
                        end
                        
                        f[q] = z
                        z = z^2 + c
                    end
                end
            end
        end
        return nothing
    end

    function bench_buddhi!(iteration::Int64, maxValues, mandelbrot, r::Int8)
        f = CUDA.zeros(ComplexF64, iteration)
        f .= 3

        nn = CUDA.size(mandelbrot,1)
        mm = CUDA.size(mandelbrot,2)
        nzoom = CUDA.size(maxValues,1)
        mzoom = CUDA.size(maxValues,2)
        numblocks = ceil(Int, length(mandelbrot)/anzahlThreads)
        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks berechnungBuddhaBrot!(nzoom, mzoom, nn, mm, iteration, maxValues, mandelbrot, f, r)
        end
    end

    function zeichnen(n::Int64, m::Int64, horizontal::Int64, vertical::Int64, Values, img, maxLanding::Int128)
        indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        strideX = blockDim().x * gridDim().x
        indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
        strideY = blockDim().y * gridDim().y

        for i = indexX:strideX:n
            for j = indexY:strideY:m
                img[i, j] = RGB{Float64}(Values[i+horizontal-1, j+vertical-1]/maxLanding, Values[i+horizontal-1, j+vertical-1]/maxLanding, Values[i+horizontal-1, j+vertical-1]/maxLanding)
            end
        end
    end

    function bench_zeich!(horizontal::Int64, vertical::Int64, Values, img, maxLanding::Int128)
        n = CUDA.size(img, 1)
        m = CUDA.size(img, 2)

        numblocks = ceil(Int, length(Values)/anzahlThreads)
        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks zeichnen(n, m, horizontal, vertical, Values, img, maxLanding)
        end
    end

    # mandelbrot
    mandelbrot = CUDA.Array([CUDA.ones(Int8,2, 2) for _ in 1:4])
    maxValues = CUDA.zeros(Int128, floor(Int, n*zoom), floor(Int, m*zoom))

    # Kontrolle ob Zoom vereinfachung moeglich
    if zoom > 2     
        if (-horizontal2 > (28.7*n)*vertical2/(4*m) - m*9.5/4) #4. Quadrant
            mandelbrot[1] = zeros(Int8, round(Int, n/2*zoom), round(Int, m*zoom/3))
        end
        if ((((vertical2-4340)^2+(horizontal2-n/2)^2>(2.5/3*m)^2) & (horizontal2 < n/2-50)) || (horizontal > 9 * n/(2*m^2)*(vertical-m)^2+n/2)) #3. Quadrant
            mandelbrot[2] = zeros(Int8, round(Int, n/2*zoom), round(Int, 2*m*zoom/3))
        end
        if ((((vertical2-4340)^2+(horizontal-n/2)^2>(2.5/3*m)^2) & (horizontal > n/2+50)) || (horizontal2 < -9 * n/(2*m^2)*(vertical-m)^2+n/2)) #2. Quadrant
            mandelbrot[3] = zeros(Int8, round(Int, n/2*zoom), round(Int, 2*m*zoom/3))
        end
        if (horizontal > (28.7*n)*vertical2/(3.8*m) - n*13/3.8) #1. Quadrant
            mandelbrot[4] = zeros(Int8, round(Int, n/2*zoom), round(Int, m*zoom/3))
        end
    else
        for i = 1:2
            mandelbrot[i^2] = zeros(Int8, round(Int, n/2*zoom), round(Int, m*zoom/3))
            mandelbrot[i+1] = zeros(Int8, round(Int, n/2*zoom), round(Int, 2*m*zoom/3))
        end
    end
    println("Startet MainProgramm")
    # Hauptprogramm
    for r::Int8 = 1:4
        if CUDA.size(mandelbrot[r]) != CUDA.size(CUDA.zeros(Int8,2,2))
            bench_mandel!(iteration, mandelbrot[r], r)
            bench_buddhi!(iteration, maxValues, mandelbrot[r], r)
        end
    end

    # loeschen und neuerstellen von Arrays aufgrund Speicher-Handling
    mandelbrot = nothing
    img = CUDA.zeros(RGB{Float64}, n, m)
    println("Startet Drawing")
    bench_zeich!(horizontal, vertical, maxValues, img, maximum(maxValues))
    println(maximum(maxValues))

    maxValues = nothing
    img_cpu = zeros(RGB{Float64}, n, m)
    img_cpu .= img
    save(string(@__DIR__) * "/Pictures/gpu/film/BuddhabrotmengeWithZoomGPU$(zoom)ToPoint$(zoomPoint)WithIteration$(iteration)withResolution$(m)x$(n).png", img_cpu)
end