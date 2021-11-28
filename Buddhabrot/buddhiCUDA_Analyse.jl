using Images, Colors, CUDA

@time begin
    #Varierende variabeln
    n = Int(2667)
    m = Int(floor(n/720*1080))

    iteration = 1000
    anzahlThreads = 256

    # zoom = 6.25  #zoom != 0
    # zoomPoint = -1.25 + 0im
    # zoomPoint = -0.5 + 0.5im
    zoomPoint = -0.5 + 0im
    zoom = 1

    #Berechnete variabeln
    zoomPointAsMatrixPoint = ((-imag(zoomPoint) + 1)*n*zoom/2 + 1, (real(zoomPoint) + 2)*m*zoom/3 + 1)
    println(zoomPointAsMatrixPoint)

    # verschiebung des Bildes im gesamt array
    horizontal = Int(floor(zoomPointAsMatrixPoint[1] - n/2))        # zuerst die auf der Komplexe-ebene rechtere
    vertical = Int(floor(zoomPointAsMatrixPoint[2] - m/2))      # zuerst die auf der Komplexe-ebene höchere
    println(string(horizontal) * " " * string(vertical))

    # erstellen der array

    mandelbrot = CUDA.zeros(Int8, round(Int, n*zoom), round(Int, m*zoom))

    # erstellen der Funktionen

    function mandelbrotmenge!(nzoom::Int64, mzoom::Int64, n::Int64, m::Int64, iteration::Int64, mandelbrot, f::CuDeviceVector{ComplexF64, 1}) #bearbeitet das mandelbrot array so das nunroch 1 und 0 gibt, 1 für drausen und 0 für in der Menge
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

    function berechnungBuddhaBrot!(nzoom::Int64, mzoom::Int64, n::Int64, m::Int64, iteration::Int64, maxValues, mandelbrot, f::CuDeviceVector{ComplexF64, 1}, starty, startx, endy, endx, q) #schaut was die Maximale Iterationzahl war, und speichert in maxValues wie oft dort ein Punkt ankam
        indexY = (blockIdx().x - 1) * blockDim().x + threadIdx().x + starty
        strideY = blockDim().x * gridDim().x
        indexX = (blockIdx().y - 1) * blockDim().y + threadIdx().y + startx
        strideX = blockDim().y * gridDim().y
        
        for i = indexY:strideY:endy
            for j = indexX:strideX:endx
                if mandelbrot[i, j] != 0
                    y = (2*i - n)/n * 1im # definitionbereich = [-1im, 1im]
                    x = (3*j - 2*m)/m # definitionbereich = [-2, 1]
                    z = x + y
                    c = x + y
                    if q==1 ⊻ (abs(z) > 1)
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

    function bench_buddhi!(nzoom::Int64, mzoom::Int64, n::Int64, m::Int64, iteration::Int64, maxValues, mandelbrot, starty, endy, startx, endx, q)
        f = CUDA.zeros(ComplexF64, iteration)
        f .= 3
        numblocks = ceil(Int, length(maxValues)/anzahlThreads)
        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks berechnungBuddhaBrot!(nzoom, mzoom, n, m, iteration, maxValues, mandelbrot, f, starty, startx, endy, endx, q)
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
                img[i,j] ^= .5
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
    big_max = CUDA.zeros(Int128, n, m)
    big_img = CUDA.zeros(RGB{Float64}, n, m)
    bench_mandel!(n*zoom, m*zoom, n, m, iteration, mandelbrot)
    for r = 1:4
        for q = 0:1
            global maxValues = CUDA.zeros(Int128, round(Int,n*zoom), round(Int,m*zoom))
            global img = CUDA.zeros(RGB{Float64}, n, m)
            bench_buddhi!(n*zoom, m*zoom, n, m, iteration, maxValues, mandelbrot, round(Int64, ((imag(floor((r-1)/2)*-1im)) + 1) * (n*zoom) / 2)+1, round(Int64, ((imag(-1im*floor((r-3)/2)) + 1) * (n*zoom) / 2)), round(Int64, ((2*abs(r-2.5)-1) * (m*zoom) / 3))+1, round(Int64, ((abs(r-2.5)+1.5) * (m*zoom) / 3)), q)
            bench_zeich!(n,m,horizontal,vertical,maxValues,img,maximum(maxValues))
            println(maximum(maxValues))
            global img_cpu = zeros(RGB{Float64}, n, m)
            img_cpu .= img
            big_max .+= maxValues
            save(string(@__DIR__) * "/Pictures/gpu/analyse/second/$(r)$(q)AnalyBuddhabrotmengeWithZoomGPU$(zoom)ToPoint$(zoomPoint)WithIteration$(iteration)withResolution$(m)x$(n).png", img_cpu)
        end
    end
    println(maximum(big_max))
    bench_zeich!(n,m,horizontal,vertical,big_max,big_img,maximum(big_max))
    big_img_cpu = zeros(RGB{Float64}, n, m)
    big_img_cpu .= big_img
    save(string(@__DIR__) * "/Pictures/gpu/analyse/AnalyBuddhabrotmengeWithZoomGPU$(zoom)ToPoint$(zoomPoint)WithIteration$(iteration)withResolution$(m)x$(n).png", big_img_cpu)
end