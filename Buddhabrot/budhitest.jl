using Images, Colors, CUDA

@time begin
    println(@__FILE__)
    # Varierende variabeln
    n = Int(2667)
    m = Int(round(n/720*1080))
    # n = 360
    # m = 640
    nm = CUDA.Array([n,m])

    iteration = 1000
    anzahlThreads = 256

    zoom = 1  # zoom != 0
    # zoomPoint = -1.25 + 0im
    # zoomPoint = -0.5 + 0.5im
    zoomPoint = -0.5 + 0im

    # Berechnete variabeln
    # verschiebung des Bildes im gesamt array
    horizontal = Int(round((-imag(zoomPoint) + 1)*n*zoom/2 + 1 - n/2))        # zuerst die auf der Komplexe-ebene rechtere
    vertical = Int(round((real(zoomPoint) + 2)*m*zoom/3 + 1 - m/2))      # zuerst die auf der Komplexe-ebene höchere
    Ho_Ve = CUDA.Array([horizontal, vertical])
    println(string(horizontal) * " " * string(vertical))
    maxLanding = CUDA.Array([0])

    # erstellen der array
    img = CUDA.zeros(RGB{Float64}, nm[1], nm[2])
    maxValues = CUDA.zeros(Int128, nm[1]*zoom, nm[2]*zoom)

    mandelbrot = []
    # mandelbrot = CUDA.CuArray(mandelbrot)
    # erstellen der Funktionen

    function mandelbrotmenge(i, j) #bearbeitet das mandelbrot array so das nunroch 1 und 0 gibt, 1 für drausen und 0 für in der Menge 
        y = (2*i - n)/n * 1im # definitionbereich = [-1im, 1im]
        x = (3*j - 2*m)/m # definitionbereich = [-2, 1]
        z = x + y
        old_z = x + y
        f = []
        for _ = 1:iteration
            if z in f
                break
            elseif abs(z) > 2
                # filter!(e->e≠[i,j], mandelbrot)
                return nothing
            end
            append!(f, z)
            z = z^2 + old_z
        end
        append!(mandelbrot, [[i, j]])
        return nothing
    end

    function berechnungBuddhaBrot(i, j) #schaut was die Maximale Iterationzahl war, und speichert in maxValues wie oft dort ein Punkt ankam
        global maxLanding
        
        y = (2*i - n)/n * 1im # definitionbereich = [-1im, 1im]
        x = (3*j - 2*m)/m # definitionbereich = [-2, 1]
        z = x + y
        old_z = x + y
        f = []
        for _ = 1:iteration
            y = PointImagToIndex(z)
            x = PointRealToIndex(z)

            if z in f
                break
            elseif abs(z) > 2
                break
            elseif true in (x > (m*zoom) , x < 1 , y > (n*zoom) , y < 1)
                break
            end
            maxValues[y, x] += 1
            append!(f, z)
            z = z^2 + old_z
        end
        return nothing
    end

    zeichnen(y, x) = RGB{Float64}(maxValues[y, x]/maxLanding[1], maxValues[y, x]/maxLanding[1], maxValues[y, x]/maxLanding[1])

    PointImagToIndex(c) = round(Int64, ((imag(c) + 1) * (n*zoom) / 2))

    PointRealToIndex(c) = round(Int64, ((real(c) + 2) * (m*zoom) / 3))

    # Hauptprogramm

    if (real(zoomPoint) - 3/(zoom*2) < -2) || (real(zoomPoint) + 3/(zoom*2) > 1) || (imag(zoomPoint) - 1/(zoom) < -1) || (imag(zoomPoint) + 1/(zoom) > 1)
        println("Zoom auserhalb der Vordefinierte Bildreichweite")
        exit()
    end

    for i = 1:(n*zoom)
        for j = 1:(m*zoom)
            mandelbrotmenge(i,j)
        end
    end

    # @cuda threads=anzahlThreads mandelbrotmenge
    # for i = 1:length(mandelbrot)
    #     berechnungBuddhaBrot(mandelbrot[i][1], mandelbrot[i][2])
    # end
    for index in mandelbrot
        berechnungBuddhaBrot(index[1], index[2])
    end
    maxLanding[1] = maximum(maxValues)
    mandelbrot = nothing
    # println(maxValues)
    for i = 1:n
        for j = 1:m
            img[i, j] = zeichnen(i + horizontal - 1, j + vertical - 1)
        end
    end
    # Bildstellung
    img_cpu = zeros(RGB{Float64}, n, m)
    img_cpu .= img
    # print(img_cpu)
    save(string(@__DIR__) * "/Pictures/ABuddhabrotmengeWithZoomGPU$(zoom)ToPoint$(zoomPoint)WithIteration$(iteration)withResolution$(m)x$(n).png", img_cpu)
end