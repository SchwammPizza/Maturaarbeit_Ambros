using Images, Colors

# variabeln deffinieren
@time begin
    #Varierende variabeln
    println(@__FILE__)
    n = Int(2667)
    m = Int(floor(n/720*1080))

    iteration = 1000

    zoom = 1 #zoom != 0
    zoomPoint = -0.5 + 0im

    #Berechnete variabeln
    zoomPointAsMatrixPoint = ((-imag(zoomPoint) + 1)*n*zoom/2 + 1, (real(zoomPoint) + 2)*m*zoom/3 + 1)
    println(zoomPointAsMatrixPoint)

    # verschiebung des Bildes im gesamt array
    horizontal = Int(floor(zoomPointAsMatrixPoint[1] - n/2))        # zuerst die auf der Komplexe-ebene rechtere
    vertical = Int(floor(zoomPointAsMatrixPoint[2] - m/2))      # zuerst die auf der Komplexe-ebene höchere
    println(string(horizontal) * " " * string(vertical))
    maxLanding = 0

    # erstellen der array
    mandelbrot = zeros(Int8, n*zoom, m*zoom)
    img = zeros(RGB{Float64}, n, m)
    bigimg = zeros(RGB{Float64}, n, m)
    maxValues = zeros(Int128, n*zoom, m*zoom)
    bigmaxValues = zeros(Int128, n*zoom, m*zoom)
    
    # erstellen der Funktionen

    function mandelbrotmenge() #bearbeitet das mandelbrot array so das nunroch 1 und 0 gibt, 1 für drausen und 0 für in der Menge
        for i = 1:(n*zoom)
            for j = 1:(m*zoom)
                y = (2*i - n)/n * 1im # definitionbereich = [-1im, 1im]
                x = (3*j - 2*m)/m # definitionbereich = [-2, 1]
                z = x + y
                old_z = x + y
                f = []
                for _ = 1:iteration
                    if z in f
                        break
                    elseif abs(z) > 2
                        mandelbrot[i, j] = 1
                        break
                    end
                    append!(f, z)
                    z = z^2 + old_z
                end
            end
        end
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
            if maxValues[y, x] > maxLanding
                maxLanding = maxValues[y, x]
            end
            append!(f, z)
            z = z^2 + old_z
        end
    end

    function zeichnen(y, x, liste = maxValues, maxL = maximum(maxValues))
        return RGB{Float64}(liste[y, x]/maxL, liste[y, x]/maxL, liste[y, x]/maxL)
    end

    PointImagToIndex(c) = round(Int64, ((imag(c) + 1) * (n*zoom) / 2))

    PointRealToIndex(c) = round(Int64, ((real(c) + 2) * (m*zoom) / 3))

    # Hauptprogramm

    if (real(zoomPoint) - 3/(zoom*2) < -2) || (real(zoomPoint) + 3/(zoom*2) > 1) || (imag(zoomPoint) - 1/(zoom) < -1) || (imag(zoomPoint) + 1/(zoom) > 1)
        println("Zoom auserhalb der Vordefinierte Bildreichweite")
        exit()
    end

    mandelbrotmenge()
    #Analyse
    for r = 1:4
        img .= 0
        maxValues .= 0

        for i = PointImagToIndex(floor((r-1)/2)*-1im)+1:PointImagToIndex(-1im*floor((r-3)/2)) #Funktionen für die Bereichen der 1:4 Quadranten
            for j = PointRealToIndex(2*abs(r-2.5)-3)+1:PointRealToIndex(abs(r-2.5)-0.5)
                if mandelbrot[i, j] != 0
                    berechnungBuddhaBrot(i, j)
                end
            end
        end
        for i = 1:n
            for j = 1:m
                img[i, j] = zeichnen(i + horizontal - 1, j + vertical - 1)
            end
        end

        println(string(maximum(maxValues)) * " " * string(r))

        bigmaxValues .+= maxValues
        bigimg .+= img

        # Bildstellung
        save(string(@__DIR__) * "/Pictures/analyse/Analyse" * string(r) * "Zoom" * string(zoom) * "ToPoint" * string(zoomPoint) * "WithIteration" * string(iteration) * "withResolution" * string(m) * "x" * string(n) * ".png", img)

    end 
    for i = 1:n
        for j = 1:m
            img[i, j] = zeichnen(i + horizontal - 1, j + vertical - 1, bigmaxValues, maximum(bigmaxValues))
        end
    end
    save(string(@__DIR__) * "/Pictures/analyse/BigAnalyse1" * "Zoom" * string(zoom) * "ToPoint" * string(zoomPoint) * "WithIteration" * string(iteration) * "withResolution" * string(m) * "x" * string(n) * ".png", img)
    save(string(@__DIR__) * "/Pictures/analyse/BigAnalyse2" * "Zoom" * string(zoom) * "ToPoint" * string(zoomPoint) * "WithIteration" * string(iteration) * "withResolution" * string(m) * "x" * string(n) * ".png", bigimg)
end