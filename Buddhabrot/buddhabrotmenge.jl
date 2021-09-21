using Images, Colors

# variabeln deffinieren
n = Int(2434)
m = Int(floor(n/438*720))

color = false # entscheidet ob man ein Mandelbrot-Bild mitmacht oder nicht

maxLanding = 0

# erstellen der aray
img = zeros(RGB{Float64}, n, m)
maxValues = zeros(Int128, n, m)

if color
    mandelbrotimg = zeros(RGB{Float64}, n, m)
else
    mandelbrotimg = zeros(Int8, n, m)
end

iteration = 1000

# erstellen der Funktionen
function farbkreis(i, max_iteration=100)
    r = max_iteration/6

    if i < 0
        i += max_iteration
    elseif i > max_iteration
        i -= max_iteration
    end
    
    if i < r
        return RGB{Float64}(1, 0, i/r)
    elseif  i == r
        return RGB{Float64}(1, 0, 1)
    elseif i < 2*r
        return RGB{Float64}(1 - (i - r)/r, 0, 1)
    elseif  i == 2*r
        return RGB{Float64}(0, 0, 1)
    elseif i < 3*r
        return RGB{Float64}(0, (i - 2*r)/r, 1)
    elseif i == 3*r
        return RGB{Float64}(0, 1, 1)
    elseif i < 4*r
        return RGB{Float64}(0, 1, 1 - (i - 3*r)/r)
    elseif i == 4*r
        return RGB{Float64}(0, 1, 0)
    elseif i < 5*r
        return RGB{Float64}((i - 4*r)/r, 1, 0)
    elseif i == 5*r
        returnRGB{Float64}(1, 1, 0)
    elseif i < max_iteration
        return RGB{Float64}(1, (i - 5*r)/r, 0)
    else
        return RGB{Float64}(1, 0, 0)
    end
end

function mandelbrotmenge()
    for i = 1:n
        for j = 1:m
            y = (2*i - n)/n * 1im
            x = (3*j - 2*m)/m
            z = x + y
            old_z = x + y
            f = []
            for l = 1:iteration
                if z in f
                    break
                elseif abs(z)-abs(old_z) > 2
                    if color
                        mandelbrotimg[i, j] = farbkreis(l, iteration)
                    else
                        mandelbrotimg[i, j] = 1
                    end
                    break
                end
                append!(f, z)
                z = z^2 + old_z
            end
        end
    end
end

function berechnungBuddhaBrot(i, j)
    global maxLanding
    
    y = (2*i - n)/n * 1im
    x = (3*j - 2*m)/m
    z = x + y
    old_z = x + y
    f = []
    for _ = 1:iteration
        y = PointImagToIndex(z)
        x = PointRealToIndex(z)

        if z in f
            break
        elseif abs(z)-abs(old_z) > 2
            break
        elseif true in (x > m , x < 1 , y > n , y < 1)
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

function zeichnen(y, x)
    return RGB{Float64}(maxValues[y, x]/maxLanding, maxValues[y, x]/maxLanding, maxValues[y, x]/maxLanding)
end

PointImagToIndex(c) = floor(Int64, ((imag(c) + 1) * n / 2))

PointRealToIndex(c) = floor(Int64, ((real(c) + 2) * m / 3))

# Hauptprogramm
mandelbrotmenge()

if color
    for i = 1:n
        for j = 1:m
            if mandelbrotimg[i, j] != RGB{Float64}(0, 0, 0)
                berechnungBuddhaBrot(i, j)
            end
        end
    end
else
    for i = 1:n
        for j = 1:m
            if mandelbrotimg[i, j] != 0
                berechnungBuddhaBrot(i, j)
            end
        end
    end
end
for i = 1:n
    for j = 1:m
        img[i, j] = zeichnen(i, j)
    end
end

# Bilderstellung
save(string(@__DIR__) * "/Pictures/Buddhabrotmenge.png", img)
if color
    save(string(@__DIR__) * "/Pictures/Mandelbrotmenge.png", img)
end