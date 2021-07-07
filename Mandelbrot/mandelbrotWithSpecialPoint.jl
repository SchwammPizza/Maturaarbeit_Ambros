using Images, Colors

@time begin
    midPoint = 0.25 + 0im
    zoom = 1.62 * 10^32 # zoom > 0
    n = Int(720*100)
    m = Int(n/72*108)

    img = zeros(RGB{Float64}, n, m)

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

    function mandelbrotmenge(farbe = true)
        zoom_y = 2im/(zoom * n)

        for i = 1:n
            for j = 1:m
                y = (n/2 - (i - 1)) * zoom_y - imag(midPoint)
                x = (-m/2 + (j - 1)) * 3/(zoom * m) + real(midPoint)
                z = x + y
                old_z = x + y
                f = []
                for l = 1:iteration
                    if z in f
                        break
                    elseif abs(z)-abs(old_z) > 2
                        if farbe
                            img[i, j] = farbkreis(l, iteration)
                        else
                            img[i, j] = 1
                        end
                        break
                    end
                    append!(f, z)
                    z = z^2 + old_z
                end
            end
        end
    end

    iteration = 200

    mandelbrotmenge()
    
    save(string(@__DIR__) * "/Pictures/Zooms/MandelbrotmengeZoomToPoint" * string(midPoint) * " withZoom" * string(zoom) * ".png", img)
end