using Images, Colors

@time begin
    n = Int(2668)
    m = Int(n/720*1080)

    img = zeros(RGB{Float64}, n, m)

    function farbkreis(i, max_iteration)
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
            return RGB{Float64}(1, 1, 0)
        elseif i < max_iteration
            return RGB{Float64}(1, (i - 5*r)/r, 0)
        else
            return RGB{Float64}(1, 0, 0)
        end
    end

    function mandelbrotmenge()
        for i = 1:n
            for j = 1:m
                y = (2*i - n)/n * 1im # definitionbereich = [-1im, 1im]
                x = (3*j - 2*m)/m # definitionbereich = [-2, 1]
                z = x + y
                old_z = x + y
                f = []
                for l = 1:iteration
                    if z in f
                        break
                    elseif abs(z) > 2
                        img[i, j] = farbkreis(l, iteration)
                        break
                    end
                    append!(f, z)
                    z = z^2 + old_z
                end
            end
        end
    end

    iteration = 1000

    mandelbrotmenge()
    
    save(string(@__DIR__) * "/Pictures/Mandelbrotmenge.png", img)
end