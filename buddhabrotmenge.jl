using Images, Colors

@time begin
    n = Int(720)
    m = Int(n/720*1080)

    img = zeros(RGB{Float64}, n, m)

    PointImagToIndex(c) = floor(Int64, ((imag(c) + 1) * n / 2))
    
    PointRealToIndex(c) = floor(Int64, ((real(c) + 2) * m / 3))

    iteration = 150

    for i = 1:n
        for j = 1:m
            y = (2*i - n)/n * 1im
            x = (3*j - 2*m)/m
            z = x + y
            old_z = x + y
            f = []
            for l = 1:iteration
                y = PointImagToIndex(z)
                x = PointRealToIndex(z)

                if z in f
                    break
                elseif abs(z)-abs(old_z) > 2
                    break
                elseif true in (x > m , x < 1 , y > n , y < 1)
                    break
                end
                
                if green(img[y, x]) < 1 - 1/iteration
                    img[y, x] = RGB{Float64}(red(img[y, x]), green(img[y, x]) + 1/iteration, blue(img[y, x]))
                end
                append!(f, z)
                z = z^2 + old_z
            end
        end
    end

    save("/Users/ambros.anrig/Documents/GitHub/Maturaarbeit_Ambros/Anti-Buddhabrotmenge.png", img)
    # save("/Users/Ambros D. Anrig/OneDrive - Kanton Glarus/Dokumente/GitHub/Maturaarbeit_Ambros/Buddhabrotmenge.png", img)
end