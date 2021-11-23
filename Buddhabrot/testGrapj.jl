using Colors, Images

n = 2667
m = 4000

img = zeros(RGB{Float64}, n, m)

for i = 1:n
    for j = 1:m
        if (-i > (28.7*n)*j/(4*m) - m*9.5/4)
            img[i, j] = RGB{Float64}(1,1,1)
        end
    end
end

save(string(@__DIR__) * "/Pictures/gpu/analyse/test.png", img)