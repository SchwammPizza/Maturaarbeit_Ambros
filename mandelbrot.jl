using Images, Colors

n = 720
m = Int(n/720*1080)

img = zeros(RGB{Float64},n,m)

function iterationBerechung(c, list)
    for i = 1:100
        i
    end
end

for i = 1:n
    for j = 1:m
        y = (2*i - n)/n * 1im
        x = (3*j - 2*m)/m
        println(x + y)

    end
end
