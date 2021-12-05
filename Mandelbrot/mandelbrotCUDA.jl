using CUDA, Images, Colors

@time begin
    n = 2668
    m = 4002

    anzahlThreads = 256

    iteration = 150

    mandelbrot = CUDA.zeros(RGB{Float64}, n, m)

    function mandelbrotberechnung!(nn::Int64, mm::Int64, iteration::Int64, mandelbrotPart, f::CuDeviceVector{ComplexF64, 1}) #bearbeitet das mandelbrot array so das nunroch 1 und 0 gibt, 1 für drausen und 0 für in der Menge
        indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        strideX = blockDim().x * gridDim().x
        indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
        strideY = blockDim().y * gridDim().y

        for i = indexX:strideX:nn
            for j = indexY:strideY:mm
                y = (2*i - nn)/nn # definitionbereich = [-1im, 1im]
                x = (3*j - 2*mm)/mm # definitionbereich = [-2, 1]
                z = x + y*1im
                c = x + y*1im
     
                for q = 1:iteration
                    if (z in f)
                        break
                    elseif (abs(z) >= 2)
                        r = iteration/6

                        if r == 0
                            @cuprintln(r)
                        end
            
                        if q < r
                            l = RGB{Float64}(1, 0, q/r)
                        elseif  q == r
                            l = RGB{Float64}(1, 0, 1)
                        elseif q < 2*r
                            l = RGB{Float64}(1 - (q - r)/r, 0, 1)
                        elseif  i == 2*r
                            l = RGB{Float64}(0, 0, 1)
                        elseif q < 3*r
                            l = RGB{Float64}(0, (q - 2*r)/r, 1)
                        elseif q == 3*r
                            l = RGB{Float64}(0, 1, 1)
                        elseif q < 4*r
                            l = RGB{Float64}(0, 1, 1 - (q - 3*r)/r)
                        elseif q == 4*r
                            l = RGB{Float64}(0, 1, 0)
                        elseif q < 5*r
                            l = RGB{Float64}((q - 4*r)/r, 1, 0)
                        elseif q == 5*r
                            l = RGB{Float64}(q, 1, 0)
                        elseif i < iteration
                            l = RGB{Float64}(1, (q - 5*r)/r, 0)
                        else
                            l = RGB{Float64}(1, 0, 0)
                        end

                        mandelbrotPart[i, j] = l
                        break
                    end
                    f[q] = z
                    z = z^2 + c
                end
            end
        end
        return nothing
    end
    
    function bench_mandel!(iteration::Int64, mandelbrotPart)
        f = CUDA.zeros(ComplexF64, iteration)
        f .= 3

        nn = CUDA.size(mandelbrotPart, 1)
        mm = CUDA.size(mandelbrotPart, 2)

        numblocks = ceil(Int, length(mandelbrotPart)/anzahlThreads)

        CUDA.@sync begin
            @cuda threads=anzahlThreads blocks=numblocks mandelbrotberechnung!(nn, mm, iteration, mandelbrotPart, f)
        end
    end

    bench_mandel!(iteration, mandelbrot)

    img = zeros(RGB{Float64}, n, m)
    img .= mandelbrot
    save(string(@__DIR__)*"/Pictures/Manddddddelbrot$(n)x$m.png", img)
end