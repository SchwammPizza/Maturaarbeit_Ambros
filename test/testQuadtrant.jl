using CUDA, Colors, Images


function mandelbrotberechnung!(nn::Int64, mm::Int64, iteration::Int64, mandelbrotPart, f::CuDeviceVector{ComplexF64, 1}, r::Int8) #bearbeitet das mandelbrot array so das nunroch 1 und 0 gibt, 1 für drausen und 0 für in der Menge
    indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    strideX = blockDim().x * gridDim().x
    indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    strideY = blockDim().y * gridDim().y

    for i = indexX:strideX:nn
        for j = indexY:strideY:mm
            y = floor((r-1)/2)*-1 + i/nn #im # definitionbereich = [-1im, 1im]
            x = 2*abs(r-2.5)-3 + 2*j/mm  # definitionbereich = [-2, 1]
            if !((y <= 0) & (x <= 0))
                mandelbrotPart[i,j] = 1
            end
        end
    end
    return nothing
end
    
function bench_mandel!(iteration::Int64, mandelbrotPart, r::Int8)
    f = CUDA.zeros(ComplexF64, iteration)
    f .= 3

    nn = CUDA.size(mandelbrotPart, 1)
    mm = CUDA.size(mandelbrotPart, 2)

    numblocks = ceil(Int, length(mandelbrotPart)/anzahlThreads)

    CUDA.@sync begin
        @cuda threads=anzahlThreads blocks=numblocks mandelbrotberechnung!(nn, mm, iteration, mandelbrotPart, f, r)
    end
end

function zeichnen(n::Int64, m::Int64, horizontal::Int64, vertical::Int64, Values, img, maxLanding::Int128)
    indexX = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    strideX = blockDim().x * gridDim().x
    indexY = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    strideY = blockDim().y * gridDim().y
        
    for i = indexX:strideX:n
        for j = indexY:strideY:m
            img[i, j] = RGB{Float64}(Values[i+horizontal-1, j+vertical-1]/maxLanding, Values[i+horizontal-1, j+vertical-1]/maxLanding, Values[i+horizontal-1, j+vertical-1]/maxLanding)
        end
    end
end
    
function bench_zeich!(horizontal::Int64, vertical::Int64, Values, img, maxLanding::Int128)
    n = CUDA.size(img, 1)
    m = CUDA.size(img, 2)

    numblocks = ceil(Int, length(Values)/anzahlThreads)
    CUDA.@sync begin
        @cuda threads=anzahlThreads blocks=numblocks zeichnen(n, m, horizontal, vertical, Values, img, maxLanding)
    end
end
anzahlThreads = 256
a = CUDA.zeros(Int128, 2668, 4002)
img = CUDA.zeros(RGB{Float64},2668,4002)
bench_mandel!(1,a,Int8(1))
bench_zeich!(1,1,a,img,Int128(1))
img_cpu = zeros(RGB{Float64},2668,4002)
img_cpu .= img
print(@__DIR__)
save(string(@__DIR__)*"/quadrant.png", img_cpu)