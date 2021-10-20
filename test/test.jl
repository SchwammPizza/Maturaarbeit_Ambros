n = Int(2667)
m = Int(round(n/720*1080))

PointImagToIndex(c) = round(Int64, ((imag(c) + 1) * (n*zoom) / 2))

PointRealToIndex(c) = round(Int64, ((real(c) + 2) * (m*zoom) / 3))

println()
println("start")
println()

for r = 1:4
    println(PointRealToIndex(2*abs(r-2.5)-3))
    println(PointRealToIndex(abs(r-2.5)-0.5))
    println()
end       
for r = 1:4
    println(PointRealToIndex(2*abs(r-2.5)-3)+1:PointRealToIndex(abs(r-2.5)-0.5))
    println()
end     