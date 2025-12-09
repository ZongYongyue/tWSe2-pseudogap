include("./tools.jl")
using Distributed
#开启多进程
spawn(7)
@everywhere begin
include("./tools.jl")
include("./L12MIM.jl")
end 

Vim = range(14, 16, 21)

musim =  [Vector(range(11.851839464882943, 11.320066889632107,11))..., Vector(range(11.320066889632107,11.581270903010033,10))...]

vcas = pmap(i->moireVCA(4.2, Vim[i], 4.7), 1:21)

sdfg = pmap(i->SDFG(vcas[i], musim[i]), 1:21)

saveData(sdfg, "L12IM_SU4p7_SDFG.jls")