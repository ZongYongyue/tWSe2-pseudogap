include("./tools.jl")
using Distributed
#开启多进程
spawn(21)
@everywhere begin
include("./tools.jl")
include("./L12MIM.jl")
end 

vcas = pmap(i->moireVCA(4.2, V[i], 4.7), 1:21)

sdfg = pmap(i->SDFG(vcas[i], mus[i]), 1:21)

saveData(sdfg, "L12MIM_SU4p7_SDFG.jls")