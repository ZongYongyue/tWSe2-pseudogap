include("./tools.jl")
using Distributed
#开启多进程
spawn(6)
@everywhere begin
include("./tools.jl")
include("./L12MIM_SO.jl")
end 
spim =
[
0.4560302734374999
0.45246394230769227
0.44889761117788457
0.44533128004807687
0.4417649489182691
0.43819861778846153
0.43463228665865383
0.43106595552884613
0.42749962439903844
0.42393329326923074
0.42036696213942304
0.4168006310096154
0.4132342998798077
0.40966796875
0.40351562499999993
0.39736328125
0.3896728515625
0.3822021484375
0
0
0
]


musim = [
12.4
12.0
11.6
11.4
11.2
11.0
10.9
10.9
10.8
10.8
10.7
10.6
10.5
10.4
10.3
10.2
10.1
10.1
10.1
10.1
10.1
]

Vim = range(14,16,21)

vcas = pmap(i->moireVCA(4.2, Vim[i], 4.7, spim[i]), 1:21)

sdfg = pmap(i->SDFG(vcas[i], musim[i]), 1:21)

saveData(sdfg, "L12IM_SO_SU4p7_SDFG.jls")


#=
ops = pmap(i->OrderParameters(:f, moireVCA(4.2, Vim[i], 4.7, spim[i]), rz, so, musim[i]), 1:length(spim))
saveData(ops, "L12IM_OP_U4p7_ops.jls")
=#
#=
vcas = pmap(i->moireVCA(4.2, V[i], 4.7, sp[i]), 1:length(sp))
gre = pmap(i->Gre(vcas[i], mus[i]), 1:length(mus))
saveData(gre, "L12IM_SO_SU4p7_Gre_4.jls")
=#



