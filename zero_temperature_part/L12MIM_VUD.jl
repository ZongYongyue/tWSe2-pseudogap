include("./tools.jl")
using Distributed
#开启多进程
spawn(12)
@everywhere begin
include("./tools.jl")
include("./L12MIM.jl")
    function VUD(vca, μ)
        unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
        rz = ReciprocalZone(reciprocals(unitcell.vectors); length=100)
        ω_range = range(-60, 60, length=600)
        GG = singleParticleGreenFunction(:f, vca, rz, ω_range; η=0.01*4.755330863071197, μ=μ)
        D = densityofstates(GG)
        return D
    end
end 

us = [4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5]
vcas = pmap(i->moireVCA(4.2, 9, us[i]), 1:12)
musu = 4.755330863071197*us/2
vud = pmap(i->VUD(vcas[i], musu[i]), 1:12)

saveData(vud, "L12MIM_VUD_2.jls")