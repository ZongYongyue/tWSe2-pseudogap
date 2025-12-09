using LinearAlgebra
BLAS.set_num_threads(1)
using Pkg
Pkg.activate(".")
include("./tools.jl")
using Distributed

spawn(21)
@everywhere begin
include("./tools.jl")
include("./L12MIM_SO.jl")
end 


vcas = pmap(i->moireVCA(4.2, V[i], 4.7, sp[i]), 1:length(sp))

sdfg = pmap(i->SDFG(vcas[i], mus[i]), 1:length(mus))

saveData(sdfg, "L12MIM_SO_SU4p7_SDFG.jls")

vcas = pmap(i->moireVCA(4.2, V[i], 4.7, sp[i]), 1:length(sp))

qpis = pmap(i->QPI(vcas[i], mus[i]), 1:length(mus))

saveData(qpis, "L12MIM_SO_SU4p7_QPI.jls")

vca4 = moireVCA(4.2, V[4], 4.7, sp[4])
vca16 = moireVCA(4.2, V[16], 4.7, sp[16])
ws = range(-10,10,21)
qpis4 = pmap(i->QPI(vca4, mus[4], ws[i]), 1:length(ws))
qpis16 = pmap(i->QPI(vca16, mus[16], ws[i]), 1:length(ws))
saveData(qpis4, "L12MIM_SO_SU4p7_QPI_WS_4.jls")
saveData(qpis16, "L12MIM_SO_SU4p7_QPI_WS_16.jls")

vca4 = moireVCA(4.2, V[4], 4.7, sp[4])
ws = range(-3,3,21)
qpis4 = pmap(i->QPI(vca4, mus[4], ws[i]), 1:length(ws))
saveData(qpis4, "L12MIM_SO_SU4p7_QPI_WS_4_2.jls")

vca4 = moireVCA(4.2, V[4], 4.7, sp[4])
vca16 = moireVCA(4.2, V[16], 4.7, sp[16])
ws = range(-10,10,21)
qpis4 = pmap(i->QPI(vca4, mus[4], ws[i], Vector(1:1)), 1:length(ws))
qpis16 = pmap(i->QPI(vca16, mus[16], ws[i],Vector(1:1)), 1:length(ws))
saveData(qpis4, "L12MIM_SO_SU4p7_QPI_WS_4_spind.jls")
saveData(qpis16, "L12MIM_SO_SU4p7_QPI_WS_16_spind.jls")


vca4 = moireVCA(4.2, V[4], 4.7, sp[4])
ws = range(-3,3,21)
qpis4 = pmap(i->QPI(vca4, mus[4], ws[i], Vector(1:1)), 1:length(ws))
saveData(qpis4, "L12MIM_SO_SU4p7_QPI_WS_4_spind_2.jls")
