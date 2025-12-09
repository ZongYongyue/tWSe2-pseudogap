using Pkg
Pkg.activate(".")
using TensorKit
using MPSKit
using DynamicalCorrelators
using QuantumLattices
using JLD2: save, load
using MoireSuperlattices
using TightBindingApproximation


gfxts = Vector(undef, 10)
βs = [0.16, 0.18, 0.2, 0.22, 0.24, 0.28, 0.32, 0.4, 0.5, 0.68]
for i in 1:10
    gfxts[i] = load("./V=20/gfs/gfu_β=$(βs[i]).jld2", "gfu")
end
ts = 0:0.05:50

ws = 4.755330863071197*range(-12, 12, length=400)

gfrws = Vector(undef, 10)
for i in 1:10
    gfrws[i] = fourier_rw(gfxts[i], ts, ws; broadentype=(maximum(ts), "P"), ifsum=true)
    save("data/moire/gf_rw_β=$(βs[i]).jld2", "gfrw", gfrws[i], "ws", ws)
end


unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
lattice  = Lattice([0,0],[√3/2, 1/2], [√3/2, 3/2], [0, 1], [-√3/2, 1/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [0, 3], [-√3/2, 5/2], [-√3, 2]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)

hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

function moireHubbard(θ, Vᶻ, u, mu)
    parameters = (a₀=3.28, m=0.45, θ=θ, Vᶻ=Vᶻ, μ=8.31, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
    recipls = bltmd.frontend.reciprocallattice.translations
    lattice = MoireTriangular(1, reciprocals(recipls))
    brillouinzone = BrillouinZone(recipls, 24)
    t = terms(bltmd, lattice, brillouinzone; tol=10^-6)[1:end-1]
    μ = Onsite(:μ, -mu)
    return (t..., μ)
end  

term = moireHubbard(4.2, 20.0, 4.7, 10.0)
origiterms = term
referterms = term

Gs = Vector(undef, 10)
for i in 1:10
    gfrw2 = zeros(ComplexF64, 24, 24, length(ws))
    gfrw2[13:24,13:24,:] .= gfrws[i]
    cpt = CPT(unitcell, lattice, hilbert, origiterms, referterms, gfrw2)
    #k_path = ReciprocalPath(reciprocals([[√3/2, 1/2], [0, 1]]), hexagon"Γ-M-K-Γ, 120°", length=300)
    rz = ReciprocalZone(reciprocals([[√3/2, 1/2], [0, 1]]); length=100)
    @time Gs[i] = singleParticleGreenFunction(cpt, rz)
end
V = 20.0
As = [[[(-1/π)*(G[i][j][2,2]).im for i in 1:length(ws)] for j in 1:10000] for G in Gs]
save("Gkw_V=$(V).jld2", "Gs", Gs)
save("Akw_V=$(V).jld2", "As", As)