using Pkg
Pkg.activate("../.")
using TensorKit
using MPSKit
using DynamicalCorrelators
using QuantumLattices
using JLD2: save, load
using MoireSuperlattices
using TightBindingApproximation


function moireHubbard(θ, Vᶻ, u, mu)
    parameters = (a₀=3.28, m=0.45, θ=θ, Vᶻ=Vᶻ, μ=8.31, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
    recipls = bltmd.frontend.reciprocallattice.translations
    lattice = MoireTriangular(3, reciprocals(recipls))
    brillouinzone = BrillouinZone(recipls, 24)
    t = terms(bltmd, lattice, brillouinzone; tol=10^-6)[[1,2,5,6]]
    U = Hubbard(:U, 4.755330863071197*u)
    μ = Onsite(:μ, -mu)
    return (t..., U, μ)
end  
mu =  10.0
term = moireHubbard(4.2, 20.0, 4.7, mu)
# lattice  = Lattice([0,0],[√3/2, 1/2], [√3/2, 3/2], [0, 1], [-√3/2, 1/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [0, 3], [-√3/2, 5/2], [-√3, 2]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)
lattice  = Lattice([0,0],[√3/2, 1/2], [√3/2, 3/2], [0, 1], [-√3/2, 1/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [0, 3], [-√3/2, 5/2], [-√3, 2])

hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

H = hamiltonian(term, lattice, hilbert)

elt = ComplexF64
L = 12
filling = (1, 1)
ψ = randFiniteMPS(elt, U1Irrep, U1Irrep, L; md=50, filling=filling)
truncs = [truncdim(d) for d in [1024, 1024, 1024, 1024, 1024, 1024]]
ϵs = [1.0 for _ in 1:L]
for i in 1:length(truncs)
    if i < 2
        pinning = nothing
        mydmrg =  DefaultDMRG2(1e-6, 3)
    elseif i < 3
        pinning = nothing
        mydmrg = DefaultDMRG2(1e-6, 3)
    elseif i < 5
        pinning = nothing
        mydmrg = DefaultDMRG2(1e-12, 8)
    else
        pinning = nothing
        mydmrg = DefaultDMRG2(1e-14, 12)
    end
    dmrg2_sweep!(i, ψ, 
                H,
                truncs[i], ϵs;
                alg = mydmrg,
                filename="tWSe2_V=20_GS.jld2", 
                verbose=2);
end