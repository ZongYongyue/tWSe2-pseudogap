using Pkg
Pkg.activate(".")
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
mu =  10
term = moireHubbard(4.2, 20.0, 4.7, mu)
# lattice  = Lattice([0,0],[√3/2, 1/2], [√3/2, 3/2], [0, 1], [-√3/2, 1/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [0, 3], [-√3/2, 5/2], [-√3, 2]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)
lattice  = Lattice([0,0],[√3/2, 1/2], [√3/2, 3/2], [0, 1], [-√3/2, 1/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [0, 3], [-√3/2, 5/2], [-√3, 2])

hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

H = hamiltonian(term, lattice, hilbert)
N = hamiltonian((Onsite(:μ, 1), ), lattice, hilbert)

rs = [0.08, 0.16, 0.18, 0.2, 0.22, 0.24, 0.28, 0.32, 0.4, 0.5, 0.68, 1]/2
ts = [0:0.04:0.08; [0.09, 0.1, 0.11, 0.12, 0.14, 0.16, 0.2, 0.25]; 0.3:0.04:0.34; 0.36:0.02:0.5]

function find_id(rs, ts)
    id = []
    for (i, t) in enumerate(ts)
        if t in rs
            push!(id, i)
        end
    end
    return id
end

rho = evolve_mps(H, -im*round.(ts, digits=2); filename="./rhos/rho_mu=$(mu).jld2", n=2, trscheme=truncdim(512), save_id=find_id(rs, ts))
rhos = load("./rhos/rho_mu=$(mu).jld2")
Ns = map(find_id(rs, ts)) do i
    rho = rhos["t=$(-im*round.(ts, digits=2)[i])"]
    dot(rho, FiniteMPO(N), rho)/dot(rho,rho)
end
println(Ns)