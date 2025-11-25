using Pkg
Pkg.activate("../.")
using TensorKit
using MPSKit
using DynamicalCorrelators
using QuantumLattices
using JLD2: save, load
using MoireSuperlattices
using TightBindingApproximation
using JLD2:save, load
using Distributed

# i = parse(Int, ARGS[1])

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

term = moireHubbard(4.2, 20.0, 4.7, 10)
# lattice  = Lattice([0,0],[√3/2, 1/2], [√3/2, 3/2], [0, 1], [-√3/2, 1/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [0, 3], [-√3/2, 5/2], [-√3, 2]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)
lattice  = Lattice([0,0],[√3/2, 1/2], [√3/2, 3/2], [0, 1], [-√3/2, 1/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [0, 3], [-√3/2, 5/2], [-√3, 2])

hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

H = hamiltonian(term, lattice, hilbert)
N = hamiltonian((Onsite(:μ, 1), ), lattice, hilbert)

elt = Float64
epu = e_plus(elt, U1Irrep, U1Irrep; spin=:up)
epd = e_plus(elt, U1Irrep, U1Irrep; spin=:down)
emu = e_min(elt, U1Irrep, U1Irrep; spin=:up)
emd = e_min(elt, U1Irrep, U1Irrep; spin=:down)

ts = 0:0.05:50
n = 2
trscheme=truncerr(1e-3)
# D = 512
# trscheme=truncdim(D)

mu = 10
bs = [0.08, 0.16, 0.18, 0.2, 0.22, 0.24, 0.28, 0.32, 0.4, 0.5, 0.68, 1]
rs = bs/2



addprocs(24)

@everywhere begin
    using Pkg
    Pkg.activate("../.")
    using TensorKit
    using MPSKit
    using DynamicalCorrelators
    using QuantumLattices
    using JLD2: save, load
    using MoireSuperlattices
    using TightBindingApproximation
    using JLD2:save, load
    using Distributed

    function moireHubbard(θ, Vᶻ, u, mu)
        parameters = (a₀=3.28, m=0.45, θ=θ, Vᶻ=Vᶻ, μ=8.31, V=-1.28, ψ=22.7, w=-12.9)
        bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
        recipls = bltmd.frontend.reciprocallattice.translations
        lattice = MoireTriangular(3, reciprocals(recipls))
        brillouinzone = BrillouinZone(recipls, 24)
        t = terms(bltmd, lattice, brillouinzone; tol=10^-6)[1:end-1]
        U = Hubbard(:U, 4.755330863071197*u)
        μ = Onsite(:μ, -mu)
        return (t..., U, μ)
    end  

end


# rho = load("./rhos/rho_mu=$(mu).jld2", "t=$(-im*round.(rs[i], digits=2))")
# gfu = dcorrelator(rho, H, (epu, emu); 
#                     trscheme=trscheme, 
#                     times=ts, 
#                     beta=bs[i], 
#                     n=n, 
#                     verbose=true, 
#                     gf_path="./gfu/beta=$(bs[i])",
#                     rho_path="./rhos/beta=$(bs[i])"
#                     );

# save("gfu_β=$(bs[i]).jld2", "gfu", gfu)

# gfu = nothing
# GC.gc()

# gfd = dcorrelator(rho, H, (epd, emd); 
#                     trscheme=trscheme, 
#                     times=ts, 
#                     beta=bs[i], 
#                     n=n, 
#                     verbose=true, 
#                     gf_path="./gfd/beta=$(bs[i])",
#                     rho_path="./rhos/beta=$(bs[i])"
#                     );

# save("gfd_β=$(bs[i]).jld2", "gfd", gfd)

gs = load("tWSe2_V=20_GS.jld2", "sweep_6_ψ")

gfu = dcorrelator(gs, H, (epu, emu);
                    verbose=true, 
                    gf_path="./gfu/beta=Inf",
                    times=ts, 
                    n=n, 
                    trscheme=trscheme
                    )

save("gfu_β=Inf.jld2", "gfu", gfu)

# @everywhere begin
#     GC.gc()
# end

# gfd = dcorrelator(gs, H, (epd, emd);
#                     verbose=true, 
#                     gf_path="./gfd/beta=Inf",
#                     times=ts, 
#                     n=n, 
#                     trscheme=trscheme
#                     )

# save("gfd_β=Inf.jld2", "gfd", gfd)