function moireVCA(θ, Vᶻ, u)
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
    cluster = Lattice([0,0],[√3/2, 1/2], [-√3/2, 1/2], [0, 1], [√3/2, 3/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [-√3, 2], [-√3/2, 5/2], [0, 3]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)
    hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
    bs = Sector(hilbert, ParticleNumber(12))
    t = moiretuning(θ, Vᶻ)
    U = Hubbard(:U, 4.755330863071197*u)
    origiterms = (t..., U)
    referterms = (t..., U)
    return VCA(:S, unitcell, cluster, hilbert, origiterms, referterms, bs; modelname="Moire-MIM", cachepath="/fsa/home/jxl_zongyy/data/VCA_Moire_MIM/")
end
#14.2

V = range(0, 20, 21)
unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
cluster = Lattice([0,0],[√3/2, 1/2], [-√3/2, 1/2], [0, 1], [√3/2, 3/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [-√3, 2], [-√3/2, 5/2], [0, 3]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)
rz = ReciprocalZone(reciprocals(cluster.vectors); length=100)


mus = [ 14.279933110367892
14.249832775919732
14.149498327759197
14.049163879598662
13.848494983277591
13.617725752508361
13.32675585284281
13.01571906354515
12.734782608695653
12.544147157190636
12.493979933110367
12.554180602006689
12.664548494983277
12.403678929765887
11.851839464882943
11.320066889632107
11.581270903010033
11.581270903010033
11.581270903010033
11.581270903010033
11.681270903010033]



function mindens(vca, arr)
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
    rz = ReciprocalZone(reciprocals(unitcell.vectors); length=200)
    dens = Vector{Float64}(undef, length(arr))
    for i in eachindex(arr)
        dens[i] = densityofstates(singleParticleGreenFunction(:f, vca, rz, 0:1:0; η=0.08*4.755330863071197, μ=arr[i]))[1]
    end
    dmin_index = argmin(dens)
    return arr[dmin_index]
end

function SDFG(vca, μ)
    function specre(gfpathv::AbstractVector; select::AbstractVector=Vector(1:size(gfpathv[1][1],2)))
        A = zeros(Float64, length(gfpathv), length(gfpathv[1]))
            for i in eachindex(gfpathv)
                for j in eachindex(gfpathv[i])
                    A[i, j] = (det(gfpathv[i][j][select, select])).re
                end
            end
            return A
    end
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
    rz = ReciprocalZone(reciprocals(unitcell.vectors); length=200)
    dz = ReciprocalZone(reciprocals([[0.5,0],[0,0.5]]); length=1000)
    k_path = ReciprocalPath(reciprocals(unitcell.vectors), hexagon"K-Γ-M₂-K, 120°", length=300)
    ω_range = 4.755330863071197*range(-10, 10, length=500)
    G = singleParticleGreenFunction(:f, vca, k_path, ω_range; η=0.08*4.755330863071197, μ=μ)
    S = spectrum(G)
    Sd = spectrum(G; select=[1])
    Su = spectrum(G; select=[2])
    GG = singleParticleGreenFunction(:f, vca, rz, ω_range; η=0.08*4.755330863071197, μ=μ)
    D = densityofstates(GG)
    GGG = singleParticleGreenFunction(:f, vca, dz, 0:1:0; η=0.08*4.755330863071197, μ=μ)
    F = reshape(spectrum(GGG), Int(√(length(dz))), Int(√(length(dz))))
    Fd = reshape(spectrum(GGG; select=[1]), Int(√(length(dz))), Int(√(length(dz))))
    Fu = reshape(spectrum(GGG; select=[2]), Int(√(length(dz))), Int(√(length(dz))))
    GSd = specre(G; select=[1])
    GSu = specre(G; select=[2])
    GFd = reshape(specre(GGG; select=[1]), Int(√(length(dz))), Int(√(length(dz))))
    GFu = reshape(specre(GGG; select=[2]), Int(√(length(dz))), Int(√(length(dz))))
    return S, Su, Sd, D, F, Fu, Fd, GSu, GSd, GFu, GFd
end

function QPI(vca,μ)
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
    rz = BrillouinZone(reciprocals(unitcell.vectors); length=100)
    return qpi = δ_QPI(:f, vca, rz, 0, 0.05; μ=μ, η=0.05)
end