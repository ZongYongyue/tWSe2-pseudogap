using Makie
using CairoMakie
using JLD2: save, load
using Colors

ts = abs.([4.6+1.8*im, 4.55+2.2*im, 4.47+2.6*im, 4.4+3*im, 4.3+3.36*im, 4.17+3.8*im])
e = 1.602176634*(10^(-19))
ħ = 1.054571817*(10^(-34))
S = (3.28*(10^(-10)))^2
N = 10000
C = (2π*e^2)/(ħ*N*S)

Cs = [1/C/e/t  for t in ts]

cols = ["#0b2e6e"
"#1f6fb2"
"#6eadd9"
"#c7daee"
"#fee0d4"
"#fbbba1"
"#f86d44"
"#af4c2f"
"#68000d"]

colors = reverse([parse(Colorant, col) for col in cols][[1,2,3,6,7,9]])

re_data = load("resistivity.jld2")
bs = re_data["betas"]
re9 = re_data["V=9"]
re11 = re_data["V=11"]
re13 = re_data["V=13"]
re15 = re_data["V=15"]
re17 = re_data["V=17"]
re19 = re_data["V=19"]

ys = [re9, re11, re13, re15, re17, re19].*Cs/100

set_theme!(Theme(fontsize = 15, rowgap = 1, colgap = 3))
#  86/ 0.352778/0.4063 = 600
# 86/ 0.352778/0.75 = 325
fig = Figure(size = (600, 400))

main_grid = fig[1:1, 1:1] 


ax = Axis(main_grid[1,1];
    xgridvisible = false,
    ygridvisible = false,
    xticksize = 5, 
    yticksize = 5,
    xtickalign = 1.0, 
    ytickalign = 1.0,
    # xminorticksvisible = true, 
    # yminorticksvisible = true, 
    yminortickalign = 1.0,
    yminorticksize = 3,
    # xminorticks = IntervalsBetween(5),
    # xticks = ([0.4, 0.6, 0.8, 1.0 ,1.2, 1.4]),
    # xticks = ([1, 2, 3, 4,5, 6, 7,8,9, 10, 11],["1","2","","","5","","","","","10",""]),
    xlabel = "T/t",
    xlabelsize = 16,
    xlabelpadding = -0, 
    ylabel = "ρ(kΩ)",
    ylabelsize = 16,
    ylabelpadding = 1, 
    )
ylims!(ax, 0.0, 8)
xlims!(ax, 0.0, 1.27)
xs = [1/b for b in bs]/5

using LsqFit

xs2 = range(0,0.55,50)

@. model(T, p) = p[1] * exp(p[2] / T) + p[3]
p0 = [1.0, 1.0, 0.0]
id = [8,9,10]
g = LsqFit.curve_fit(model, xs[id], ys[1][id], p0)
fit1 = lines!(ax, xs2, model(xs2, g.param),color =colors[1], linestyle=:dash)

# @. model(T, p) = p[1] * exp(p[2] / T) + p[3]
# p0 = [1.0, 1.0, 0.0]
# id = [8,9,10]
# g = LsqFit.curve_fit(model, xs[id], ys[2][id], p0)
# fit2 = lines!(ax, xs2, model(xs2, g.param),color =colors[2], linestyle=:dash)

@. model(T, p) = p[1] * exp(p[2] / T) + p[3]
p0 = [1.0, 1.0, 0.0]
id = [8,9,10]
g = LsqFit.curve_fit(model, xs[id], ys[3][id], p0)
fit3 = lines!(ax, xs2, model(xs2, g.param),color =colors[3], linestyle=:dash)

xs2 = range(0,1.0,50)

@. model(T, p) = p[1] * T + p[2]
p0 = [1.0, 0.0]
id = [6,7,8,9,10]
g = LsqFit.curve_fit(model, xs[id], ys[4][id], p0)
fit4 = lines!(ax, xs2, model(xs2, g.param),color =colors[4], linestyle=:dash)

# @. model(T, p) = p[1] * T^2 + p[2]
# p0 = [1.0, 0.0]
# id = [3,4,5,6,7,8,9,10]
# g = LsqFit.curve_fit(model, xs[id], ys[5][id], p0)
# fit5 = lines!(ax, xs2, model(xs2, g.param),color =colors[5], linestyle=:dash)

@. model(T, p) = p[1] * T^2 + p[2]
p0 = [1.0, 0.0]
id = [3,4,5,6,7,8,9,10]
g = LsqFit.curve_fit(model, xs[id], ys[6][id], p0)
fit5 = lines!(ax, xs2, model(xs2, g.param),color =colors[6], linestyle=:dash)


lins = []
scas = []

for i in 1:6
    lin = lines!(ax, xs, ys[i], linewidth = 2, color = colors[i])
    push!(lins, lin)
    sca = scatter!(ax, xs, ys[i], color = colors[i], strokecolor=colors[i],markersize = 6, strokewidth=1)
end

x = 0.92
y = [1.8,1.5,1.2,0.9,0.6,0.3]*0.1

Legend(main_grid[1,1],
[lins[1]],
[""], 
    tellwidth = false,
    tellheight = false,
    halign = x,
    valign = y[1],
    labelsize = 12,
    framevisible = false)

Legend(main_grid[1,1],
[lins[2]],
[""], 
    tellwidth = false,
    tellheight = false,
    halign = x,
    valign = y[2],
    labelsize = 12,
    framevisible = false)

Legend(main_grid[1,1],
[lins[3]],
[""], 
    tellwidth = false,
    tellheight = false,
    halign = x,
    valign = y[3],
    labelsize = 12,
    framevisible = false)

Legend(main_grid[1,1],
[lins[4]],
[""], 
    tellwidth = false,
    tellheight = false,
    halign = x,
    valign = y[4],
    labelsize = 12,
    framevisible = false)

Legend(main_grid[1,1],
[lins[5]],
[""], 
    tellwidth = false,
    tellheight = false,
    halign = x,
    valign = y[5],
    labelsize = 12,
    framevisible = false)

Legend(main_grid[1,1],
[lins[6]],
[""], 
    tellwidth = false,
    tellheight = false,
    halign = x,
    valign = y[6],
    labelsize = 12,
    framevisible = false)

text!(ax, 1.08, 1.7, text = "9meV", fontsize = 15)
text!(ax, 1.08, 0.2, text = "19meV", fontsize = 15)

display(fig)


# save("result/resistivity_tn.pdf", fig,pt_per_unit=0.4063)
