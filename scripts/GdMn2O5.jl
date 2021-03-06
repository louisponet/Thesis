#%%
using Revise
using DrWatson
@quickactivate()
using GdMn2O5
const Gd = GdMn2O5
using LinearAlgebra
using LaTeXStrings
using ProgressMeter
using FileIO
using PStdLib.Images
pyplot()
#%%

l  = Ledger(best_model())
rad2deg(l[L][2].ϕ)
optimize!(l, Spin, L, Pb)
const start = 3.0
const stop  = 5.5
const len   = 120
const Hr    = H_sweep_range(start, stop, len)
const hsweep = generate_Hsweep_states(l, 10, Hr)
const N = 152
const margin = 2.1 * Plots.Measures.mm * 2
const BLACK = RGB(0, 0, 0)
const DARK_GREEN = RGB(0, 0.5, 0)
const DARK_BLUE = RGB(0, 0, 0.5)
const LIGHT_GREEN = RGB(0.1, 1, 0.1) 
const LIGHT_BLUE = RGB(0.1, 0.1, 1.0)
const RED = RGB(1, 0, 0)
const PATH_COLOR = RGB(240 / 255, 92 / 255, 7 / 255)
const PATH_COLOR = RGB(1.0, 1.0, 1.0)
function plot_3angle_Pb(l; a1 = 9.80, a2 = 10, a3 = 10.35)
    optimize(l, Spin, L)
    Hr = H_sweep_range(0, 10, 100)[1:end - 99]
    h_sweep_state0 = generate_Hsweep_states(l, a1, Hr)

    h_sweep_state11 = generate_Hsweep_states(l, a2, Hr)
    h_sweep_state20 = generate_Hsweep_states(l, a3, Hr)
    colors = [fill(BLACK, 99);fill(DARK_GREEN, 99); fill(DARK_BLUE, 99); fill(LIGHT_GREEN, 99)]
    p3 = plot(Hr[1:201], Pb(h_sweep_state0[1:201]),  color = [fill(BLACK, 99);fill(LIGHT_GREEN, 101)])
    p4 = plot(Hr, Pb(h_sweep_state11), color = colors);
    p5 = plot(Hr[1:201], Pb(h_sweep_state20[1:201]), color = [fill(BLACK, 99);fill(LIGHT_GREEN, 101)]);
    return p3, p4, p5
end

function calculate_Etots(hsweep...)
    p = Progress(length(hsweep))
    L1_ϕs = range(0, 2π, length = N)
    L2_ϕs = range(0, 2π, length = N)
    E_tots = [zeros(N, N) for i = 1:length(hsweep)]
    for ih = 1:length(hsweep)
        thread_ls = [deepcopy(hsweep[ih]) for t = 1:Threads.nthreads()]
        Threads.@threads for i1 = 1:N
            for (i2, ϕ2) in enumerate(L2_ϕs)
                ϕ1 = L1_ϕs[i1]
                tl = thread_ls[Threads.threadid()]
                tl[L].data[1] = GdMn2O5.L(ϕ = ϕ1) 
                tl[L].data[2] = GdMn2O5.L(ϕ = ϕ2)
                optimize!(tl, Spin, Pb)
                E_tots[ih][i1,i2] = calc_etot(tl).e
            end
        end
        next!(p)
    end
    return E_tots
end
function plot_heatmap(h_ledger, hsweep, L1_traj_shift, L2_traj_shift, offset = 1.0)
    ϕ1_r = range(0, 2π, length = N)
    ϕ2_r = range(0, 2π, length = N)
    E_tots =calculate_Etots(h_ledger)[1]
    E_tots = log.(E_tots.-minimum(E_tots).+offset)
    cs = :rainbow
    p = heatmap(ϕ1_r, ϕ2_r, E_tots, colorbar_title = "E", color = cs, aspect_ratio = 1, legend = false)
    plot!(p, map(x->x[L][2].ϕ + L2_traj_shift, hsweep), map(x->x[L][1].ϕ + L1_traj_shift, hsweep), color = PATH_COLOR)
    scatter!(p, [h_ledger[L][2].ϕ + L2_traj_shift], [h_ledger[L][1].ϕ + L1_traj_shift], color = PATH_COLOR, markersize = 4)
    return p
end
#%%
hms = [hsweep[1]; hsweep[60:10:120]]
Etots = calculate_Etots(hms...)
ϕ1_r = range(0, 2π, length = 152)
ϕ2_r = range(0, 2π, length = 152)

fontsize = 25 * 1
plts = []
for (i, E) in enumerate(Etots)
    # push!(plts, heatmap(ϕ1_r, ϕ2_r, E,
    #         colorbar_title = "E",
    #         color          = :rainbow,
    #         aspect_ratio   = 1,
    #         linewidth      = 5,
    #         markersize     = 15,
    #         framestyle     = :box,
    #         legend=false))
    push!(plts, heatmap(ϕ1_r, ϕ2_r, E,
            title = "|H| = $(round(norm(singleton(hms[i], H).v), digits=2))", 
            color          = :rainbow,
            xtickfontsize  = fontsize,
            ytickfontsize  = fontsize,
            yguidefontsize = fontsize,
            xguidefontsize = fontsize,
            titlefontsize  = fontsize,
            yticks         = ([ϕ1_r[76], ϕ1_r[end]], [L"\pi", L"2\pi"]),
            ylabel         = L"ϕ_{L_1}",
            xlabel         = L"ϕ_{L_2}",
            ylims          = [0,2π],
            xticks         = ([ϕ1_r[1],ϕ1_r[76], ϕ1_r[end]], [" 0 ",L"\pi", L"2\pi"]),
            xlims          = [0,2π],
            aspect_ratio   = 1,
            framestyle     = :box,
            legend=false))
    # savefig(plts[i], "GdMn2O5/Images/field_heatmap$i.png")
    # clip_image("GdMn2O5/Images/field_heatmap$i.png")
end
plts[1]
plot(plts[1], colorbar=true,colorbarticks=nothing)
scatter!(plts[1], [ϕ1_r[67], ϕ1_r[67+75], ϕ1_r[67+75], ϕ1_r[67]],[ϕ2_r[81], ϕ2_r[81], ϕ2_r[81-75], ϕ2_r[81-75]],
         color=:white,
         markersize=20,
        yticks         = ([π, 2π], ["π", "2π"]),
        ylabel         = L"ϕ_{L_1}",
        xlabel         = L"ϕ_{L_2}",
        ylims          = [0,2π],
        xticks         = ([ϕ1_r[1],ϕ1_r[76], ϕ1_r[end]], [" 0 "," pi ", " 2pi "]),
        xlims          = [0,2π],
        aspect_ratio   = 1,
         framestyle=:box, legend=false)
savefig(plts[1], "GdMn2O5/Images/field_heatmap_1.png")
clip_image("GdMn2O5/Images/field_heatmap_1.png")
plts[2]
plot(plts...)
findfirst(x->x==π, ϕ1_r)
#%%
Entity(l, GdMn2O5.VisualizationSettings(Gd_color=GdMn2O5.Gl.GREEN, spin_arrow_thickness=0.08f0, spin_arrow_length=0.8f0, Gd_sphere_radius=0.4f0, Mn_text_offset=Vec3f0(0.4,0.4,0.0), Gd_text_offset=Vec3f0(1.4,1.4,0.0)))
l[H][1].ϕ = deg2rad(10)
dio = Diorama(l)
dio[GdMn2O5.OptimizationSettings][1].continuous = true
expose(dio)
dio.loop

#%%

e_Gd_H = []
e_L_H = []
e_LL = []
e_L2H2 = []
e_LK = []
e_Gd_L= []

for h in hsweep
    h[Etot][1].e = 0.0
    update(E_Gd_H(), h)
    push!(e_Gd_H, h[Etot][1].e)
    h[Etot][1].e = 0.0
    update(E_L_H(), h)
    push!(e_L_H, h[Etot][1].e)
    h[Etot][1].e = 0.0
    update(Gd.E_Gd_L(), h)
    push!(e_Gd_L, h[Etot][1].e)
    h[Etot][1].e = 0.0
    update(E_L_easy(), h)
    push!(e_LK, h[Etot][1].e)
    h[Etot][1].e = 0.0
    update(Gd.E_LL(), h)
    push!(e_LL, h[Etot][1].e)
    h[Etot][1].e = 0.0
    update(Gd.E_L2H2(), h)
    push!(e_L2H2, h[Etot][1].e)
end

plot(Hr, e_Gd_H, yguide="E (meV)", dpi=200, xguide="|H| (T)", label="Gd Zeeman", yguidefontsize=15,xguidefontsize=15, xtickfontsize=15, ytickfontsize=15, legendfontsize=10)
plot!(Hr, e_L_H, label="L Zeeman")
plot!(Hr, e_Gd_L.+60, label = L"\left[{\rm Gd} \leftrightarrow {\rm L}\right] + {\rm 60\,meV}")
plot!(Hr, e_LL, label = L"{\rm L_1} \leftrightarrow {\rm L_2}")
plot!(Hr, e_LK, label = "L easy-axis anisotropy")
plot!(Hr, e_L2H2, label = "L2H2")
plot(e_LL.+e_LK.+e_Gd_L .+ e_L2H2)


# savefig(papersdir("GdMn2O5/Images/energy_contributions.png"))
#%%





