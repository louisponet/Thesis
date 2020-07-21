#%%
using DrWatson
quickactivate(@__DIR__)
using Revise
using DFWannier
const DFW = DFWannier
using LinearAlgebra
using Plots
using LaTeXStrings
#%%
## GeTe
cd(datadir("Rashba/GeTe"))
job = DFJob("../../GeTe/NSOC")
job_soc = DFJob("../../GeTe/SOC")

bands = readbands(job)
bands_soc = readbands(job_soc)
fermi = readfermi(job)
fermi_soc = readfermi(job_soc)
bands[20].eigvals.-fermi
plot(bands[15].eigvals[90:110] .- fermi , ylims=[-0.8,-0.4], ylabel = L"E - E_f (eV)", xticks=([6, 12, 18], [L"A \leftarrow Z", L"Z", L"Z \rightarrow U"]), label="No SOC", color=:blue, xtickfontsize=15, ytickfontsize=15, yguidefontsize=20, linewidth=2, legendfontsize=15, legend=:bottom, linestyle=:dash, dpi=150)
plot!([(bands_soc[30].eigvals[90:110] .- fermi_soc) (bands_soc[29].eigvals[90:110] .- fermi_soc)], ylims=[-0.8,-0.4], ylabel = L"E - E_f (eV)", xticks=([6, 12, 18], [L"A \leftarrow Z", L"Z", L"Z \rightarrow U"]), label=["SOC" ""], color=:red, linewidth=2)
savefig("../papers/Rashba/Images/intro_dispersion.png")
wgrid = DFW.read_points_from_xsf(joinpath(job, "wan_00001.xsf"))
wfuncs = [DFW.WannierFunction(joinpath(job, "wan_0000$i.xsf"), wgrid) for i = 1:8]


L_up = DFW.WannierFunction(wgrid, wfuncs[6]+1.0im*wfuncs[8])


#all shifting is done in the negative direction
function find_shift_indices(wfunction, cell)
    origin = wfunction.points[1,1,1]
    id1 = Tuple(findmin(map(x -> norm(x - (origin .+ ustrip.(cell[:,1]))), wfunction.points))[2]).-1
    id2 = Tuple(findmin(map(x -> norm(x - (origin .+ ustrip.(cell[:,2]))), wfunction.points))[2]).-1
    id3 = Tuple(findmin(map(x -> norm(x - (origin .+ ustrip.(cell[:,3]))), wfunction.points))[2]).-1
    return (id1, id2, id3)
end


function create_shifted_neg(wfunction, axis_id, cell)
    origin  = wfunction.points[1,1,1]
    id      = find_shift_indices(wfunction, cell)[axis_id]
    shift   = ustrip.(cell[:,axis_id])
    shifted = zeros(wfunction)

    shift_norm = norm(shift)
    for i in CartesianIndices(wfunction.points)
        if i[axis_id] <= id[axis_id]
            shifted.values[(Tuple(i) .+ 2 .* id)...] = wfunction.values[i]
        else
            shifted.values[(Tuple(i) .- id)...] = wfunction.values[i]
        end
    end
    return shifted
end

function create_shifted_pos(wfunction, axis_id, cell)
    origin  = wfunction.points[end,end,end]
    id      = find_shift_indices(wfunction, cell)[axis_id]
    shift   = ustrip.(cell[:,axis_id])
    shifted = similar(wfunction)

    shift_norm = norm(shift)
    for i in CartesianIndices(wfunction.points)
        if i[axis_id] > 2*id[axis_id]
            shifted.values[(Tuple(i) .- 2 .* id)...] = wfunction.values[i]
        else
            shifted.values[(Tuple(i) .+ id)...] = wfunction.values[i]
        end
    end
    return shifted
end

shifted_neg = create_shifted_neg(L_up, 2, cell(job))
shifted_neg = DFW.WannierFunction(L_up.points, exp(-2π*1.0im * 0.6) * create_shifted_neg(wfuncs[8], 2, cell(job)))
shifted_pos = DFW.WannierFunction(L_up.points, exp(2π*1.0im * 0.6) * create_shifted_pos(wfuncs[8], 2, cell(job)))

DFW.calc_dip(shifted_neg, L_up)+DFW.calc_dip(shifted_pos, L_up)

sum(shifted_neg.values .- )
sum(L_up.values)
dot(L_up, L_up)
dot(shifted,L_up)
shifted_pos = create_shifted_pos(L_up, 2, cell(job))
sum(shifted_pos)
DFW.write_xsf_file(joinpath(job, "testwan2.xsf"), shifted_pos, job.structure)
DFW.calc_dip(shifted, L_up)

norm(DFW.calc_dip(shifted_neg, shifted_neg) .- DFW.calc_dip(L_up, L_up))
norm(cell(job)[:,2])
DFW.calc_dip(shifted_pos, L_up) + DFW.calc_dip(shifted_neg, L_up)

#%%
## HfO2
using Plots
pyplot()
cd(datadir("Rashba/HfO2"))
job = DFJob("SOC")

bands = readbands(job)
fermi = readfermi(job)
hami = readhami(job)
wbands = wannierbands(hami, bands)

plot(wbands, bands, fermi=fermi, ylims=[-5,10])
bands[1].k_points_cryst[1]
ops = load(joinpath(job, "Operators.jld2"))
Smat, Lmat = ops["S"], ops["L"]

findlast(x->maximum(x.eigvals)<fermi+1, wbands)
plot(wbands[49])
plot(wbands[48])
plot!(wbands[47])
plot(wbands[49])
plot!(wbands[52])
wbands[1].kpoints_cryst[155]
#%%

tot_S = zero(DFW.Vec3{ComplexF64})
for (i1,v1) in enumerate(wbands[52].eigvec[161]), (i2,v2) in enumerate(wbands[52].eigvec[161])
    global tot_S += v1' * v2 * Smat[i1, i2]
end

tot_S = zero(DFW.Vec3{ComplexF64})
for (i1,v1) in enumerate(wbands[51].eigvec[161]), (i2,v2) in enumerate(wbands[51].eigvec[161])
    global tot_S += v1' * v2 * Smat[i1, i2]
end

tot_S = zero(DFW.Vec3{ComplexF64})
for (i1,v1) in enumerate(wbands[49].eigvec[155]), (i2,v2) in enumerate(wbands[49].eigvec[155])
    global tot_S += v1' * v2 * Smat[i1, i2]
end

tot_S


#%%
using Plots
wbands = wannierbands(hami, bands)
length(bands[1].k_points_cryst)
plot(wbands[4].eigvals)

plot(abs.(getindex.(wbands[5].eigvec[90:110],1)))
plot!(abs.(getindex.(wbands[6].eigvec[90:110],1)))






