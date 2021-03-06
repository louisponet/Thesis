#%%
using DrWatson
quickactivate(@__DIR__)
using Revise
using DFControl
using DFWannier
const DFW = DFWannier
using LinearAlgebra
using Plots
using LaTeXStrings
using FileIO
pyplot()
#%%
## GeTe
# cd(datadir("Rashba/GeTe"))
job_nsoc = DFJob(datadir("Rashba/GeTe/NSOC"))
bands_nsoc = readbands(job_nsoc)
fermi_nsoc = readfermi(job_nsoc)
hami_nsoc = readhami(job_nsoc)
wfuncs_nsoc = load(joinpath(job_nsoc, "wfuncs.jld2"))["wfuncs"]
wbands_nsoc = wannierbands(hami_nsoc, bands_nsoc)
# plot(wbands_nsoc, bands_nsoc, ylims=[-2.1, 16])
p = plot(job_nsoc, -5.3, 5, 0.1, dpi=200, layout = grid(1, 2, widths=(0.7,0.3)))
plot!(p[1], xticks = ([1, 101, 201], ["A", "Z", "U"]))
plot!(p[2], legendfontsize = 10)
plot(bands_nsoc[11:20], fermi=fermi_nsoc, ylims=[-12.5,6], dpi=200, color=hcat([[:red for i=1:2];[:green for i =1:3];[:blue for i =1:5]]...), xticks = ([1, 101, 201], ["A", "Z", "U"]), title=nothing, xtickfontsize=15, ytickfontsize=15, yguidefontsize=15)
plot!([-10 for i = 1:length(bands_nsoc[1].k_points_cryst)])
savefig(papersdir("Rashba", "Images", "band_windows.png"))


DFWannier.write_xsf("test.xsf", wfuncs_nsoc[1], job_nsoc.structure, value_func = x -> sign(cos(x[1]))*abs(x[1]))
map(x -> abs(x[1])<1e-3 ? SVector((0.0im)) : x, wfuncs_nsoc[1].values)
# wfuncs_nsoc = WannierFunction[]
# for i = 1:8
#     push!(wfuncs1_soc, WannierFunction(wfuncs1[i].points, map(x-> SVector(x[1], 0.0), wfuncs1[i].values)))
#     push!(wfuncs1_soc, WannierFunction(wfuncs1[i].points, map(x-> SVector(0.0, x[1]), wfuncs1[i].values)))
# end

job_soc = DFJob(datadir("Rashba/GeTe/SOC"))
bands_soc = readbands(job_soc)
fermi_soc = readfermi(job_soc)
hami_soc = readhami(job_soc)
wfuncs_soc = load(joinpath(job_soc, "wfuncs.jld2"))["wfuncs"]
Sx_soc, Sy_soc, Sz_soc = DFW.readspin(job_soc)
wbands_soc = wannierbands(hami_soc, bands_soc)
p = plot(job_soc, -5.3, 5, 0.1, dpi=200, layout = grid(1, 2, widths=(0.7,0.3)))
plot!(p[1], xticks = ([1, 101, 201], ["A", "Z", "U"]))
plot!(p[2], legendfontsize = 10)
savefig(papersdir("Rashba", "Images", "SOC_dos.png"))
plot(bands_soc, fermi=fermi_soc, ylims=[-12.5,5])
plot(plot(wbands_nsoc, bands_nsoc, fermi=fermi_nsoc, ylims=[-5.3, 5], legend=:topright), plot(wbands_soc, bands_soc, fermi=fermi_soc, ylims=[-5.3, 5]), dpi=200, title = "",xticks = ([1, 101, 201], ["A", "Z", "U"]))
savefig(papersdir("Rashba", "Images", "wanvsdft.png"))


L = zeros(Vec3{ComplexF64}, 16, 16)
for (i1, w1) in enumerate(wfuncs_soc)
    for (i2,w2) in enumerate(wfuncs_soc)
        L[i1, i2] = DFW.calc_angmom(w1, w2, ustrip(job_soc.structure.atoms[1].position_cart), 2)
    end
end
evec=wbands_soc[10].eigvec[95]
bfunc = similar(wfuncs_soc[1])
for (i, e) in enumerate(evec)
    bfunc += e * wfuncs_soc[i]
end
DFW.calc_angmom(bfunc,bfunc, ustrip(job_soc.structure.atoms[1].position_cart), 3)
evec' * map(x->x[3], L) * evec



const smat = DFWannier.σx(2)
for ib = 1:10
evec = wbands_soc[ib].eigvec[95]
bfunc = similar(wfuncs_soc[1])
for (i, e) in enumerate(evec)
    bfunc += e * wfuncs_soc[i]
end
write_xsf(datadir("Rashba","GeTe","SOC","bf_95_$ib.xsf"), bfunc, job_soc.structure) 
end


wbands_soc[4].eigvals[101]
for i = 1:16
    DFW.write_xsf(joinpath(job1, "$i.xsf"), wfuncs1_soc[i], job1.structure, value_func = x -> abs(x[1]) > abs(x[2]) ? sign(real(x[1])) * norm(x) : sign(real(x[2]))*norm(x))
end
wbands[5].eigvec[101]
t = similar(wfuncs[1])
Lx, Ly, Lz = zeros(ComplexF64, 8, 8), zeros(ComplexF64, 8, 8), zeros(ComplexF64, 8, 8)
Sx, Sy, Sz = zeros(ComplexF64, 8, 8), zeros(ComplexF64, 8, 8), zeros(ComplexF64, 8, 8)
for i = 1:4, j = 1:4
    Sx[i,j], Sy[i,j], Sz[i,j] = DFW.calc_spin(wfuncs_soc[i], wfuncs_soc[j])
    Sx[i+4,j+4], Sy[i+4,j+4], Sz[i+4,j+4] = DFW.calc_spin(wfuncs_soc[i+4], wfuncs_soc[j+4])
end
save(joinpath(job_soc, "operators.jld2"), "Lx", Lx, "Ly", Ly, "Lz", Lz,"Sx", Sx, "Sy", Sy, "Sz", Sz, "hash", wan_hash(job_soc))
Lx + Sx
DFW.write_xsf(joinpath(job, "Z.xsf"), t, job.structure, value_func = x -> sign(real(x[1])) * norm(x))
#%%
bands[20].eigvals.-fermi
plot(bands[15].eigvals[90:110] .- fermi , ylims=[-0.8,-0.4], ylabel = L"E - E_f (eV)", xticks=([6, 12, 18], [L"A \leftarrow Z", L"Z", L"Z \rightarrow U"]), label="No SOC", color=:blue, xtickfontsize=15, ytickfontsize=15, yguidefontsize=20, linewidth=2, legendfontsize=15, legend=:bottom, linestyle=:dash, dpi=150)
plot!([(bands_soc[30].eigvals[90:110] .- fermi_soc) (bands_soc[29].eigvals[90:110] .- fermi_soc)], ylims=[-0.8,-0.4], ylabel = L"E - E_f (eV)", xticks=([6, 12, 18], [L"A \leftarrow Z", L"Z", L"Z \rightarrow U"]), label=["SOC" ""], color=:red, linewidth=2)
savefig("../papers/Rashba/Images/intro_dispersion.png")
#%%


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


## Potentials figure
#%%
using DrWatson
quickactivate(@__DIR__)
using Revise
using Glimpse

dio = Diorama()
grid = Glimpse.Grid([Point3f0(x, y, z) for x = -2:0.01:2, y=-2:0.01:2, z=-2:0.01:2])
density = map(x->1/sqrt(x[1]^2 + x[2]^2) + x[3], grid.points)
grid2 = Glimpse.Grid([Point3f0(x, y, z) for x = -2:0.01:6, y=-2:0.01:6, z=-2:0.01:2])
density2 = map(x-> -exp(-sqrt(x[1]^2 + x[2]^2)^2/0.50) + x[3], grid.points)
density3 = map(x->cos(x[1]/3) - x[3], grid2.points)
#%%
esph1 = Entity(dio, Glimpse.assemble_sphere(Point3f0(0.0,0.0,0.0), radius=0.5f0)...)
esph2 = Entity(dio, Glimpse.assemble_sphere(Point3f0(4.0,0.0,0.0), radius=0.5f0)...)
e1 = Entity(dio, grid, Glimpse.DensityGeometry(density, 0.4), Glimpse.UniformColor(RGB{Float32}(1.0,0.0,0.0)),Glimpse.Alpha(0.8f0), Glimpse.Material(), Glimpse.Spatial())  
e2 = Entity(dio, grid, Glimpse.DensityGeometry(density, 0.4), Glimpse.UniformColor(RGB{Float32}(1.0,0.0,0.0)),Glimpse.Alpha(0.8f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(4,0,0)))  
esph3 = Entity(dio, Glimpse.assemble_sphere(Point3f0(0.0,4.0,0.0), radius=0.5f0)...)
esph4 = Entity(dio, Glimpse.assemble_sphere(Point3f0(4.0,4.0,0.0), radius=0.5f0)...)
e3 = Entity(dio, grid, Glimpse.DensityGeometry(density, 0.4), Glimpse.UniformColor(RGB{Float32}(1.0,0.0,0.0)),Glimpse.Alpha(0.8f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(0,4,0)))  
e4 = Entity(dio, grid, Glimpse.DensityGeometry(density, 0.4), Glimpse.UniformColor(RGB{Float32}(1.0,0.0,0.0)),Glimpse.Alpha(0.8f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(4,4,0)))
e5 = Entity(dio, grid, Glimpse.DensityGeometry(density2, 0.4), Glimpse.UniformColor(RGB{Float32}(0.0,1.0,0.0)),Glimpse.Alpha(0.5f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(0,0,0.0)))
e6 = Entity(dio, grid, Glimpse.DensityGeometry(density2, 0.4), Glimpse.UniformColor(RGB{Float32}(0.0,1.0,0.0)),Glimpse.Alpha(0.5f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(4,0,0.0)))
e7 = Entity(dio, grid, Glimpse.DensityGeometry(density2, 0.4), Glimpse.UniformColor(RGB{Float32}(0.0,1.0,0.0)),Glimpse.Alpha(0.5f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(4,4,0.0)))
e8 = Entity(dio, grid, Glimpse.DensityGeometry(density2, 0.4), Glimpse.UniformColor(RGB{Float32}(0.0,1.0,0.0)),Glimpse.Alpha(0.5f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(0,4,0.0)))
e9 = Entity(dio, Glimpse.assemble_box(Point3f0(-4.0, -4.0, 2.0), Point3f0(20, 20, 2.01), color=RGB{Float32}(0.5,0.2,0.4))..., Glimpse.Rotation(Vec3(1.0,0.0,0.0), deg2rad(20)), Glimpse.Alpha(0.4f0))
e10 = Entity(dio, grid2, Glimpse.DensityGeometry(density3, 0.4), Glimpse.UniformColor(RGB{Float32}(1.0,1.0,0.0)),Glimpse.Alpha(0.5f0), Glimpse.Material(), Glimpse.Spatial(position=Point3f0(0,0,1.0)))
dio[Entity(2)] = Glimpse.Spatial(position=Point3f0(0,-10,10.0))
e11 = Entity(dio, Glimpse.assemble_wire_axis_box(position=Point3f0(-2,-2,-2), x=Vec3f0(4,0,0), y=Vec3f0(0,4,0), z=Vec3f0(0,0,4),color=Glimpse.BLACK)...)
e12 = Entity(dio, Glimpse.assemble_wire_axis_box(position=Point3f0(2,-2,-2), x=Vec3f0(4,0,0), y=Vec3f0(0,4,0), z=Vec3f0(0,0,4),color=Glimpse.BLACK)...)
e13 = Entity(dio, Glimpse.assemble_wire_axis_box(position=Point3f0(2,2,-2), x=Vec3f0(4,0,0), y=Vec3f0(0,4,0), z=Vec3f0(0,0,4),color=Glimpse.BLACK)...)
e14 = Entity(dio, Glimpse.assemble_wire_axis_box(position=Point3f0(-2,2,-2), x=Vec3f0(4,0,0), y=Vec3f0(0,4,0), z=Vec3f0(0,0,4),color=Glimpse.BLACK)...)
arrow_entities = map(x->Entity(dio, x...), Glimpse.assemble_axis_arrows(Point3f0(-2.5,-2.5,0), thickness=0.05f0, axis_length=2f0)[[1,3,4]])
dio[Glimpse.Camera3D][1].camerakind = Glimpse.Orthographic
expose(dio)
#%%
delete!(dio, e9)
dio[e9] = Glimpse.Alpha(0.8f0)
dio[e9] = Glimpse.Spatial(position=Point3f0(0,0,0))

dio[Glimpse.DioEntity]


all(x->x∈valid_entities(dio), arrow_entities)
dio[Glimpse.UniformColor]
for c in components(dio)
    if !isempty(c)
        @show eltype(c)
        @show c.indices
    end
end
dio[Mesh]
delete!(dio, arrow_entities
dio.ledger.free_entities
dio.ledger.components



using Glimpse
using ColorSchemes

N = 10

px(x, y, x0 = 0.0) = (x-x0) * exp(-sqrt((x-x0)^2 + y^2))
py(x, y, x0 = 0.0) = y * exp(-sqrt((x-x0)^2 + y^2))
px3(x, y, z, x0 = 0.0) = (x - x0) * exp(-sqrt((x-x0)^2 + y^2 +z^2))
py3(x, y, z, x0 = 0.0) = z * exp(-sqrt((x-x0)^2 + y^2 + z^2))

xrange = range(-2N,2N,length=600)
yrange = range(-N,N, length=300)
zrange = range(-N,N, length=300)
px_orb1 = [px(x, y, 7) for x in xrange, y in yrange]'
py_orb = [py(x, y) for x in xrange, y in yrange]'
px_orb2 = [px(x, y, -7) for x in xrange, y in yrange]'

p(K, r_part, args...) = exp(-1im*2pi * K)*(r_part * px(args...) + 1im*(1-r_part)*py(args...))
p3(K, args...) = exp(-1im*2pi * K)* (px3(args...) + 1im*py3(args...))
# p3(K, args...) = exp(-1im*2pi * K)* (px3(args...) +py3(args...))


const K = 0.1
const x0 = 4
wfc  = [p3(0.0,  x, y, z) for x in xrange, y in yrange, z in zrange]


wfc1 = [p3(K,  x, y, z,-x0) for x in xrange, y in yrange, z in zrange]
wfc2 = [p3(-K,  x, y, z,x0) for x in xrange, y in yrange, z in zrange]
sum(abs.(wfc1))
sum(abs.(wfc))
# dens = [(w = p3(K, x, y, z, -5) + p3(0.0, x, y, z) + p3(-K, x, y, z, 5); w' * w) * z for x in xrange, y in yrange, z in zrange]
dens = [(w = p3(2K, x, y, z, -2x0) + p3(K, x, y, z, -x0) + p3(0.0, x, y, z) + p3(-K, x, y, z, x0) + p3(-2K, x, y, z, 2x0); abs(w)/sqrt(5)) * z for x in xrange, y in yrange, z in zrange]
dio = Diorama()
col = (x) ->Glimpse.RGBf0(get(ColorSchemes.rainbow, (angle(p3(0.0, x[1], x[2], x[3]))+ pi)/2pi))
col1 = (x) ->Glimpse.RGBf0(get(ColorSchemes.rainbow, (angle(p3(K, x[1], x[2], x[3], -x0))+ pi)/2pi))
col2 = (x) ->Glimpse.RGBf0(get(ColorSchemes.rainbow, (angle(p3(-K, x[1], x[2], x[3], x0))+ pi)/2pi))
m = minimum(real.(dens))
M = maximum(real.(dens))
dipcol = (x) ->Glimpse.RGBf0(get(ColorSchemes.rainbow, (angle(p3(K, x[1], x[2], x[3], -x0) + p3(0.0, x[1], x[2], x[3]) + p3(-K, x[1], x[2], x[3], x0))+ pi)/2pi))
# dipcol = (x) ->p3(K, x, y, z, -5) + p3(0.0, x, y, z) + p3(-K, x, y, z, 5) > 0 ? Glimpse.RGBf0(1.0,0.0,0.0) : Glimpse.RGBf0(0.0,0.0,1.0)
grid = Glimpse.Grid([Point3f0(x, y, z) for x in xrange, y in yrange, z in zrange])

Entity(dio, Glimpse.DensityGeometry(Float32.(abs.(wfc)), 0.30),Glimpse.Spatial(), grid, Glimpse.Alpha(1.0f0), Glimpse.Material(), Glimpse.FunctionColor(col))
Entity(dio, Glimpse.DensityGeometry(Float32.(abs.(wfc1)), 0.30),Glimpse.Spatial(), grid, Glimpse.Alpha(1.0f0), Glimpse.Material(), Glimpse.FunctionColor(col1))
Entity(dio, Glimpse.DensityGeometry(Float32.(abs.(wfc2)), 0.30),Glimpse.Spatial(), grid, Glimpse.Alpha(1.0f0), Glimpse.Material(), Glimpse.FunctionColor(col2))
Entity(dio, Glimpse.DensityGeometry(Float32.(abs.(real.(dens))), 0.26),Glimpse.Spatial(), grid, Glimpse.Alpha(1.0f0), Glimpse.Material(), Glimpse.UniformColor(Glimpse.DEFAULT_COLOR), Glimpse.Alpha(0.6f0))

dio[Glimpse.Camera3D].data[1].camerakind = Glimpse.Orthographic
dio[Glimpse.Spatial][Entity(2)] = Glimpse.Spatial(position=Glimpse.Point3f0(200,-200,200))
expose(dio)

maximum(Float32.(abs.(wfc1)))

dio.loop=nothing
empty!(dio)
empty!(dio.ledger.free_entities)

k = 0.3
t1 = [p(0, sqrt(2)/2, x, y) * p(0, sqrt(2)/2, x, y)' for x in xrange, y in yrange] 
d1 = 30*[(p(k, sqrt(2)/2, x, y, -7)' * p(0, sqrt(2)/2, x, y) + p(-k, sqrt(2)/2, x, y,7)' * p(0, sqrt(2)/2, x, y)) *y for x in xrange, y in yrange]'
sum(d1)
k = -K
d2 = 30*[(p(k, sqrt(2)/2, x, y, -7)' * p(0, sqrt(2)/2, x, y) + p(-k, sqrt(2)/2, x, y,7)' * p(0, sqrt(2)/2, x, y)) *y for x in xrange, y in yrange]'

d1.+d2

d1

heatmap(xrange,yrange, real.(t1) )


K = 0.1
k = K
d1 = 30*[((exp(-k*2π*1im/7) * (px(x, y, 7) + 1im*py(x, y, 7))+ exp(k*2π*1im/7) * (px(x, y, -7) + 1im*py(x, y, -7)) *(1im*py(x, y)+px(x, y))) + (exp(-k*2π*1im/7)*(1im*py(x, y, 7) + px(x, y, 7))  + exp(k*2π*1im/7)*(1im* py(x, y, -7) + px(x, y, -7))*(px(x,y) + 1im*py(x,y))))*y for x in xrange, y in yrange]'
k = -K
d2 = 30*[((exp(-k*2π*1im/7) * (px(x, y, 7) + 1im*py(x, y, 7))+ exp(k*2π*1im/7) * (px(x, y, -7) + 1im*py(x, y, -7)) *(1im*py(x, y)+px(x, y))) + (exp(-k*2π*1im/7)*(1im*py(x, y, 7) + px(x, y, 7))  + exp(k*2π*1im/7)*(1im* py(x, y, -7) + px(x, y, -7))*(px(x,y) + 1im*py(x,y))))*y for x in xrange, y in yrange]'

heatmap(xrange, yrange, imag.(d1.+d2))

K = 0.7
k = K
d3 = 30*[((exp(-k*2π*1im/7) * (px(x, y, 7) + 1im*py(x, y, 7))+ exp(k*2π*1im/7) * (px(x, y, -7) + 1im*py(x, y, -7)) *(1im*py(x, y)+px(x, y))) + (exp(-k*2π*1im/7)*(1im*py(x, y, 7) + px(x, y, 7))  + exp(k*2π*1im/7)*(1im* py(x, y, -7) + px(x, y, -7))*(px(x,y) + 1im*py(x,y))))*y for x in xrange, y in yrange]'
k = -K
d4 = 30*[((exp(-k*2π*1im/7) * (px(x, y, 7) + 1im*py(x, y, 7))+ exp(k*2π*1im/7) * (px(x, y, -7) + 1im*py(x, y, -7)) *(1im*py(x, y)+px(x, y))) + (exp(-k*2π*1im/7)*(1im*py(x, y, 7) + px(x, y, 7))  + exp(k*2π*1im/7)*(1im* py(x, y, -7) + px(x, y, -7))*(px(x,y) + 1im*py(x,y))))*y for x in xrange, y in yrange]'

sum(d3 .+ d4 .- d2 .-d1)

