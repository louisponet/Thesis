#%%
using DrWatson
@quickactivate()
using Revise
using BTO
using Plots
pyplot()
using LaTeXStrings
#%%
f=0.5
xsize = 10e-9
ysize = 10e-9
zsize = 10e-9
elsize = 0.5e-9
tipsize = 1.0e-9
cs = [1, 2, 5, 10]
name = "noflexo_noq44_force_$(f)_$(xsize)x$(ysize)x$(zsize)x$(elsize)_$(tipsize)"
params    = ModelParams(xsize       = xsize,
                        ysize       = ysize,
                        zsize       = zsize,
                        elementsize = elsize,
                        q44=0.0,
                        F11=0.0,
                        F12=0.0,
                        F44=0.0)

simparams = SimParams(workdir       = datadir()*"$name/",
                      maxsteps      = 10000,
                      tipradius     = tipsize,
                      minimizer=BTO.ConjugateGradient,
                      centers       = cs,
                      forceconstant = f*0.5*1e1,
                      timelimit     = 150*60,
                      decayrate     = 2e-9,
                      force_potential = BTO.element_potential_uforce)

m = simmodel(params,simparams)
using BenchmarkTools
@btime F($(m[1].dofs), $(m[1]))

@btime BTO.element_potential($(rand(24)), $(m[1].threadcaches[1].cellvalues), $(m[1].threadcaches[1].extradata), $(BTO.Ffunc(params)), 1.0) 

simulation(params, simparams)


### Peierls Nabarro barriers
all_energies = Vector{Float64}[]
meshsizes = [round(i*0.1, digits=3) for i in 24:2:36]
for j = meshsizes
    mp = ModelParams(xsize = 10e-9, ysize=1e-9, zsize=5e-9, elementsize = j*1e-10, F11 =0.0, F12 = 0.0, F44 =0.0)
    sim = SimParams(tipradius = 1e-10, forceconstant = 2.0e29,  gtol=1e-5, centers = [0.0, 1.0e-9, 2.0e-9, 5.0e-9],minimizer=BTO.ConjugateGradient, workdir=datadir("test2"), timelimit=1000000000000, maxsteps=1000)
    m, dofs  = simmodel(mp, sim);
    res = optimize(m)
    m.dofs.=res.minimizer
    BTO.save_converged_model(mp, m, m.dofs)
    tree = BTO.Landau.BallTree(m)
    new_dofs = zeros(length(m.dofs))
    energies = Float64[]
    for i = 0.1:0.1:10
        Threads.@threads for n in m.dofnodes
            v = Vec((n.coord .+ Vec(i*1e-10,0.0,0.0))...)
            # if v[1] > 5e-9
            #     new_dofs[n.dofs[:u]] .= m.dofs[n.dofs[:u]]
            #     new_dofs[n.dofs[:P]] .= m.dofs[n.dofs[:P]]
            # else
            # end
            d = BTO.Landau.extract_data(m, [v], m.dofs, tree)
            new_dofs[n.dofs[:u]] .= d[1].u
            new_dofs[n.dofs[:P]] .= d[1].P
        end
        push!(energies, F(new_dofs, m))
    end
    push!(all_energies, energies)
end
for i in all_energies
    i .-= minimum(i)
    i .*= 1e-20
    i .*= 1000
    i ./= (1e-9 * 5e-9)
end
plot(plot([i for i = 0.1:0.1:10], all_energies, yguide = L"E- E_0\,{\rm (mJ/m^2)}", xguide = L"$\delta_w$ (Å)", label=hcat(["$i Å" for i in meshsizes]...), legendtitle = "Elementsize"), plot(meshsizes, maximum.(all_energies) .- minimum.(all_energies), legend=false, yguide = L"E_{PN}\,{\rm (mJ/m^2)}", xguide="Elementsize (Å)"), layout=(2,1), dpi=200, size=(1260, 1260), margin=30Plots.mm)
savefig(papersdir("Softening","Images","peierls_nabarro.png"))

## Breathing mode
mp = ModelParams(xsize = 10e-9, ysize=1e-9, zsize=5e-9, elementsize = 3e-10, F11 =0.0, F12 = 0.0, F44 =0.0)
sim = SimParams(tipradius = 10e-9, forceconstant = 2.0e29,  gtol=1e-5, centers = [0.0, 1.0e-9, 2.0e-9, 5.0e-9],minimizer=BTO.ConjugateGradient, workdir=datadir("test2"), timelimit=1000000000000, maxsteps=1000)
m, dofs  = simmodel(mp, sim);
res = optimize(m)
m.dofs .= res.minimizer

tree = BTO.Landau.BallTree(m)
force = zeros(length(m.dofs))
minfunc = generate_minizing_function(sim,
                                     (x, y) -> F(x, y, force),
                                     (x, y, z) -> ∇F!(x, y, z, force),
                                     (x, y, z) -> ∇²F!(x, y, z, force))
force_model = LandauModel(m, BTO.fill_sim_params(sim, BTO.fill_model_params(mp, sim.force_potential)))
centers = BTO.find_force_center.((force_model,), sim.centers)
left, right  = BTO.left_right(mp)
xyz = range(Vec(left[1], 0.0, right[3]), Vec(right[1], 0.0, right[3]), length=Int(ceil(mp.xsize/mp.elementsize)))
center=xyz[18]
P_res = Vector{Float64}[]
for f in 1:8:80
    sim = SimParams(tipradius = 1e-9, forceconstant = f*1e28,  gtol=1e-5, centers = [0.0, 1.0e-9, 2.0e-9, 5.0e-9],minimizer=BTO.ConjugateGradient, workdir=datadir("test2"), timelimit=1000000000000, maxsteps=1000)
    force_model.dofs = copy(m.dofs)
    BTO.calculateforce!(force,
    				force_model.dofhandler.grid.nodes,
    				center,
    				sim.forceconstant,
    				sim.tipradius,
    				sim.decayrate)
	res = minfunc(force_model)
	vtk_save("test$f.vtu", force_model, res.minimizer)
	push!(P_res, map(x->x.P[3], BTO.extract_data(force_model, xyz, res.minimizer, tree)))
end

#%%
plot(map(x->abs.(x), P_res))
