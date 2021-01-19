using DrWatson
@quickactivate()
using Revise
using BTO

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
mp = ModelParams(xsize = 10e-9, ysize=1e-9, zsize=5e-9, elementsize = 3e-10, F11 =0.0, F12 = 0.0, F44 =0.0)
sim = SimParams(tipradius = 1e-10, forceconstant = 2.0e29,  gtol=1e-5, centers = [0.0, 1.0e-9, 2.0e-9, 5.0e-9],minimizer=BTO.ConjugateGradient, workdir=datadir("test2"), timelimit=1000000000000, maxsteps=1000)
m, dofs  = simmodel(mp, sim);
res = optimize(m)
m.dofs.=res.minimize
tree = BTO.Landau.BallTree(m)
new_dofs = zeros(length(m.dofs))
energies = Float64[]
for i = 0.1:0.1:10
    for n in m.dofnodes
        v = Vec((n.coord .+ Vec(i*1e-10,0.0,0.0))...)
        d = BTO.Landau.extract_data(m, [v], m.dofs, tree)
        new_dofs[n.dofs[:u]] .= d[1].u
        new_dofs[n.dofs[:P]] .= d[1].P
    end
    push!(energies, F(new_dofs, m))
end
plot(energies)
r
