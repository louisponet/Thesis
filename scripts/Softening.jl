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
