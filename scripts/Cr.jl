#%%
using DrWatson
@quickactivate()
using Revise
using CrSDW
using Plots
sets=load_experimental_data()
using LaTeXStrings
pyplot()
using Plots.Measures: mm
#%%


hm = CrSDW.cdw_heatmap(sets[11])
#%%
plot(
     plot(plot(sets[11], label=["Experiment" "Theory"], xticks=[], bottom_margin = 10mm, legendfontsize=5),
          plot(sets[6], legend=false, xlabel = "t (ps)", bottom_margin = 18mm),
          ylabel = L"y/y_0",
          left_margin=20mm,
          layout=(2,1)),
     plot(heatmap(hm..., ylabel="pump-pump delay (ps)", xlabel="pump-probe delay (ps)", right_margin=20mm),
          plot(load("/home/ponet/Documents/PhD/Thesis/papers/Cr/Images/exp_heatmap.jpg")),
          left_margin=28mm,
          layout=(2,1)),
     xguidefontsize=5,
     yguidefontsize=5,
     xtickfontsize=5,
     ytickfontsize=5,
     dpi=200)

savefig("/home/ponet/Documents/PhD/Thesis/papers/Cr/Images/Theory_Fit.pdf")
#%%
mp = get_model_params(Cb = 1e6, ξ=0.04)
envelope1 = x -> x < 2π ? 0.5*sin(x) : ( x < 3π ? 0.5/π*(x-2π) : -0.5/π*(x-3π) + 0.5)
fitrange = -1.0:0.01:4π
t1 = fit_signal(envelope1, fitrange, params=mp, optim_params=get_optim_params(iterations=500),threaded=true)
fitrange = -1.0:0.002:4π
# plot(fitrange, t1, mp, envelope)
# savefig(papersdir("Cr/Images/fit1.pdf"))
# mp = get_model_params(Cb = 1e6, ξ=0.04)
envelope2 = x -> x < 2π ? 0.5*CrSDW.gaussian(x, π, π/3) : ( x < 3π ? 0.5/2π*(exp((x-2π)/(π/log(2π+1)))-1) : -0.5/2π*(exp((x-3π)/(π/log(2π+1)))-1) + 0.5)
fitrange = -1.0:0.01:4π
t2 = fit_signal(envelope2, fitrange, params=mp, optim_params=get_optim_params(iterations=500), threaded=true)
fitrange = -1.0:0.002:4π
# plot(fitrange, t2, mp, envelope)

plot(plot(fitrange, t1, mp, envelope1), plot(fitrange, t2, mp, envelope2, legend=nothing), ylims=[-1.5,2.5])
savefig(papersdir("Cr/Images/fits.pdf"))

#%%
## For maximum amplitude
envelope = x -> 4.5
fitrange = -1.0:0.01:4π
t1 = fit_signal(envelope, fitrange, params=mp, optim_params=get_optim_params(iterations=500),threaded=true)
t2 = fit_signal(envelope, fitrange, params=mp, optim_params=get_optim_params(iterations=500),threaded=true, max_amplitude=true)
fitrange = -1.0:0.002:4π
plot(plot(fitrange, t1, mp, envelope), plot(fitrange, t2, mp, envelope, legend=nothing))
savefig(papersdir("Cr/Images/maximumamplitude.pdf"))
#%%

plot(plot(sets[3], dpi=200, xticks = ([0, 2, 4, 6, 8], ["", "","","",""]), yguide = "PLD Amplitude (a.u.)"), plot(sets[1], dpi=200, legend=nothing, xticks = ([0, 2, 4, 6, 8], ["", "","","",""]), yticks = ([-0.5, 0.0, 0.5, 1.0, 1.5], ["","","","",""])), plot(sets[6], yguide = "PLD Amplitude (a.u.)", dpi=200, legend=nothing, xguide = "Time (ps)",), plot(sets[7], dpi=200, legend=nothing, xguide = "Time (ps)", yticks = ([-0.5, 0.0, 0.5, 1.0, 1.5], ["","","","",""])), ylims=[-0.5,1.5])
savefig(papersdir("Cr/Images/exp_fits.pdf"))
#%%
