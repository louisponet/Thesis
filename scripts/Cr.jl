#%%
using DrWatson
cd("/home/ponet/Documents/PhD/CrSDW")
using Pkg
Pkg.activate(pwd())
using Revise
using CrSDW
using Plots
sets=load_sets()
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
