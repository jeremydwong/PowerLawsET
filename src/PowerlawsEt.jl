module PowerlawsEt
using Statistics 

include("reach_et.jl")

using .reach_et #bring in the functions from the reach_et module
export plot_vt_surfaces, plot_vt_surfaces2 #export for now only the plot_vt_surfaces function

end