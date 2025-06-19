module PowerLawsET
using Statistics 

include("reach_et.jl")

using .reach_et #bring in the functions from the reach_et module
export main_sequence, replicate_pub_whk_sim, plot_vt_surfaces #export for now only the plot_vt_surfaces function

end