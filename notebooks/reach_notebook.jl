### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ a1b2c3d4-e5f6-7890-abcd-ef1234567890
begin
    using Pkg
	if !isdefined(Main, :PowerLawsET)
        Pkg.activate(joinpath(@__DIR__, ".."))
        Pkg.instantiate()
    end
	using GLMakie
    using GLM 
    using DataFrames
    using MAT
	using PlutoUI
    # Configure WGLMakie for better Pluto performance
    WGLMakie.activate!(inline=true, format="svg", px_per_unit=2)
	using PowerLawsET
end

# ╔═╡ b2c3d4e5-f607-4901-bcde-f23456789012
begin
    f, mainseq, fit = replicate_pub_whk_sim(fname="reachgrid_ctmech.mat")
    f
end

# ╔═╡ c3d4e5f6-a708-9012-cdef-345678901234
function overlay_subject_data(fig, mainseq::main_sequence; mrkr=:circle, color=:red, size=10)
    """
    Overlays subject data points on the existing surface plots.
    
    Parameters:
    -----------
    fig : GLMakie.Figure
        The figure with surfaces created by plot_vt_surfaces
    mainseq : main_sequence
        
    ms : Symbol
        Marker symbol to use for data points
    color : Symbol or RGB
        Color for the subject data points
    size : Number
        Size of the markers
        
    Returns:
    --------
    fig : Figure
        The updated figure with overlaid data points
    """
        
    # Extract relevant fields
    subject_distances = mainseq.distance
    subject_durations = mainseq.duration  
    subject_peak_speeds = mainseq.peak_speed  
    # ct is one number per subject, so we need to repeat it for each distance
    subject_ct = ones(length(mainseq.distance))*mainseq.time_valuation
    
    # Get the existing axes from the figure
    axes = fig.content
    duration_axis = axes[1]
    peak_speed_axis = axes[3] 
    
    # Add data points to the duration surface
    GLMakie.scatter!(duration_axis, 
                    subject_ct, 
                    subject_distances, 
                    subject_durations,
                    markersize=size, 
                    color=color, 
                    marker=mrkr)
    
    # Add data points to the peak speed surface
    GLMakie.scatter!(peak_speed_axis, 
                    subject_ct, 
                    subject_distances, 
                    subject_peak_speeds,
                    markersize=size, 
                    color=color, 
                    marker=mrkr)
    
    return fig
end

# ╔═╡ d4e5f6a7-b8c9-0123-defa-456789012345
md"""
### Function: PowerLaw fit to subject data.

"""

# ╔═╡ e5f6a7b8-c9d0-1234-efab-567890123456
function fit_ct_for_vt(k_v,k_t,ms::main_sequence;norm_v=1.0,norm_t=1.0,only_v1_t2=0,M=1)
  """
  Fit the main sequence data (V, T) simultaneously via the two powerlaws for V and T, each of which depend on Ct and D. 
  #Arguments:
    k_v : float64
      the constant for the peak_speeds powerlaw
    k_t : float64
      the constant for the durations powerlaw
    ms : main_sequence object
      the main sequence data
    Optional:
    norm_v : float64
      value to normalize the peak_speeds
    norm_t : float64
      value to normalize the durations
    only_v1_t2 : int
      if 1, fit only the peak_speeds, if 2, fit only the durations
    
  #Returns:
      ms: main_sequence
        the main sequence data with the time_valuation, 
        predicted peak_speeds and durations, and R² values 
    
  """
  # normalize the peak_speeds and durations
  V  = ms.peak_speed
  T = ms.duration
  L = ms.distance

  # # Define power law functions
  # V, peak_speed.
    fn_v = (c_t, L) -> k_v*c_t.^(1/4) .* ((1/M)^(1/4)) .* L.^(3/4)
    
    A_v = (k_v/norm_v) * L.^(3/4)
    b_v = (V ./ norm_v)

    # T, duration. note: c_t is in the denominator.
    fn_t = (c_t, L) -> k_t*(1 ./ c_t).^(1/4) .* ((M)^(1/4)) .* L.^(1/4)

    A_t = T/norm_t
    b_t = (k_t/norm_t) * L.^(1/4)
    
    # stack T and V. if we want to fit only one, do it here.
    A = [A_t; A_v]
    b = [b_t; b_v]

  if only_v1_t2 == 1
      A = [A_v]
      b = [b_v]
  elseif only_v1_t2 == 2
      A = [A_t]
      b = [b_t]
  end

  # Convert to DataFrame (equivalent to MATLAB's table)
  valid_indices = .!isnan.(vec(A)) .& .!isnan.(vec(b))
  df = DataFrame(A=vec(A)[valid_indices], b=vec(b)[valid_indices])
  #display A to screen to debug
  print("A: ")
  println(A)
  print("b: ")
  println(b)

  # Use GLM.jl for linear modeling 
  # In Julia GLM, this is done with the @formula syntax
  
  model = lm(@formula(b ~ A - 1), df)

  # Extract coefficient and raise to power 4 
  # In MATLAB: stats.Coefficients.Estimate(1)
  # In Julia: coef(model)[1]
  CTstar = (coef(model)[1])^4  # we solve for c_t^(1/4)
  #
  # Calculate R² value
  # In MATLAB: stats.Rsquared.Ordinary
  # In Julia: r2(model)
  r2_val = r2(model)
  # generate predictions using fn_v and fn_t
  peak_speed_star = fn_v(CTstar, L)
  duration_star = fn_t(CTstar, L)
  time_valuation = CTstar
  ms_star = main_sequence(peak_speed = peak_speed_star, duration = duration_star, distance = L, time_valuation = time_valuation)
  return (ms_star, r2_val)
end

# ╔═╡ 6530ee37-78c7-45cc-8fb5-80237728ea4b
md"""
## Step 2: fit experiment condition data
Everyday, Fast, Slow and Preferred movements from each subject get fit with a condition-specific $c_t$.  
"""

# ╔═╡ f6a7b8c9-d0e1-2345-fabc-678901234567
begin
    k_v = fit.k_v
    k_t = fit.k_t
    
	# loop from subject # 1 to 9, fit_ct_for_vt for each
    
	# pre-allocate main_sequence vector of length 9 
    main_seq_star_vec = Vector{main_sequence}(undef, 9)
    main_seq_vec = Vector{main_sequence}(undef, 9)
    R2_vec = zeros(9)
    conditions = ["eday","fast","slow","infp"]
    for i in 1:9
        #construct filename
        
        fname = "ms_"*conditions[1]*"_subject_"*string(i)*".mat"
        dir   = dirname(dirname(@__FILE__))
        fpath = joinpath(dir,"data","reaching","subjects",fname)
        data  = matread(fpath)
        ms    = main_sequence(peak_speed = data["peak_speed"], duration = data["duration"], distance = data["distance"],time_valuation = 0.0)
        
        ms_fit = fit_ct_for_vt(k_v,k_t,ms,norm_v=.5,norm_t=1.0,only_v1_t2=0)
        main_seq_star_vec[i] = ms_fit[1]

		#add the estimated time_valuation to the main_sequence
		main_seq_vec[i]     = main_sequence(peak_speed = data["peak_speed"], duration = data["duration"], distance = data["distance"], time_valuation = ms_fit[1].time_valuation)
		
        R2_vec[i] = ms_fit[2]
    end;
	
end

# ╔═╡ a7b8c9d0-e1f2-3456-abcd-789012345678
md"""
## Subject Data Overlay - All Subjects
Shows experimental data (red crosses) and model fits (green circles) for all 9 subjects overlaid on the power law surfaces.
"""

# ╔═╡ f2a3b4c5-d6e7-8901-bcde-234567890123
begin
    # Create a copy of the base figure and overlay ALL subject data
    fig_all = f
    
    # Loop through all subjects and overlay their data
    for i in 1:9
        # Overlay model predictions (green circles) 
        overlay_subject_data(fig_all, main_seq_star_vec[i], mrkr=:circle, color=:green, size=8)
        # Overlay actual data (red crosses)
        overlay_subject_data(fig_all, main_seq_vec[i], mrkr=:xcross, color=:red, size=6)
    end
    
    println("Overlaid data for all 9 subjects:")
    for i in 1:9
        println("Subject $i: R² = $(round(R2_vec[i], digits=3))")
    end
    
    fig_all
end

# ╔═╡ b8c9d0e1-f2a3-4567-bcde-890123456789
fit.allfit

# ╔═╡ c9d0e1f2-a3b4-5678-cdef-901234567890
fit.k_t

# ╔═╡ 3d90dfdc-32ba-41c7-aea2-061e97490785
main_seq_vec[1][1:5]

# ╔═╡ Cell order:
# ╠═a1b2c3d4-e5f6-7890-abcd-ef1234567890
# ╠═b2c3d4e5-f607-4901-bcde-f23456789012
# ╠═c3d4e5f6-a708-9012-cdef-345678901234
# ╟─d4e5f6a7-b8c9-0123-defa-456789012345
# ╠═e5f6a7b8-c9d0-1234-efab-567890123456
# ╟─6530ee37-78c7-45cc-8fb5-80237728ea4b
# ╠═f6a7b8c9-d0e1-2345-fabc-678901234567
# ╟─a7b8c9d0-e1f2-3456-abcd-789012345678
# ╠═f2a3b4c5-d6e7-8901-bcde-234567890123
# ╠═b8c9d0e1-f2a3-4567-bcde-890123456789
# ╠═c9d0e1f2-a3b4-5678-cdef-901234567890
# ╠═3d90dfdc-32ba-41c7-aea2-061e97490785
