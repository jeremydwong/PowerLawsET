### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ a1b2c3d4-e5f6-7890-abcd-ef1234567890
begin
    using PowerLawsE
    using Revise
    using GLMakie
    using GLM 
    using DataFrames
    using MAT
    using Pkg
    GLMakie.activate!(inline=true)
    # activate the package. We are in /src.
    if !isdefined(Main, :PowerLawsET)
        Pkg.activate(joinpath(@__DIR__, ".."))
        Pkg.instantiate()
    end
    Pkg.activate("..")
end

# ╔═╡ b2c3d4e5-f6g7-8901-bcde-f23456789012
begin
    f, mainseq, fit = replicate_pub_whk_sim("reachgrid_ctmech.mat")
    f
end

# ╔═╡ c3d4e5f6-g7h8-9012-cdef-345678901234
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

# ╔═╡ d4e5f6g7-h8i9-0123-defg-456789012345
md"""
### Subject PowerLaw fit 
fn_v = fit.
"""

# ╔═╡ e5f6g7h8-i9j0-1234-efgh-567890123456
function fit_ct_2vandt(k_v,k_t,ms::main_sequence;norm_v=1.0,norm_t=1.0,only_v1_t2=0,M=1)
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

  # use GLM.jl for linear modeling instead of MATLAB's fitlm
  # Fit linear model without intercept (the -1 in MATLAB notation)
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

# ╔═╡ f6g7h8i9-j0k1-2345-fghi-678901234567
begin
    k_v = fit.k_v
    k_t = fit.k_t
    # loop from subject # 1 to 9, fit_ct_2vandt
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
        data = matread(fpath)
        ms = main_sequence(peak_speed = data["peak_speed"], duration = data["duration"], distance = data["distance"],time_valuation = 0.0)
        main_seq_vec[i] = ms
        ms_fit = fit_ct_2vandt(k_v,k_t,ms,norm_v=.5,norm_t=1.0,only_v1_t2=0)
        main_seq_star_vec[i] = ms_fit[1]
        R2_vec[i] = ms_fit[2]
    end
end

# ╔═╡ g7h8i9j0-k1l2-3456-ghij-789012345678
begin
    # loop across the two ms vectors and plot the human data and model fits
    for i in 1:9
        fig = overlay_subject_data(f, main_seq_star_vec[i], mrkr=:circle, color=:green, size=10)
        display(f)
    end
end

# ╔═╡ h8i9j0k1-l2m3-4567-hijk-890123456789
fit.allfit

# ╔═╡ i9j0k1l2-m3n4-5678-ijkl-901234567890
fit.k_t

# ╔═╡ Cell order:
# ╠═a1b2c3d4-e5f6-7890-abcd-ef1234567890
# ╠═b2c3d4e5-f6g7-8901-bcde-f23456789012
# ╠═c3d4e5f6-g7h8-9012-cdef-345678901234
# ╟─d4e5f6g7-h8i9-0123-defg-456789012345
# ╠═e5f6g7h8-i9j0-1234-efgh-567890123456
# ╠═f6g7h8i9-j0k1-2345-fghi-678901234567
# ╠═g7h8i9j0-k1l2-3456-ghij-789012345678
# ╠═h8i9j0k1-l2m3-4567-hijk-890123456789
# ╠═i9j0k1l2-m3n4-5678-ijkl-901234567890