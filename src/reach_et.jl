module reach_et
"""
reach_et.jl: Functions for simulating reaching movements with energy optimization.

"""

using Interpolations
using Statistics
using ColorSchemes
using MAT
using GLMakie

#define main_sequence struct to hold
# duration
# distance
# peak speed
# time valuation
struct main_sequence
    distance::Array{Float64,2}
    duration::Array{Float64,2}
    peak_speed::Array{Float64,2}
    #time_valuation can either be an array, or a float64
    time_valuation::Union{Array{Float64,2}, AbstractFloat}
  function main_sequence(;distance, duration, peak_speed, time_valuation)
    new(distance, duration, peak_speed, time_valuation)
  end
end

mutable struct fit_
    fn_vk::Function
    fn_tk::Function
    k_v::Float64
    k_t::Float64
    allfit::Dict{String, Any}
end

export main_sequence, replicate_pub_whk_sim, plot_vt_surfaces, fit_powerlaws_to_oc

"""
function replicate_whk()
  loads default data, plots surface, and returns main sequence object

"""
function replicate_pub_whk_sim(;fname = "reachgrid_ctmech.mat")
  # make simfile location relative to this file. 
  base = dirname(dirname(@__FILE__))
  # concatenate
  simfile = joinpath(base, "data","reaching",fname)
  data      = matread(simfile)
  distance  = data["distance"]
  duration  = data["duration"]
  peak_speed = data["peak_speed"]
  timeValuation = data["time_valuation"]
  mainseq = main_sequence(distance=distance, duration=duration, peak_speed=peak_speed,time_valuation=timeValuation)

  f=plot_vt_surfaces(mainseq)
  ret_fit = fit_powerlaws_to_oc(mainseq)
  return (f, mainseq, ret_fit)
end

function plot_vt_surfaces(ms::main_sequence, az=27, el=33.6,plotstyle="interactive")
"""
function plot_vt_surfaces(ms::main_sequence; simfile="data/reachgrid.mat", az=27, el=33.6,plotstyle="GLMakie")
  
    Plot the surfaces of distance, duration, and peak speed as functions of time valuation and distance.

    Parameters
    ----------
    simfile : str
        Path to the MAT file containing the simulation results.
    f : int
        Figure number for the plot.
    az : float
        Azimuth angle for the 3D plot.
    el : float
        Elevation angle for the 3D plot.
    plotstyle : str
        Plotting style: "Plots" or "GLMakie".

    Returns
    -------
    fig : Figure
        The generated figure object.
"""
  # Handle default parameters:

  # view
  vclps = [az, el]

  # Colors in RGB format (0-1)
  # cs = range(RGB(239/255, 237/255, 245/255), RGB(117/255, 107/255, 177/255), length=100)
  
  # Fill missing values in the grid
  (ct, distance, duration)  = simple_grid_fillmissing(ms.time_valuation, ms.distance, ms.duration)
  (_, _, peakspeed)         = simple_grid_fillmissing(ms.time_valuation, ms.distance, ms.peak_speed)
  
  # Add zero columns
  distance  = hcat(zeros(size(distance, 1)), distance)
  duration  = hcat(zeros(size(distance, 1)), duration)
  peakspeed = hcat(zeros(size(distance, 1)), peakspeed)
  ct = hcat(ct[:, 1], ct)
  
  ct_mech = ct

  # Plotting
      # Create a new figure
    fig = Figure(;size=(1000, 700), fontsize=16)
    
    # Add a 3D axis for the duration surface
    ax1 = Axis3(fig[1, 1], 
               xlabel="Time valuation Cₜ (W)",
               ylabel="Distance L (m)",
               zlabel="Duration T (s)",
               title="Duration")
    
    # Create the surface
    surf1 = GLMakie.surface!(ax1, ct_mech, distance, duration,
               colormap=:plasma,
               transparency=true,
               alpha=0.85,
               colorrange=(0, 2.5))
    
    # Set limits
    GLMakie.xlims!(ax1, 0, 15)
    GLMakie.ylims!(ax1, 0, 0.55)
    GLMakie.zlims!(ax1, 0, 2.5)
    
    # Add colorbar
    Colorbar(fig[1, 2], surf1, label="Duration (s)")
    
    # Add a 3D axis for the peak speed surface
    ax2 = Axis3(fig[2, 1], 
               xlabel="Time valuation Cₜ (W)",
               ylabel="Distance L (m)",
               zlabel="Peak Speed V (m/s)",
               title="Peak Speed")
    
    # Create the surface
    surf2 = GLMakie.surface!(ax2, ct_mech, distance, peakspeed,
               colormap=:plasma,
               transparency=true,
               alpha=0.85,
               colorrange=(0, 1.1))
    
    # Set limits
    GLMakie.xlims!(ax2, 0, 15)
    GLMakie.ylims!(ax2, 0, 0.55)
    GLMakie.zlims!(ax2, 0, 1.1)
    
    # Add colorbar
    Colorbar(fig[2, 2], surf2, label="Peak Speed (m/s)")
    
    # Construct 2 cameras based on the ax1 and ax2 scenes
    cam1 = Camera3D(ax1.scene)
    cam2 = Camera3D(ax2.scene)

    rotate_cam!(ax1.scene, cam1, (deg2rad(el), deg2rad(az), 0))
    rotate_cam!(ax2.scene, cam2, (deg2rad(el), deg2rad(az), 0))
    
    if plotstyle == "interactive"
        
    
    # Add camera control sliders
    azimuth_slider = Slider(fig[3, 1], range=0:1:360, startvalue=az)
    elevation_slider = Slider(fig[4, 1], range=0:1:90, startvalue=el)
    
    # Labels for sliders
    Label(fig[3, 1, Top()], "Azimuth")
    Label(fig[4, 1, Top()], "Elevation")
    
    # Connect sliders to camera rotation
    on(azimuth_slider.value) do val
        rotate_cam!(ax1.scene, cam1, deg2rad(elevation_slider.value[]), deg2rad(val), 0)
        rotate_cam!(ax1.scene, cam2, deg2rad(elevation_slider.value[]), deg2rad(val), 0)
    end
    
    on(elevation_slider.value) do val
        rotate_cam!(ax1.scene, cam1, deg2rad(val), deg2rad(azimuth_slider.value[]), 0)
        rotate_cam!(ax1.scene, cam2, deg2rad(val), deg2rad(azimuth_slider.value[]), 0)
    end
    end
    return fig
end

function simple_grid_fillmissing(xs, ys, vals)
  """
  A straightforward fill of missing values in a 2D grid surface.
  This function assumes data is organized in a grid structure and handles
  missing coordinates and values (either can be 0s or nans).
  
  Parameters:
  -----------
  xs : Array
      x-coordinates of the data points
  ys : Array
      y-coordinates of the data points
  vals : Array
      Values at each (x,y) point
      
  Returns:
  --------
  grid_xs : Array
      Filled x-coordinates in grid structure
  grid_ys : Array
      Filled y-coordinates in grid structure
  grid_vals : Array
      Filled values in grid structure
  """
  # Create copies to avoid modifying originals
  grid_xs = copy(xs)
  grid_ys = copy(ys)
  grid_vals = copy(vals)
  
  # Step 1: Find all unique valid x and y values to define our grid
  flat_xs = grid_xs[:]
  flat_ys = grid_ys[:]
  
  valid_xs = unique(filter(x -> x != 0.0 && !isnan(x), flat_xs))
  valid_ys = unique(filter(y -> y != 0.0 && !isnan(y), flat_ys))
  
  sort!(valid_xs)
  sort!(valid_ys)
  
  # Step 2: If we have a complete, uniform grid, then each row should have the same x values
  # and each column should have the same y values
  
  # Get dimensions of our grid
  nx = length(valid_xs)
  ny = length(valid_ys)
  
  # Create a new, perfectly regular grid
  perfect_grid_xs = zeros(nx, ny)
  perfect_grid_ys = zeros(nx, ny)
  perfect_grid_vals = zeros(nx, ny)
  
  # Fill the perfect grid with coordinate values
  for i in 1:nx
      for j in 1:ny
          perfect_grid_xs[i, j] = valid_xs[i]
          perfect_grid_ys[i, j] = valid_ys[j]
          perfect_grid_vals[i, j] = NaN  # Initially all values are NaN
      end
  end
  
  # Step 3: Transfer known values from original grid to perfect grid
  for i in 1:size(grid_xs, 1)
      for j in 1:size(grid_xs, 2)
          # Skip if either coordinate is missing
          if (grid_xs[i, j] == 0.0 || isnan(grid_xs[i, j]) || 
              grid_ys[i, j] == 0.0 || isnan(grid_ys[i, j]) ||
              isnan(grid_vals[i, j]))
              continue
          end
          
          # Find the position in our perfect grid
          x_index = findfirst(x -> isapprox(x, grid_xs[i, j], atol=1e-10), valid_xs)
          y_index = findfirst(y -> isapprox(y, grid_ys[i, j], atol=1e-10), valid_ys)
          
          if !isnothing(x_index) && !isnothing(y_index)
              perfect_grid_vals[x_index, y_index] = grid_vals[i, j]
          end
      end
  end
  
  # Step 4: Interpolate missing values in the perfect grid
  # We'll use linear interpolation along rows and columns
  
  # First interpolate along rows (x direction)
  for j in 1:ny
      row_values = perfect_grid_vals[:, j]
      valid_indices = findall(!isnan, row_values)
      
      if length(valid_indices) >= 2  # Need at least 2 points for interpolation
          # For each NaN value, linearly interpolate
          for i in 1:nx
              if isnan(row_values[i])
                  # Find nearest valid points to the left and right
                  left_indices = valid_indices[valid_indices .< i]
                  right_indices = valid_indices[valid_indices .> i]
                  
                  if !isempty(left_indices) && !isempty(right_indices)
                      left_idx = maximum(left_indices)
                      right_idx = minimum(right_indices)
                      
                      left_x = valid_xs[left_idx]
                      right_x = valid_xs[right_idx]
                      left_val = perfect_grid_vals[left_idx, j]
                      right_val = perfect_grid_vals[right_idx, j]
                      
                      # Linear interpolation
                      x = valid_xs[i]
                      t = (x - left_x) / (right_x - left_x)
                      perfect_grid_vals[i, j] = left_val * (1 - t) + right_val * t
                  end
              end
          end
      end
  end
  
  # Then interpolate along columns (y direction)
  for i in 1:nx
      col_values = perfect_grid_vals[i, :]
      valid_indices = findall(!isnan, col_values)
      
      if length(valid_indices) >= 2  # Need at least 2 points for interpolation
          # For each NaN value, linearly interpolate
          for j in 1:ny
              if isnan(col_values[j])
                  # Find nearest valid points above and below
                  below_indices = valid_indices[valid_indices .< j]
                  above_indices = valid_indices[valid_indices .> j]
                  
                  if !isempty(below_indices) && !isempty(above_indices)
                      below_idx = maximum(below_indices)
                      above_idx = minimum(above_indices)
                      
                      below_y = valid_ys[below_idx]
                      above_y = valid_ys[above_idx]
                      below_val = perfect_grid_vals[i, below_idx]
                      above_val = perfect_grid_vals[i, above_idx]
                      
                      # Linear interpolation
                      y = valid_ys[j]
                      t = (y - below_y) / (above_y - below_y)
                      perfect_grid_vals[i, j] = below_val * (1 - t) + above_val * t
                  end
              end
          end
      end
  end
  
  # Step 5: For any remaining NaN values, use nearest neighbor
  # This handles extrapolation for edge cases
  for i in 1:nx
      for j in 1:ny
          if isnan(perfect_grid_vals[i, j])
              # Find all valid values and their distances to this point
              valid_mask = .!isnan.(perfect_grid_vals)
              
              if any(valid_mask)  # If there are any valid values
                  valid_indices = findall(valid_mask)
                  
                  # Calculate distances to all valid points
                  distances = [
                      sqrt((i - idx[1])^2 + (j - idx[2])^2) 
                      for idx in valid_indices
                  ]
                  
                  # Find nearest valid point
                  nearest_idx = valid_indices[argmin(distances)]
                  perfect_grid_vals[i, j] = perfect_grid_vals[nearest_idx]
              else
                  # If no valid values exist (shouldn't happen), set to 0
                  @warn "No valid values exist for interpolation, setting to 0"
                  perfect_grid_vals[i, j] = 0.0
              end
          end
      end
  end
  
  return perfect_grid_xs, perfect_grid_ys, perfect_grid_vals
end

"""
fit_powerlaws_to_oc(L_dist, V_ps, T_dur, timeValuation; verbose=0)

Fits power law functions to velocity and duration data based on distance and time valuation.

# Arguments
- `L_dist`: Matrix of distances
- `V_ps`: Matrix of peak speeds
- `T_dur`: Matrix of durations
- `timeValuation`: Matrix of time valuation costs
- `verbose`: (Optional) Set to 1 to enable plotting and additional output
- `M`: (Optional) Mass parameter for power law functions

# Returns
- `fn_vk`: Function for predicting velocity given time valuation and distance
- `fn_tk`: Function for predicting duration given time valuation and distance
- `k_v`: Parameter for velocity function
- `k_t`: Parameter for duration function
- `allfit`: Dictionary with statistics on fit quality
"""
function fit_powerlaws_to_oc(ms::main_sequence;M=1, verbose = 0)
    # time valuation changes across dim 1.
    # distance changes across dim 2.
    # have two surfaces, ct vs dist -> spd, ct vs dist -> dur
    # wanna replace two functions of V and T with k^1/4 * ct^1/4 * L^3/4 and k^1/4 * ct^-1/4 * L^1/4
    L_dist          = ms.distance
    T_dur           = ms.duration
    V_ps            = ms.peak_speed
    time_valuation  = ms.time_valuation

    cts_x = time_valuation
    col = findfirst(sum(cts_x .== 0, dims=1) .== 0) # don't use a sim row that has empty spots
    col = col[2]
    ct = cts_x[:, col]

    # Call the fillmissingsurf function from the reach_et module
    (cta, L_dist, T_dur)  = reach_et.simple_grid_fillmissing(time_valuation, L_dist, T_dur)
    (_, _, V_ps)          = reach_et.simple_grid_fillmissing(time_valuation, L_dist, V_ps)

    # Define power law functions
    fn_v = (c_t, L) -> c_t.^(1/4) .* ((1/M)^(1/4)) .* L.^(3/4)
    Av = fn_v(cta, L_dist)
    bv = V_ps
    
    fn_t = (c_t, L) -> (1 ./ c_t).^(1/4) .* ((M)^(1/4)) .* L.^(1/4)
    At = fn_t(cta, L_dist)
    bt = T_dur
    
    #

    # Solve for parameters using least squares
    k_v = Av[:] \ bv[:]
    k_t = At[:] \ bt[:]

    # Create the scaled functions
    fn_vk = (c_t, L) -> k_v * fn_v(c_t, L)
    fn_tk = (c_t, L) -> k_t * fn_t(c_t, L)
    
    # Calculate predicted values
    V_star = fn_vk(cta, L_dist)
    T_star = fn_tk(cta, L_dist)

    # Calculate fit statistics
    ssr_v = sum((V_ps .- V_star).^2)      
    sst_v = sum((V_ps .- mean(V_ps)).^2)
    r2v = 1 - ssr_v/sst_v
    
    ssr_t = sum((T_dur .- T_star).^2)
    sst_t = sum((T_dur .- mean(T_dur)).^2)
    r2t = 1 - ssr_t/sst_t

    RMSE_t = sqrt(mean((T_dur .- T_star).^2))
    rmse_v = sqrt(mean((V_ps .- V_star).^2))
    
    # Compile all fit statistics into a dictionary
    allfit = Dict(
        "f_ctd2t" => fn_tk,
        "f_ctd2v" => fn_vk,
        "k_t" => k_t,
        "k_v" => k_v,
        "R2_t" => r2t,
        "R2_v" => r2v,
        "RMSE_t" => RMSE_t,
        "rmse_v" => rmse_v,
        "rmsep_t" => RMSE_t/mean(T_dur),
        "rmsep_v" => rmse_v/mean(V_ps)
    )

    # Visualization section would go here if verbose=1
    # I'm excluding it for now as it would require equivalent Julia plotting functions
    
    if verbose > 0
        # This would be implemented using Plots.jl or GLMakie.jl
        # Similar to the plot_vt_surfaces functions in the reach_et module
        println("Verbose output enabled, but visualization is not implemented in this conversion.")
        println("Fit statistics:")
        println("R² for velocity: $(r2v)")
        println("R² for duration: $(r2t)")
        println("RMSE for velocity: $(rmse_v)")
        println("RMSE for duration: $(RMSE_t)")
    end

    return fit_(fn_vk, fn_tk, k_v, k_t, allfit)
end
end
