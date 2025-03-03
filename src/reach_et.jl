module reach_et
"""
reach_et.jl: Functions for simulating reaching movements with energy optimization.

"""

using Interpolations
using Statistics
using ColorSchemes
using MAT
using Plots
using GLMakie

export plot_vt_surfaces, plot_vt_surfaces2

function plot_vt_surfaces(; simfile="data/reachgrid.mat", cf=4.2, az=27, el=33.6,plotstyle="GLMakie")
"""
function plot_vt_surfaces(; simfile="data/reachgrid.mat", cf=4.2, az=27, el=33.6,plotstyle="GLMakie")
  
    Plot the surfaces of distance, duration, and peak speed as functions of time valuation and distance.

    Parameters
    ----------
    simfile : str
        Path to the MAT file containing the simulation results.
    f : int
        Figure number for the plot.
    cf : float
        Conversion factor for time valuation.
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
  cs = range(RGB(239/255, 237/255, 245/255), RGB(117/255, 107/255, 177/255), length=100)
  
  # Load MAT file
  data = matread(simfile)
  distance = data["distance"]
  duration = data["duration"]
  peakspeed = data["peakspeed"]
  timeValuation = data["timeValuation"]
  
  (ct, distance, duration)  = simple_grid_fillmissing(timeValuation, distance, duration)
  (_, _, peakspeed)         = simple_grid_fillmissing(timeValuation, distance, peakspeed)
  
  # Add zero columns
  distance  = hcat(zeros(size(distance, 1)), distance)
  duration  = hcat(zeros(size(distance, 1)), duration)
  peakspeed = hcat(zeros(size(distance, 1)), peakspeed)
  ct = hcat(ct[:, 1], ct)
  
  # Convert ct to ct_mech
  eff_mech = 0.25
  ct_mech = ct * eff_mech

  # Plotting
  if plotstyle == "Plots"
      
    # Create a new figure
    fig = Plots.plot(layout=(2, 3), size=(575, 350), background_color=:white)
    
    # First subplot (equivalent to subplot(231) in MATLAB)
    subplot = Plots.plot!(fig[1], ct_mech, distance, duration, 
        st=:surface, 
        color=cs, 
        alpha=0.5, 
        showaxis=true, 
        xlabel="Time valuation Cₜ (W)", 
        ylabel="Distance L (m)", 
        zlabel="Duration T (s)",
        xlim=(0, 50/cf),
        ylim=(0, 0.55),
        zlim=(0, 2.5),
        camera=(vclps[1], vclps[2])
    )
    
    # Second subplot (equivalent to subplot(234) in MATLAB)
    subplot = Plots.plot!(fig[4], ct_mech, distance, peakspeed, 
        st=:surface, 
        color=cs, 
        alpha=0.5, 
        showaxis=true, 
        xlabel="Cₜ (W)", 
        ylabel="Distance L (m)", 
        zlabel="Peak speed V (m/s)",
        xlim=(0, 50/cf),
        ylim=(0, 0.55),
        zlim=(0, 1),
        camera=(vclps[1], vclps[2])
    )
    
    display(fig)
    return fig
  elseif plotstyle == "GLMakie"
    
    
    # Create a new figure
    fig = Figure(resolution=(1000, 700), fontsize=16)
    
    # Add a 3D axis for the duration surface
    ax1 = Axis3(fig[1, 1], 
               xlabel="Time valuation Cₜ (W)",
               ylabel="Distance L (m)",
               zlabel="Duration T (s)",
               title="Duration Surface")
    
    # Create the surface
    surf1 = GLMakie.surface!(ax1, ct_mech, distance, duration,
               colormap=:plasma,
               transparency=true,
               alpha=0.85,
               shading=true)
    
    # Set limits
    GLMakie.xlims!(ax1, 0, 50/cf)
    GLMakie.ylims!(ax1, 0, 0.55)
    GLMakie.zlims!(ax1, 0, 2.5)
    
    # Add colorbar
    Colorbar(fig[1, 2], surf1, label="Duration (s)")
    
    # Add a 3D axis for the peak speed surface
    ax2 = Axis3(fig[2, 1], 
               xlabel="Time valuation Cₜ (W)",
               ylabel="Distance L (m)",
               zlabel="Peak Speed V (m/s)",
               title="Peak Speed Surface")
    
    # Create the surface
    surf2 = GLMakie.surface!(ax2, ct_mech, distance, peakspeed,
               colormap=:plasma,
               transparency=true,
               alpha=0.85,
               shading=true)
    
    # Set limits
    GLMakie.xlims!(ax2, 0, 50/cf)
    GLMakie.ylims!(ax2, 0, 0.55)
    GLMakie.zlims!(ax2, 0, 1.0)
    
    # Add colorbar
    Colorbar(fig[2, 2], surf2, label="Peak Speed (m/s)")
    
    # Construct 2 cameras based on the ax1 and ax2 scenes
    cam1 = Camera3D(ax1.scene)
    cam2 = Camera3D(ax2.scene)

    rotate_cam!(ax1.scene, cam1, (deg2rad(el), deg2rad(az), 0))
    rotate_cam!(ax2.scene, cam2, (deg2rad(el), deg2rad(az), 0))
    
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
    
    return fig
  else
    error("Invalid plotstyle: $plotstyle")
  end
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
end
