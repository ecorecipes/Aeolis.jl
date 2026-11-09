# Processes

These functions advance one process of a time step; [`AeolisProblem`](@ref)'s integrator calls
them in AeoLiS's order. They are not exported.

```@docs
Aeolis.setup_model
Aeolis.aeolis_step!
Aeolis.set_timestep!
Aeolis.initialize_grid!
Aeolis.grid_rotate!
Aeolis.initialize_wind!
Aeolis.interpolate_wind!
Aeolis.roughness_length
Aeolis.velocity_stress!
Aeolis.stress_velocity!
Aeolis.compute_shear!
Aeolis.WindShear
Aeolis.wind_shear!
Aeolis.RotationGrid
Aeolis.rotated_shear!
Aeolis.okin_shear!
Aeolis.compute_threshold!
Aeolis.threshold_grainsize!
Aeolis.threshold_moisture!
Aeolis.threshold_salt!
Aeolis.threshold_sheltering!
Aeolis.threshold_nelayer!
Aeolis.grainspeed_duran
Aeolis.grainspeed_constant!
Aeolis.equilibrium!
Aeolis.compute_weights
Aeolis.renormalize_weights!
Aeolis.solve_trunk
Aeolis.solve_euler_forward
Aeolis.solve_pieter
Aeolis.solve_steady_state
Aeolis.sweep!
Aeolis.initialize_bed!
Aeolis.mixtoplayer!
Aeolis.wet_bed_reset!
Aeolis.update_bed!
Aeolis.average_change!
Aeolis.interpolate_hydro!
Aeolis.update_hydro!
Aeolis.stockdon_runup
Aeolis.angle_of_repose!
Aeolis.avalanche!
Aeolis.initialize_vegetation!
Aeolis.vegetation_shear!
Aeolis.germinate!
Aeolis.grow_vegetation!
Aeolis.initialize_fences!
Aeolis.update_fences!
Aeolis.dune_erosion!
```

## Numerical helpers

```@docs
Aeolis.InputArray
Aeolis.AeolisModel
Aeolis.state_dims
Aeolis.copyto_state!
Aeolis.parse_value
Aeolis.rotate
Aeolis.rotate!
Aeolis.normalize_along
Aeolis.prevent_tiny_negatives!
Aeolis.apply_mask!
Aeolis.CircularSeries
Aeolis.CircularStepSeries
Aeolis.clamped_interp
Aeolis.diags
Aeolis.np_arange
Aeolis.exact_grid
Aeolis.fftfreq
Aeolis.pad_odd_reflect
Aeolis.pad_edge
Aeolis.regular_grid_interp
Aeolis.linear_fit
Aeolis.max_diff
Aeolis.gaussian_filter1d
Aeolis.gaussian_filter
Aeolis.grain_size_percentile
Aeolis.mean_grain_size
Aeolis.pyrange
```
