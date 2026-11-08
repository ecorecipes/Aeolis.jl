# Model parameters. One table lists every parameter with its type, default, unit and
# description; the `AeolisParameters` struct, its keyword constructor, `parameter_info` and the
# configuration-file reader are all generated from it. Names are those of AeoLiS configuration
# files (`aeolis.txt`), so a file and a keyword call use the same words.

"""
    InputArray

Type of the parameters that hold input data (grids, time series, masks): `nothing` or an array.
In a configuration file they name a file, which [`read_config`](@ref) loads.
"""
const InputArray = Union{Nothing,AbstractArray}

# (name, type, default, unit, description)
const PARAMETER_TABLE = Tuple{Symbol,Type,Any,String,String}[
                                                             # processes
                                                             (:process_wind, Bool, true,
                                                              "-",
                                                              "Enable the process of wind"),
                                                             (:process_transport, Bool,
                                                              true, "-",
                                                              "Enable the process of transport"),
                                                             (:process_bedupdate, Bool,
                                                              true, "-",
                                                              "Enable the process of bed updating"),
                                                             (:process_threshold, Bool,
                                                              true, "-",
                                                              "Enable the process of threshold"),
                                                             (:th_grainsize, Bool, true,
                                                              "-",
                                                              "Enable wind velocity threshold based on grain size"),
                                                             (:th_bedslope, Bool, false,
                                                              "-",
                                                              "Enable wind velocity threshold based on bed slope (no effect, as in AeoLiS)"),
                                                             (:th_moisture, Bool, false,
                                                              "-",
                                                              "Enable wind velocity threshold based on moisture"),
                                                             (:th_drylayer, Bool, false,
                                                              "-",
                                                              "Enable threshold based on drying of layer (not implemented in AeoLiS)"),
                                                             (:th_humidity, Bool, false,
                                                              "-",
                                                              "Enable wind velocity threshold based on humidity (not implemented in AeoLiS)"),
                                                             (:th_salt, Bool, false, "-",
                                                              "Enable wind velocity threshold based on salt"),
                                                             (:th_sheltering, Bool, false,
                                                              "-",
                                                              "Enable wind velocity threshold based on sheltering by roughness elements"),
                                                             (:th_nelayer, Bool, false, "-",
                                                              "Enable wind velocity threshold based on a non-erodible layer"),
                                                             (:process_avalanche, Bool,
                                                              false, "-",
                                                              "Enable the process of avalanching"),
                                                             (:process_shear, Bool, false,
                                                              "-",
                                                              "Enable the process of wind shear"),
                                                             (:process_tide, Bool, false,
                                                              "-",
                                                              "Enable the process of tides"),
                                                             (:process_wave, Bool, false,
                                                              "-",
                                                              "Enable the process of waves"),
                                                             (:process_runup, Bool, false,
                                                              "-",
                                                              "Enable the process of wave runup"),
                                                             (:process_moist, Bool, false,
                                                              "-",
                                                              "Enable the process of moisture"),
                                                             (:process_mixtoplayer, Bool,
                                                              false, "-",
                                                              "Enable mixing of the top layer by waves"),
                                                             (:process_wet_bed_reset, Bool,
                                                              false, "-",
                                                              "Enable resetting the bed in the intertidal zone"),
                                                             (:process_meteo, Bool, false,
                                                              "-",
                                                              "Enable the process of meteo"),
                                                             (:process_salt, Bool, false,
                                                              "-",
                                                              "Enable the process of salt"),
                                                             (:process_humidity, Bool,
                                                              false, "-",
                                                              "Enable the process of humidity"),
                                                             (:process_groundwater, Bool,
                                                              false, "-",
                                                              "Enable the process of groundwater"),
                                                             (:process_scanning, Bool,
                                                              false, "-",
                                                              "Enable scanning curves in soil water retention"),
                                                             (:process_inertia, Bool, false,
                                                              "-",
                                                              "Inertia (unused in AeoLiS)"),
                                                             (:process_separation, Bool,
                                                              false, "-",
                                                              "Enable the separation bubble"),
                                                             (:process_vegetation, Bool,
                                                              false, "-",
                                                              "Enable the process of vegetation"),
                                                             (:process_fences, Bool, false,
                                                              "-",
                                                              "Enable the process of sand fencing"),
                                                             (:process_dune_erosion, Bool,
                                                              false, "-",
                                                              "Enable wave-driven dune erosion"),
                                                             (:process_seepage_face, Bool,
                                                              false, "-",
                                                              "Enable groundwater seepage (only for positive beach slopes)"),
                                                             (:visualization, Bool, false,
                                                              "-",
                                                              "Plot the model setup (not used in Julia)"),
                                                             (:output_sedtrails, Bool,
                                                              false, "-",
                                                              "Generate extra output for SedTRAILS"),
                                                             (:nfraction_sedtrails, Int, 0,
                                                              "-",
                                                              "Index of the fraction for SedTRAILS"),
                                                             # input data
                                                             (:xgrid_file, InputArray,
                                                              nothing, "m",
                                                              "x-coordinates of grid cells"),
                                                             (:ygrid_file, InputArray,
                                                              nothing, "m",
                                                              "y-coordinates of grid cells"),
                                                             (:bed_file, InputArray,
                                                              nothing, "m",
                                                              "Bed level of grid cells"),
                                                             (:wind_file, InputArray,
                                                              nothing, "-",
                                                              "Wind time series: time [s], speed [m/s], direction [deg]"),
                                                             (:tide_file, InputArray,
                                                              nothing, "-",
                                                              "Water level time series: time [s], level [m]"),
                                                             (:wave_file, InputArray,
                                                              nothing, "-",
                                                              "Wave time series: time [s], height [m], period [s]"),
                                                             (:meteo_file, InputArray,
                                                              nothing, "-",
                                                              "Meteorological time series: time [s], T, Q, RH, P, U"),
                                                             (:bedcomp_file, InputArray,
                                                              nothing, "kg/m^2",
                                                              "Initial bed composition"),
                                                             (:threshold_file, InputArray,
                                                              nothing, "m/s",
                                                              "Shear velocity threshold"),
                                                             (:fence_file, InputArray,
                                                              nothing, "m",
                                                              "Sand fence height above the bed"),
                                                             (:ne_file, InputArray, nothing,
                                                              "m", "Non-erodible layer"),
                                                             (:veg_file, InputArray,
                                                              nothing, "-",
                                                              "Initial vegetation density"),
                                                             (:supply_file, InputArray,
                                                              nothing, "kg/m^2/s",
                                                              "Prescribed sediment supply"),
                                                             (:wave_mask, InputArray,
                                                              nothing, "-",
                                                              "Mask for wave height"),
                                                             (:tide_mask, InputArray,
                                                              nothing, "-",
                                                              "Mask for tidal elevation"),
                                                             (:runup_mask, InputArray,
                                                              nothing, "-",
                                                              "Mask for run-up"),
                                                             (:threshold_mask, InputArray,
                                                              nothing, "-",
                                                              "Mask for the shear velocity threshold"),
                                                             (:gw_mask, InputArray, nothing,
                                                              "-",
                                                              "Mask for the groundwater level"),
                                                             (:vver_mask, InputArray,
                                                              nothing, "-",
                                                              "Mask for vertical vegetation growth"),
                                                             # grid and time
                                                             (:nx, Int, 0, "-",
                                                              "Number of grid cells in x (set from the grid)"),
                                                             (:ny, Int, 0, "-",
                                                              "Number of grid cells in y (set from the grid)"),
                                                             (:dt, Float64, 60.0, "s",
                                                              "Time step"),
                                                             (:dx, Float64, 1.0, "m",
                                                              "Cell size of the computational grids for shear in x"),
                                                             (:dy, Float64, 1.0, "m",
                                                              "Cell size of the computational grids for shear in y"),
                                                             (:CFL, Float64, 1.0, "-",
                                                              "CFL number for the explicit scheme"),
                                                             (:accfac, Float64, 1.0, "-",
                                                              "Numerical acceleration factor"),
                                                             (:max_bedlevel_change, Float64,
                                                              999.0, "m",
                                                              "Maximum bed level change in one step (999 to disable)"),
                                                             (:tstart, Float64, 0.0, "s",
                                                              "Start time of the simulation"),
                                                             (:tstop, Float64, 3600.0, "s",
                                                              "End time of the simulation"),
                                                             (:restart,
                                                              Union{Nothing,Float64},
                                                              nothing, "s",
                                                              "Interval for writing restart files"),
                                                             (:dzb_interval, Float64,
                                                              86400.0, "s",
                                                              "Interval for vegetation growth"),
                                                             (:output_times, Float64, 60.0,
                                                              "s", "Output interval"),
                                                             (:output_file,
                                                              Union{Nothing,String},
                                                              nothing, "-",
                                                              "Name of the netCDF output file"),
                                                             (:output_vars, Vector{Symbol},
                                                              [:zb, :zs, :Ct, :Cu, :uw,
                                                               :udir, :uth, :mass, :pickup,
                                                               :w], "-",
                                                              "Variables written to the output"),
                                                             (:output_types, Vector{Symbol},
                                                              Symbol[], "-",
                                                              "Statistics written to the output (avg, sum, var, min, max)"),
                                                             (:external_vars,
                                                              Vector{Symbol}, Symbol[], "-",
                                                              "Variables set by a coupled model"),
                                                             # sediment
                                                             (:grain_size, Vector{Float64},
                                                              [225e-6], "m",
                                                              "Grain size of each fraction"),
                                                             (:grain_dist,
                                                              VecOrMat{Float64}, [1.0], "-",
                                                              "Initial distribution over fractions (one row per layer, or one row)"),
                                                             (:nlayers, Int, 3, "-",
                                                              "Number of bed layers"),
                                                             (:layer_thickness, Float64,
                                                              0.01, "m",
                                                              "Thickness of bed layers"),
                                                             # physical constants
                                                             (:g, Float64, 9.81, "m/s^2",
                                                              "Gravitational acceleration"),
                                                             (:v, Float64, 0.000015,
                                                              "m^2/s", "Air viscosity"),
                                                             (:rhoa, Float64, 1.225,
                                                              "kg/m^3", "Air density"),
                                                             (:rhog, Float64, 2650.0,
                                                              "kg/m^3", "Grain density"),
                                                             (:rhow, Float64, 1025.0,
                                                              "kg/m^3", "Water density"),
                                                             (:porosity, Float64, 0.4, "-",
                                                              "Sediment porosity"),
                                                             (:Aa, Float64, 0.085, "-",
                                                              "Constant in the grain size threshold"),
                                                             (:z, Float64, 10.0, "m",
                                                              "Measurement height of wind velocity"),
                                                             (:h, Union{Nothing,Float64},
                                                              nothing, "m",
                                                              "Height of the saltation layer (defaults to z)"),
                                                             (:k, Float64, 0.001, "m",
                                                              "Bed roughness"),
                                                             (:L, Float64, 100.0, "m",
                                                              "Typical length scale of dunes (shear perturbation)"),
                                                             (:l, Float64, 10.0, "m",
                                                              "Inner layer height (shear perturbation)"),
                                                             (:c_b, Float64, 0.2, "-",
                                                              "Slope at the lee side of the separation bubble"),
                                                             (:mu_b, Float64, 30.0, "deg",
                                                              "Minimum slope for flow separation"),
                                                             (:buffer_width, Float64, 10.0,
                                                              "m",
                                                              "Buffer around the rotated grid for shear"),
                                                             (:sep_filter_iterations, Int,
                                                              0, "-",
                                                              "Filter iterations on the separation bubble"),
                                                             (:zsep_y_filter, Bool, false,
                                                              "-",
                                                              "Filter the separation bubble in y"),
                                                             (:Cb, Float64, 1.5, "-",
                                                              "Constant in the Bagnold formulation"),
                                                             (:Ck, Float64, 2.78, "-",
                                                              "Constant in the Kawamura formulation"),
                                                             (:Cl, Float64, 6.7, "-",
                                                              "Constant in the Lettau formulation"),
                                                             (:Cdk, Float64, 5.0, "-",
                                                              "Constant in the DK formulation"),
                                                             (:m, Float64, 0.5, "-",
                                                              "Ratio of average to maximum shear stress (sheltering; absent from AeoLiS defaults)"),
                                                             (:kappa, Float64, 0.41, "-",
                                                              "Von Karman constant"),
                                                             (:sigma, Float64, 4.2, "-",
                                                              "Ratio of basal to frontal area of roughness elements"),
                                                             (:beta, Float64, 130.0, "-",
                                                              "Ratio of drag coefficients of roughness elements and bare surface"),
                                                             (:bi, Float64, 1.0, "-",
                                                              "Bed interaction factor"),
                                                             (:T, Float64, 1.0, "s",
                                                              "Adaptation time scale in the advection equation"),
                                                             (:Tdry, Float64, 3600.0 * 1.5,
                                                              "s",
                                                              "Adaptation time scale for soil drying"),
                                                             (:Tsalt, Float64,
                                                              3600.0 * 24.0 * 30.0, "s",
                                                              "Adaptation time scale for salination"),
                                                             (:Tbedreset, Float64, 86400.0,
                                                              "s",
                                                              "Time scale of the intertidal bed reset"),
                                                             (:eps, Float64, 1e-3, "m",
                                                              "Minimum water depth for a flooded cell"),
                                                             (:gamma, Float64, 0.5, "-",
                                                              "Maximum ratio of wave height to depth"),
                                                             (:xi, Float64, 0.3, "-",
                                                              "Surf similarity parameter"),
                                                             (:facDOD, Float64, 0.1, "-",
                                                              "Ratio of depth of disturbance to wave height"),
                                                             (:csalt, Float64, 35e-3, "-",
                                                              "Maximum salt concentration in the surface layer"),
                                                             (:cpair, Float64, 1.0035e-3,
                                                              "MJ/kg/oC",
                                                              "Specific heat capacity of air"),
                                                             # moisture and groundwater
                                                             (:fc, Float64, 0.11, "-",
                                                              "Moisture content at field capacity (volumetric)"),
                                                             (:w1_5, Float64, 0.02, "-",
                                                              "Moisture content at wilting point (gravimetric)"),
                                                             (:resw_moist, Float64, 0.01,
                                                              "-",
                                                              "Residual moisture content, wetting"),
                                                             (:satw_moist, Float64, 0.35,
                                                              "-",
                                                              "Saturated moisture content, wetting"),
                                                             (:resd_moist, Float64, 0.01,
                                                              "-",
                                                              "Residual moisture content, drying"),
                                                             (:satd_moist, Float64, 0.5,
                                                              "-",
                                                              "Saturated moisture content, drying"),
                                                             (:nw_moist, Float64, 2.3, "-",
                                                              "Pore-size distribution index, wetting"),
                                                             (:nd_moist, Float64, 4.5, "-",
                                                              "Pore-size distribution index, drying"),
                                                             (:mw_moist, Float64, 0.57, "-",
                                                              "van Genuchten m, wetting"),
                                                             (:md_moist, Float64, 0.42, "-",
                                                              "van Genuchten m, drying"),
                                                             (:alfaw_moist, Float64, -0.070,
                                                              "cm^-1",
                                                              "Inverse air-entry value, wetting"),
                                                             (:alfad_moist, Float64, -0.035,
                                                              "cm^-1",
                                                              "Inverse air-entry value, drying"),
                                                             (:thick_moist, Float64, 0.002,
                                                              "m",
                                                              "Thickness of the surface moisture layer"),
                                                             (:K_gw, Float64, 0.00078,
                                                              "m/s",
                                                              "Hydraulic conductivity"),
                                                             (:ne_gw, Float64, 0.3, "-",
                                                              "Effective porosity"),
                                                             (:D_gw, Float64, 12.0, "m",
                                                              "Aquifer depth"),
                                                             (:tfac_gw, Float64, 10.0, "-",
                                                              "Reduction factor of the groundwater time step"),
                                                             (:Cl_gw, Float64, 0.7, "m",
                                                              "Groundwater overheight due to runup"),
                                                             (:in_gw, Float64, 0.0, "m",
                                                              "Initial groundwater level"),
                                                             (:GW_stat, Float64, 1.0, "m",
                                                              "Landward static groundwater boundary"),
                                                             (:max_moist, Float64, 10.0,
                                                              "%",
                                                              "Volumetric moisture above which the threshold is infinite"),
                                                             # avalanching and vegetation
                                                             (:theta_dyn, Float64, 33.0,
                                                              "deg",
                                                              "Dynamic angle of repose"),
                                                             (:theta_stat, Float64, 34.0,
                                                              "deg",
                                                              "Static angle of repose"),
                                                             (:avg_time, Float64, 86400.0,
                                                              "s",
                                                              "Averaging period of bed level change"),
                                                             (:gamma_vegshear, Float64,
                                                              16.0, "-",
                                                              "Roughness factor of vegetation shear reduction"),
                                                             (:hveg_max, Float64, 1.0, "m",
                                                              "Maximum vegetation height"),
                                                             (:dzb_opt, Float64, 0.0,
                                                              "m/year",
                                                              "Burial for optimal growth"),
                                                             (:V_ver, Float64, 0.0,
                                                              "m/year",
                                                              "Vertical growth potential"),
                                                             (:germinate, Float64, 0.0,
                                                              "1/year",
                                                              "Probability of germination per year"),
                                                             (:lateral, Float64, 0.0,
                                                              "1/year",
                                                              "Probability of lateral expansion per year"),
                                                             (:veg_gamma, Float64, 1.0, "-",
                                                              "Influence of sediment burial"),
                                                             (:veg_sigma, Float64, 0.0, "-",
                                                              "Width of the Gaussian filter on vegetation cover"),
                                                             (:sedimentinput, Float64, 0.0,
                                                              "-",
                                                              "Constant boundary influx (unused)"),
                                                             # numerics and boundaries
                                                             (:scheme, Symbol,
                                                              :euler_backward, "-",
                                                              "Numerical scheme: euler_forward, euler_backward or crank_nicolson"),
                                                             (:solver, Symbol, :trunk, "-",
                                                              "Transport solver: trunk, pieter, steadystate or steadystatepieter"),
                                                             (:boundary_lateral, Symbol,
                                                              :circular, "-",
                                                              "Lateral boundary: circular, constant, gradient or flux"),
                                                             (:boundary_offshore, Symbol,
                                                              :constant, "-",
                                                              "Offshore boundary: flux, constant, uniform or gradient"),
                                                             (:boundary_onshore, Symbol,
                                                              :gradient, "-",
                                                              "Onshore boundary: flux, constant, uniform or gradient"),
                                                             (:boundary_gw, Symbol,
                                                              :no_flow, "-",
                                                              "Landward groundwater boundary: no_flow or static"),
                                                             (:method_moist_threshold,
                                                              Symbol, :belly_johnson, "-",
                                                              "Moisture threshold formulation"),
                                                             (:method_moist_process, Symbol,
                                                              :infiltration, "-",
                                                              "Moisture process: infiltration or surf_moisture"),
                                                             (:offshore_flux, Float64, 0.0,
                                                              "-",
                                                              "Offshore boundary flux as a fraction of saturation"),
                                                             (:constant_offshore_flux,
                                                              Float64, 0.0, "kg/m/s",
                                                              "Constant offshore input flux"),
                                                             (:onshore_flux, Float64, 0.0,
                                                              "-",
                                                              "Onshore boundary flux as a fraction of saturation"),
                                                             (:constant_onshore_flux,
                                                              Float64, 0.0, "kg/m/s",
                                                              "Constant onshore input flux"),
                                                             (:lateral_flux, Float64, 0.0,
                                                              "-",
                                                              "Lateral boundary flux as a fraction of saturation"),
                                                             (:method_transport, Symbol,
                                                              :bagnold, "-",
                                                              "Equilibrium transport formulation"),
                                                             (:method_roughness, Symbol,
                                                              :constant, "-",
                                                              "Roughness length formulation"),
                                                             (:method_grainspeed, Symbol,
                                                              :windspeed, "-",
                                                              "Grain speed: windspeed, duran, duran_full or constant"),
                                                             (:method_shear, Symbol, :fft,
                                                              "-",
                                                              "Shear perturbation: fft, quasi2d or duna2d"),
                                                             (:max_error, Float64, 1e-6,
                                                              "-",
                                                              "Tolerance of the implicit iterations"),
                                                             (:max_iter, Int, 1000, "-",
                                                              "Maximum iterations of the implicit schemes"),
                                                             (:max_iter_ava, Int, 1000, "-",
                                                              "Maximum iterations of avalanching"),
                                                             (:refdate, String,
                                                              "2020-01-01 00:00", "-",
                                                              "Reference date of the output"),
                                                             (:callback,
                                                              Union{Nothing,String},
                                                              nothing, "-",
                                                              "Callback reference (Python only; use SciMLBase callbacks in Julia)"),
                                                             (:wind_convention, Symbol,
                                                              :nautical, "-",
                                                              "Wind direction convention: nautical or cartesian"),
                                                             (:alfa, Float64, 0.0, "deg",
                                                              "Real-world grid orientation with respect to North"),
                                                             (:dune_toe_elevation, Float64,
                                                              3.0, "m",
                                                              "Dune toe elevation (PH12 dune erosion)"),
                                                             (:beach_slope, Float64, 0.1,
                                                              "-",
                                                              "Beach slope (runup and PH12 dune erosion)"),
                                                             (:Cs, Float64, 2.5e-3, "-",
                                                              "Dune erosion coefficient (PH12)"),
                                                             (:veg_min_elevation, Float64,
                                                              -10.0, "m",
                                                              "Minimum elevation for vegetation"),
                                                             (:vegshear_type, Symbol,
                                                              :raupach, "-",
                                                              "Vegetation shear: raupach or okin"),
                                                             (:okin_c1_veg, Float64, 0.48,
                                                              "-",
                                                              "Okin spatial reduction factor for vegetation"),
                                                             (:okin_c1_fence, Float64, 0.48,
                                                              "-",
                                                              "Okin spatial reduction factor for fences"),
                                                             (:okin_initialred_veg, Float64,
                                                              0.32, "-",
                                                              "Okin initial reduction for vegetation"),
                                                             (:okin_initialred_fence,
                                                              Float64, 0.32, "-",
                                                              "Okin initial reduction for fences"),
                                                             (:veggrowth_type, Symbol,
                                                              :orig, "-",
                                                              "Vegetation growth: orig or duranmoore14"),
                                                             (:rhoveg_max, Float64, 0.5,
                                                              "-",
                                                              "Maximum vegetation density (duranmoore14)"),
                                                             (:t_veg, Float64, 3.0, "day",
                                                              "Time scale of vegetation growth (duranmoore14)"),
                                                             (:v_gam, Float64, 1.0, "-",
                                                              "Burial response (duranmoore14)"),
                                                             # marram grass (Aeolis.jl only)
                                                             (:veg_model, Symbol, :legacy,
                                                              "-",
                                                              "Vegetation model: legacy (AeoLiS) or marram"),
                                                             (:veg_traits_file,
                                                              Union{Nothing,String},
                                                              nothing, "-",
                                                              "TOML file of marram grass traits (defaults when None)"),
                                                             (:veg_water_file, InputArray,
                                                              nothing, "-",
                                                              "Root-zone water forcing: time [s], rainfall [mm/day], potential evapotranspiration [mm/day]"),
                                                             (:veg_interval, Float64,
                                                              86400.0, "s",
                                                              "Interval of the marram grass update in accelerated time")]

const PARAMETER_NAMES = Tuple(first.(PARAMETER_TABLE))
const PARAMETER_INDEX = Dict(row[1] => row for row in PARAMETER_TABLE)

let fields = [:($(row[1])::$(row[2])) for row in PARAMETER_TABLE]
    @eval begin
        """
            AeolisParameters(; kwargs...)

        The parameters of an AeoLiS model: process switches, input data, grid and time settings,
        physical constants and numerical options. Every keyword has the name it has in an
        AeoLiS configuration file; [`parameter_info`](@ref) lists them with defaults, units and
        descriptions, and [`read_config`](@ref) builds parameters from a file.

        Values are converted to the parameter's type: numbers to `Float64` or `Int`, strings to
        `Symbol` for method names, scalars to one-element vectors for list parameters. An
        unknown keyword raises an `ArgumentError`.

        # Examples

        ```julia
        p = AeolisParameters(; xgrid_file=collect(0.0:1.0:100.0), bed_file=zeros(101),
                             wind_file=[0.0 10.0 270.0; 3600.0 10.0 270.0], tstop=3600.0)
        p.grain_size   # [0.000225]
        ```
        """
        struct AeolisParameters
            $(fields...)
        end
    end
end

_convert_param(::Type{T}, v) where {T} = convert(T, v)
_convert_param(::Type{Symbol}, v::AbstractString) = Symbol(lowercase(strip(v)))
_convert_param(::Type{Symbol}, v::Symbol) = Symbol(lowercase(String(v)))
function _convert_param(::Type{Vector{Symbol}}, v)
    v === nothing && return Symbol[]
    v isa Union{Symbol,AbstractString} && return [Symbol(v)]
    return Symbol[Symbol(x) for x in v]
end
_convert_param(::Type{Vector{Float64}}, v::Real) = [Float64(v)]
_convert_param(::Type{Vector{Float64}}, v::AbstractArray) = Float64.(vec(v))
_convert_param(::Type{VecOrMat{Float64}}, v::Real) = [Float64(v)]
_convert_param(::Type{VecOrMat{Float64}}, v::AbstractVector) = Float64.(v)
_convert_param(::Type{VecOrMat{Float64}}, v::AbstractMatrix) = Float64.(v)
_convert_param(::Type{Union{Nothing,Float64}}, v::Real) = Float64(v)
_convert_param(::Type{Union{Nothing,String}}, v::Symbol) = String(v)
_convert_param(::Type{String}, v::Symbol) = String(v)
_convert_param(::Type{Bool}, v::AbstractString) = uppercase(strip(v)) in ("T", "TRUE", "1")

function _default_value(row)
    v = row[3]
    return v isa AbstractArray ? copy(v) : v
end

function _closest_names(name::Symbol, n::Int=3)
    s = lowercase(String(name))
    score(k) = _edit_distance(s, lowercase(String(k)))
    return sort(collect(PARAMETER_NAMES); by=score)[1:min(n, length(PARAMETER_NAMES))]
end

function _edit_distance(a::AbstractString, b::AbstractString)
    a, b = collect(a), collect(b)
    d = collect(0:length(b))
    for (i, ca) in enumerate(a)
        prev, d[1] = d[1], i
        for (j, cb) in enumerate(b)
            cur = min(d[j + 1] + 1, d[j] + 1, prev + (ca == cb ? 0 : 1))
            prev, d[j + 1] = d[j + 1], cur
        end
    end
    return d[end]
end

AeolisParameters(; kwargs...) = _make_parameters(Dict{Symbol,Any}(kwargs))

function _make_parameters(kw::Dict{Symbol,Any})
    for k in keys(kw)
        haskey(PARAMETER_INDEX, k) ||
            throw(ArgumentError("unknown AeoLiS parameter :$k; closest names are " *
                                join((":" * String(c) for c in _closest_names(k)), ", ")))
    end
    vals = Vector{Any}(undef, length(PARAMETER_TABLE))
    for (i, row) in enumerate(PARAMETER_TABLE)
        name, T = row[1], row[2]
        v = haskey(kw, name) ? kw[name] : _default_value(row)
        vals[i] = try
            _convert_param(T, v)
        catch err
            err isa Union{MethodError,InexactError,ArgumentError} || rethrow()
            throw(ArgumentError("parameter :$name expects a value of type $T, got $(repr(v))"))
        end
    end
    return AeolisParameters(vals...)
end

"""
    remake_parameters(p::AeolisParameters; kwargs...)

A copy of `p` with the given parameters replaced.
"""
function remake_parameters(p::AeolisParameters; kwargs...)
    current = Dict{Symbol,Any}(n => getfield(p, n) for n in PARAMETER_NAMES)
    merge!(current, Dict{Symbol,Any}(kwargs))
    return _make_parameters(current)
end

"""
    parameter_info([name])

A table of every parameter as `(name, type, default, unit, description)` named tuples, or the
entry of one parameter.
"""
function parameter_info()
    return [(; name=r[1], type=r[2], default=r[3], unit=r[4], description=r[5])
            for r in PARAMETER_TABLE]
end
function parameter_info(name::Symbol)
    haskey(PARAMETER_INDEX, name) ||
        throw(ArgumentError("unknown AeoLiS parameter :$name; closest names are " *
                            join((":" * String(c) for c in _closest_names(name)), ", ")))
    r = PARAMETER_INDEX[name]
    return (; name=r[1], type=r[2], default=r[3], unit=r[4], description=r[5])
end

function Base.show(io::IO, ::MIME"text/plain", p::AeolisParameters)
    on = [String(n)
          for n in PARAMETER_NAMES
          if (startswith(String(n), "process_") || startswith(String(n), "th_")) &&
        getfield(p, n) === true]
    println(io, "AeolisParameters")
    println(io, "  time: ", p.tstart, " to ", p.tstop, " s, dt = ", p.dt, " s")
    println(io, "  scheme: ", p.scheme, ", solver: ", p.solver)
    println(io, "  fractions: ", length(p.grain_size), ", layers: ", p.nlayers)
    return print(io, "  enabled: ", join(on, ", "))
end
Base.show(io::IO, p::AeolisParameters) = print(io, "AeolisParameters(...)")

function Base.:(==)(a::AeolisParameters, b::AeolisParameters)
    return all(n -> isequal(getfield(a, n), getfield(b, n)), PARAMETER_NAMES)
end
