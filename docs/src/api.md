# Problems and solutions

```@docs
Aeolis
AeolisProblem
init(::AeolisProblem, ::AbstractAeolisAlgorithm)
solve(::AeolisProblem, ::AbstractAeolisAlgorithm)
step!(::AeolisIntegrator, ::Integer)
solve!(::AeolisIntegrator)
remake(::AeolisProblem)
reinit!(::AeolisIntegrator)
AeolisIntegrator
AeolisSolution
transect
statistic
save_restart
AbstractAeolisAlgorithm
EulerBackward
CrankNicolson
EulerForward
default_algorithm
AbstractTransportSolver
Trunk
Pieter
SteadyState
SteadyStatePieter
SciMLSteadyState
read_config
write_config
run_config
write_netcdf
load_array
```

## Exceptions

```@docs
AeolisError
ConfigurationError
UnknownOptionError
UnsupportedOptionError
```
