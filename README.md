# Shocktube Setup

This script gives a framework to set up arbitrary shocktubes for standard hydro tests in Gadget.

## Install

The files `Project.toml` and `Manifest.toml` contain all dependencies with the correct version numbers. The first line of the script installs dependencies and activates the virtual environment used for the setup:

```julia
using Pkg; Pkg.instantiate(); Pkg.activate()
```

## General setup

You can set up a shocktube by calling:

```julia
setup_shocktube(param)
```

where `param` is of type `ShockParameters`.
You can define a `ShockParameters` type with the custom constructor:

```julia
ShockParameters(glass_file::String="",                  # path to a suitable glass file
                output_file::String="",                 # path for the output file
                U::Vector{Float64}  = zeros(2),         # Internal energy at the left and right end of the tube
                B::Array{Float64,2} = zeros(2,3),       # magnetic field vectors at left and right of the tube
                v::Array{Float64,2} = zeros(2,3);       # initial velocities at left and right of the tube
                turb::Bool=false,                       # flag if turbulent magnetic field is desired
                B0::Float64=0.0,                        # normalisation for a turbulent magnetic field
                n_blocks::Int64=70)                     # number of stacks of the glass files in each direction
```

You can set `U` by hand, or get them by fixing the left value (where the high density region is) and solving for the right value in such a way that the resulting shock will have a target Mach number.

For that we use the package [AnalyticMHDTestSolutions.jl](https://github.com/LudwigBoess/AnalyticMHDTestSolutions.jl).
To get the left and right internal energies for an arbitrary Mach number use:

```julia
Mach_number = 100.0
riemann_par = RiemannParameters(Ul=100.0, Mach=Mach_number, t=1.5)

riemann_par.Ul  # internal energy at the left of the tube
riemann_par.Ur  # internal energy at the right
```

## Sod Shocks

Sod shocks, after [Sod (1978)](https://doi.org/10.1016%2F0021-9991%2878%2990023-2) are defined by a density gradient of 8/1  and an arbitrary temperature / internal energy gradient over the contact discontinuity.
You can set them up by using:

```julia
glass_file = "/path/to/glass.file"

out_file = "path/to/out.file"

Mach_number = 100.0 # change accordingly

riemann_par = RiemannParameters(Ul=100.0, Mach=Mach_number, t=1.5)


param = ShockParameters(glass_file,
                        out_file,
                        [riemann_par.Ul, riemann_par.Ur])      # U

setup_shocktube(param)
```

## Ryu&Jones MHD shocks

[Ryu&Jones (1995)](https://ui.adsabs.harvard.edu/link_gateway/1995ApJ...442..228R/doi:10.1086/175437) proposed a series of MHD shocks with analytic solutions.

As an example, you can set up their shock from figure 5a with:

```julia
glass_file = "/path/to/glass.file"

out_file = "path/to/out.file"

param = ShockParameters(glass_file,
                        out_file,
                        [1.0, 0.1],                    # U
                        [0.0 1.0 0.0; 0.0 -1.0 0.0])   # B

setup_shocktube(param)
```

Be cautions with the 2a, that requires no density gradient and is still to be implemented!


## Shocks with an angle between magnetic field vector and shock normal

For shocks with cosmic ray acceleration the angle between magnetic field vector and shock normal is relevant.
To set up such shocks for a given angle `theta` you can use the helper function:

```julia
B = get_bfield_from_angle(theta, reduction_scale=1.e-10)
```

which returns the magnetic field vector `B`. The arbitrary argument `reduction_scale` rescales the magnetic field vector to avoid a kinetic impact of the magnetic field on the shock and defaults to `reduction_scale=1.e-10`.

You can the set up a shock with:

```julia
glass_file = "/path/to/glass.file"

out_file = "path/to/out.file"

Mach_number = 10.0 # change accordingly

riemann_par = RiemannParameters(Ul=100.0, Mach=Mach_number, t=1.5)

theta = 45.0  # degrees

B = get_bfield_from_angle(theta)

param = ShockParameters(glass_file,
                        out_file,
                        [riemann_par.Ul, riemann_par.Ur],
                        [B';  B'] )

setup_shocktube(param)
```

## Shocks with a turbulent magnetic field

To set up a shock with a random, but divergence-free magnetic field you can use:

```julia
glass_file = "/path/to/glass.file"

out_file = "path/to/out.file"

Mach_number = 10.0 # change accordingly

riemann_par = RiemannParameters(Ul=100.0, Mach=Mach_number, t=1.5)


param = ShockParameters(glass_file,
                        out_file,
                        [riemann_par.Ul, riemann_par.Ur],    # U
                        B0   = 1.e-20,
                        turb = true )

setup_shocktube(param)
```