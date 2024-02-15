# Dust modelled as a low-pressure fluid

## Solver and general options

Dust can also be modelled in the code as an additional low-pressure fluid which is coupled to the gas via gas drag. This requires setting `DustFluid` to `Yes` in the input parameter file. Contrary to the Lagrangian approach, a low-pressure fluid can only simulate dust of a given size (parameter `Sizepart`). Unless the short-friction time approximation is used (by setting `SFTApprox` to `Yes`), the same momentum equation as for the gas is solved for the dust fluid on the polar grid, with the addition of the gas drag term. In particular, if `DustFeelSG` is set to `Yes` in the parameter file, the self-gravitating acceleration of the gas is included in the dust’s momentum equation (although, just like for the Lagrangian approach, the dust’s self-gravity is currently ignored). The viscous term can be included for the dust fluid upon specification of the dust’s alpha turbulent viscosity (parameter `DAlphaViscosity`). Alternatively, and probably preferentially, dust turbulence can be modelled as a diffusion term in the dust’s continuity equation (see Zhu et al. 2012, ApJ, 755 ; [ADS](https://ui.adsabs.harvard.edu/abs/2012ApJ...755....6Z/abstract) – their equation 8), where, just like for Lagrangian particles, the dust’s turbulent diffusivity, $D$, is taken to be $D = \nu(1 + 4St^2 )/(1 + St^2 )^2$, with $\nu$ the local kinematic viscosity of the gas, and $St$ the Stokes number of the dust. This requires the keyword `DustDiffusion` set to `Yes` in the input parameter file.

## Internal density and code units

Just like when dust is modelled as Lagrangian particles, the dust’s internal density needs to be specified in $\rm{g.cm^{−3}}$, which implies that the code’s units of length and mass also need to be specified in the input parameter file (see the [](dust_particles.md#internal-density-and-code-units) subsection).

## Initial conditions

Like for the gas, the initial conditions for the dust fluid are related to its surface density (which is set by the parameter `DustToGasDensityRatio`) and its aspect ratio (parameters `DAspectRatio` and `DFlaringIndex`). The latter points to the main caveat of modelling dust as a low-pressure fluid, which is
to determine how low the dust’s pressure (or sound speed) can be, which is problem- and resolution-dependent (see tests in Zhu et al. 2012, ApJ, 755 ; [ADS](https://ui.adsabs.harvard.edu/abs/2012ApJ...755....6Z/abstract)). As a first guess inspired by many test cases, setting the dust’s aspect ratio to $1/10$th that of the gas gives good agreement with the Lagrangian approach.

## Boundary conditions

Wave-killing zones with the same radial extent as for the gas can be used for the dust fluid, where the density and velocity components are damped towards their instantaneous, axisymmetric profiles. Alternatively, an open boundary condition for the dust fluid can be used by setting `InnerBoundaryDust` to `O` in the parameter file.

## Outputs

The code will write `dustX.dat` binary files with 2D arrays of the dust surface density and velocity field, just like for the gas.
