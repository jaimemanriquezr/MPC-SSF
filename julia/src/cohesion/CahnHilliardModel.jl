# Port of src/cohesion/@CahnHilliardModel/CahnHilliardModel.m
#
# Cohesion submodel for the biofilm. Supplies:
#   kappa               - interfacial/diffusion coefficient (used by the
#                         Cahn-Hilliard diffusion matrix)
#   zeta_0, zeta_1      - cohesion coefficients
#   mobility(u)         - mobility function
#   potential_gradient(u) - chemical-potential gradient dψ/du
#
# In MATLAB `mobility`/`potential_gradient` are function_handle fields; here they
# are callables. `kappa`/`zeta_0` are required (presets always set them);
# `zeta_1` defaults to 0. Functions are written with broadcasting so they accept
# scalars or vectors.

"""
    CahnHilliardModel(; kappa, zeta_0, zeta_1=0.0,
                        mobility=(u -> u .* (1 .- u)),
                        potential_gradient=(u -> 0.25 .* (u .^ 2 .* (1 .- u) .^ 2)))

Biofilm cohesion submodel. Defaults for `mobility` and `potential_gradient`
mirror `CahnHilliardModel.m`.
"""
Base.@kwdef struct CahnHilliardModel
    kappa::Float64
    zeta_0::Float64
    zeta_1::Float64 = 0.0
    mobility::Function = u -> u .* (1 .- u)
    potential_gradient::Function = u -> 0.25 .* (u .^ 2 .* (1 .- u) .^ 2)
end
