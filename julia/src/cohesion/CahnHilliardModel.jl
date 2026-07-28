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
#
# `potential_gradient` defaults to the published form. The manuscript defines
# psi(u) = u^3 (u/2 - zeta_1) / 2, hence dpsi/du = u^2 (u - 3 zeta_1 / 2), and the
# legacy models store exactly that as a baked handle (`SDparameters.mat` and
# `SDparameters_pathogen.mat`: `@(phi) phi.^2.*(phi-3*(1/200)/2)` with
# zeta_1 = 1/200; `thesis_model.mat`: the same with zeta_1 = 1/100). Because the
# default closes over `zeta_1`, a caller changing `zeta_1` must reconstruct
# WITHOUT passing `potential_gradient`, so the handle regenerates -- otherwise the
# field changes and the physics does not.

"""
    CahnHilliardModel(; kappa, zeta_0, zeta_1=0.0,
                        mobility=(u -> u .* (1 .- u)),
                        potential_gradient=(u -> u .^ 2 .* (u .- 1.5 * zeta_1)))

Biofilm cohesion submodel. `potential_gradient` defaults to the published
dpsi/du = u^2 (u - 3 zeta_1 / 2), matching the handles stored in the legacy
`.mat` models. With `zeta_1 = 0` it degenerates to `u^3`.
"""
Base.@kwdef struct CahnHilliardModel
    kappa::Float64
    zeta_0::Float64
    zeta_1::Float64 = 0.0
    mobility::Function = u -> u .* (1 .- u)
    potential_gradient::Function = u -> u .^ 2 .* (u .- 1.5 * zeta_1)
end
