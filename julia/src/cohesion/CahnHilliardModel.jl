# Port of src/cohesion/@CahnHilliardModel/CahnHilliardModel.m
#
# Cohesion submodel for the biofilm: supplies the chemical-potential gradient,
# the mobility function, and the parameter zeta_0 used to assemble the
# Cahn-Hilliard system ("SOLVER A") in simulate.
#
# In MATLAB these are function_handle fields; in Julia they are stored as
# callables (functions / closures).
# TODO: confirm defaults and signatures against the MATLAB class + presets.

"""
    CahnHilliardModel(; potential_gradient, mobility, zeta_0)

Cohesion submodel. `potential_gradient(u)` returns dψ/du, `mobility(u)` the
mobility, and `zeta_0` the cohesion coefficient.
"""
Base.@kwdef struct CahnHilliardModel
    potential_gradient::Function
    mobility::Function
    zeta_0::Float64
end
