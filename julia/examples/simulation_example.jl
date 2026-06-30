# Julia analogue of examples/simulation_example.m (target API once ported).
# This is the intended end-state usage — it will not run until the solver is
# implemented.

using MPCSSF

# Build a filter and discretize it.
filter = SandFilter(temperature = 25 + 273)
filter = addgridpoints(filter, 200)              # TODO: port addgridpoints

# Load / build a model (Phase 1: mirror presets/modelLund.m).
# model = modelLund()                            # TODO: port preset

# clean_filter = State(filter, model; preset = :clean)
# inflow = [1e-4, 2e-4, 0, 0, 1e-5, 1e-6, 1e-9, 0, 1e-5]
# results = simulate(clean_filter; inflow_concentrations = inflow,
#                    simulation_time = 10.0)
