filter = SandFilter();
filter = filter.addGridPoints(200);

model = load("./data/SimpleModel.mat").model;
model.CohesionSubModel.Kappa = 1e-2;
model.CohesionSubModel.Zeta0 = 1.0;

model.DetachmentFunction = @(v) sqrt(abs(v));
inflow = dictionary("Microorganism", 1e-2, "Nutrient", 3.0);

filterState = State(filter, model);
results = filterState.simulate(InflowConcentrations=inflow, SimulationTime=2.0);
