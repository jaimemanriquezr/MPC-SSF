# MPC-SSF
MultiPhase Continuum simulations of Slow Sand Filtration

## Example

```matlab
initpath;
import MPC.*

filter = SandFilter();
filter.Temperature = 25; % in Celsius

% add grid points
filter = filter.addGridPoints(200);

% load a preset
model = load("./data/LundMPCModel.mat").model;

cleanFilter = State(filter, model, preset="clean");

inflow = [1e-4, 2e-4, 0, 0, 1e-5, 1e-6, 1e-9, 0, 1e-5];
results = simulate(cleanFilter, InflowConcentrations=inflow, SimulationTime=10.0);
```
## Creating an ecological model
```matlab
initpath;
import MPC.*

model = Model();

%% add components
HET = ecological.Particle(Name="HET", Density=1.0)
POM = ecological.Particle(Name="POM", Density=2.0)
model.Components = [HET, POM]


%% add reactions
mu = 1.0;
sigma = dictionary("HET", -1.0, "POM", 0.5)
order = dictionary("HET", 1)
% 2HET -> POM
% with a rate mu*HET
inactivationRx = ecological.Reaction(...
        NominalRate=1.0, ...
        StoichiometricCoefficients=sigma, ..
        Order = order)

model.Reactions = [inactivationRx];
```
