model = Model();
model.WaterDensity = 1.0;

P = Particle(Name="Microorganism");
P.AttachmentMatrix = 10.0;
P.AttachmentSand = 1.0;
P.Density = 1.1;

P0 = P;
P0.Name = "POM";

L = Liquid(Name="Nutrient");
L.Density = 1.0;

P_death = Reaction(Name="Inactivation", Order = dictionary(P, 1), ...
                   StoichiometricCoefficients=dictionary(P, -1.0, P0, 1.0));
P_growth = Reaction(Name="Growth", Order = dictionary(P, 1), ...
                   HalfSaturationConstants = dictionary(L, 1e-03), ...
                   StoichiometricCoefficients=dictionary(P, 1.0, L, -1.0));
Hydrolysis = Reaction(Name="Hydrolysis", Order = dictionary(P, 1), ...
                    HalfSaturationConstants=dictionary("Microorganism/POM", 2e-05), ...
                    StoichiometricCoefficients=dictionary(P0, -1.0, L, 1.0));

model.Components = [P, P0, L];
model.Reactions = [P_death, P_growth, Hydrolysis];
plot(model);