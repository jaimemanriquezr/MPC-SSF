function probeLimitDiag(stack, opts)
% Report which substrate binds the Liebig min, per reaction and region.
% stack: "manuscript" | "corrected"
arguments
    stack (1,1) string {mustBeMember(stack, ["manuscript","corrected"])}
    opts.NCells (1,1) double = 30
    opts.Tsim   (1,1) double = 2.0
end
here = fileparts(mfilename("fullpath"));
W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis"));

m = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
rx = m.Reactions; nm = [rx.Name];
if stack == "corrected"
    iH = nm == "Heterotroph growth";  iP = nm == "Phototroph growth";
    rx(iH).NominalRate = 0.042*24;    rx(iP).NominalRate = 0.125*24;
    rx(iH).HalfSaturationConstants = dictionary(["O2","DOM","NH4","HPO4"], [2.0e-4, 4.0e-3, 1.0e-6, 2.0e-5]);
    rx(iP).HalfSaturationConstants = dictionary(["IC","NH4","HPO4"],       [1.2e-3, 2.0e-5, 2.0e-5]);
end
m = Model(m.Components, rx, Kappa=1e-7, Zeta0=1e2, Zeta1=1e-2, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
    BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);

f = SandFilter(Temperature=19); f = f.addGridPoints(opts.NCells);
f.LightIrradiation = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

r = simulate(State(f,m), InflowConcentrations=infl, SimulationTime=opts.Tsim, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, ...
    FrameNumber=24, ImplicitOsmosis=true, Quiet=true, RecordLimitation=true);

L = r.Frames.Limitation;  z = f.GridPoints.Centers(:);
sup = z < 0;  top = z >= 0 & z <= 0.02;
fprintf('\n===== %s stack (%d cells, %.1f d) flag=%s =====\n', stack, opts.NCells, opts.Tsim, r.Flag);
for reac = ["Heterotroph growth","Phototroph growth"]
    j = find(L.ReactionNames == reac);
    for reg = ["supernatant","top 0-2cm"]
        mask = sup; if reg == "top 0-2cm", mask = top; end
        for phase = ["Biofilm","Flowing"]
            v = L.(phase)(mask, 2:end, j); v = v(:); v = v(v > 0);
            mo = L.("Monod"+phase)(mask, 2:end, j); mo = mo(:);
            if isempty(v), continue; end
            u = unique(v); cnt = arrayfun(@(x) sum(v==x), u);
            [~,o] = sort(cnt,'descend');
            parts = arrayfun(@(k) sprintf('%s %.0f%%', L.Names(u(o(k))), 100*cnt(o(k))/numel(v)), ...
                             1:min(3,numel(u)), UniformOutput=false);
            fprintf('  %-18s %-11s %-8s : %s  | mean min-Monod %.4g\n', reac, reg, phase, strjoin(parts,', '), mean(mo,'omitnan'));
        end
    end
end
la = L.LightAttenuated;
fprintf('  light: surface Ihat max %.3g | at z=0 max %.3g | first bed cell max %.3g\n', ...
    max(la(1,:)), max(la(find(z>=0,1),:)), max(la(find(z>0,1)+1,:)));
end
