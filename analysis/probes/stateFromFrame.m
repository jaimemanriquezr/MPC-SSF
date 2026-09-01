function st = stateFromFrame(r, k)
% STATEFROMFRAME  Snapshot of frame k as a re-homeable struct.
%
% Extracted verbatim (behaviour-preserving) from the local `frameState` in
% analysis/manuscriptExperiments.m so that long runs can be CHAINED: run 10 d,
% snapshot, restart from the snapshot, run the next 10 d. Long jobs then yield
% partial results instead of all-or-nothing.
%
% WARNING carried over from the original: `grab` replaces a scalar entry with
% ZEROS. That is the silent-data-loss path -- if any (component, volume) pair is
% ever stored as a scalar, the restart starts from zero for that field with no
% error. probeRestart.m exists to detect exactly that, by differencing a chained
% run against a single continuous one.
arguments
    r            % results struct; not type-constrained (STRUCT-on-object warning)
    k (1,1) double
end
C = r.Frames.Concentrations;
pNames = [r.Model.Particles.Name];  lNames = [r.Model.Liquids.Name];
st.pNames = pNames;  st.lNames = lNames;
st.Matrix            = grab(C, pNames, "Matrix", k);
st.EnclosedParticles = grab(C, pNames, "Enclosed", k);
st.FlowingParticles  = grab(C, pNames, "Flowing", k);
st.EnclosedLiquids   = grab(C, lNames, "Enclosed", k);
st.FlowingLiquids    = grab(C, lNames, "Flowing", k);
densityL = mean([r.Model.Liquids.Density]);
st.EnclosedWaterVolume = C{"Water", "Enclosed"}{1}(:, k)/densityL;
st.VelocityBiofilm = r.Frames.Velocity.Biofilm(:, k);
st.Time = r.Frames.Time(k);
% Record which fields came back all-zero, so a caller can assert on it rather
% than discover it three runs later.
st.allZero = string.empty;
for fn = ["Matrix","EnclosedParticles","FlowingParticles","EnclosedLiquids","FlowingLiquids","EnclosedWaterVolume"]
    if all(st.(fn)(:) == 0), st.allZero(end+1) = fn; end
end
end

function M = grab(C, names, vol, k)
M = [];
for nm = names
    v = C{nm, vol}{1};
    if isscalar(v)
        warning("stateFromFrame:scalarField", ...
            "component %s in volume %s is stored as a scalar; restart will zero it.", nm, vol);
        v = zeros(size(C{names(1), "Enclosed"}{1}, 1), size(C{names(1), "Enclosed"}{1}, 2));
    end
    M = [M, v(:, k)]; %#ok<AGROW>
end
end
