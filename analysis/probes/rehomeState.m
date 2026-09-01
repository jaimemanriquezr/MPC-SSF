function s = rehomeState(f, m, st, t0)
% REHOMESTATE  Rebuild a State from a stateFromFrame snapshot.
%
% Extracted verbatim from the local `rehome` in analysis/manuscriptExperiments.m.
%
% Velocity.Flowing is deliberately NOT restored: simulate.m recomputes it before
% the time loop from velBiofilm and phiBiofilm (see the velFlowing assignment
% just after `velBiofilm = startingConditions.Velocity.Biofilm`), so it is a
% derived quantity, not state. Restoring it would be harmless but redundant.
arguments
    f SandFilter
    m Model
    st struct
    t0 (1,1) double
end
s = State(f, m);
s.GlobalConcentration.Matrix            = st.Matrix;
s.GlobalConcentration.EnclosedParticles = st.EnclosedParticles;
s.GlobalConcentration.FlowingParticles  = st.FlowingParticles;
s.GlobalConcentration.EnclosedLiquids   = st.EnclosedLiquids;
s.GlobalConcentration.FlowingLiquids    = st.FlowingLiquids;
s.EnclosedWaterVolume = st.EnclosedWaterVolume;
s.Velocity.Biofilm    = st.VelocityBiofilm;
s.Time = t0;
end
