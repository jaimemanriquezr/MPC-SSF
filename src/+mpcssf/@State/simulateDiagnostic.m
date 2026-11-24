function results = simulateDiagnostic(obj, varargin)
% SIMULATEDIAGNOSTIC  simulate with the research paths unlocked.
%
%   results = simulateDiagnostic(state, ...)   % identical call to simulate
%
% Until 2026-09-22 this was a frozen copy of simulate carrying ~650 duplicated
% lines. The two are
% now one implementation; this forwards with Research=true, which unlocks the
% "bailo" and "matched" cohesion schemes, TrackCFL and PositiveMobility.
%
% Kept as a named entry point because six probes and the four diag_* regression
% configs call it, and because "which runs are research?" is worth being able to
% grep for.
results = simulate(obj, varargin{:}, "Research", true);
end
