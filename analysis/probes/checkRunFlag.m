function checkRunFlag(r, tag)
% CHECKRUNFLAG  Refuse to write results from a run that did not complete.
%
% simulate() returns Flag = "OK" only when the time loop ran to completion. On
% "CLOGGED" / "BIOFILM" / "FLOWING" it aborts early, leaves the remaining frames
% as ZEROS, and still returns a results struct that looks entirely normal.
%
% That is not hypothetical. On 2026-08-24 the 20 d zeta_1 arms at MaxDt = 1e-5
% aborted with Flag = "BIOFILM" at t = 15.5 d of 20; every probe metric was then
% computed from an all-zero final frame (bedPeak 0.0000, bedDecline 1.6e306) and
% written over a good 3e-6 result under the same filename. Nothing in the output
% said the run had failed.
%
% So: check the flag, write nothing, and fail loudly. simulate() has already
% printed the time and cell of the failure to stdout, so the diagnosis is in the
% log -- what must not happen is a partial run silently replacing a complete one.
arguments
    r        % results struct from simulate(); not type-constrained, since
             % `arguments r struct` provokes a STRUCT-on-object warning
    tag (1,1) string
end
if r.Flag == "OK", return, end
error("probe:runFailed", ...
    "simulate returned Flag = ""%s"" for %s -- run did NOT complete. " + ...
    "No result written (a partial run must never overwrite a good one). " + ...
    "See the abort time and cell printed above.", r.Flag, tag);
end
