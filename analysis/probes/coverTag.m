function tag = coverTag(cfg, opts)
% COVERTAG  Deterministic filename-safe tag for one cover-sweep configuration.
%
%   tag = coverTag(cfg)                  full tag, unique per configuration
%   tag = coverTag(cfg, OmitCover=true)  the pair key: same tag minus the cover
%                                        field, so the two arms of a paired
%                                        contrast share it
%
% Every varied field appears in the tag. That is the whole point: a stale .mat
% from an edited registry can then never be silently reused, because the tag it
% was saved under no longer exists. collectCover additionally asserts that the
% cfg stored inside the .mat is isequal to the registry entry.
%
% Non-varying fields (those still at their registry default) are omitted, so the
% common case stays readable: "A_cov1p0_30d" rather than a 12-field string.
%
% See also COVERCASES, PROBECOVER, COLLECTCOVER.

arguments
    cfg (1,1) struct
    opts.OmitCover (1,1) logical = false
end

% Defaults; a field equal to its default is left out of the tag.
D = struct("phoIn", 1.00e-2, "nh4In", 2.00e-5, "hpo4In", 0.0, "icIn", 6.23e-3, ...
           "ncells", 100, "maxdt", 3e-6, "kinetics", "corrected", "etaSand", 1500);

% NOTE: `stage` is deliberately NOT part of the tag. A configuration is
% identified by its physics, not by which stage asked for it -- otherwise
% A_cov1_30d and B_cov1_30d are the same run computed twice (2 x 25 min). With
% stage excluded they collapse to one tag and coverCases dedups them.
% Consequence: collectCover's integrity assert compares the physics fields only,
% via coverConfigFields(), not the whole struct.
parts = strings(1,0);

if ~opts.OmitCover
    parts(end+1) = "cov" + num2sci(cfg.cover);
end

parts(end+1) = num2sci(cfg.tsim) + "d";

if cfg.phoIn  ~= D.phoIn,  parts(end+1) = "pho"  + num2sci(cfg.phoIn);  end
if cfg.nh4In  ~= D.nh4In,  parts(end+1) = "nh4"  + num2sci(cfg.nh4In);  end
if cfg.hpo4In ~= D.hpo4In, parts(end+1) = "hpo4" + num2sci(cfg.hpo4In); end
if cfg.icIn   ~= D.icIn,   parts(end+1) = "ic"   + num2sci(cfg.icIn);   end
if cfg.ncells ~= D.ncells, parts(end+1) = "n"    + string(cfg.ncells);  end
if cfg.maxdt  ~= D.maxdt,  parts(end+1) = "dt"   + num2sci(cfg.maxdt);  end
if cfg.etaSand ~= D.etaSand, parts(end+1) = "eta" + num2sci(cfg.etaSand); end
if cfg.kinetics ~= D.kinetics, parts(end+1) = cfg.kinetics;             end

tag = strjoin(parts, "_");
end

function s = num2sci(v)
% Compact, filename-safe, and injective over the values the registry uses.
% "." -> "p", "-" -> "m", so 1e-2 -> "1em2", 0.001 -> "1em3", 1.0 -> "1p0".
if v == 0
    s = "0"; return
end
if v >= 1e-3 && v < 1e4 && abs(v - round(v, 4)) < eps(v)*8
    s = string(sprintf("%g", v));
else
    s = string(sprintf("%.3g", v));
end
s = replace(s, ".", "p");
s = replace(s, "-", "m");
s = replace(s, "+", "");
end

