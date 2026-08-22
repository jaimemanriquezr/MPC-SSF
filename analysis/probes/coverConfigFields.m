function names = coverConfigFields()
% COVERCONFIGFIELDS  Fields that define a cover-sweep configuration physically.
%
% `stage`, `tag` and `pairKey` are bookkeeping and are excluded, so the same
% physics requested by two different stages compares equal. collectCover uses
% this for its stale-artefact assertion.
%
% See also COVERCASES, COVERTAG, COLLECTCOVER.
names = ["cover","phoIn","nh4In","hpo4In","icIn","tsim","ncells","maxdt", ...
         "kinetics","etaSand"];
end
