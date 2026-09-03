% Winter phi_b profiles: presets x influent phosphate, at 60 and 90 d.
% Gate + figure. Local plotting only.
R = fileparts(fileparts(mfilename("fullpath")));
R = fileparts(R);
addpath(genpath(fullfile(R,'src'))); addpath(fullfile(R,'analysis')); addpath(fullfile(R,'analysis','probes'));
D = fullfile(R,'analysis','probes','data','chain');

sets = { % file60, file90, label, style
  'chain_hpo4zero_winter_leg6.mat','',                        'OLD presets, HPO_4 = 0',     'k--'
  'chain_hpo4base_winter_leg6.mat','chain_fld2x_winter_leg9.mat','OLD presets, HPO_4 = 5e-6','k-'
  'chain_prod_win_p0_leg6.mat',   'chain_prod_win_p0_leg9.mat','NEW presets, HPO_4 = 0',     'r--'
  'chain_prod_win_p5_leg6.mat',   'chain_prod_win_p5_leg9.mat','NEW presets, HPO_4 = 5e-6',  'r-'
};

fprintf('=== GATE ===\n');
f = figure(Visible="off", Position=[0 0 1250 500], Color="w");
tl = tiledlayout(f,1,2,TileSpacing="compact",Padding="compact");
for col = 1:2
  ax = nexttile; hold(ax,'on'); grid(ax,'on');
  for s = 1:size(sets,1)
    fn = sets{s,col};
    if isempty(fn), continue; end
    Q = load(fullfile(D,fn),'results'); r = Q.results;
    z = r.SandFilter.GridPoints.Centers(:);
    pb = r.getVolumeFractions().Biofilm(:,end);
    if col==1
      fprintf('  %-28s Flag %-3s  t=%5.1f d  rows=%d  T=%g C\n', erase(fn,'.mat'), ...
        string(r.Flag), max(r.Frames.Time), numel(z), r.SandFilter.Temperature);
    end
    w = z>=0 & z<=0.30;
    plot(ax, z(w), pb(w), sets{s,4}, LineWidth=1.8, DisplayName=sets{s,3});
    [~,im] = max(pb(w)); zz = z(w);
    plot(ax, zz(im), max(pb(w)), 'o', MarkerSize=7, LineWidth=1.4, ...
         Color=get(gca().Children(1),'Color'), HandleVisibility='off');
  end
  xlabel(ax,'Depth z [m]'); ylabel(ax,'\phi_b [-]');
  title(ax, sprintf('winter, t = %d d', 60 + 30*(col-1)));
  xlim(ax,[0 0.30]); set(ax,FontSize=12);
  if col==1, legend(ax,Location='northeast',FontSize=10,Box='off'); end
end
title(tl,'Winter biofilm profile: presets \times influent phosphate (circles mark the maximum)','FontSize',13);
exportgraphics(f, fullfile(R,'analysis','results','figures','winter_bulge_2x2.png'), Resolution=150);
close(f);
fprintf('\nwrote analysis/results/figures/winter_bulge_2x2.png\n');
