function [r,G,x1] = calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv)
%UNTITLED Take an initial value of tracers, return change at end of a year
%   Detailed explanation goes here

tName = tracer_names(0);    % no CISO tracers
tendStr   = strjoin(tName(sim.selection));
gStr = sprintf('G( %s )', tendStr);

persistent gFileCnt x0_prev
if isempty(gFileCnt)
    gFileCnt = 1;
    fprintf('\ncall #%d to %s\n', gFileCnt, gStr);
    fprintf('norm(x0         ) = %-#15.7g\n', norm(x0         ));
else
    gFileCnt = gFileCnt +1;
    fprintf('\ncall #%d to %s\n', gFileCnt, gStr);
    fprintf('norm(x0_prev    ) = %-#15.7g\n', norm(x0_prev    ));
    fprintf('norm(x0         ) = %-#15.7g\n', norm(x0         ));
    dx0 = x0 -x0_prev;
    fprintf('norm(x0 -x0_prev) = %-#15.7g\n', norm(dx0) );
    figure (503); plot(dx0)    ; title('dx0'); xlabel('idx FP'); ylabel(strjoin(tName(sim.selection)));
end
Npt = -gFileCnt;

if sim.debug_disable_phi
    %         bgc.tracer = -bgc.tracer;
    x0 = 0 *x0;
end

% Combine current tracers being changed by Nsoli, aka "sim.selection",
% with ones not being fiddled with to create full set of tracer to input
% to phi() aka MARBL.

numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);
sz = [numWaterParcels, numTracers];

% "_bgc" in tracer name means it has all 32 tracers

x0_bgc = replaceSelectedTracers(sim, c0, x0, sim.selection);
bgc.tracer = nsoli2bgc(sim, bgc, x0_bgc);   % marbl format x0


% initial_moles = global_moles(bgc.tracer, sim);  % DEBUG
[sim, bgc, ~] = phi(sim, bgc, time_series, forcing, MTM);
% final_moles = global_moles(bgc.tracer, sim);    % DEBUG


x1_bgc = bgc2nsoli(sim, bgc.tracer); % unitless end of year values
% checkNegAndHisto(sim, selectedTracers(sim, x1_bgc, sim.selection), 100.0, 'x', 900+gFileCnt);

G = reshape(x1_bgc -x0_bgc, sz);
G = G(:,sim.selection);             % just selected cols
G = G(:);                           % nsoli format

% % x1 -x0 = phi(x0) -x0
% x1 = reshape(x1_bgc, sz);
% x1 = x1(:,sim.selection);
x1 = x0 +G;                         % x1 only of selection

% Precondition residual

if ~sim.debug_PQ_inv
    r = mfactor(PQ_inv, G) - G;
else
    fprintf('\n\n ********Disabling mfactor in calc_G*******\n\n');
    r = G;
end


fprintf('%s.m: Npt = %d %s norm(r,2)= %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), norm(r));
fprintf('%s.m: Npt = %d %s norm(G,2)= %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), norm(G));
if (sim.verbose_debug)
    disp('  ')
    fprintf('%s.m: Npt = %d %s max(x0)  = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), max(x0));
    fprintf('%s.m: Npt = %d %s min(x0)  = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), min(x0));
    fprintf('%s.m: Npt = %d %s mean(x0) = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), mean(x0));
    fprintf('%s.m: Npt = %d %s std(x0)  = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), std(x0));
    disp('  ')
    fprintf('%s.m: Npt = %d %s max(G)   = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), max(G));
    fprintf('%s.m: Npt = %d %s min(G)   = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), min(G));
    fprintf('%s.m: Npt = %d %s mean(G)  = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), mean(G));
    fprintf('%s.m: Npt = %d %s std(G)   = %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), std(G));
    disp('  ')
    fprintf('%s.m: Npt = %d %s median(G)= %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), median(G));
    fprintf('%s.m: Npt = %d %s madMed(G)= %1.7g\n',  mfilename, Npt, strjoin(tName(sim.selection)), max((G)));
    fprintf('%s.m: Npt = %d %s madAvg(G)= %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), mad(G,1));
    disp('  ')
    fprintf('||G(x)|| = (max(abs(%s))) = %g \n', gStr, max(abs(G)));
    fprintf('||Precon( %s )|| = (max(abs(r))) = %g \n', gStr, max(abs(r)));
    % % DEBUG
    tmp = replaceSelectedTracers(sim, c0, G, sim.selection);
    res_moles = global_moles(nsoli2bgc(sim, bgc, tmp), sim);
    res_moles = res_moles(sim.selection);
    fprintf('%s.m: Npt = %d %s moles/y = %1.7g\n',   mfilename, Npt, gStr, res_moles);

    figure (900); scatter(x0,r); title(strjoin(["scatter( r(",gStr,"), ",strjoin(tName(sim.selection)),")"]));    xlabel(strjoin(tName(sim.selection)));  ylabel(strjoin(["r(",gStr,")"]))
    figure (901); plot(r);       title(strjoin(["r(",gStr,")"]));          xlabel('idx FP');   ylabel(strjoin(["r(",gStr,")"]))
    figure (902); qqplot(r);     title(strjoin(["qqplot( r(",gStr,"))"]))
    figure (602); histogram(r);  title(strjoin(["histogram( r(",gStr,"))"]));     xlabel(strjoin(["r(",gStr,")"]));   ylabel('Count')

    myRng = 1:20;
    [maxR,idxMaxR]     = sort(   (G),"descend",'MissingPlacement','last');
    % fprintf("max(%s) %g\n", gStr, maxR(1));
    [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(idxMaxR(myRng), sim, 996);  title('Most postive')

    [minR,idxMinR]     = sort(   (G),"ascend",'MissingPlacement','last');
    % fprintf("min(%s) %g\n", gStr, minR(1));
    [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(idxMinR(myRng), sim, 997);  title('Most Negative')

    [maxAbsR,idxAbsR]  = sort(abs(G),"descend",'MissingPlacement','last');
    [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(idxAbsR(myRng), sim, 998); title('Largest Abs')
end

x0_prev = x0;

end % G()
