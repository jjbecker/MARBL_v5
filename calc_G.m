function [r,G,x1] = calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv)
%UNTITLED Take an initial value of tracers, return change at end of a year
%   Detailed explanation goes here

Npt = -123;
tName = tracer_names(0);    % no CISO tracers
% selection = [ ...
%     find( strcmp(tName,'SiO3') ) ];     % #3
tendStr   = strjoin(tName(sim.selection));
gStr = sprintf('G( %s )', tendStr);
% fprintf('%s.m: %s...\n',mfilename, gStr);

persistent gFileCnt x0_prev
if isempty(gFileCnt)
    gFileCnt = 1;
    fprintf('\ncall #%d to %s\n', gFileCnt, gStr);
    fprintf('norm(x0         ) = %f\n', norm(x0         ));
    % checkNegAndHisto(sim, x0, 100.0, 'x0', 900);
    %     figure (500); plot(x0); title('x0')
else
    gFileCnt = gFileCnt +1;
    fprintf('\ncall #%d to %s\n', gFileCnt, gStr);
    fprintf('norm(x0_prev    ) = %f\n', norm(x0_prev    ));
    fprintf('norm(x0         ) = %f\n', norm(x0         ));
    dx0 = x0 -x0_prev;
    fprintf('norm(x0 -x0_prev) = %.6f\n', norm(dx0) );
    %     figure (501); plot(x0_prev); title('x0_prev', 'Interpreter', 'none')
    %     figure (502); plot(x0)     ; title('x0')
    figure (503); plot(dx0)    ; title('dx0'); xlabel('idx FP'); ylabel(strjoin(tName(sim.selection)));
end

% Check for negative tracers

% accept = 1;
% c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
% [~,c0] = ck_constraints(c0,-2);     % corrected any negatives
%
% tmp = x0;
% [accept,x0] = ck_constraints(x0,Npt);
% if (accept == false)
%     disp('accept = false; histo our input in fig(999)')
%     figure(999); histogram(tmp)
%     keyboard
%     r = NaN(size(x0));
%     return
% end % if not within constraints


% Combine current tracers being changed by Nsoli, aka "sim.selection",
% with ones not being fiddled with to create full set of tracer to input
% to phi() aka MARBL.

numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);
sz = [numWaterParcels, numTracers];

% "_bgc" in tracer name means it has all 32 tracers

x0_bgc = replaceSelectedTracers(sim, c0, x0, sim.selection);
bgc.tracer = nsoli2bgc(sim, bgc, x0_bgc);   % marbl format

initial_moles = global_moles(bgc.tracer, sim);  % DEBUG
[sim, bgc, ~] = phi(sim, bgc, time_series, forcing, MTM);
final_moles = global_moles(bgc.tracer, sim);    % DEBUG

x1_bgc = bgc2nsoli(sim, bgc.tracer); % unitless end of year values
% checkNegAndHisto(sim, selectedTracers(sim, x1_bgc, sim.selection), 100.0, 'x', 900+gFileCnt);
% x1 = reshape(x1_bgc, sz);
% x1 = x1(:,sim.selection);

G = reshape(x1_bgc -x0_bgc, sz);    % x1 -x0 = phi(x0) -x0
G = G(:,sim.selection);             % just selected cols
G = G(:);                           % nsoli format
fprintf('||G(x)|| = (max(abs(%s))) = %g \n', gStr, max(abs(G)));
x1 = x0 +G;


% Precondition the residual

r = mfactor(PQ_inv, G) - G;

fprintf('||Precon( %s )|| = (max(abs(r))) = %g \n', gStr, max(abs(r)));



% DEBUG
disp([mfilename,'.m: Moles  start of phi() = ',num2str(initial_moles,7)])
disp([mfilename,'.m: Moles  end of phi()   = ',num2str(final_moles,7)])
disp([mfilename,'.m: Moles  delta          = ',num2str(final_moles-initial_moles,7)])
ppm = ((final_moles-initial_moles)./ final_moles *1e6);
disp([mfilename,'.m: Moles  delta (ppm)    = ',num2str(ppm,7)])
tmp = replaceSelectedTracers(sim, c0, G, sim.selection);
res_moles = global_moles(nsoli2bgc(sim, bgc, tmp), sim);
res_moles = res_moles(sim.selection);
% res_moles ./ final_moles(sim.selection) *1e6

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
fprintf('%s.m: Npt = %d %s norm(G,2)= %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), norm(G));
disp('  ')
fprintf('%s.m: Npt = %d %s median(G)= %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), median(G));
fprintf('%s.m: Npt = %d %s madMed(G)= %1.7g\n',  mfilename, Npt, strjoin(tName(sim.selection)), max((G)));
fprintf('%s.m: Npt = %d %s madAvg(G)= %1.10g\n', mfilename, Npt, strjoin(tName(sim.selection)), mad(G,1));
fprintf('%s.m: Npt = %d %s moles/y = %1.7g\n',   mfilename, Npt, gStr, res_moles);

% fprintf('%s.m: Npt = %d %s min     = %1.7g\n',  mfilename, Npt, gStr, min((G)));
% fprintf('%s.m: Npt = %d %s median  = %1.7g\n',  mfilename, Npt, gStr, median(G));
% fprintf('%s.m: Npt = %d %s mean    = %1.7g\n',  mfilename, Npt, gStr, mean(G));
% fprintf('%s.m: Npt = %d %s std     = %1.7g\n',  mfilename, Npt, gStr, std(G));
% fprintf('%s.m: Npt = %d %s mad(avg)= %1.7g\n',  mfilename, Npt, gStr, mad(G,0));
% fprintf('%s.m: Npt = %d %s mad(med)= %1.7g\n',  mfilename, Npt, gStr, mad(G,1));

if (0)
    % Matlab load() has trouble with filenames that space and so on.
    % KISS
    myGfile = sprintf('%s/G_%s_%d.mat', sim.outputRestartDir, strjoin(tName(sim.selection)), round(gFileCnt));
    fprintf('%s.m: Saving "%s"...\n', mfilename,myGfile);
    variables = who;
    toexclude = {'MTM','PQ_inv'};
    variables = variables(~ismember(variables, toexclude));
    save(myGfile, variables{:}, '-v7.3');
end

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


%     tracer = bgc.tracer;
%     copyfile( sim.inputRestartFile, myGfile);
%     save( myGfile, 'tracer', '-append' ); % overwrites tracer from input
%     save( myGfile, 'r', '-append' );      % overwrites tracer from input
%     save( myGfile, 'x0', '-append' );      % overwrites tracer from input
%     save( myGfile, 'x0_prev', '-append' );      % overwrites tracer from input
%     save( myGfile, 'x0_prev', '-append' );      % overwrites tracer from input
% Save everything in G, which is not everything in the whole sim

x0_prev = x0;


end % G()
