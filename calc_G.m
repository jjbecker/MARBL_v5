function [r,G,x1] = calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv)
%UNTITLED Take an initial value of tracers, return change at end of a year
%   Detailed explanation goes here

Npt = -123;

persistent gFileCnt x0_prev
if isempty(gFileCnt)
    gFileCnt = 1;
    fprintf('call #%d to G\n', gFileCnt);
    fprintf('norm(x0         ) = %f\n', norm(x0         ));
    % checkNegAndHisto(sim, x0, 100.0, 'x0', 900);
%     figure (500); plot(x0); title('x0')
else
    gFileCnt = gFileCnt +1;
    fprintf('call #%d to G\n', gFileCnt);
    fprintf('norm(x0_prev    ) = %f\n', norm(x0_prev    ));
    fprintf('norm(x0         ) = %f\n', norm(x0         ));
    dx0 = x0 -x0_prev;
    fprintf('norm(x0 -x0_prev) = %.6f\n', norm(dx0) );
%     figure (501); plot(x0_prev); title('x0_prev', 'Interpreter', 'none')
%     figure (502); plot(x0)     ; title('x0')
    figure (503); plot(dx0)    ; title('dx0')
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
fprintf('||G(x)|| = (max(abs(G))) = %g \n', max(abs(G)));
x1 = x0 +G;


% Precondition the residual

r = mfactor(PQ_inv, G) - G;

fprintf('||Precon( G(x) )|| = (max(abs(r))) = %g \n', max(abs(r)));



% DEBUG
disp([mfilename,'.m: Moles  start of phi() = ',num2str(initial_moles,7)])
disp([mfilename,'.m: Moles  end of phi()   = ',num2str(final_moles,7)])
disp([mfilename,'.m: Moles  delta          = ',num2str(final_moles-initial_moles,7)])
ppm = ((final_moles-initial_moles)./ final_moles *1e6);
disp([mfilename,'.m: Moles  delta (ppm)    = ',num2str(ppm,7)])

fprintf('%s.m: Npt = %d G norm    = %1.10g\n', mfilename, Npt, norm(G));
fprintf('%s.m: Npt = %d G max     = %1.7g\n',  mfilename, Npt, max((G)));
fprintf('%s.m: Npt = %d G min     = %1.7g\n',  mfilename, Npt, min((G)));
fprintf('%s.m: Npt = %d G median  = %1.7g\n',  mfilename, Npt, median(G));
fprintf('%s.m: Npt = %d G mean    = %1.7g\n',  mfilename, Npt, mean(G));
fprintf('%s.m: Npt = %d G std     = %1.7g\n',  mfilename, Npt, std(G));
fprintf('%s.m: Npt = %d G mad(avg)= %1.7g\n',  mfilename, Npt, mad(G,0));
fprintf('%s.m: Npt = %d G mad(med)= %1.7g\n',  mfilename, Npt, mad(G,1));
tmp = replaceSelectedTracers(sim, c0, G, sim.selection);
res_moles = global_moles(nsoli2bgc(sim, bgc, tmp), sim);
res_moles = res_moles(sim.selection);
% res_moles ./ final_moles(sim.selection) *1e6
fprintf('%s.m: Npt = %d G moles   = %1.7g\n',  mfilename, Npt, res_moles);

if (0)
    myGfile = sprintf('%s/G_%d.mat', sim.outputRestartDir, round(gFileCnt));
    fprintf('%s.m: Saving "%s"...\n', mfilename,myGfile);
    variables = who;
    toexclude = {'MTM','PQ_inv'};
    variables = variables(~ismember(variables, toexclude));
    save(myGfile, variables{:}, '-v7.3');
end

    figure (900); scatter(x0,r); title("scatter(x0,r)");    xlabel('x0');   ylabel('r')
    figure (901); plot(r);       title("plot(r)");          xlabel('idx FP');ylabel('r')
    figure (902); qqplot(r);     title("qqplot(r)")
    figure (602); histogram(r);  title("histogram(r)");     xlabel('r');   ylabel('Count')

    myRng = 1:20;
    [maxR,idxMaxR]     = sort(   (G),"descend",'MissingPlacement','last');
    fprintf("max(G) %g\n", maxR(1));
    [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(idxMaxR(myRng), sim, 996);  title('Most postive')

    [minR,idxMinR]     = sort(   (G),"ascend",'MissingPlacement','last');
    fprintf("min(G) %g\n", minR(1));
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
