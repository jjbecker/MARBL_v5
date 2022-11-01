function [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% DEBUG
if sim.debug_disable_phi
    fprintf('\n\n\n%s.m: ********* phi() is short circuited to return x0 or -x0 as x1 for debugging  *********\n\n\n',mfilename)
    %     return
end

tName = tracer_names(0);    % no CISO tracers
tendStr   = strjoin(tName(sim.selection));
phiStr = sprintf('phi( %s )', tendStr);
persistent phiFileCnt %x0_prev
if isempty(phiFileCnt)
    phiFileCnt = 1;
else
    phiFileCnt = phiFileCnt +1;
end
fprintf('\ncall #%d to %s', phiFileCnt, phiStr);

timer_total = tic;

% Adams-Bashforth code from Ann requires "monthly" units of time.
% Go ahead and make that multuples of 12; aka years.

total_months = 12* round(sim.num_time_steps *sim.dt /sim.const.sec_y);
if (total_months ~= 12)
    keyboard
end
fprintf('\n%s.m: Start  integration of %s years: ',mfilename, int2str(total_months/12))
fprintf('%s\n', datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')))

% save initial state

bgc_0 = bgc;

numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);
sz = [numWaterParcels, numTracers];

x0_bgc = bgc2nsoli(sim, bgc_0.tracer);
x0 = reshape(x0_bgc, sz);
% x0 = x0(:,sim.selection);

%%%% Set up many agro loop indexes, initialize C1, and run

% allow for case of total_months==0, etc etc
n = 0;
current_month = 0;
years_gone_by = -1; 

initial_moles = global_moles(bgc_0.tracer, sim);
while current_month < total_months
    current_month = current_month+1;
    myMonth = mod(current_month-1,12)+1;        % 1 to 12 if multi year
    years_gone_by = floor((current_month-1)/12);
    current_yr = round(sim.start_yr+years_gone_by);

    if sim.debug_disable_phi
        if ( phiFileCnt<=3 || phiFileCnt>=5 )
            bgc.tracer = (1. -   phiFileCnt/10) +0*bgc.tracer;  % normal case phi decrease from previous value
        else
            bgc.tracer = (1. - 0*phiFileCnt/10) +0*bgc.tracer;  % trigger Armijo steps
        end
        fprintf('%s.m: Finish integration of %s years: %s\n',mfilename,num2str(1+years_gone_by,2),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
        return
    else

        [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, n, forcing(myMonth), MTM(myMonth), myMonth);

    end

    % IMPORTANT!
    % after calling time_step_ann(), "bgc" has "x1" not "x0"!

    if (sim.logTracers)
        small_plots(sim, time_series, n, sim.time_series_loc, sim.time_series_lvl);
    end


    if mod(current_month, 12) == 0    % This runs after last time step of every y

        final_moles = global_moles(bgc.tracer, sim);

        x1_bgc = bgc2nsoli(sim, bgc.tracer);    % unitless end of year values
        x1 = reshape(x1_bgc, sz);
        % x1 = x1(:,sim.selection);

        G = reshape(x1_bgc -x0_bgc, sz); % needed G size is sz; aka 32 col

        [x0, G] = calcStats(x0, G, initial_moles, final_moles, sim.selection, current_yr, mfilename);

        % save annual restart in a way that doesnt fill up disk...
        %%%

        % save "x0" or initial state file...
        myRestartFile = sprintf('%s/restart_x0.mat', sim.outputRestartDir);
        [sim, bgc_0] = saveRestartFiles(sim, bgc_0, bgc_0.tracer, myRestartFile);

        % save "x1" or final state file...
        myRestartFile = sprintf('%s/restart_x1.mat', sim.outputRestartDir);
        [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile);

        %%%

        % update for next year, if multiple year run
        x0 = x1;
        x0_bgc = x1_bgc;
        bgc_0 = bgc;

    end
end % loop over time steps

fprintf('%s.m: Finish integration of %s years: %s\n',mfilename,num2str(1+years_gone_by,2),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
fprintf('%s.m: Total Runtime: %s min, rate = %s hr/sim_y\n', mfilename, ...
    num2str(toc(timer_total)/60,'%.2f'), ...
    num2str((toc(timer_total)/3600/(current_month/12)),'%.2f'));

end % of function