function [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

timer_total = tic;

% Adams-Bashforth code from Ann requires "monthly" units of time.
% Go ahead and make that multuples of 12; aka years.

total_months = 12* round(sim.num_time_steps *sim.dt /sim.const.sec_y);

fprintf('\n%s.m: Start integration of %s years: ',mfilename, int2str(total_months/12))
fprintf('%s\n', datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')))

% save initial state

tracer_0 = bgc.tracer;  % [7881, 60, 32]
x0_bgc = bgc2nsoli(sim, tracer_0);    % unitless start of year values
numWaterParcels = numel(sim.domain.iwet_JJ);
numTracers = sim.bgc_struct_base.size.tracer(2);
sz = [numWaterParcels, numTracers];
x0 = reshape(x0_bgc, sz);
% x0 = x0(:,sim.selection);



% % % Note: mean of tendency is -NOT- tendency of mean tracer, but it
% % % is not far off. Both mean of tracer and mean tendency are captured so
% % % error from this assumption can be determined.
% %
% % sim.average_tracer   = zeros(size(bgc.tracer  ));
% % sim.average_tendency = zeros(size(bgc.tendency));

%%%% Set up many agro loop indexes, initialize C1, and run

n = 0;
current_month = 0;
years_gone_by = -1; % allow for case of total_months==0, etc etc


% DEBUG
if sim.debug_disable_phi
    fprintf('\n\n\n%s.m: ********* phi() is short circuited to return x0 as x1 for debugging  *********\n\n\n',mfilename)
    %     return
end


initial_moles = global_moles(tracer_0, sim);
while current_month < total_months
    current_month = current_month+1;
    years_gone_by = floor((current_month-1)/12);
    myMonth = mod(current_month-1,12)+1;
    current_yr = round(sim.start_yr+years_gone_by);

    if sim.debug_disable_phi
        % avoid millions of lines of text
        %         fprintf('\n\n\n%s.m: ********* phi() is short circuited to return x0 as x1 for debugging  *********\n\n\n',mfilename)
        %         return
    else
        [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, n, forcing(myMonth), MTM(myMonth), myMonth);
    end

    % IMPORTANT!
    % after calling time_step_ann(), "bgc" has "x1" not "x0"!

    if (sim.logTracers)
        %         toc(timer_total);
        %         fprintf('%s.m: Starting plots...\n',mfilename);
        small_plots(sim, time_series, n, sim.time_series_loc, sim.time_series_lvl);
        %         toc(timer_total);
    end


    if mod(current_month, 12) == 0    % This runs after last time step of every y
        final_moles = global_moles(bgc.tracer, sim);

        x1_bgc = bgc2nsoli(sim, bgc.tracer);    % unitless end of year values
        % x1 = reshape(x1_bgc, sz);
        % x1 = x1(:,sim.selection);
        G = reshape(x1_bgc -x0_bgc, sz); % needed G size is sz; aka 32 col
        [x0, G] = calcStats(x0, G, initial_moles, final_moles, sim.selection, current_yr, mfilename);

        %%%
        % save annual restart in a way that doesnt fill up disk...
        %
        % save restart file(s) every few years for DEBUG and recovery from crashes and so on
        %         if mod(current_yr, sim.yearsBetweenRestartFiles) == 0    % This runs after last time step of every 10 y
        % save "x0" or initial state file...
        %%%


        % save "x0" or initial state file...
        myRestartFile = sprintf('%s/restart_x0.mat', sim.outputRestartDir);
        [sim, bgc] = saveRestartFiles(sim, bgc, tracer_0, myRestartFile);

        % save "x1" or final state file...
        myRestartFile = sprintf('%s/restart_x1.mat', sim.outputRestartDir);
        [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile);

        %%%

        %         x0_bgc = x1_bgc;        % update for next year

        %%%
        % fprintf('%s.m: Finished integration of year %s: ',mfilename, int2str(sim.start_yr+years_gone_by))
        % fprintf('%s, rate = %s hr/sim_y, norm(%s,2) = %1.10g\n', ...
        %     datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')), ...
        %     num2str((toc(timer_total)/3600/(current_month/12)),'%.2f'), ...
        %     gStr, norm(tmpG));
    end


    % %     %     tic; tend_log.tendency(n,1:prod(size(bgc.tendency))) = bgc.tendency(:)'; toc
    % %     %     tic; A_log.A          (n,1:prod(size(bgc.A)))        = bgc.A(:)'       ; toc
    % %     sim.average_tracer   = sim.average_tracer   + bgc.tracer   /nsteps;
    % %     sim.average_tendency = sim.average_tendency + bgc.tendency /nsteps;

end % loop over time steps

fprintf('%s.m: Finish integration of %s years: %s\n',mfilename,num2str(1+years_gone_by,2),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
fprintf('%s.m: Total Runtime: %s min, rate = %s hr/sim_y\n', mfilename, ...
    num2str(toc(timer_total)/60,'%.2f'), ...
    num2str((toc(timer_total)/3600/(current_month/12)),'%.2f'));

end % of function