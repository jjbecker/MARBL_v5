function [sim, bgc, time_series, tracer_0] = phi(sim, bgc, time_series, forcing, MTM)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% The Adams-Bashforth code from Ann requires "monthly" units of time.
% Go ahead and make that multuples of 12; aka years.

total_months = 12* round(sim.num_time_steps *sim.dt /sim.const.sec_y);

fprintf('\n%s.m: Start integration of %s years: ',mfilename, int2str(total_months/12))
fprintf('%s\n', datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')))

timer_total = tic;

% % % Note: The mean of tendency is -NOT- tendency of the mean tracer, but it
% % % is not far off. Both mean of tracer and mean tendency are captured so
% % % error from this assumption can be determined.
% %
% % sim.average_tracer   = zeros(size(bgc.tracer  ));
% % sim.average_tendency = zeros(size(bgc.tendency));

%%%% Set up many agro loop indexes, initialize C1, and run

n = 0;
current_month = 0;
years_gone_by = -1; % allow for case of total_months==0, etc etc
tName = tracer_names(0);    % no CISO tracers

tracer_0 = bgc.tracer;
x0_bgc = bgc2nsoli(sim, bgc.tracer);    % unitless start of year values

% DEBUG
initial_moles = global_moles(bgc.tracer, sim);

while current_month < total_months
    current_month = current_month+1;
    years_gone_by = floor((current_month-1)/12);
    myMonth = mod(current_month-1,12)+1;
    %     fakeMonth = mod(6+current_month-1,12)+1;

    [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, n, forcing(myMonth), MTM(myMonth), myMonth);

    % IMPORTANT!
    % at this point "bgc" is "x1" not "x0"!

    if (sim.logTracers)
        %         toc(timer_total);
        %         fprintf('%s.m: Starting plots...\n',mfilename);
        small_plots(sim, time_series, n, sim.time_series_loc, sim.time_series_lvl);
        %         toc(timer_total);
    end

    %         fprintf('%s.m: Runtime to end of myMonth %d: %s min, rate = %s hr/sim_y\n', mfilename, current_month, ...
    %             num2str(toc(timer_total)/60,'%.2f'), ...
    %             num2str((toc(timer_total)/3600/(current_month/12)),'%.2f'));

    if mod(current_month, 12) == 0    % This runs after last time step of every y

        final_moles = global_moles(bgc.tracer, sim);
        x1_bgc = bgc2nsoli(sim, bgc.tracer);    % unitless end of year values

        numWaterParcels = numel(sim.domain.iwet_JJ);
        numTracers = sim.bgc_struct_base.size.tracer(2);
        sz = [numWaterParcels, numTracers];

        x0 = reshape(x0_bgc, sz);
        x0 = x0(:,sim.selection);               % just selected cols

        tmpG_all = reshape(x1_bgc -x0_bgc, sz); % needed G size is sz; aka 32 col
        tmpG = tmpG_all(:,sim.selection);       % just selected cols
        tmpG = tmpG(:);                         % nsoli format

        % DEBUG
        fprintf(        '%s.m: Year %d               %s\n',mfilename,round(sim.start_yr+years_gone_by),strjoin(pad(tName,14)));
        disp([mfilename,  '.m: Year ',num2str(round(sim.start_yr+years_gone_by)),' Moles start = ',num2str(initial_moles,'%-#15.7g')])
        disp([mfilename  ,'.m: Year ',num2str(round(sim.start_yr+years_gone_by)),' Moles end   = ',num2str(final_moles,'%-#15.7g')])
        disp([mfilename  ,'.m: Year ',num2str(round(sim.start_yr+years_gone_by)),' Moles delta = ',num2str(final_moles-initial_moles,'%-#15.7g')])
        ppm = ((final_moles-initial_moles)./ final_moles *1e6);
        disp([mfilename  ,'.m: Year ',num2str(round(sim.start_yr+years_gone_by)),' Moles (ppm) = ',num2str(ppm,'%-#15.7g')])
        disp([mfilename  ,'.m: Year ',num2str(round(sim.start_yr+years_gone_by)),' norm G      = ', num2str(max(abs(tmpG_all)),'%-#15.7g')])
        fprintf(        '%s.m: Year %d               %s\n',mfilename,round(sim.start_yr+years_gone_by),strjoin(pad(tName,14)));

        tName = tracer_names(0);    % no CISO tracers
        % selection = [ ...
        %     find( strcmp(tName,'SiO3') ) ];     % #3
        tendStr   = strjoin(tName(sim.selection));
        gStr = sprintf('G( %s )', tendStr);
        figure (700); scatter(x0,tmpG); title(strjoin(["scatter(",gStr,", ",strjoin(tName(sim.selection)),")"]));    xlabel(strjoin(tName(sim.selection)));   ylabel(gStr)
        figure (701); plot(tmpG);       title(strjoin(["plot(",gStr,")"]));         xlabel('idx FP');                        ylabel(gStr)
        figure (702); qqplot(tmpG);     title(strjoin(["qqplot(",gStr,")"]))
        figure (601); histogram(tmpG);  title(strjoin(["histogram(",gStr,")"]));    xlabel(gStr);                            ylabel('Count')


        x0_bgc = x1_bgc;        % update for next year

        fprintf('%s.m: Finished integration of year %s: ',mfilename, int2str(sim.start_yr+years_gone_by))
        fprintf('%s, rate = %s hr/sim_y, norm(%s,2) = %1.10g\n', ...
            datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')), ...
            num2str((toc(timer_total)/3600/(current_month/12)),'%.2f'), ...
            gStr, norm(tmpG));

    end

    if mod(current_month, 12*sim.yearsBetweenRestartFiles) == 0    % This runs after last time step of every 10 y
        [sim, bgc] = saveRestartFiles(sim, bgc, tracer_0, years_gone_by);

        %         myRestartFile = sprintf('%s/restart_%d_%s_x1.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));
        %         fprintf('%s.m: Saving "%s"...\n', mfilename,myRestartFile);
        %         % copy original restart file, then replace original "tracer" with
        %         % the current bgc.tracer. Surprisingly fast!
        %         copyfile( sim.inputRestartFile, myRestartFile);
        %
        %         tracer = bgc.tracer;                            % --- these are x1 everywhere  ---
        %         save( myRestartFile, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init
        %
        %
        %
        %         myRestartFile = sprintf('%s/restart_%d_%s_x0.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));
        %         fprintf('%s.m: Saving "%s"...\n', mfilename,myRestartFile);
        %         % copy original restart file, then replace original "tracer" with
        %         % the current bgc.tracer. Surprisingly fast!
        %         copyfile( sim.inputRestartFile, myRestartFile);
        %
        %         tracer = tracer_0;                            % --- these are x0 everywhere  ---
        %         tracer(:,:,sim.selection) = bgc.tracer(:,:,sim.selection);   % only selected x1
        %         save( myRestartFile, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init
    end

    % %     %     tic; tend_log.tendency(n,1:prod(size(bgc.tendency))) = bgc.tendency(:)'; toc
    % %     %     tic; A_log.A          (n,1:prod(size(bgc.A)))        = bgc.A(:)'       ; toc
    % %     sim.average_tracer   = sim.average_tracer   + bgc.tracer   /nsteps;
    % %     sim.average_tendency = sim.average_tendency + bgc.tendency /nsteps;

end % time steps

fprintf('%s.m: Finish integration of %s years: %s\n',mfilename,num2str(1+years_gone_by,2),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
fprintf('%s.m: Total Runtime: %s min, rate = %s hr/sim_y\n', mfilename, ...
    num2str(toc(timer_total)/60,'%.2f'), ...
    num2str((toc(timer_total)/3600/(current_month/12)),'%.2f'));

end