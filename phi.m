function [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM)

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

tracer_0 = bgc.tracer;
x0_bgc = bgc2nsoli(sim, bgc.tracer);    % unitless start of year values

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

        x1_bgc = bgc2nsoli(sim, bgc.tracer);    % unitless end of year values

        numWaterParcels = numel(sim.domain.iwet_JJ);
        numTracers = sim.bgc_struct_base.size.tracer(2);
        sz = [numWaterParcels, numTracers];

        x0 = reshape(x0_bgc, sz);
        x0 = x0(:,sim.selection);               % just selected cols

        tmpG = reshape(x1_bgc -x0_bgc, sz);     % needed G size is sz; aka 32 col
        tmpG = tmpG(:,sim.selection);           % just selected cols
        tmpG = tmpG(:);                         % nsoli format

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
        %         toc(timer_total);
        %         save the entire workspace. Surprisingly slow. Perhaps v7.3 compress of so much data is slow...
        %         allFile = sprintf('%s/all_%d.mat', sim.outputRestartDir, round(1+years_gone_by));
        %         fprintf('%s.m: Saving "%s"...\n', mfilename, allFile);
        %         save(allFile,'-v7.3');
        %
        % Matlab load() has trouble with filenames that space and so on.
        % KISS
        myRestartFile = sprintf('%s/restart_%d_%s_x1.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));
        fprintf('%s.m: Saving "%s"...\n', mfilename,myRestartFile);
        % copy original restart file, then replace original "tracer" with
        % the current bgc.tracer. Surprisingly fast!
        copyfile( sim.inputRestartFile, myRestartFile);

        tracer = bgc.tracer;                            % --- these are x1 everywhere  ---
        save( myRestartFile, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init



        myRestartFile = sprintf('%s/restart_%d_%s_x0.mat', sim.outputRestartDir, round(sim.start_yr+years_gone_by),strjoin(tName(sim.selection)));
        fprintf('%s.m: Saving "%s"...\n', mfilename,myRestartFile);
        % copy original restart file, then replace original "tracer" with
        % the current bgc.tracer. Surprisingly fast!
        copyfile( sim.inputRestartFile, myRestartFile);

        tracer = tracer_0;                            % --- these are x0 everywhere  ---
        tracer(:,:,sim.selection) = bgc.tracer(:,:,sim.selection);   % only selected x1
        save( myRestartFile, 'tracer',  '-append' );    % overwrites tracer ONLY, keep forcing from init
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