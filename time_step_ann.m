function [sim, bgc, time_series, n] = time_step_ann (sim, bgc, time_series, n, forcing, TR, myMonth)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%

% persistent tracerCapture tendencyCapture
% if (sim.captureAllSelectedTracers == 1)
%     if isempty(tendencyCapture)
%         tmp = packMarbl( bgc.tracer, sim.domain.iwet_JJ );
%         tmp = tmp(:,sim.selection);
%         myCaptureSz = [size(tmp(:),1), sim.num_time_steps];
%         tracerCapture = zeros(myCaptureSz);
%         fprintf('\n%s.m: HUGE tendency log created: %g bytes\n\n',mfilename, getMemSize(tracerCapture));
%         tendencyCapture = zeros(myCaptureSz);
%         fprintf('\n%s.m: HUGE tracer log created: %g bytes\n\n',mfilename, getMemSize(tendencyCapture));
%     end
% end

%%
dt               = sim.dt;
tsperDay         = round(sim.const.sec_d/dt);
days_in_month    = [31 28 31 30 31 30 31 31 30 31 30 31];
steps_per_period = tsperDay*days_in_month(myMonth);

steps_between_plots = 1e+10;
%    steps_between_plots = 8*7;
%    steps_between_plots = 8*1;
%%

% DEBUG: n = n+steps_per_period; return; % quick return to debug G

% forcing, river, and transport created monthly. We get only this myMonth's
%  in this routine

bgc.forcing      = forcing.interior;
bgc.surf_forcing = forcing.surf_forcing;

% river flux is not a tendency until divided by surface layer thickness
% *** IN CM ***
% Yes! Units are mixed m and cm because MARBL depths are cm.

surfaceLayerThickness_MARBL = (sim.domain.dzt(1) *sim.domain.MARBL_depth_per_m);

river_tendency = forcing.river_flux /surfaceLayerThickness_MARBL;

%%

C_1 = packMarbl( bgc.tracer, sim.domain.iwet_JJ );
[C_2, C_3, C_4] = deal( zeros( size(C_1)));

%%
% Advection Adams-Bashforth method
%
% source/sink term
% S = bgc.tendency + bgc.river_flux/thickness;
%
% LHS is implicit diffusion
% RHS is explitc advection
%
%     LHS  = d0(1+xi(:,2)) + dt.*(-TR.D);   % -CONSTANT- same LHS on all tracers
%     FLHS = mfactor(LHS);
%     RHS = C(:,k-1) + dt.*(A*w1(:,1) + H*C(:,k-1)+ S);
%
% Solve for C: (I - dt.*D)*C(:,k) = RHS
%
%     C(:,k) = mfactor(FLHS,RHS);           % aka time stepped tracers

if n<0
    n = 0;
    returnAfterOneTimeStep = true;
else
    returnAfterOneTimeStep = false;
    A    = TR.A;
    H    = TR.H;
    % tau  = 1/sim.const.sec_y/1e+6;            % 1/1e+6 years (s)
    % msk = sim.domain.M3d;
    % R    = d0( msk(sim.domain.iwet_JJ ) / tau );   % (1/sec)
    LHS  = d0(1+dt*TR.dxidt) + dt.*(-TR.D);     % LHS on all tracers
    FLHS = mfactor(LHS);    % run time = 0.44 sec for 3x3 any amount of tracers
end


for it = 1:steps_per_period

    % This "short circuit is to calc Jacobian withOUT getting transport
    % involved...

    n = n+1;
    [~,bgc] = calculate_forcing(sim, bgc, n);

    k = 4; if (it <= 2), k = it+1; end

    % calculate tendency
    tic  % comment this tic out to get cumulative time

% keyboard
% save('MARBL_IC_3d.mat', 'dt', 'bgc', 'sim', 'forcing', 'time_series', 'surfaceLayerThickness_MARBL','river_tendency',  '-v7.3','-nocompression');
    if (sim.runInParallel)
        [bgc, time_series] = MARBL_loop_parallel (n, sim, bgc,time_series);
    else
        [bgc, time_series] = MARBL_loop          (n, sim, bgc,time_series);
    end

    % This "short circuit is to calc Jacobian withOUT getting transport
    % involved...
    %
    % we have MARBL output, and that is ALL we want. No transport, etc.
    if returnAfterOneTimeStep
        break
    end

%     if (sim.captureAllSelectedTracers == 1)
% 
%         tmp = packMarbl( bgc.tracer, sim.domain.iwet_JJ );
%         tracerCapture  (:,n) = tmp(:,sim.selection);
% 
%         tmp = packMarbl( bgc.tendency, sim.domain.iwet_JJ );
%         tendencyCapture(:,n) = tmp(:,sim.selection);
% 
%     end
    % add river flux, as a tendency, to MARBL tendency

    bgc.tendency(:,1,:) = bgc.tendency(:,1,:) +river_tendency; % S = MARBL +river
    % bgc.tendency = 0 *bgc.tendency;                          % S = 0
    % bgc.tendency(:,1,:) = river_tendency;                    % S = 0 +river

    S  = packMarbl( bgc.tendency, sim.domain.iwet_JJ );

    %     toc
    %%
    % Captured initial value of tracers and tendency in MARBL_LOOP
    % Capture integrated values here
    if sim.logTracers
        time_series.moles (:, n) = global_moles(bgc.tracer          , sim);
        time_series.Dmoles(:, n) = 1E-3 *sum( S .* sim.domain.dVt_FP);
    end

    switch (it)
        case 1
            w1 = C_1;
            RHS = C_1 + dt.*(A*w1 + H*C_1+ S);
            C_2 = mfactor(FLHS,RHS);
            bgc.tracer = unpackMarbl( C_2, sim.domain.iwet_JJ, size(bgc.tracer));
            %             bgc.accumulate = bgc.accumulate + ( S ./C_1 ); % FIXME: time stepped tracer or input tracer

        case 2
            w1 = (1/2).*(3*C_2 - C_1 );
            RHS = C_2 + dt.*(A*w1 + H*C_2+ S);
            C_3 = mfactor(FLHS,RHS);
            bgc.tracer = unpackMarbl( C_3, sim.domain.iwet_JJ, size(bgc.tracer));
            %             bgc.accumulate = bgc.accumulate + ( S ./C_2 ); % FIXME: time stepped tracer or input tracer

        otherwise
            w1 = ( (1/12).*( 23*C_3 - 16*C_2 + 5*C_1));
            RHS = C_3 + dt.*(A*w1 + H*C_3+ S);
            C_4 = mfactor(FLHS,RHS);
            bgc.tracer = unpackMarbl( C_4, sim.domain.iwet_JJ,size(bgc.tracer));
            %             bgc.accumulate = bgc.accumulate + ( S ./C_3 ); % FIXME: time stepped tracer or input tracer
    end

    if sim.logTracers && ( n+1 == size(time_series.moles,2))    % capture value expected at start of next year.

        time_series.tracer(:, :, n+1) = squeeze(bgc.tracer(sim.time_series_loc,:,:));

        time_series.moles (:, n+1)    = global_moles(bgc.tracer , sim);

        % FIXME: need to call MARBL_loop to get SFO and tendency and diags EXPECTED at n+1
        % time_series.sfo           (:, n+1)    = 0; its already zero
        %         if (sim.logDiags)
        %             time_series.diag      (:, :, n+1) = interior.diag';
        %             time_series.surf_diag (:, n+1)    = surface.diag';
        %         end
        %         time_series.Dmoles(:, n) = 1E-3 *sum( S .* sim.domain.dVt_FP);

    end

    % Update temps storing previous time steps
    if (it > 2)
        C_1 = C_2; C_2 = C_3; C_3 = C_4;
    end

    %     disp(['Volume Integral: ',mfilename,'/test() end of Period #',num2str(n),' = ',num2str(time_series.moles(:, n)',12)])

    %     toc

    %%
    % Make occasional plots (slow) to monitor progress of sim and debug

    if mod(n, steps_between_plots) == 0
        if (sim.logTracers)
            small_plots(sim, time_series, n, sim.time_series_loc, sim.time_series_lvl);
            %         autoArrangeFigures(0,0);
        end
    end
    %     if returnAfterOneTimeStep
    %         break
    %     end
end % steps_per_period

% mpiprofile viewer

if k <4 % error
    %     keyboard
end

% if (sim.logTracers)
% %     small_plots(sim, time_series, n, sim.time_series_loc, sim.time_series_lvl);
% %     disp(['Volume Integral: ',mfilename,' end of Month #',num2str(myMonth),' = ',num2str(time_series.moles(:, n)',12)])
% end

% if (sim.captureAllSelectedTracers == 1 && myMonth == 12)
% 
%     myFileName = sprintf('%s/all_tracers_all_times.mat', sim.outputRestartDir);
%     fprintf('%s.m: Saving "%s"...\n', mfilename, myFileName);
%     tracer = tracerCapture;
%     % tmp_time_series = squeeze(time_series.tracer(sim.time_series_lvl,7,:));
%     % tmp = unpackMarbl(tracerCapture,sim.domain.iwet_JJ,[7881,60,730]);
%     % tmp_log = squeeze(tmp(sim.time_series_loc, sim.time_series_lvl,:));
%     % tst = [tmp_time_series tmp_log];
%     save( myFileName, 'tracer' ,'-v7.3','-nocompression');
% 
%     myFileName = sprintf('%s/all_tendency_all_times.mat', sim.outputRestartDir);
%     fprintf('%s.m: Saving "%s"...\n', mfilename, myFileName);
%     tendency = tendencyCapture;
%     save( myFileName, 'tendency' ,'-v7.3','-nocompression');
% 
% end

end % integration function
